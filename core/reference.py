"""Build and locate the transcriptome reference the off-target gate searches.

Fetches the RefSeq H_sapiens mRNA_Prot RNA FASTA set, concatenates it, runs
bowtie-build, and writes a manifest recording exactly what went in (RefSeq
release, source URLs, sha256, bowtie version, build date) so a build is
auditable. core.offtarget.BowtieSearcher consumes what this produces.

Fetching the reference is a one-time, network-heavy setup step, deliberately
separate from designing: the download is multi-gigabyte and the index is not
redistributable in an image. HUSH is installed first and the reference fetched
afterwards, which is why this is a user-facing command (hush reference build)
rather than an ops-side script. Only the standard library is used.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

from core.errors import ReferenceBuildError

# RefSeq's current human mRNA/ncRNA transcripts. These files track the live
# RefSeq release; the release number is read separately and recorded in the
# manifest. Point base_url at an archived release directory to pin exactly.
DEFAULT_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot"
RELEASE_NUMBER_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"

# Where the index lives unless told otherwise. The container sets
# HUSH_REFERENCE_DIR to its mounted volume so the command needs no arguments.
REFERENCE_DIR_ENV = "HUSH_REFERENCE_DIR"
DEFAULT_REFERENCE_DIR = Path("data/reference")
INDEX_NAME = "refseq_human_mrna"


def reference_dir() -> Path:
    """The directory holding the bowtie index, from the environment or default."""
    configured = os.environ.get(REFERENCE_DIR_ENV)
    return Path(configured) if configured else DEFAULT_REFERENCE_DIR


def default_index_prefix() -> Path:
    """The bowtie index prefix HUSH uses when none is given explicitly."""
    return reference_dir() / INDEX_NAME


def index_exists(prefix: str | Path) -> bool:
    """Whether a usable bowtie index sits at this prefix."""
    return Path(f"{prefix}.1.ebwt").is_file()


def manifest_path(prefix: str | Path) -> Path:
    return Path(f"{prefix}.manifest.json")


def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def _read_release_number() -> str:
    try:
        with urllib.request.urlopen(RELEASE_NUMBER_URL, timeout=30) as resp:
            return resp.read().decode().strip()
    except Exception:  # pragma: no cover - network best-effort
        return "unknown"


def _download(url: str, dest: Path) -> None:
    _log(f"  downloading {url}")
    with urllib.request.urlopen(url, timeout=120) as resp, dest.open("wb") as out:
        shutil.copyfileobj(resp, out)


def _fetch_transcripts(base_url: str, workdir: Path) -> list[str]:
    """Download human.N.rna.fna.gz until one 404s; return the source URLs used."""
    urls: list[str] = []
    index = 1
    while True:
        url = f"{base_url}/human.{index}.rna.fna.gz"
        dest = workdir / f"human.{index}.rna.fna.gz"
        try:
            _download(url, dest)
        except urllib.error.HTTPError as exc:
            if exc.code == 404 and index > 1:
                break  # ran off the end of the numbered set
            raise
        urls.append(url)
        index += 1
    if not urls:
        raise ReferenceBuildError(f"No transcript files found under {base_url}")
    return urls


def _concatenate(workdir: Path, fasta_out: Path) -> None:
    """Gunzip every downloaded shard into one plain FASTA."""
    _log(f"  concatenating into {fasta_out}")
    shards = sorted(workdir.glob("human.*.rna.fna.gz"))
    with fasta_out.open("wb") as out:
        for shard in shards:
            with gzip.open(shard, "rb") as gz:
                shutil.copyfileobj(gz, out)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _bowtie_build(fasta: Path, index_prefix: Path) -> None:
    _log(f"  bowtie-build -> {index_prefix}")
    subprocess.run(["bowtie-build", str(fasta), str(index_prefix)], check=True)


def _bowtie_version() -> str:
    try:
        out = subprocess.run(
            ["bowtie-build", "--version"], capture_output=True, text=True
        )
        return out.stdout.splitlines()[0].strip() if out.stdout else "unknown"
    except Exception:  # pragma: no cover
        return "unknown"


def build_index(
    outdir: Path | None = None,
    *,
    base_url: str = DEFAULT_BASE_URL,
    index_name: str = INDEX_NAME,
    keep_fasta: bool = False,
) -> Path:
    """Fetch the transcriptome and build the bowtie index. Returns its prefix."""
    if shutil.which("bowtie-build") is None:
        raise ReferenceBuildError(
            "bowtie-build not found on PATH. The HUSH container ships bowtie1, so "
            "running this inside it needs no further setup; on a host install, "
            "install bowtie1 from your package manager or bioconda first."
        )

    outdir = Path(outdir) if outdir is not None else reference_dir()
    outdir.mkdir(parents=True, exist_ok=True)
    workdir = outdir / "_download"
    workdir.mkdir(exist_ok=True)
    fasta = outdir / f"{index_name}.fasta"
    index_prefix = outdir / index_name

    release = _read_release_number()
    _log(f"RefSeq release: {release}")
    urls = _fetch_transcripts(base_url, workdir)
    _concatenate(workdir, fasta)
    _bowtie_build(fasta, index_prefix)

    manifest = {
        "reference": "RefSeq H_sapiens mRNA_Prot (RNA)",
        "refseq_release": release,
        "base_url": base_url,
        "source_urls": urls,
        "index_prefix": str(index_prefix),
        "fasta_sha256": _sha256(fasta),
        "bowtie_build_version": _bowtie_version(),
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "built_date": date.today().isoformat(),
    }
    path = manifest_path(index_prefix)
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    _log(f"Wrote manifest {path}")

    shutil.rmtree(workdir, ignore_errors=True)
    if not keep_fasta:
        fasta.unlink(missing_ok=True)
    return index_prefix


def describe_index(prefix: str | Path) -> dict | None:
    """The manifest recorded for this index, or None if there is no manifest."""
    path = manifest_path(prefix)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="ascii"))
    except (OSError, ValueError):
        return None
