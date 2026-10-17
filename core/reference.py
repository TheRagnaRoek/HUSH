"""Build and locate the bowtie index of the human transcriptome.

build_index() downloads the RefSeq H_sapiens mRNA_Prot RNA FASTA files, runs
bowtie-build and writes a manifest of what went in. Standard library only.
See docs/design-rules.md (Transcriptome off-target) and docs/cli.md.
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

# The live RefSeq release; use an archived release directory to pin one.
DEFAULT_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot"
RELEASE_NUMBER_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"

REFERENCE_DIR_ENV = "HUSH_REFERENCE_DIR"
DEFAULT_REFERENCE_DIR = Path("data/reference")
INDEX_NAME = "refseq_human_mrna"


def reference_dir() -> Path:
    """Index directory: HUSH_REFERENCE_DIR if set, else data/reference."""
    configured = os.environ.get(REFERENCE_DIR_ENV)
    return Path(configured) if configured else DEFAULT_REFERENCE_DIR


def default_index_prefix() -> Path:
    """Index prefix used when none is given."""
    return reference_dir() / INDEX_NAME


def index_exists(prefix: str | Path) -> bool:
    """Whether a bowtie index exists at the prefix."""
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
    """Download human.N.rna.fna.gz until one 404s; return the URLs used."""
    urls: list[str] = []
    index = 1
    while True:
        url = f"{base_url}/human.{index}.rna.fna.gz"
        dest = workdir / f"human.{index}.rna.fna.gz"
        try:
            _download(url, dest)
        except urllib.error.HTTPError as exc:
            if exc.code == 404 and index > 1:
                break  # end of the numbered set
            raise
        urls.append(url)
        index += 1
    if not urls:
        raise ReferenceBuildError(f"No transcript files found under {base_url}")
    return urls


def _concatenate(workdir: Path, fasta_out: Path) -> None:
    """Gunzip all downloaded files into one FASTA."""
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
    """Download the transcriptome and build the bowtie index; return its prefix."""
    if shutil.which("bowtie-build") is None:
        raise ReferenceBuildError(
            "bowtie-build not found on PATH. The HUSH container ships bowtie1, so "
            "running this inside it needs no further setup; on a host install, "
            "install bowtie1 from your package manager or bioconda first."
        )

    outdir = Path(outdir) if outdir is not None else reference_dir()
    outdir.mkdir(parents=True, exist_ok=True)
    if not os.access(outdir, os.W_OK):
        # A volume from an older image keeps that image's uid. Reading still
        # works; only a rebuild fails, with a bare EACCES partway through.
        uid = getattr(os, "geteuid", lambda: None)()
        who = f" (uid {uid})" if uid is not None else ""
        raise ReferenceBuildError(
            f"{outdir} is not writable by this user{who}. A reference volume "
            "created by an older HUSH image is owned by a different uid. Take "
            "ownership of it once:\n"
            "  docker run --rm --user root -v <volume>:/data/reference <image> "
            "chown -R 10001:10001 /data/reference\n"
            "or point --outdir at a fresh volume."
        )
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
    """The index manifest, or None when it is missing or unreadable."""
    path = manifest_path(prefix)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="ascii"))
    except (OSError, ValueError):
        return None
