#!/usr/bin/env python3
"""Fetch a RefSeq human transcriptome and build a bowtie index for off-target search.

Downloads the RefSeq H_sapiens mRNA_Prot RNA FASTA set, concatenates it into one
file, runs `bowtie-build`, and writes a manifest recording exactly what went in
(RefSeq release, source URLs, sizes, sha256, bowtie version, build date) so a
build is reproducible. The index this produces is what `hush design --offtarget
--offtarget-db <prefix>` consumes; core.offtarget.BowtieSearcher reads it.

This is a deliberate, one-off, network-heavy build step - not part of the design
pipeline or the test suite. It is an ops-side job run once by us, not by the end
user: HUSH ships browser-based with the index prebuilt and hosted server-side.
The reference is large; nothing here is committed (data/reference is gitignored).

Preferred: run it in the HUSH container, which already carries bowtie1, so
nothing needs installing on the host (see docker/README.md):
    docker compose run --rm build-index

Direct usage (needs bowtie1's `bowtie-build` on PATH, e.g. inside that image):
    python tools/build_offtarget_db.py            # into data/reference/
    python tools/build_offtarget_db.py --outdir data/reference --keep-fasta

Only the Python standard library is used otherwise.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

# RefSeq's current human mRNA/ncRNA transcripts. These files track the live
# RefSeq release; the actual release number is read from RELEASE_NUMBER and
# recorded in the manifest so the build is auditable. Point --base-url at an
# archived release directory for an exact, permanent pin.
DEFAULT_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot"
RELEASE_NUMBER_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"

DEFAULT_OUTDIR = Path("data/reference")
INDEX_NAME = "refseq_human_mrna"


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
        raise SystemExit(f"No transcript files found under {base_url}")
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
    subprocess.run(
        ["bowtie-build", str(fasta), str(index_prefix)],
        check=True,
    )


def _bowtie_version() -> str:
    try:
        out = subprocess.run(
            ["bowtie-build", "--version"], capture_output=True, text=True
        )
        return out.stdout.splitlines()[0].strip() if out.stdout else "unknown"
    except Exception:  # pragma: no cover
        return "unknown"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                    help="Directory for the index and manifest (default data/reference).")
    ap.add_argument("--base-url", default=DEFAULT_BASE_URL,
                    help="RefSeq mRNA_Prot directory URL (point at an archived "
                    "release to pin exactly).")
    ap.add_argument("--index-name", default=INDEX_NAME,
                    help="Basename of the bowtie index prefix.")
    ap.add_argument("--keep-fasta", action="store_true",
                    help="Keep the concatenated FASTA after building the index.")
    args = ap.parse_args(argv)

    if shutil.which("bowtie-build") is None:
        raise SystemExit(
            "bowtie-build not found on PATH. Run this in the HUSH container, which "
            "ships bowtie1: 'docker compose run --rm build-index' (see "
            "docker/README.md). Installing bowtie1 on the host is not required."
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    workdir = args.outdir / "_download"
    workdir.mkdir(exist_ok=True)
    fasta = args.outdir / f"{args.index_name}.fasta"
    index_prefix = args.outdir / args.index_name

    release = _read_release_number()
    _log(f"RefSeq release: {release}")
    urls = _fetch_transcripts(args.base_url, workdir)
    _concatenate(workdir, fasta)
    _bowtie_build(fasta, index_prefix)

    manifest = {
        "reference": "RefSeq H_sapiens mRNA_Prot (RNA)",
        "refseq_release": release,
        "base_url": args.base_url,
        "source_urls": urls,
        "index_prefix": str(index_prefix),
        "fasta_sha256": _sha256(fasta),
        "bowtie_build_version": _bowtie_version(),
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "built_date": date.today().isoformat(),
    }
    manifest_path = args.outdir / f"{args.index_name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="ascii")
    _log(f"Wrote manifest {manifest_path}")

    shutil.rmtree(workdir, ignore_errors=True)
    if not args.keep_fasta:
        fasta.unlink(missing_ok=True)
    _log(f"Done. Use: hush design ... --offtarget --offtarget-db {index_prefix}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
