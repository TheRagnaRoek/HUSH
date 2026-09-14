#!/usr/bin/env sh
# HUSH launcher (Docker). Runs the containerized HUSH CLI - Python, bowtie1 and
# the RefSeq off-target index - so nothing needs installing on the host. The
# current directory is mounted at /work (HUSH reads and writes there); the index
# comes from the reference volume at /data/reference. Arguments pass through to
# hush, e.g.:
#   ./docker/hush.sh design mytarget.fasta --offtarget-db /data/reference/refseq_human_mrna
# Override the image or volume name with the HUSH_IMAGE / HUSH_REF_VOLUME env vars.
image="${HUSH_IMAGE:-hush:latest}"
ref="${HUSH_REF_VOLUME:-hush_hush_ref}"
exec docker run --rm -v "$ref:/data/reference" -v "$PWD:/work" -w /work "$image" hush "$@"
