#!/usr/bin/env pwsh
# HUSH launcher (Docker). Runs the containerized HUSH CLI - Python, bowtie1 and
# the RefSeq off-target index - so nothing needs installing on the host. Your
# current directory is mounted at /work (HUSH reads and writes there); the index
# comes from the reference volume at /data/reference. Arguments pass through to
# hush, e.g.:
#   .\docker\hush.ps1 design mytarget.fasta --offtarget-db /data/reference/refseq_human_mrna
# Override the image or volume name with the HUSH_IMAGE / HUSH_REF_VOLUME env vars.
$ErrorActionPreference = 'Stop'
$image = if ($env:HUSH_IMAGE) { $env:HUSH_IMAGE } else { 'hush:latest' }
$ref   = if ($env:HUSH_REF_VOLUME) { $env:HUSH_REF_VOLUME } else { 'hush_hush_ref' }
docker run --rm -v "${ref}:/data/reference" -v "$($PWD.Path):/work" -w /work $image hush @args
exit $LASTEXITCODE
