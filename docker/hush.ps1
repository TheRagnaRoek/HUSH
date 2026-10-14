#!/usr/bin/env pwsh
# HUSH launcher (Docker). Runs the containerized HUSH CLI - Python and bowtie1
# included - so nothing needs installing on the host. Your current directory is
# mounted at /work (HUSH reads and writes there), and the off-target reference
# comes from the volume at /data/reference. Arguments pass through to hush:
#   .\docker\hush.ps1 reference build          # once, before designing
#   .\docker\hush.ps1 design mytarget.fasta
# The reference is found automatically, so no --offtarget-db is needed.
# Override the image or volume name with the HUSH_IMAGE / HUSH_REF_VOLUME env vars.
$ErrorActionPreference = 'Stop'
$image = if ($env:HUSH_IMAGE) { $env:HUSH_IMAGE } else { "ghcr.io/theragnaroek/hush:1.0.0" }
$ref   = if ($env:HUSH_REF_VOLUME) { $env:HUSH_REF_VOLUME } else { 'hush_hush_ref' }

# Attach an interactive TTY only at a real console, so 'hush shell' (a REPL that
# reads stdin) stays open while piped/redirected runs still work. Without it the
# container's stdin is closed and the shell hits EOF and exits immediately.
$tty = @()
if (-not [Console]::IsInputRedirected -and -not [Console]::IsOutputRedirected) {
    $tty = @('-it')
}
docker run --rm @tty -v "${ref}:/data/reference" -v "$($PWD.Path):/work" -w /work $image hush @args
exit $LASTEXITCODE
