#!/usr/bin/env sh
# HUSH launcher (Docker). Runs the containerized HUSH CLI - Python and bowtie1
# included - so nothing needs installing on the host. The current directory is
# mounted at /work (HUSH reads and writes there), and the off-target reference
# comes from the volume at /data/reference. Arguments pass through to hush:
#   ./docker/hush.sh reference build          # once, before designing
#   ./docker/hush.sh design mytarget.fasta
# The reference is found automatically, so no --offtarget-db is needed.
# Override the image or volume name with the HUSH_IMAGE / HUSH_REF_VOLUME env vars.
image="${HUSH_IMAGE:-ghcr.io/theragnaroek/hush:1.0.0}"
ref="${HUSH_REF_VOLUME:-hush_hush_ref}"

# Attach an interactive TTY only when stdin and stdout are terminals, so
# 'hush shell' (a REPL that reads stdin) stays open while pipes and CI still
# work. Without it the container's stdin is closed and the shell exits at once.
tty=""
if [ -t 0 ] && [ -t 1 ]; then
  tty="-it"
fi
exec docker run --rm $tty -v "$ref:/data/reference" -v "$PWD:/work" -w /work "$image" hush "$@"
