# Running HUSH in Docker

HUSH is distributed to its users rather than hosted, so the container is the
unit you run. One image carries everything: the web GUI, the HTTP API, the CLI,
and the `bowtie1` binary the off-target check needs. Nothing else goes on your
machine.

The one thing the image does **not** carry is the transcriptome reference. It is
several gigabytes and tracks a live RefSeq release, so you fetch it once into a
Docker volume and every later run finds it automatically.

## Quickstart

```bash
docker run --rm -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:latest hush reference build
docker run -p 8000:8000 -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:latest
```

Then open <http://localhost:8000>. The first command is the one-time reference
fetch; the second serves the GUI and the API together on port 8000.

With Compose, from the repository root:

```bash
docker compose run --rm reference-build
docker compose up hush
```

## The reference, and why it matters

The transcriptome off-target check is HUSH's safety floor: it is what tells you
a guide does not also silence something else. It is on by default and needs no
flag once the reference exists.

HUSH still runs without it, because you may want to swap in a different
reference genome, but every such run is labelled: the CLI prints an
`UNREFERENCED RUN` warning and the design TSV carries an `offtarget_checked`
column, so an empty hit count can never be mistaken for "no off-targets found".

Check what you have:

```bash
docker run --rm -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:latest hush reference status
```

That reports the RefSeq release, the bowtie version and the build date recorded
in the index's manifest. To pin an exact reference rather than the live release,
pass `--base-url` an archived RefSeq release directory. To use a different
reference entirely, build it yourself and point `--offtarget-db` at the prefix.

## Upgrading from an older image

A reference volume built by an older HUSH image is owned by that image's user,
not this one. Designs are unaffected, because reading the index does not need
write access, but **rebuilding the reference into that volume fails**. HUSH
reports this up front rather than part-way through a multi-gigabyte download.

Take ownership of the existing volume once:

```bash
docker run --rm --user root -v hush_hush_ref:/data/reference   ghcr.io/theragnaroek/hush:latest chown -R 10001:10001 /data/reference
```

Or simply build into a fresh volume and drop the old one once you are happy.

## The CLI

The image's default command serves the GUI. Pass a command to use the CLI
instead:

```bash
docker run --rm -v "$PWD:/work" -w /work -v hush_ref:/data/reference \
  ghcr.io/theragnaroek/hush:latest hush design mytarget.fasta -o selected.tsv
```

Mounting the current directory at `/work` is what lets HUSH read your real files
and write results back out.

Two wrappers do that for you, and pass everything through to `hush`:

```powershell
.\docker\hush.ps1 design mytarget.fasta -o selected.tsv
```

```bash
./docker/hush.sh design mytarget.fasta -o selected.tsv   # chmod +x once
```

They attach an interactive terminal when run from a console, so the REPL
(`hush shell`) stays open, and drop it when input or output is redirected, so
pipes and CI still work. `HUSH_IMAGE` and `HUSH_REF_VOLUME` override the image
and volume names.

### Making it feel native

A one-line shell function makes the wrapper look like a plain `hush`:

```powershell
# PowerShell: add to your profile (notepad $PROFILE), using your checkout path
function hush { & "C:\path\to\HUSH\docker\hush.ps1" @args }
```

```bash
# bash/zsh: add to ~/.bashrc or ~/.zshrc, using your checkout path
hush() { /path/to/HUSH/docker/hush.sh "$@"; }
```

### An interactive shell inside the image

```bash
docker run --rm -it -v hush_ref:/data/reference -v "$PWD:/work" -w /work \
  ghcr.io/theragnaroek/hush:latest bash
```

`hush`, `bowtie` and `bowtie-build` are all on `PATH`.

## Building the image yourself

The Dockerfile ships with the source, so you never depend on the registry:

```bash
docker build -f docker/Dockerfile -t hush:local .
HUSH_IMAGE=hush:local docker compose up hush
```

Two things in the build are deliberate:

- The web client is built in a `node` stage pinned to `BUILDPLATFORM`. The
  bundle is identical on every architecture, so building it natively keeps the
  arm64 image off QEMU for that whole stage.
- The runtime base is pinned to `bookworm`. That is the Debian suite where
  bowtie1 is packaged for both amd64 and arm64, which is what makes a multi-arch
  image possible without compiling the aligner.

## Configuration

| Variable | Default | What it does |
|---|---|---|
| `HUSH_API_HOST` | `0.0.0.0` in the image, `127.0.0.1` otherwise | Interface the API binds |
| `HUSH_API_PORT` | `8000` | Port the API binds |
| `HUSH_REFERENCE_DIR` | `/data/reference` in the image | Where the off-target index lives |
| `HUSH_WEB_DIR` | `/app/web` in the image | Built web client to serve; unset serves the API alone |

The API binds loopback by default outside the container, which is the right
default for a local tool; the image overrides it because a container's loopback
is not reachable from a published port.

## A note on exposing it

The API has no authentication, and a session accepts a filesystem path for an
alternative reference index. It is built to run on your own machine. Do not
publish the port beyond localhost without putting your own authentication and a
reverse proxy in front of it.
