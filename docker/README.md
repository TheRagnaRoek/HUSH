# Running HUSH in Docker

HUSH is distributed to users and not hosted, so the container is the unit that is
run. One image contains the web GUI, the HTTP API, the CLI and the `bowtie1`
binary that the off-target check needs. Nothing else has to be installed.

The image does **not** contain the transcriptome reference. It is several
gigabytes and tracks a live RefSeq release, so it is fetched once into a Docker
volume, and every later run finds it automatically.

## Quickstart

From the repository root, building from source:

```bash
docker compose run --rm reference-build
docker compose up --build hush
```

Then open <http://localhost:8000>. The first command is the one-time reference
download; the second serves the GUI and the API together on port 8000. The image
is built from `docker/Dockerfile` and tagged `hush:local`, so this needs no
registry access. Compose reuses an existing `hush:local` without rebuilding, so
keep `--build` after pulling new code, or the container will run the old source.

The same two steps with the published image, pinned to an exact version:

```bash
docker run --rm -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:1.0.0 hush reference build
docker run -p 8000:8000 -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:1.0.0
```

## The reference

The transcriptome off-target check tests whether a guide also silences other
transcripts. It is on by default and needs no flag once the reference exists.

HUSH also runs without the built reference, for example to use a different
reference genome, but it then labels the run: the CLI prints an
`UNREFERENCED RUN` warning and the design TSV has an `offtarget_checked` column,
so an empty hit count cannot be mistaken for "no off-targets found".

Check the installed reference:

```bash
docker run --rm -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:1.0.0 hush reference status
```

This reports the RefSeq release, the bowtie version and the build date recorded in
the manifest of the index. To fix the reference to an exact release instead of the
live one, pass an archived RefSeq release directory as `--base-url`. To use a
different reference entirely, build it separately and set `--offtarget-db` to its
prefix.

## Upgrading from an older image

A reference volume built by an older HUSH image is owned by the user of that image,
not by the current one. Designs are not affected, because reading the index does
not need write access, but **rebuilding the reference into that volume fails**.
HUSH reports this before the download starts and not part-way through it.

Take ownership of the existing volume once:

```bash
docker run --rm --user root -v hush_hush_ref:/data/reference   ghcr.io/theragnaroek/hush:1.0.0 chown -R 10001:10001 /data/reference
```

Alternatively, build into a fresh volume and delete the old one afterwards.

## The CLI

The default command of the image serves the GUI. Pass a command to use the CLI
instead:

```bash
docker run --rm -v "$PWD:/work" -w /work -v hush_ref:/data/reference \
  ghcr.io/theragnaroek/hush:1.0.0 hush design mytarget.fasta -o selected.tsv
```

Mounting the current directory at `/work` lets HUSH read input files and write
results.

Two wrapper scripts do this and pass all arguments to `hush`:

```powershell
.\docker\hush.ps1 design mytarget.fasta -o selected.tsv
```

```bash
./docker/hush.sh design mytarget.fasta -o selected.tsv   # chmod +x once
```

The wrappers attach an interactive terminal when run from a console, so
`hush shell` works, and omit it when input or output is redirected, so pipes and
CI work. `HUSH_IMAGE` and `HUSH_REF_VOLUME` override the image and volume names.

### Calling the wrapper as `hush`

A one-line shell function lets the wrapper be called as `hush`:

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
  ghcr.io/theragnaroek/hush:1.0.0 bash
```

`hush`, `bowtie` and `bowtie-build` are all on `PATH`.

## Building the image

The Dockerfile is part of the source, so building does not depend on a registry:

```bash
docker build -f docker/Dockerfile -t hush:local .
HUSH_IMAGE=hush:local docker compose up hush
```

The build has two fixed choices:

- The web client is built in a `node` stage pinned to `BUILDPLATFORM`. The bundle
  is identical on every architecture, so this stage runs natively and not under
  QEMU emulation in the arm64 build.
- The runtime base is pinned to `bookworm`. This is the Debian suite in which
  bowtie1 is packaged for both amd64 and arm64, so a multi-arch image can be built
  without compiling the aligner.

## Configuration

| Variable | Default | Description |
|---|---|---|
| `HUSH_API_HOST` | `0.0.0.0` in the image, `127.0.0.1` otherwise | Interface the API binds |
| `HUSH_API_PORT` | `8000` | Port the API binds |
| `HUSH_REFERENCE_DIR` | `/data/reference` in the image | Location of the off-target index |
| `HUSH_WEB_DIR` | `/app/web` in the image | Built web client to serve; if unset, only the API is served |

Outside the container, the API binds to the loopback interface by default, which
suits a local tool. The image overrides this because the loopback interface of a
container is not reachable through a published port.

## Security

The API has no authentication, and a session accepts a file system path for an
alternative reference index. It is intended to run on the local machine. Do not
publish the port beyond localhost without adding authentication and a reverse
proxy.
