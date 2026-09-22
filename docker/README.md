# HUSH container image (bowtie1 + off-target index)

The transcriptome off-target stage needs the `bowtie1` binary and a prebuilt
RefSeq index. Neither is an end-user dependency: HUSH ships as a hosted,
browser-based tool, so bowtie and the index live on **our** side. This image is
where they live.

One image, two jobs:

1. **Build the index once** with `tools/build_offtarget_db.py` (downloads RefSeq,
   runs `bowtie-build`, writes a reproducibility manifest).
2. **Serve designs** - the future web backend runs in this same image and shells
   out to the bowtie already installed here.

Because it is a Linux container, it also removes the pain of installing bowtie1
via bioconda on native Windows: everything below runs the same on any host with
Docker.

## Build the index (one time)

From the repository root:

```bash
docker compose run --rm build-index
```

That builds the image if needed, downloads the RefSeq human transcriptome, runs
`bowtie-build`, and writes the index plus `refseq_human_mrna.manifest.json` into
the `hush_ref` Docker volume. It is network-heavy and runs once; the volume is
reused afterwards. Point the builder at an archived RefSeq release for a
permanent pin by passing `--base-url` (see the script's `--help`).

## Run the HUSH CLI without installing anything

The repo's root launchers (`./hush`, `hush.bat`) run the **pip-only** CLI from a
local venv - no bowtie, so transcriptome off-target is skipped. To get the
**full** pipeline (bowtie + the RefSeq index) without installing anything on the
host, run the CLI inside the container. Two ready-made wrappers do this; both
mount your current directory at `/work` (so HUSH reads and writes your real
files) and mount the reference volume at `/data/reference`.

Prerequisites: the image is built (`docker build -f docker/Dockerfile -t
hush:latest .`) and the index exists (`docker compose run --rm build-index`).

**Windows (PowerShell)** - [`docker/hush.ps1`](hush.ps1):

```powershell
.\docker\hush.ps1 design mytarget.fasta --offtarget-db /data/reference/refseq_human_mrna
```

**macOS / Linux** - [`docker/hush.sh`](hush.sh) (once, `chmod +x docker/hush.sh`):

```bash
./docker/hush.sh design mytarget.fasta --offtarget-db /data/reference/refseq_human_mrna
```

Everything after the script name passes straight through to `hush`, so `--help`,
`-o out.tsv`, `inspect`, `--panel`, etc. all work. Outputs land in the folder you
ran from. The image name and reference-volume name are overridable with the
`HUSH_IMAGE` and `HUSH_REF_VOLUME` environment variables (the volume defaults to
`hush_hush_ref`, the name Docker Compose gives it - see the "Volume ... Created"
line from the build).

### Make it feel native

There is no `activate` step like a venv, but a one-line shell function makes the
wrapper look like a plain `hush` command from anywhere. Point it at your checkout:

```powershell
# PowerShell: add to your profile (notepad $PROFILE), using your repo path
function hush { & "D:\Dateien\Dokumente\iGEM2026\HUSH\docker\hush.ps1" @args }
```

```bash
# bash/zsh: add to ~/.bashrc or ~/.zshrc, using your repo path
hush() { /path/to/HUSH/docker/hush.sh "$@"; }
```

### Or enter the environment interactively

For a session where you run several commands, open a shell inside the container -
the closest thing to "activating" it. `hush` and `bowtie` are on `PATH`; `exit`
leaves:

```bash
docker run --rm -it -v hush_hush_ref:/data/reference -v "${PWD}:/work" -w /work hush:latest bash
```

For the HUSH REPL specifically (`hush shell` - load a target once, then tune and
re-run `design`/`convert`), the wrappers work directly:

```powershell
.\docker\hush.ps1 shell
```

They pass an interactive terminal through when run from a console, so the REPL
stays open. When their input or output is redirected (a pipe, or CI) they drop
the interactive flags automatically, so one-shot `design`/`convert` runs are
unaffected.

## Hosting the index

Build the index once, then publish the volume's contents (the `.ebwt` shards and
the manifest) to wherever the backend reads them - a mounted disk, an object
store, or a release asset the runtime downloads on startup. End users only ever
hit the web UI; they never install bowtie and never build or download the index
themselves.
