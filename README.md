<p align="center">
  <img src="web/public/hush_logo.svg" alt="HUSH" width="420">
</p>

# HUSH

**H**uman genome **U**niversal **S**ilencing **H**airpin-designer - an
allele-aware siRNA/shRNA design pipeline for polymorphic gene families.

Give HUSH a target transcript and, optionally, a panel of related sequences. It
enumerates siRNA candidates, screens them against the classical efficacy rules
and off-target filters, reasons about conserved versus polymorphic regions for
broad or allele-specific knockdown, and converts what survives into cloneable
stem-loop shRNA constructs.

HUSH was inspired by the **HLA** system, among the most polymorphic loci in the
human genome; it is that allelic diversity which motivates the allele-aware
stage. HLA is not baked in - it is the proof-of-concept panel, and any panel of
related sequences works. HUSH is developed as part of **Project HIDE** (HLA
Interference through Disguise and Erasure).

## Quickstart

From a clone of this repository, with nothing else installed but Docker:

```bash
docker compose run --rm reference-build   # one-time, see below
docker compose up hush
```

Then open <http://localhost:8000> for the web GUI. Compose builds the image
from this repository's own Dockerfile, so these commands depend on nothing but
the source you already have. The container carries everything else, including
the `bowtie1` binary the off-target check needs.

A prebuilt image is published for convenience. It is pinned to an exact
version, never a moving tag, so a given revision of this repository always
names the image built from it:

```bash
docker run --rm -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:1.0.0 hush reference build
docker run -p 8000:8000 -v hush_ref:/data/reference ghcr.io/theragnaroek/hush:1.0.0
```

For the CLI rather than the GUI:

```bash
docker run --rm -v "$PWD:/work" -w /work -v hush_ref:/data/reference \
  ghcr.io/theragnaroek/hush:1.0.0 hush design mytarget.fasta -o selected.tsv
```

See [`docker/README.md`](docker/README.md) for Compose, the launcher wrappers and
the configuration variables.

## The reference check

The transcriptome off-target check is HUSH's safety floor: it is what tells you
a guide does not also silence something else. It is **on by default** and needs
no flag once the reference exists.

The reference is not shipped inside the image. It is several gigabytes and
tracks a live RefSeq release, so you fetch it once with `hush reference build`
and every later run finds it automatically. `hush reference status` reports what
you have - the RefSeq release, the bowtie version and the build date recorded in
the index manifest.

HUSH still runs without a reference, because you may want to supply a different
one, but it never does so quietly: the CLI prints an `UNREFERENCED RUN` warning
and the design TSV carries an `offtarget_checked` column, so an empty off-target
hit count can never be mistaken for "no off-targets found".

## Install from source

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# POSIX:    source .venv/bin/activate
pip install -e ".[dev,api]"
```

To install the exact versions HUSH was built and tested with, use the lockfile
instead - this is what the container does, and what to use when reproducing
results:

```bash
pip install -r requirements.txt -r requirements-dev.txt
pip install --no-deps -e .
```

Pure pip, with one external dependency: **bowtie1**, used only for the
transcriptome off-target search. Install it from your package manager or
bioconda, or use the container, which already has it. Everything else - the
design engine, the MSA stage, the conversion stage, the tests - runs with no
binary at all.

## Usage

### Interactive shell

`hush shell` keeps one instance alive: load the target and panel once, tune the
criteria, then re-run `design`/`convert` without re-parsing. The panel alignment
is computed once and reused.

```
hush> load data/examples/target_example.fasta
hush> panel data/examples/hla_panel_example.fasta
hush> set --no-seed --allele-specific "HLA-DEMO*01:01"
hush> design
hush> convert --loop jensen
hush> save run1/
```

`help` inside the shell lists every command. The one-shot commands below do the
same work per invocation and are the way to script HUSH.

### Inspect an input

```bash
hush inspect data/examples/target_example.fasta \
  --panel data/examples/hla_panel_example.fasta
```

A **protein** target (for example parsed from a PDB structure) is flagged:
design needs the coding transcript, so a protein must first be resolved to its
mRNA/CDS.

### Design siRNA candidates

```bash
hush design data/examples/demo_lowgc_cds.fasta -o selected.tsv
```

HUSH enumerates every candidate site (19-nt core, 2-nt 3' overhangs) and reports
those clearing the enabled gates. It is a **filter, not a score**: every gate is
explicit and each candidate keeps its full per-criterion breakdown, so any
verdict is legible. The run prints the criteria, a funnel of how many candidates
clear each gate, and the selected table.

> **GC-rich families:** the seed filter is the siDirect industry standard but is
> deliberately strict - on a very GC-rich CDS such as HLA it can reject every
> candidate. Relax it there with `--no-seed`.

Every flag, its default and what it does: [`docs/cli.md`](docs/cli.md), or
`hush design --help`.

### Allele-aware design

`--panel` engages the MSA stage: the panel is aligned to the target and each
candidate is additionally gated on whether its window **hits** the alleles you
want silenced and **discriminates** against the ones you want spared.

```bash
# Conserved sites across the whole family (the default intent)
hush design target.fasta --no-seed --panel panel.fasta

# Hit one panel member, discriminate against the rest
hush design target.fasta --no-seed --panel panel.fasta \
  --allele-specific "HLA-DEMO*01:01"

# Silence the loaded target, spare the whole reference panel
hush design target.fasta --no-seed --target-specific --panel panel.fasta --explain
```

A request that spares something is ranked by discrimination by default, and the
table gains `minDiff` and `disc` columns describing the **weakest** spared
allele - the one the window distinguishes least well. `--explain` prints the
per-allele breakdown. See [`docs/design-rules.md`](docs/design-rules.md).

### Convert to shRNA

```bash
hush convert selected.tsv --loop jensen -o shrna.tsv
```

Each selected candidate's **full strands** - paired core plus 3' overhang - are
grafted onto a stem-loop scaffold in the sense-loop-antisense orientation, so
the hairpin expresses exactly the guide the design stage validated. Output is
the bare insert as top/bottom DNA oligos plus a Pol III poly-T terminator,
vector-agnostic: add your vector's sticky-end tails before cloning.

Constructs whose body carries a premature Pol III terminator are withheld by
default and counted in the summary; `--keep-premature-terminators` shows them.
Details in [`docs/shrna-conversion.md`](docs/shrna-conversion.md).

## HTTP API and web GUI

The container serves both on port 8000. To run them from a source checkout:

```bash
python -m api          # API on :8000, or: uvicorn api.app:app --reload
cd web && npm run dev  # GUI on :5173, proxying the API
```

Interactive API docs are at `/docs`, the OpenAPI schema at `/openapi.json`. See
[`docs/api.md`](docs/api.md) and [`web/README.md`](web/README.md).

## Reproducing the results

Three claims, each reproducible from a clean checkout. Install from the
lockfile first so the versions match.

**1. The rule pre-screen reproduces siDirect.** siDirect is HUSH's
gold-standard reference; on matched settings the pre-screen reproduces its
picks candidate-for-candidate.

```bash
python tools/compare_sidirect.py \
  --sidirect data/testing/siDirect_PDCD1-CDS.ods \
  --hush data/testing/pdcd1_hush_relaxed.tsv
```

**2. Allele-aware design on a real HLA target.** A GC-rich class I CDS, which
is why the seed filter is relaxed here:

```bash
hush design data/hla_a01.fasta --no-seed -o selected.tsv
hush convert selected.tsv -o shrna.tsv
```

Both runs print a funnel showing how many candidates clear each gate, so the
numbers are auditable rather than asserted. Designs run without a reference
index are labelled `UNREFERENCED RUN`; build the reference first for the
off-target gate to participate.

**3. The test suite**, including the 46-siRNA Amarzguioui Table 3
exact-reproduction fixture:

```bash
pytest -q
```

The bowtie parity test skips when `bowtie` is absent, so run it in the
container to exercise that path - see [`docs/development.md`](docs/development.md).

## Pipeline

| Stage | Module | Status |
|---|---|---|
| 1. Input | `core.inputs`, `core.sequence` | implemented |
| 2. siRNA design | `core.sirna` / `rules` / `thermo` / `seed` / `offtarget` / `asymmetry` / `hairpin` / `selection` | implemented, validated against siDirect |
| 3. MSA / allele | `msa` | implemented, via `hush design --panel` |
| 4. shRNA conversion | `shrna` | implemented, via `hush convert` |
| 5. Output | `app.cli` (tables, TSV), `api.serialize` (JSON), `web` (GUI) | implemented |

## Documentation

| Document | Covers |
|---|---|
| [`docs/architecture.md`](docs/architecture.md) | How the stages and layers fit together, and where a change belongs |
| [`docs/design-rules.md`](docs/design-rules.md) | The science: every gate, why it exists, and its sources |
| [`docs/shrna-conversion.md`](docs/shrna-conversion.md) | Stage 4: hairpin geometry, loops, QC read-outs |
| [`docs/cli.md`](docs/cli.md) | Every command and flag |
| [`docs/api.md`](docs/api.md) | HTTP endpoints and request/response shapes |
| [`docs/development.md`](docs/development.md) | Tests, layout, conventions |
| [`docker/README.md`](docker/README.md) | Running and building the container |
| [`docs/igem-mirror.md`](docs/igem-mirror.md) | The frozen iGEM GitLab copy |

## Project layout

```
core/     Domain model and the siRNA design engine
msa/      MSA / allele stage (alignment, variation, windows, targeting)
shrna/    shRNA conversion (loops, convert)
service/  HushSession: the orchestration seam the CLI and API both drive
app/      CLI (cli.py) and interactive shell (shell.py)
api/      FastAPI backend over the same session layer
web/      React + Vite browser client
data/     Small example and panel FASTAs
tests/    pytest suite
docker/   Dockerfile and launcher wrappers
```

## Development

```bash
pytest
cd web && npm run typecheck && npm run lint
```

The bowtie parity test skips when `bowtie` is not on PATH. See
[`docs/development.md`](docs/development.md).

## Data, models and disclosure

**No AI or machine-learning model is used.** HUSH is a deterministic,
rule-based pipeline: every verdict comes from a published criterion or a
thermodynamic calculation, and the same input always gives the same output. No
model is trained, fine-tuned or called at runtime, and no commercial model API
is contacted. An ML efficacy layer was considered and deliberately left out of
scope; the reasoning is in [`docs/design-rules.md`](docs/design-rules.md).

**No large files live in this repository.** The only bulky artefact is the
transcriptome off-target index, which is built on the user's machine by
`hush reference build` directly from NCBI's RefSeq FTP. It is not redistributed
here, so the reference always matches a RefSeq release the user can name, and
the index manifest records exactly which one.

**Dependencies are pinned.** `requirements.txt`, `requirements-dev.txt` and
`web/package-lock.json` record the exact versions used; `pyproject.toml` holds
the supported ranges.

## License

MIT. See [`LICENSE`](LICENSE).
