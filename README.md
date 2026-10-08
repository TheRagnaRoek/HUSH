# HUSH

**H**uman genome **U**niversal **S**ilencing **H**airpin-designer - a modern,
allele-aware siRNA/shRNA design pipeline for polymorphic gene families.

HUSH is an open-source, from-sequence RNAi design tool: give it a target gene
sequence (and, optionally, a panel of related sequences) and it enumerates siRNA
candidates, screens them against the classical efficacy rules and off-target
filters, and reasons about conserved vs. polymorphic regions for broad or
allele-specific knockdown. The pipeline is sequence-agnostic - any panel of
related sequences works.

HUSH was inspired by the **HLA** system, among the most polymorphic loci in the
human genome: it is that allelic diversity which motivates the allele-aware
design stage. HLA is **not** baked in - it is the proof-of-concept panel. HUSH is
developed as part of **Project HIDE** (HLA Interference through Disguise and
Erasure). See [`docs/HUSH_concept_dossier.md`](docs/HUSH_concept_dossier.md) for
the full rationale and prior-art landscape.

## Pipeline

| Stage | Module | Status |
|---|---|---|
| 1. Input Layer | `core.inputs` + `core.sequence` | implemented |
| 2. siRNA Design | `core.sirna` / `rules` / `thermo` / `seed` / `offtarget` / `asymmetry` / `hairpin` / `selection` | implemented, validated against siDirect |
| 3. MSA / Allele | `msa` | implemented, wired into `hush design --panel` |
| 4. shRNA Conversion | `shrna` | implemented (stem-loop scaffold), wired into `hush convert` |
| 5. Output Layer | `core` / `app` | planned |

## Install

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# POSIX:    source .venv/bin/activate
pip install -e ".[dev]"
```

Pure-pip; the core needs no external binaries. The one non-pip piece is
**bowtie1**, used only for transcriptome-scale off-target search
(`hush design --offtarget-db`). It is **not** an end-user dependency: HUSH ships
as a hosted, browser-based tool, so bowtie and the RefSeq index live server-side
in the container image under [`docker/`](docker/README.md) - built and hosted
once by us, never installed or built by the user. For local development against
the transcriptome, that same image is the supported way to get bowtie and the
index (see [`docker/README.md`](docker/README.md)); the pure-Python core, tests,
and small-panel path all run with no binary at all.

Once installed, the thin launchers call the venv's interpreter directly, so
`hush` works from the project folder without activating the venv:

```bash
# Windows:  .\hush design data/examples/hla_a01.fasta
# POSIX:    ./hush  design data/examples/hla_a01.fasta
```

## Usage

### Interactive shell

For iterative work, `hush shell` keeps one instance alive: load the target (and
panel) once, tune the criteria, then re-run `design`/`convert` without
re-loading or re-parsing. The panel alignment is computed once and reused across
runs.

```bash
hush shell
```

Inside the shell:

| Command | Does |
|---|---|
| `load <path>` | load the target mRNA/CDS |
| `panel <path>` | load a related-sequence panel (turns on allele targeting) |
| `set [flags]` | set persistent design criteria (same flags as `hush design`); no flags prints them |
| `design [flags]` | design with the current criteria (inline flags apply to that run only) |
| `convert [flags]` | convert the last design's selected candidates to shRNA |
| `save <dir>` | write `design.tsv` (and `shrna.tsv`) into `<dir>` |
| `show` / `inspect` | review the criteria / preview the design input |
| `limit <n>`, `reset`, `help`, `quit` | display cap, restore defaults, help, exit |

```
hush> load data/examples/target_example.fasta
hush> panel data/examples/hla_panel_example.fasta
hush> set --no-seed --allele-specific "HLA-DEMO*01:01"
hush> design
hush> convert --loop jensen
hush> save run1/
```

The one-shot commands below do the same work per invocation and remain the way
to script HUSH.

### Inspect an input

List every sequence in a FASTA/PDB/mmCIF file, or preview the design input with a
panel (target/panel roles, molecule types, readiness):

```bash
hush inspect data/examples/hla_panel_example.fasta
hush inspect data/examples/target_example.fasta \
  --panel data/examples/hla_panel_example.fasta
```

A **protein** target (e.g. parsed from a PDB structure) is flagged: siRNA/shRNA
design needs the coding transcript, so a protein must first be resolved to its
mRNA/CDS.

### Design siRNA candidates

```bash
hush design data/examples/demo_lowgc_cds.fasta
```

`hush design` enumerates every candidate site (19-nt core by default, with 2-nt
3' overhangs) and reports the ones that clear the enabled gates, in
target-position order. It is a **filter, not a ranking**: every gate is explicit
and each candidate keeps its full per-criterion breakdown, so any verdict is
legible. The run prints the criteria, a funnel of how many candidates clear each
gate, and the selected table; `-o file.tsv` writes them all.

The **default paradigm** (best-practice, not the strictest or the loosest):

| Gate | Default | Flag |
|---|---|---|
| Classical rules (Ui-Tei / Reynolds / Amarzguioui), combined by **union** | **on** | `--rules`, `--mode` |
| Seed-duplex off-target Tm (siDirect) | **on** | `--seed` / `--no-seed` |
| Homopolymer-run gate (rejects e.g. `GGGG`) | **on** | `--forbid-homopolymers` / `--allow-homopolymers` |
| Transcriptome off-target (needs an index) | **on**\* | `--offtarget-db`, `--no-offtarget` |
| GC-content window | off | `--gc` |
| Internal-hairpin gate | off | `--hairpin` |
| Duplex-end loading asymmetry | off | `--asymmetry` |

\* On by default, but it needs a reference index; with none configured it is
skipped with a warning rather than erroring. In the shipped tool the index is
prebuilt and served server-side. To build one for local development, use the
container (`docker compose run --rm build-index`, see
[`docker/README.md`](docker/README.md)), then pass `--offtarget-db`.

> **GC-rich families:** the seed filter is the siDirect industry standard but is
> deliberately strict - on a very GC-rich CDS (e.g. HLA) it can reject every
> candidate. Relax it there with `--no-seed`.

```bash
# GC-rich target: drop the seed filter, add the GC and hairpin gates
hush design data/examples/hla_a01.fasta --no-seed --gc --hairpin -o selected.tsv
```

See [`docs/design_options.md`](docs/design_options.md) for a one-page reference to
every option, its default, and what it does (or `hush design --help`).

### Allele-aware design

Add `--panel` to engage the MSA/allele stage: the panel is aligned to the target
and each candidate is additionally gated on whether its window **hits** the
alleles you want silenced and **discriminates** against the ones you want spared.
Four targeting intents:

- `--pan-allele` - conserved sites across the whole family (the default);
- `--allele-specific ID` - single one panel member out, spare the rest (use when
  the allele to silence is itself in the panel);
- `--target-specific` - spare the **whole** panel and hit none of it (use when
  the target is loaded apart from its reference panel: it is silenced by
  construction, so specificity means discriminating against every panel member);
- `--hit` (+ `--spare`) - explicit subsets.

```bash
# Pan-allele (conserved-region) knockdown across the family
hush design data/examples/target_example.fasta --no-seed \
  --panel data/examples/hla_panel_example.fasta

# Allele-specific: hit one panel member, discriminate against the rest
hush design data/examples/target_example.fasta --no-seed \
  --panel data/examples/hla_panel_example.fasta \
  --allele-specific "HLA-DEMO*01:01"

# Target-specific: silence the loaded target, spare the whole reference panel
hush design data/examples/my_target.fasta --no-seed --target-specific \
  --panel data/examples/hla_reference_panel.fasta --explain
```

The criteria block and funnel gain an `allele:` line. A request that spares
something is ranked by discrimination by default (`--sort`), and the table adds
`minDiff` and `disc` columns describing the **weakest** spared allele - the one
the window distinguishes least well. `--explain` (or `assess <region>` in the
shell) prints the full per-allele breakdown for every spared allele; `-o` adds
the allele columns to the TSV.

### Convert siRNA candidates to shRNA

`hush convert` reads a `hush design -o` TSV and grafts each selected candidate's
19-bp core (unchanged, so the validated guide is preserved) onto an established
stem-loop scaffold in the sense-loop-antisense orientation. Output is the bare
hairpin insert as top/bottom DNA oligos plus a Pol III poly-T terminator - it is
vector-agnostic, so add your vector's restriction sticky-end tails before
cloning. Each construct is folded to confirm it forms the intended stem-loop (not
a competing structure), with QC flags for Pol III 5'-A/G initiation, loop/stem
terminator-safety, and homopolymer runs.

```bash
hush design data/examples/hla_a01.fasta --no-seed -o selected.tsv
hush convert selected.tsv --loop jensen -o shrna.tsv
```

The loop is selectable (`--loop jensen|psuper|trc`, default the least
stem-dependent Jensen UGUGCUU); `--no-terminator` and `--terminator-length` tune
the poly-T terminator. Scaffold rationale and references are in
[`docs/shrna_design_note.md`](docs/shrna_design_note.md).

## HTTP API

An optional FastAPI backend exposes the same pipeline as JSON, for the browser
client. Install the extra and run it:

```bash
pip install -e ".[api]"
python -m api                       # or: uvicorn api.app:app --reload
```

Interactive docs are at `/docs` and the OpenAPI schema at `/openapi.json` (from
which a TypeScript client can be generated). The target and panel are sent as
FASTA text in the JSON body; GC bounds and results mirror the CLI.

Two ways to call it. Stateless, one request carries everything:

```bash
curl -s localhost:8000/design -H 'content-type: application/json' -d '{
  "target_fasta": ">t\nACGT...",
  "params": {"use_seed_filter": false}
}'
```

Or a session: post the inputs once, then design/convert by id (the server caches
the panel alignment, so re-designs with new parameters do not re-align):

```
POST   /sessions                 {target_fasta, panel_fasta?}  -> {id, ...}
POST   /sessions/{id}/design     {params}                      -> results
POST   /sessions/{id}/convert    {options}                     -> constructs
DELETE /sessions/{id}
```

Sessions are an in-memory warm cache (lost on restart), not a database.

## Web GUI

A browser front-end (React + Vite + TypeScript) renders this API as a guided,
allele-aware design workflow with a construct viewer. It lives in
[`web/`](web/README.md) and runs against the FastAPI backend above (it also opens
on sample data with no backend). See [`web/README.md`](web/README.md) to develop it.

## Benchmarking against siDirect

siDirect is HUSH's gold-standard reference; on matched settings the rule
pre-screen reproduces its picks candidate-for-candidate. `tools/compare_sidirect.py`
diffs a HUSH design TSV against a siDirect export:

```bash
python tools/compare_sidirect.py \
  --sidirect data/siDirect_PDCD1-CDS.ods --hush data/pdcd1_hush_relaxed.tsv
```

## Project layout

```
core/    Domain model + the siRNA design engine (sequence, inputs, sirna, rules,
         thermo, seed, offtarget, asymmetry, hairpin, selection, errors)
msa/     MSA / allele stage (alignment, align_star, variation, windows, targeting)
shrna/   shRNA conversion (stage 4: loops, convert)
service/ Orchestration layer (HushSession) the CLI and API render over
app/     CLI: one-shot commands (app.cli) + interactive shell (app.shell)
api/     HTTP API (FastAPI): schemas, serialize, store, app (optional extra)
web/     Browser GUI (React + Vite + TypeScript) over the API - see web/README.md
data/    Example inputs and benchmarking fixtures
tools/   Developer utilities (build_offtarget_db.py, compare_sidirect.py)
docs/    Concept dossier and design notes
tests/   pytest suite
```

## Documentation

- [`docs/HUSH_concept_dossier.md`](docs/HUSH_concept_dossier.md) - project vision and prior art.
- [`docs/design_notes.md`](docs/design_notes.md) - long-form rationale, formulae and validation records for every stage.
- [`docs/design_options.md`](docs/design_options.md) - one-page `hush design` option reference.
- [`docs/shrna_design_note.md`](docs/shrna_design_note.md) - the stem-loop shRNA conversion: scope, parameters, and references.
- [`docs/offtarget_design_note.md`](docs/offtarget_design_note.md) and [`docs/seed_offtarget_design_note.md`](docs/seed_offtarget_design_note.md) - the two off-target filters in depth.
- Module docstrings carry the per-unit contract and scope.

## Development

```bash
pytest
```
