# HUSH

**H**uman genome **U**niversal **S**ilencing **H**airpin-designer - a modern,
allele-aware siRNA/shRNA design pipeline for polymorphic gene families.

HUSH is an open-source, from-sequence RNAi design tool: it takes a target gene
sequence (and, optionally, a panel of related sequences), designs siRNA
candidates, reasons about conserved vs. polymorphic regions for broad or
allele-specific knockdown, and converts validated candidates into cloneable
shRNA constructs. The pipeline itself is sequence-agnostic - any panel of
related sequences works. See
[`docs/HUSH_concept_dossier.md`](docs/HUSH_concept_dossier.md) for the full
rationale and prior-art landscape.

### Why "allele-aware", and where HLA comes in

HUSH was inspired by the HLA system - among the most polymorphic loci in the
human genome. It is that allelic diversity which motivates the MSA/allele
module: a design tool cannot simply be pointed at a highly polymorphic gene and
be trusted to hit (or spare) the right alleles. HLA is **not** baked into the
tool. It serves as the **proof-of-concept panel**, and HUSH is developed as part
of the wider **Project HIDE** (HLA Interference through Disguise and Erasure).

## Pipeline

| Stage | Module | Status |
|---|---|---|
| 1. Input Layer | `core.inputs` + `core.sequence` | implemented |
| 2. siRNA Design | `core.sirna` + `core.rules` | in progress - candidate enumeration and the classical rule pre-screen are implemented; thermodynamics next |
| 3. MSA / Allele | `msa` | planned |
| 4. shRNA Conversion | `shrna` | planned |
| 5. Output Layer | `core` / `app` | planned |

## Project layout

```
core/          Shared domain model + the siRNA design engine
  sequence.py    BioSequence / SequencePanel / DesignInput (pure stdlib)
  inputs.py      Input Layer: FASTA & PDB parsing -> domain model
  sirna.py       siRNA candidate model + window enumeration
  rules.py       Classical rule pre-screen (Ui-Tei / Reynolds / Amarzguioui)
  errors.py      Shared exception types
msa/           MSA / allele module (stage 3, placeholder)
shrna/         shRNA conversion module (stage 4, placeholder)
app/           CLI (app.cli) - the CLI-first interface
data/examples/ Runnable example inputs (synthetic; not real genes)
tests/         pytest suite
docs/          Concept dossier
```

## Getting started

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# POSIX:    source .venv/bin/activate
pip install -e ".[dev]"
```

## Usage (Input Layer)

Load and inspect a target, with an optional related-sequence panel:

```bash
hush input data/examples/target_example.fasta \
  --panel data/examples/hla_panel_example.fasta
```

List every sequence in any FASTA/PDB/mmCIF file:

```bash
hush inspect data/examples/hla_panel_example.fasta
```

The Input Layer parses FASTA and PDB/mmCIF, validates residues against the
molecule-type alphabet, auto-detects DNA/RNA/protein, and hands downstream
stages a `DesignInput` (target + optional panel). A **protein** target (e.g.
parsed from a PDB structure) is flagged: siRNA/shRNA design needs the coding
transcript, so the protein must first be resolved to its mRNA/CDS.

### Python API

```python
from core.inputs import build_design_input

design = build_design_input(
    "data/examples/target_example.fasta",
    panel_path="data/examples/hla_panel_example.fasta",
)
target = design.require_nucleotide_target()   # guards against protein targets
mrna = target.as_rna()                         # target mRNA for siRNA design
```

## siRNA design (stage 2)

`core.sirna` slides a window along the target mRNA and yields every candidate
site as an `SiRNACandidate` (sense/passenger and antisense/guide strands as
validated RNA sequences, plus its coordinates). The duplex core length is
tunable - default 19 nt, canonical range 19-24 - and the 2-nt 3' overhangs are
treated as a later design choice rather than baked into candidates.

`core.rules` then applies the three classical 2004 rule sets, each reported
separately so no verdict is a black box:

| Rule set | Shape | Verdict |
|---|---|---|
| Ui-Tei | 4-rule pass/fail filter | `satisfied` (all four), `score` 0-4 |
| Reynolds | 8-point additive checklist | `score`, `meets_threshold()` (cutoff 6) |
| Amarzguioui | additive criterion score | `motif_score`, `score`, `meets_threshold()` (cutoff 3) |

```python
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates
from core.rules import ui_tei, reynolds, amarzguioui

target = BioSequence("ATGGCTACGTTAGCCATTGCACGTAAC", id="DEMO", moltype=MolType.DNA)

for candidate in enumerate_candidates(target, length=19):
    u = ui_tei(candidate)
    r = reynolds(candidate)
    a = amarzguioui(candidate)
    if u.satisfied and r.meets_threshold():
        print(candidate.region, candidate.guide.seq, r.score, a.motif_score)
```

### Known deferrals

Two criteria depend on the not-yet-built thermodynamics stage. Rather than
assume or invent a value, both are reported as `pending` and excluded from
their score, so a checklist is never silently under- or over-counted:

- **Reynolds criterion III** (low internal secondary structure) - supply via
  `reynolds(candidate, low_internal_structure=...)`; needs a folding/Tm
  calculation.
- **Amarzguioui duplex-end asymmetry** - supply via
  `amarzguioui(candidate, end_asymmetry=...)`. The paper tabulates this term as
  an integer in `[-2, 3]`, but the mapping from the terminal-triplet A/U
  differential to those points is not confirmed, so no mapping is guessed.

## Development

```bash
pytest            # run the test suite (107 tests)
```
