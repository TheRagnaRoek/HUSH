# HUSH

**H**uman genome **U**niversal **S**ilencing **H**airpin-designer - a modern,
HLA-aware siRNA/shRNA design pipeline.

HUSH is an open-source, from-sequence RNAi design tool: it takes a target gene
sequence (and, optionally, a panel of related sequences such as an HLA allele
family), designs siRNA candidates, reasons about conserved vs. polymorphic
regions for broad or allele-specific knockdown, and converts validated
candidates into cloneable shRNA constructs. See
[`docs/HUSH_concept_dossier.md`](docs/HUSH_concept_dossier.md) for the full
rationale and prior-art landscape.

## Pipeline

| Stage | Module | Status |
|---|---|---|
| 1. Input Layer | `core.inputs` + `core.sequence` | implemented |
| 2. siRNA Design | `core` (siRNA engine) | next |
| 3. MSA / Allele | `msa` | planned |
| 4. shRNA Conversion | `shrna` | planned |
| 5. Output Layer | `core` / `app` | planned |

## Project layout

```
core/          Shared domain model + the siRNA design engine
  sequence.py    BioSequence / SequencePanel / DesignInput (pure stdlib)
  inputs.py      Input Layer: FASTA & PDB parsing -> domain model
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

## Development

```bash
pytest            # run the test suite
```
