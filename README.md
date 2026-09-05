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
| 2. siRNA Design | `core.sirna` + `core.rules` + `core.thermo` / `core.seed` / `core.offtarget` / `core.asymmetry` / `core.hairpin` + `core.selection` | complete - enumeration (with 3' overhangs), the classical rule pre-screen, the thermodynamics stage, the seed and transcriptome off-target filters, the option-rich selection filter (incl. a homopolymer gate) and a CLI (`hush design`) are implemented and validated against siDirect |
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
  thermo.py      Nearest-neighbor duplex thermodynamics (dH/dS -> dG, Tm)
  seed.py        siDirect seed-duplex Tm off-target filter
  offtarget.py   Transcriptome homology off-target filter (bowtie1 / brute force)
  asymmetry.py   Duplex-end thermodynamic asymmetry (RISC strand loading)
  hairpin.py     Intramolecular hairpin Tm (Reynolds internal-structure rule)
  selection.py   Option-rich candidate filter (which rules + filters to apply)
  errors.py      Shared exception types
msa/           MSA / allele module (stage 3, placeholder)
shrna/         shRNA conversion module (stage 4, placeholder)
app/           CLI (app.cli) - the CLI-first interface
data/examples/ Runnable example inputs (synthetic panels plus real gene fixtures
               such as PDCD1/PD-1 for siDirect benchmarking)
tools/         Developer utilities (e.g. compare_sidirect.py, HUSH-vs-siDirect diff)
tests/         pytest suite
docs/          Concept dossier and design notes
```

## Getting started

```bash
python -m venv .venv
# Windows:  .venv\Scripts\activate
# POSIX:    source .venv/bin/activate
pip install -e ".[dev]"
```

HUSH is pure-pip and needs no external binaries for its core. The one optional
dependency is **bowtie1** (a separate binary, e.g. via bioconda or your package
manager), needed only for transcriptome-scale off-target search (`hush design
--offtarget`); everything else, including small-panel off-target search, runs
without it.

## Usage (Input Layer)
### Running without activating the venv

Once installed, you do not have to activate the venv every time. The repo ships
thin launchers that call the venv's interpreter directly, so `hush` works from
the project folder (or by full path from anywhere):

```bash
# Windows:  .\hush inspect data/examples/hla_panel_example.fasta
# POSIX:    ./hush  inspect data/examples/hla_panel_example.fasta
```

## Usage (Input Layer)

List every sequence in any FASTA/PDB/mmCIF file:

```bash
hush inspect data/examples/hla_panel_example.fasta
```

Add `--panel` to treat the file as the design target and preview what the
pipeline will receive - target/panel roles, molecule types, and a readiness
check:

```bash
hush inspect data/examples/target_example.fasta \
  --panel data/examples/hla_panel_example.fasta
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
tunable - default 19 nt, canonical range 19-24 - and each candidate carries a
2-nt 3' overhang on each strand, tracked apart from the paired core so no filter
sees it. Overhangs default to the target's natural flanks (`--overhang natural`;
the offset Tuschl duplex) and can be switched to `dtdt` or `uu`.

`core.rules` then applies the three classical 2004 rule sets, each reported
separately so no verdict is a black box:

| Rule set | Shape | Verdict |
|---|---|---|
| Ui-Tei | 4-rule pass/fail filter | `satisfied` (all four), `score` 0-4 |
| Reynolds | additive score out of 10 | `score`, `positive_score`, `penalty`, `meets_threshold()` (cutoff 6) |
| Amarzguioui | additive criterion score | `motif_score`, `end_au_differential`, `score`, `meets_threshold()` (cutoff 3) |

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

## Thermodynamics (stage 2)

`core.thermo` is the nearest-neighbor engine every thermodynamic step builds on:
sequence to `(dH, dS)`, and from there dG at any temperature and Tm. Two entry
points are kept deliberately separate - `duplex_thermo()` for a complete,
isolated duplex (initiation, terminal A-U penalties and symmetry included), and
`stacking_terms()` for a window *inside* a longer helix, where the edges are not
real helix ends. Two parameter sets are exposed by name, Xia 1998 (default) and
Freier 1986, because a published threshold is only valid against the set it was
calibrated on.

`core.seed` applies the siDirect 2.0 seed-duplex Tm off-target filter to both
strands (guide and passenger positions 2-8), rejecting candidates whose seed
pairs too stably with its target. It reproduces siDirect's published figures
exactly: 4488 of the 16384 possible heptamers (27.4%) below the 21.5 C
threshold.

```python
from core.seed import seed_tm

result = seed_tm(candidate)          # both strands, default threshold 21.5 C
result.passes                        # True only if BOTH strands clear it
result.guide_tm, result.passenger_tm
```

`core.asymmetry` measures duplex-end thermodynamic asymmetry directly, in
kcal/mol, rather than through the A/U proxies the 2004 rule sets use. It
compares the stacking free energy of the terminal window at each end of the
duplex; a **positive** differential means the guide's 5' end is the looser one,
which is the orientation that loads the intended strand.

```python
from core.asymmetry import end_asymmetry

a = end_asymmetry(candidate)         # default 4-bp windows, 37 C
a.differential                       # kcal/mol; positive favours the guide
a.preferred_strand                   # "guide" | "passenger" | "none"
```

The window length is a parameter (default 4 bp) because no single canonical
value exists in the literature, and **no pass/fail threshold is defined** -
turning kcal/mol into a cut-off needs calibration data HUSH does not yet have.
The sign alone already carries the mechanistic claim.

`core.hairpin` estimates the most stable intramolecular hairpin an siRNA can
fold into, and its **intramolecular** melting temperature (`Tm = dH/dS`, no
concentration term). It supplies Reynolds' "absence of internal repeat"
criterion (`no_internal_repeat = Tm <= 20 C`). The stem energy comes from the
Xia NN engine; loop destabilisation comes from Groebe & Uhlenbeck (1988). This
is a **transparent approximation**, not a reproduction of Reynolds' unspecified
folder: perfect Watson-Crick stems only (no G-U or bulges), and the one modelling
assumption the ΔG₃₇-only loop data forces - loop destabilisation treated as
entropic - is documented in the module.

```python
from core.hairpin import internal_structure_tm, has_no_internal_repeat

internal_structure_tm(candidate.sense.seq)   # degrees C, or None if nothing folds
has_no_internal_repeat(candidate.sense.seq)   # True when Tm <= 20 C
```

> **The seed filter is a benchmark, not the intended final answer.** It is
> believed to be over-conservative, and on a GC-rich polymorphic target it can
> reject every candidate - it returns nothing for the HLA-A\*01:01 CDS in
> `data/examples`. It is kept siDirect-faithful so improvements can be measured
> against it. The critique, the supporting measurements, the proposed
> replacement and the experiment that would justify building it are recorded in
> [`docs/seed_offtarget_design_note.md`](docs/seed_offtarget_design_note.md).

## Candidate selection (stage 2)

`core.selection` applies a user's design choices to the enumerated candidates.
It is a **filter, not a ranking**: there is no single combined score, because
the signals live on incommensurable scales and any weighting would need efficacy
data HUSH does not have. Instead - in the spirit of siDirect - every criterion
is an explicit, optional gate, and each candidate keeps its full per-criterion
breakdown so the reason for any verdict is visible.

- **Rule sets** (Ui-Tei, Reynolds, Amarzguioui) toggle independently and combine
  by `RuleMode.ALL` (must satisfy every enabled set) or `RuleMode.ANY` (kept if
  any enabled set passes). No rule enabled means the rule step adds no constraint.
- **Independent filters**, ANDed on top, each optional with tunable values: a GC
  window, the seed-duplex off-target Tm, a standalone internal-structure
  (self-hairpin) gate, a homopolymer-run gate (reject a single-nucleotide run
  over `max_homopolymer_run`, e.g. `GGGG`), and the loading asymmetry. The GC
  filter is separate from Reynolds' internal criterion I (fixed at the published
  30-52%), and the hairpin gate is separate from Reynolds' internal-repeat
  sub-criterion (Xia-pinned at 20 C) - both here are the tunable versions.
- **Advanced thermodynamic knobs**, correct by default and safe to ignore: the
  asymmetry window (default 4 bp) and temperature, and the nearest-neighbor
  parameter set (Xia 1998 / Freier 1986) for the asymmetry and hairpin readouts.
  The parameter set never touches the seed filter, whose 21.5 C threshold is only
  valid against its own pinned Freier internals.

```python
from core.sirna import enumerate_candidates
from core.selection import SelectionCriteria, RuleMode, select
from core.thermo import ThermoParams

criteria = SelectionCriteria(
    use_ui_tei=True, use_reynolds=True, use_amarzguioui=False,
    rule_mode=RuleMode.ANY,          # union of the enabled rule sets
    use_seed_filter=True, seed_tm_threshold=21.5,
    use_gc_filter=True, gc_min=0.30, gc_max=0.52,
    use_hairpin_filter=True, hairpin_tm_max=20.0,
    use_asymmetry_filter=True,       # advanced knobs below keep their defaults
    asymmetry_window=4, thermo_params=ThermoParams.XIA_1998,
)
for ev in select(enumerate_candidates(target), criteria):
    print(ev.candidate.region, ev.candidate.guide.seq)
```

Selected candidates come back in target-position order - nothing is ranked or
reordered. On the GC-rich HLA-A\*01:01 CDS the strict default (all three rules,
seed on) selects nothing, while relaxing to `ANY` and dropping the seed filter
recovers candidates: the option-richness is the point.

The same options drive the `hush design` command:

```bash
hush design data/examples/hla_a01.fasta --rules ui-tei,reynolds --mode any \
  --seed-tm 21.5 --gc --gc-min 30 --gc-max 52 --hairpin --homopolymer \
  --overhang natural -o selected.tsv
```

It prints the chosen criteria, a funnel of how many candidates clear each enabled
gate, and the selected candidates (guide + passenger 5'->3', GC%, which rules
passed, seed Tm, internal-structure Tm); `-o` writes them all to a TSV. The
expert thermodynamic knobs (`--hairpin-tm`, `--asym-window`, `--asym-temp`,
`--params`) live in an "Advanced" group; run `hush design --help` for the full
option list, or see [`docs/design_options.md`](docs/design_options.md) for a
one-page reference to every option, its default, and what it does.

## Transcriptome off-target filter (stage 2)

`core.offtarget` is HUSH's homology off-target search - the open, offline
counterpart of siDirect's transcriptome step. For each candidate it searches a
reference transcriptome for near-complementary sites (0-3 mismatches, both
strands, the intended transcript excluded) and gates out candidates whose nearest
off-target is too close. This is **distinct** from `core.seed`: the seed filter is
a sequence-intrinsic proxy that never looks at a transcript, while this is a
database search. The two are orthogonal and compose.

The search runs on **bowtie1** in `-v` mode - the open equivalent of the GGGenome
k-mismatch search siDirect relies on, and the right tool for exact full-length
mismatch counting (not BLAST, a local heuristic). A pure-Python `BruteForceSearcher`
covers small references (a curated allele panel, and the tests) with no external
tool; `BowtieSearcher` scales to the whole transcriptome. Build a reference index
first (a large, one-off download):

```bash
python tools/build_offtarget_db.py          # RefSeq human mRNA -> data/reference/
```

```python
from core.offtarget import BowtieSearcher, evaluate_offtarget

searcher = BowtieSearcher("data/reference/refseq_human_mrna", max_mismatch=3)
profile = evaluate_offtarget(candidate, searcher, threshold=2,
                             self_accession="NM_005018.3")
profile.min_mismatch      # nearest off-target, or None if none within depth
profile.passes            # True unless a nearer-than-threshold off-target exists
```

Rationale (why bowtie1 and not BLAST, the RefSeq pinning, strand semantics, the
allele-aware caveat) is in
[`docs/offtarget_design_note.md`](docs/offtarget_design_note.md).

## Benchmarking against siDirect

siDirect is HUSH's gold-standard reference. On matched settings the functional
rule pre-screen reproduces siDirect's picks exactly (validated candidate-for-
candidate on the siDirect sample sequence). `tools/compare_sidirect.py` aligns a
HUSH design TSV against a siDirect result export and reports the difference:

```bash
python tools/compare_sidirect.py \
  --sidirect data/siDirect_PDCD1-CDS.ods --hush data/pdcd1_hush_relaxed.tsv
```

Without the off-target filter, HUSH selects a **superset** of siDirect's
candidates - the extras are exactly those siDirect removed for near-perfect
off-targets. Enabling `hush design --offtarget` (with a built index) closes that
gap, dropping the extra candidates so the two tools agree.

## Development

```bash
pytest            # run the test suite (321 tests)
```
