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
| 2. siRNA Design | `core.sirna` + `core.rules` + `core.thermo` / `core.seed` / `core.asymmetry` / `core.hairpin` + `core.selection` | nearly complete - enumeration (with 3' overhangs), the classical rule pre-screen, the thermodynamics stage, the option-rich selection filter (incl. a homopolymer gate) and a CLI (`hush design`) are implemented and validated against siDirect. Remaining: a transcriptome off-target filter (WIP - see below) |
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

Amarzguioui is complete, including its duplex-end A/U differential (the paper's
`DT3`), which is measured **3' end minus 5' end** over the three terminal base
pairs of the sense strand and enters the score directly. All 46 siRNAs of the
paper's own Table 3 reproduce exactly - every `DT3`, every motif and every total
score - which is what pins the sign convention down.

Reynolds is scored **out of 10**, following the original algorithm rather than
the common "8 equal criteria" simplification: five +1 binary criteria, an
**additive** 3'-terminal A/U count (0 to +5 across positions 15-19), and two
-1 penalties ('G' at 13, 'G'/'C' at 19). Position 19 carries three of these at
once, giving a clean A(+2) / U(+1) / G,C(-1) gradient. A score of 6 predicts a
functional siRNA. Its "absence of internal repeat" criterion is computed from
the sense strand by `core.hairpin` (see below); pass
`reynolds(candidate, no_internal_repeat=...)` to override it with an external
folder's verdict. All three 2004 rule sets are now fully evaluated - nothing is
deferred.

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

## Benchmarking against siDirect, and the off-target gap (WIP)

siDirect is HUSH's gold-standard reference. On matched settings the functional
rule pre-screen reproduces siDirect's picks exactly (validated candidate-for-
candidate on the siDirect sample sequence). `tools/compare_sidirect.py` aligns a
HUSH design TSV against a siDirect result export and reports the difference:

```bash
python tools/compare_sidirect.py \
  --sidirect data/siDirect_PDCD1-CDS.ods --hush data/pdcd1_hush_relaxed.tsv
```

> **Transcriptome off-target filtering is not yet implemented (WIP).** siDirect
> BLASTs each guide against the human transcriptome and drops candidates with a
> near-perfect (0-1 mismatch) off-target. HUSH has no transcriptome search yet,
> so on real targets it selects a superset of siDirect's candidates - the extras
> are exactly those siDirect removed for off-target reasons. This is the next
> task after the siRNA design core is finalised. Note this is distinct from the
> seed-duplex Tm filter (`core.seed`), which is a *seed-mediated* off-target
> proxy and is already implemented.

## Development

```bash
pytest            # run the test suite (295 tests)
```
