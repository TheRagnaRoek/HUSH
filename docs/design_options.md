# HUSH `design` options - quick reference

A one-page orientation to the optional parameters of the selection step, as
exposed by the `hush design` command and the `SelectionCriteria` dataclass. Each
row lists the CLI flag, the matching `SelectionCriteria` field, its default, and
what it does. The *why* behind each choice lives in
[`design_notes.md`](design_notes.md) (selection, seed, asymmetry and hairpin
sections); this file is just the map.

Selection is a **filter, not a ranking**: every criterion below is an independent,
optional gate, and a candidate is selected only if it clears every enabled one.
Nothing is scored into a single number, and selected candidates come back in
target-position order.

## Standard options

The everyday design surface - what a normal run turns.

| CLI flag | `SelectionCriteria` field | Default | What it does |
|---|---|---|---|
| `--length` / `-l` | (enumeration arg) | `19` | siRNA **core** length, 19-24 nt (the paired region; overhangs are added on top). |
| `--overhang` | (enumeration arg) | `natural` | 3' overhang on each strand: `natural` (mRNA-derived flanks), `dtdt` (2'-deoxythymidine, shown lowercase `tt`), or `uu`. |
| `--rules` | `use_ui_tei`, `use_reynolds`, `use_amarzguioui` | `all` | Which of the three 2004 rule sets to apply (comma list, or `all`/`none`). |
| `--mode` | `rule_mode` | `all` | Combine enabled rule sets by `all` (intersection) or `any` (union). |
| `--reynolds-threshold` | `reynolds_threshold` | `6` | Reynolds score cut-off (out of 10). |
| `--amarzguioui-threshold` | `amarzguioui_threshold` | `3` | Amarzguioui score cut-off. |
| `--seed` / `--no-seed` | `use_seed_filter` | on | Apply the siDirect seed-duplex off-target Tm filter (both strands). |
| `--seed-tm` | `seed_tm_threshold` | `21.5` | Seed-duplex Tm cut-off in C. |
| `--offtarget` / `--no-offtarget` | `use_offtarget_filter` | off | Transcriptome homology off-target filter (needs a bowtie index; see `--offtarget-db`). Distinct from `--seed`: a database search, not a sequence-intrinsic proxy. |
| `--offtarget-db` | `offtarget_searcher` | - | bowtie index prefix for the reference transcriptome (build with `tools/build_offtarget_db.py`). |
| `--offtarget-min-mismatch` | `offtarget_min_mismatch` | `2` | Reject a candidate whose nearest off-target has fewer than this many mismatches (2 matches siDirect on the PDCD1 benchmark). |
| `--offtarget-max-mismatch` | (searcher depth) | `3` | Mismatch depth the search enumerates (bowtie `-v` allows 0-3). |
| `--offtarget-self-accession` | `offtarget_self_accession` | auto | Accession of the target's own transcript to exclude from hits (auto-detected from the FASTA header). |
| `--bowtie-path` | (backend) | `bowtie` | Path to the bowtie executable. |
| `--gc` / `--no-gc` | `use_gc_filter` | off | Apply a standalone GC-content window (separate from Reynolds' internal GC criterion). |
| `--gc-min` / `--gc-max` | `gc_min` / `gc_max` | `30` / `52` | GC window bounds, in percent. |
| `--homopolymer` / `--no-homopolymer` | `use_homopolymer_filter` | off | Reject candidates whose core has a single-nucleotide run longer than `--max-run` (e.g. `GGGG`). |
| `--max-run` | `max_homopolymer_run` | `3` | Longest tolerated single-nucleotide run in the core; a run over this is rejected. |
| `--asymmetry` / `--no-asymmetry` | `use_asymmetry_filter` | off | Require duplex-end asymmetry to favour loading the guide (sign only, no magnitude threshold). |

### Duplex length and overhangs

HUSH's canonical duplex is a **19-nt paired core with symmetric 2-nt 3' overhangs**,
i.e. 21-nt strands (Elbashir et al. 2001, *EMBO J* 20:6877). `--length` sets the
**core**, not the strand: the three rule sets (Ui-Tei, Reynolds, Amarzguioui) are all
defined on the 19-nt core, so scoring stays faithful to the source papers regardless of
overhang. Overhangs are tracked apart from the core and never reach a filter.

The `natural` default draws each 3' overhang from the target's own flanks (the guide's
overhang is the reverse-complement of the 2 nt upstream of the core - the true offset
Tuschl duplex). A candidate whose flank runs past a transcript end falls back to `dtdt`
there. In the TSV, `guide_5to3` / `passenger_5to3` are the full ordered strands and the
`*_overhang` columns name the overhang explicitly; `start`/`end` are core coordinates.

Note for benchmarking against siDirect: in `natural` mode the **guide** strand still
equals siDirect's listed guide (core + upstream overhang = reverse-complement of the
full window), but the **passenger** is offset (its overhang is the downstream flank), so
compare on the guide or the target position, not the passenger.

## Allele-aware options (MSA panel targeting)

Off unless `--panel` is given. Supplying a panel engages the MSA/allele stage
(stage 3): the panel is star-aligned to the target, and each candidate is
additionally gated on whether its window hits the requested alleles and
discriminates against the spared ones. The allele gate is ANDed with every
selection filter above. The window length is the siRNA core length (`--length`),
so windows and candidates line up one-for-one.

| CLI flag | Maps to | Default | What it does |
|---|---|---|---|
| `--panel` / `-p` | (loads the panel) | - | Allele-family FASTA (e.g. an HLA locus). Its presence turns the allele gate on. |
| `--panel-format` | (format override) | auto | Override the panel format: `fasta` / `pdb` / `cif`. |
| `--pan-allele` | `AlleleTargeting.pan_allele` | on\* | Hit every panel member, spare none (conserved-region targeting). |
| `--allele-specific` | `AlleleTargeting.allele_specific` | - | Hit exactly this panel member id, spare the rest. |
| `--hit` | `AlleleTargeting.hit` | - | Comma list of member ids to hit. With `--spare` both sets are explicit; without it the rest are spared (subset). |
| `--spare` | `AlleleTargeting.spare` | - | Comma list of member ids to discriminate against; only refines an explicit `--hit`. |
| `--min-spare-mismatches` | `WindowCriteria.min_spare_mismatches` | `1` | Differences a spared allele must carry inside a window before it counts as discriminated. |

\* `--pan-allele` is the default intent when `--panel` is given without
`--allele-specific` or `--hit`. The three intents (`--pan-allele`,
`--allele-specific`, `--hit`) are mutually exclusive.

The allele gate is a **filter, not a ranking**, and coverage is required: a window
where a targeted allele is uncovered, or where a spared allele is uncovered,
cannot satisfy the request - missing data is never read as evidence. See the
[`docs/design_notes.md`](design_notes.md) MSA/allele section for the rationale.

## Advanced options (thermodynamic internals)

Correct by default and safe to ignore. These sit in the "Advanced" group of
`hush design --help`.

| CLI flag | `SelectionCriteria` field | Default | What it does |
|---|---|---|---|
| `--hairpin` / `--no-hairpin` | `use_hairpin_filter` | off | Reject candidates whose sense strand folds into a stable internal hairpin. Standalone version of Reynolds' internal-repeat criterion. |
| `--hairpin-tm` | `hairpin_tm_max` | `20.0` | Internal-structure Tm cut-off in C; reject if the most stable self-hairpin melts above this. |
| `--asym-window` | `asymmetry_window` | `4` | Terminal base pairs compared at each duplex end for the asymmetry measurement. |
| `--asym-temp` | `asymmetry_temp_c` | `37.0` | Temperature in C for the asymmetry free energy. |
| `--params` | `thermo_params` | `xia` | Nearest-neighbor parameter set (`xia` / `freier`) for the asymmetry and hairpin readouts **only**. |

## What is deliberately not tunable

Some values are pinned so a published calibration stays valid - exposing them
would let a run silently invalidate a threshold:

- **The seed filter's internals** (Freier parameters, initiation constant,
  concentration, salt). Only the seed *threshold* (`--seed-tm`) is tunable; its
  machinery is fixed, and `--params` never reaches it. siDirect's 21.5 C cut-off
  is only meaningful against exactly these numbers.
- **Reynolds' internal GC criterion** (30-52%) and **internal-repeat criterion**
  (Xia 1998, Tm <= 20 C). The standalone `--gc` and `--hairpin` filters are the
  tunable counterparts; the rule keeps its published form.
- **The asymmetry magnitude threshold** - there is none. The filter uses the sign
  of the differential only; turning kcal/mol into a cut-off needs calibration data
  HUSH does not have.

## Two example runs

Strict default (all three rules, seed filter on) - very selective on a GC-rich
target:

```bash
hush design data/examples/hla_a01.fasta
```

Relaxed, with the standalone structure and GC gates added and the union of rule
sets:

```bash
hush design data/examples/hla_a01.fasta --rules all --mode any \
  --no-seed --gc --hairpin -o selected.tsv
```

With the transcriptome off-target filter (needs a bowtie index built by
`tools/build_offtarget_db.py`):

```bash
hush design data/examples/PDCD1_CDS.fasta --rules ui-tei,reynolds --mode any \
  --offtarget --offtarget-db data/reference/refseq_human_mrna
```

Allele-aware, singling one allele out of a panel (the allele gate is ANDed on
top of the rule and filter gates):

```bash
hush design data/examples/target_example.fasta \
  --panel data/examples/hla_panel_example.fasta \
  --allele-specific "HLA-DEMO*01:01" --mode any --no-seed
```

Every run prints the chosen criteria and a funnel of how many candidates clear
each enabled gate, so an empty result shows exactly which gate removed them.

## Planned / deferred

- **Sequence-motif constraint** (e.g. force a `G` at core position 3 via an
  `NNGNN...`-style pattern). This is primarily a loop-grafting requirement for
  shRNA, where a specific base at the graft junction improves hairpin folding, so
  it will land in the shRNA module rather than the siRNA `design` step. Tracked
  here so the requirement is not lost.
