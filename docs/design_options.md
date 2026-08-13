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
| `--length` / `-l` | (enumeration arg) | `19` | siRNA core length, 19-24 nt. |
| `--rules` | `use_ui_tei`, `use_reynolds`, `use_amarzguioui` | `all` | Which of the three 2004 rule sets to apply (comma list, or `all`/`none`). |
| `--mode` | `rule_mode` | `all` | Combine enabled rule sets by `all` (intersection) or `any` (union). |
| `--reynolds-threshold` | `reynolds_threshold` | `6` | Reynolds score cut-off (out of 10). |
| `--amarzguioui-threshold` | `amarzguioui_threshold` | `3` | Amarzguioui score cut-off. |
| `--seed` / `--no-seed` | `use_seed_filter` | on | Apply the siDirect seed-duplex off-target Tm filter (both strands). |
| `--seed-tm` | `seed_tm_threshold` | `21.5` | Seed-duplex Tm cut-off in C. |
| `--gc` / `--no-gc` | `use_gc_filter` | off | Apply a standalone GC-content window (separate from Reynolds' internal GC criterion). |
| `--gc-min` / `--gc-max` | `gc_min` / `gc_max` | `30` / `52` | GC window bounds, in percent. |
| `--asymmetry` / `--no-asymmetry` | `use_asymmetry_filter` | off | Require duplex-end asymmetry to favour loading the guide (sign only, no magnitude threshold). |

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

Every run prints the chosen criteria and a funnel of how many candidates clear
each enabled gate, so an empty result shows exactly which gate removed them.
