# CLI reference

All commands and flags with their defaults. The rationale for each gate is in
[design-rules.md](design-rules.md).

| Command | Purpose |
|---|---|
| `hush inspect` | Preview a FASTA/PDB/mmCIF input, or a target plus panel |
| `hush design` | Enumerate and filter siRNA candidates |
| `hush convert` | Graft selected candidates onto a stem-loop shRNA scaffold |
| `hush shell` | Interactive session: load the inputs once and re-run the design |
| `hush reference build` | Download the transcriptome and build the off-target index |
| `hush reference status` | Report whether the index exists, and its provenance |

## `hush reference`

The transcriptome off-target check is on by default and needs this index. Run the
build once; designs then find the index automatically.

| Flag | Default | Description |
|---|---|---|
| `--outdir` | `$HUSH_REFERENCE_DIR`, else `data/reference` | Directory for the index and manifest |
| `--base-url` | RefSeq `H_sapiens/mRNA_Prot` | Use an archived release to fix the reference version |
| `--index-name` | `refseq_human_mrna` | Basename of the bowtie index prefix |
| `--keep-fasta` | off | Keep the concatenated FASTA after building |

The download is several gigabytes and needs `bowtie-build` on PATH; the container
already contains it. `hush reference status` prints the RefSeq release, bowtie
version and build date from the manifest, and exits with a non-zero status when no
index exists, so it can be used as a readiness check.

## `hush design`

```bash
hush design TARGET.fasta [-o selected.tsv]
```

Enumerates every candidate site and reports those that pass the enabled gates. It
is a **filter, not a ranking**: each candidate keeps its per-criterion breakdown,
so every verdict can be traced. The run prints the criteria, a funnel (the number
of candidates that pass each gate) and the selected table. For an empty result,
the funnel shows which gate removed the candidates.

**Defaults.** `hush design TARGET` applies all three rule sets combined by union,
the seed off-target Tm filter, the homopolymer gate and the transcriptome
off-target search. The GC, hairpin and asymmetry gates are opt-in.

### Standard options

| Flag | `SelectionCriteria` field | Default | Description |
|---|---|---|---|
| `--length` / `-l` | (enumeration arg) | `19` | siRNA **core** length, 19-24 nt. Overhangs are added on top. |
| `--overhang` | (enumeration arg) | `natural` | 3' overhang per strand: `natural` (mRNA-derived flanks), `dtdt` (shown lowercase `tt`), or `uu`. |
| `--rules` | `use_ui_tei`, `use_reynolds`, `use_amarzguioui` | `all` | Which 2004 rule sets to apply (comma list, or `all`/`none`). |
| `--mode` | `rule_mode` | `any` | Combine them by `any` (union) or `all` (intersection). |
| `--reynolds-threshold` | `reynolds_threshold` | `6` | Reynolds score cut-off (out of 10). |
| `--amarzguioui-threshold` | `amarzguioui_threshold` | `3` | Amarzguioui score cut-off. |
| `--seed` / `--no-seed` | `use_seed_filter` | **on** | Seed-duplex off-target Tm filter (Naito et al. 2009), both strands. |
| `--seed-tm` | `seed_tm_threshold` | `21.5` | Seed-duplex Tm cut-off in C. |
| `--offtarget` / `--no-offtarget` | `use_offtarget_filter` | **on** | Transcriptome homology filter. Uses the built reference automatically. |
| `--offtarget-db` | (searcher) | auto | bowtie index prefix. Defaults to the one `hush reference build` produced; set it to use a different reference genome. |
| `--offtarget-min-mismatch` | `offtarget_min_mismatch` | `2` | Reject when the nearest off-target has fewer than this many mismatches. |
| `--offtarget-max-mismatch` | (searcher depth) | `3` | Mismatch depth searched (bowtie `-v` allows 0-3). |
| `--offtarget-self-accession` | `offtarget_self_accessions` | - | Accessions of the design's own transcripts, excluded from hits (comma-separated, version-insensitive). Added to the header accession and to the transcripts found by homology. |
| `--offtarget-self-fasta` | `offtarget_self_fasta` | - | FASTA whose entries count as the design's own sequence and are resolved to transcripts by homology. For an allele family that the target and panel do not cover. |
| `--offtarget-resolve-self` / `--no-offtarget-resolve-self` | `resolve_offtarget_self` | **on** | Resolve the owned transcripts by homology to the target, the panel and `--offtarget-self-fasta`. This is the only option that covers a whole gene family. |
| `--offtarget-self-coverage` | `offtarget_self_coverage` | `10` | Percent of the k-mers of a design sequence that a transcript must match exactly to count as owned. A lower value includes more of the family; a higher value restricts the set to the target gene. |
| `--bowtie-path` | (backend) | `bowtie` | Path to the bowtie executable. CLI only; not available over HTTP. |
| `--gc` / `--no-gc` | `use_gc_filter` | off | Standalone GC-content window. |
| `--gc-min` / `--gc-max` | `gc_min` / `gc_max` | `30` / `52` | GC window bounds, in percent. |
| `--forbid-homopolymers` / `--allow-homopolymers` | `use_homopolymer_filter` | **forbid** | Reject a core with a single-nucleotide run over `--max-run`. |
| `--max-run` | `max_homopolymer_run` | `3` | Longest tolerated run. |
| `--asymmetry` / `--no-asymmetry` | `use_asymmetry_filter` | off | Require duplex-end asymmetry to favour the guide (sign only). |

### Allele-aware options

Off unless `--panel` is given. The allele gate is combined with every gate above
by AND. The window length equals `--length`, so each candidate has exactly one
window.

| Flag | Default | Description |
|---|---|---|
| `--panel` / `-p` | - | Allele-family FASTA. Its presence turns the allele gate on. |
| `--panel-format` | auto | Override: `fasta` / `pdb` / `cif`. |
| `--pan-allele` | on\* | Hit every panel member, spare none (conserved-region targeting). |
| `--allele-specific ID` | - | Hit exactly this member, spare the rest. Use when the allele to silence is itself in the panel. |
| `--target-specific` | - | Spare **every** member and hit none. Use when the target is loaded apart from its reference panel. |
| `--hit` | - | Comma list of member ids to hit. |
| `--spare` | - | Comma list to discriminate against; refines an explicit `--hit`. |
| `--min-spare-mismatches` | `1` | Differences a spared allele needs inside a window to count as discriminated. |
| `--sort` | `auto` | `auto` (rank by discrimination when something is spared, else position), `position`, or `discrimination` (errors when nothing is spared). |
| `--explain` | off | Print the full per-allele discrimination breakdown, not just the weakest. |

\* The default intent when `--panel` is given without another targeting flag. The
four intents are mutually exclusive.

The `disc` readout (for example `P16 post H+ S pp`) is explained in
[design-rules.md](design-rules.md#discrimination-ranking).

### Advanced options

The defaults are suitable for normal use.

| Flag | Default | Description |
|---|---|---|
| `--hairpin` / `--no-hairpin` | off | Reject a sense strand that folds into a stable internal hairpin. |
| `--hairpin-tm` | `20.0` | Internal-structure Tm cut-off in C. |
| `--asym-window` | `4` | Terminal base pairs compared at each duplex end. |
| `--asym-temp` | `37.0` | Temperature in C for the asymmetry free energy. |
| `--params` | `xia` | Nearest-neighbor set (`xia` / `freier`) for the asymmetry and hairpin read-outs **only**. |

### Fixed values

Some values cannot be changed because they belong to a published calibration:

- **The internals of the seed filter** (Freier parameters, initiation constant,
  concentration, salt). Only `--seed-tm` can be changed, and `--params` does not
  affect the seed filter: the 21.5 C cut-off of Naito et al. 2009 is only valid
  with those values.
- **The internal GC (30-52%) and internal-repeat (Tm <= 20 C) criteria of
  Reynolds.** The `--gc` and `--hairpin` gates are the tunable counterparts; the
  rule keeps its published form.
- **The asymmetry magnitude threshold.** There is no magnitude threshold. Only the
  sign is used, because a kcal/mol cut-off would need calibration data that HUSH
  does not have.

## `hush convert`

```bash
hush convert selected.tsv [-o shrna.tsv]
```

Reads a `hush design -o` TSV and grafts each candidate's **full strands** onto a
stem-loop scaffold. See [shrna-conversion.md](shrna-conversion.md).

| Flag | Default | Description |
|---|---|---|
| `--loop` | `jensen` | Loop scaffold: `jensen` (UGUGCUU), `psuper` (UUCAAGAGA), `trc` (CUCGAG). |
| `--terminator` / `--no-terminator` | on | Append the Pol III poly-T terminator. |
| `--terminator-length` | `6` | Length of that run (>= 4 when included). |
| `--drop-premature-terminators` / `--keep-premature-terminators` | **drop** | Withhold constructs whose body carries a premature Pol III terminator (>= 4 U). |
| `--limit` / `-n` | `20` | Max constructs to print (0 = all). |

The input TSV must carry `region`, `start`, `end`, `passenger_5to3`,
`passenger_overhang` and `guide_overhang`. The guide overhang has its own column
because under the `natural` policy it comes from the opposite (upstream) flank. It
is therefore not the reverse complement of the passenger overhang and cannot be
derived from it.

## `hush shell`

Keeps one session open: the target, the panel and the panel alignment are loaded
once and reused.

| Command | Purpose |
|---|---|
| `load <path>` | Load the target mRNA/CDS |
| `panel <path>` | Load a panel (turns on allele targeting) |
| `set [flags]` | Set persistent criteria (same flags as `design`); no flags prints them |
| `design [flags]` | Design with the current criteria (inline flags apply to that run only) |
| `convert [flags]` | Convert the last design's selection |
| `assess [region]` | Per-allele discrimination for the last design |
| `save <dir>` | Write `design.tsv` and `shrna.tsv` |
| `show` / `inspect` | Review criteria / preview the input |
| `limit <n>`, `reset`, `help`, `quit` | Display cap, restore defaults, help, exit |

## Examples

```bash
# Default settings. On a GC-rich CDS the seed filter may select nothing.
hush design data/hla_a01.fasta

# Relaxed for a GC-rich target, with the standalone structure gates
hush design data/hla_a01.fasta --no-seed --gc --hairpin -o selected.tsv

# Allele-specific: hit one panel member, discriminate against the rest
hush design data/examples/target_example.fasta --no-seed \
  --panel data/examples/hla_panel_example.fasta \
  --allele-specific "HLA-DEMO*01:01"

# Target-specific, with the full per-allele breakdown
hush design my_target.fasta --no-seed --target-specific \
  --panel reference_panel.fasta --explain

# Convert, keeping constructs that carry a premature terminator
hush convert selected.tsv --loop psuper --keep-premature-terminators
```
