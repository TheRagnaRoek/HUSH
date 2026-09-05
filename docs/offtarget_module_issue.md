# Issue: Transcriptome off-target filter (`core.offtarget`)

> **RESOLVED (2026-08-15).** Implemented in `core/offtarget.py` (bowtie1 `-v`
> backend + a pure-Python brute-force backend behind one searcher protocol),
> wired into `core.selection` and `hush design` (`--offtarget*`), with a RefSeq
> index builder (`tools/build_offtarget_db.py`), `tools/compare_sidirect.py`
> extended to check the off-target column, and tests in `tests/test_offtarget.py`.
> Design record: [`offtarget_design_note.md`](offtarget_design_note.md) and the
> `core.offtarget` section of [`design_notes.md`](design_notes.md). The original
> hand-off spec is kept below for provenance.

> Hand-off spec for a dedicated implementation session. Self-contained: it
> assumes only the current `main`/feature-branch state of HUSH. Paste into a
> GitHub issue and/or point a fresh agent session at this file.

## Context

HUSH is a from-sequence siRNA/shRNA design pipeline (see `README.md`). The
siRNA design core is otherwise complete: window enumeration with 3' overhangs
(`core.sirna`), the classical rule pre-screen (`core.rules`), thermodynamics
(`core.thermo`/`seed`/`asymmetry`/`hairpin`), and the option-rich selection
filter (`core.selection`), all driven by `hush design` (`app/cli.py`).

The **one missing piece** is a transcriptome homology off-target filter. Note
this is **distinct** from `core.seed`, which already implements siDirect's
*seed-mediated* off-target proxy (a seed-duplex Tm cut-off). That is a
sequence-intrinsic property of the guide's own seed; it is **not** a search for
where the guide could bind elsewhere in the transcriptome. This issue is the
latter: direct sequence-homology off-target detection.

## The validated gap (motivation)

Benchmarking HUSH against siDirect on the PDCD1 CDS (`data/examples/PDCD1_CDS.fasta`,
RefSeq NM_005018.3) with matched functional-rule settings gives a **perfect
subset**: siDirect's picks are all in HUSH's set, and HUSH keeps 5 extra
candidates. Reproduce with:

```bash
python tools/compare_sidirect.py \
  --sidirect data/siDirect_PDCD1-CDS.ods --hush data/pdcd1_hush_relaxed.tsv
```

The 5 HUSH-only candidates (by core region on the CDS):

| region  | HUSH rules | note |
|---------|-----------|------|
| 160-178 | `-R-`     | Reynolds-only |
| 391-409 | `-R-`     | Reynolds-only |
| **650-668** | **`URA`** | **passes all three rules + GC + homopolymer** |
| 651-669 | `-R-`     | Reynolds-only |
| 744-762 | `-R-`     | Reynolds-only |

`650-668` is the decisive clue: it clears every functional rule and every
sequence-composition gate, so no rule difference can explain siDirect dropping
it. siDirect's own output columns show why - every candidate siDirect **kept**
has "min. number of mismatches against off-targets" >= 2, i.e. siDirect BLASTs
each guide against the human transcriptome and rejects candidates with a
near-perfect (0-1 mismatch) off-target. HUSH cannot see those, so it keeps them.
Closing this gap should make the two tools agree.

## Goal

Implement `core.offtarget`: for each candidate, search a reference human
transcriptome for near-complementary sites, report the off-target profile
(min mismatch, hit counts by mismatch level, both strands, on-target excluded),
and expose an optional selection gate that rejects candidates whose nearest
off-target is below a configurable mismatch threshold.

## Requirements

Functional:
- Input: a guide sequence (and optionally the passenger - siDirect reports both).
- Search a reference transcriptome for hits at 0/1/2/3 mismatches on both strands.
- Exclude the intended on-target transcript from the hit set.
- Report per candidate: `min_mismatch` to the nearest off-target and hit counts
  by mismatch level (mirror siDirect's columns so `tools/compare_sidirect.py`
  can check them directly).
- Provide a gate: reject if `min_mismatch < threshold` (configurable).

Non-functional:
- Deterministic, offline, reproducible - a local reference DB, pinned build.
- Tractable on ~850 candidates (an 867-nt CDS yields 849 19-mers). Needs a
  seed-and-extend index or an external aligner; brute force over the whole
  transcriptome per candidate will not scale.

## Architecture / integration (follow the existing patterns)

- **New module `core/offtarget.py`**, a sibling to `core/seed.py`. Read
  `core/seed.py` and `docs/seed_offtarget_design_note.md` first for the house
  style and the off-target philosophy already written down.
- **Selection gate** in `core.selection.SelectionCriteria` - copy the
  `use_seed_filter` pattern exactly: add `use_offtarget_filter`,
  `offtarget_min_mismatch`, and a DB handle/path; validate in `__post_init__`;
  compute in `evaluate()`; add `offtarget_*` fields to `CandidateEvaluation`;
  extend `_reasons()`; add a funnel row.
- **CLI** in `app/cli.py` `design` command - mirror how `--seed`/`--seed-tm`
  are wired: `--offtarget/--no-offtarget`, `--offtarget-db`,
  `--offtarget-min-mismatch`; surface it in `_print_criteria`, `_print_funnel`,
  the table, and `_write_design_tsv` (add min-mismatch + hit-count columns).
- **Reference data** under `data/reference/` (already `.gitignore`d for large
  downloads - see the tail of `.gitignore`). Provide a fetch/build helper; do
  not commit the transcriptome.

## Design decisions to make (raise with the maintainer)

1. **Search method**: external aligner (bowtie/BLAST + local index) vs a
   pure-Python seed-and-extend k-mer index vs hybrid. Trade-off: dependency
   weight vs speed vs fidelity. HUSH is currently light (biopython/numpy/pandas/
   sklearn/typer) - adding a binary aligner is a real dependency shift.
2. **Reference transcriptome**: RefSeq human mRNA vs Ensembl cDNA vs GENCODE;
   pin the exact build/version for reproducibility.
3. **Guide only, or guide + passenger** (siDirect scores both).
4. **Threshold semantics**: a flat `min_mismatch >= N`, or position-weighted so
   seed-region mismatches count more (siDirect weights the seed). Decide whether
   and how this composes with the existing `core.seed` Tm filter.
5. **On-target exclusion**: how to identify and drop the intended transcript
   (accession match? self-alignment at 0 mismatch over full length?).

## Acceptance criteria

- On the PDCD1 benchmark, enabling the off-target gate with siDirect-matched
  settings drops the 5 HUSH-only candidates, so `tools/compare_sidirect.py`
  reports `HUSH-ONLY = 0` (or a documented, explained remainder).
- The reported `min_mismatch`/hit-count values match siDirect's columns for the
  7 shared candidates (ground truth is in `data/siDirect_PDCD1-CDS.ods`).
- Unit tests: a guide with a planted off-target in a tiny synthetic reference;
  on-target exclusion; threshold boundary behaviour; both-strand hits.
- Docs: flip the WIP notes to "implemented" in `README.md` (benchmarking
  section) and `docs/design_options.md` (Planned/deferred -> options table);
  add a `core.offtarget` section to `docs/design_notes.md`.

## Out of scope

- Genomic (non-transcript) off-targets.
- Seed-mediated off-target scoring beyond siDirect's method - `core.seed`
  already covers that angle.

## Pointers

- siDirect 2.0: Naito et al. 2009 (cited in `docs/design_notes.md`); the
  off-target columns are defined in its result export.
- Sibling module and its design note: `core/seed.py`,
  `docs/seed_offtarget_design_note.md`.
- Validation tool + benchmark fixtures: `tools/compare_sidirect.py`,
  `data/siDirect_PDCD1-CDS.ods`, `data/pdcd1_hush_relaxed.tsv`.
