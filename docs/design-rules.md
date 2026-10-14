# Design rules

What each gate tests, why it exists, and what it rests on. Implementation lives
in `core/` and `msa/`; the flags are in [cli.md](cli.md).

Selection is a **filter, not a score**. Every criterion below is an independent,
optional gate, and a candidate is selected only if it clears every enabled one.
Nothing is combined into a single number. This is a deliberate choice: a
composite score hides which property actually disqualified a candidate, and the
three rule sets below were each published with their own thresholds, which a
weighted blend would silently invalidate.

## The duplex

HUSH's canonical duplex is a **19-nt paired core with symmetric 2-nt 3'
overhangs**, i.e. 21-nt synthesised strands (Elbashir et al. 2001, *Nature*
411:494; *EMBO J* 20:6877).

`--length` sets the **core**, not the strand. All three rule sets are defined on
the core, so scoring stays faithful to the source papers regardless of overhang
policy, and overhangs are tracked apart from the core so no gate scores them.

The `natural` default draws each overhang from the target's own flanks: the
passenger's from the 2 nt 3' of the core, the guide's from the reverse-complement
of the 2 nt 5' of it - the true offset Tuschl duplex. A candidate whose flank
runs past a transcript end falls back to `dtdt` there. `dtdt` is 2'-deoxy­thymidine,
tolerated at the 3' end and the common synthetic default; `uu` is ribo-uridine.

Because the two overhangs come from opposite flanks, **the guide's overhang is
not the reverse-complement of the passenger's**. Anything reconstructing a duplex
must read both, not derive one from the other.

## The three rule sets

All three are from 2004 and are reimplemented faithfully in `core.rules`. They
are combined by **union** (`--mode any`) by default, each weighted equally.

**Ui-Tei** (Ui-Tei et al. 2004, *NAR* 32:936) - a strict four-rule pass/fail
filter. It is also the source of a tension HUSH documents rather than resolves:
its rule (ii) keeps the passenger seed GC-rich while rule (iii) pulls the guide
seed Tm down.

**Reynolds** (Reynolds et al. 2004, *Nat Biotechnol* 22:326) - a score out of 10,
implemented following the original algorithm rather than the widespread "eight
equal binary criteria" simplification. Default cut-off 6.

**Amarzguioui** (Amarzguioui & Prydz 2004, *BBRC* 316:1050) - an additive motif
score plus the raw A/U differential. Default cut-off 3. Its Table 3 (46 siRNAs)
is HUSH's exact-reproduction validation fixture: `tests/test_asymmetry.py`
reproduces the authors' own scores.

## Seed-duplex off-target Tm

**On by default.** `core.seed`, after siDirect (Naito et al. 2009, *BMC
Bioinformatics* 10:392).

Seed-mediated off-target silencing tracks the thermodynamic stability of the
seed duplex (Ui-Tei et al. 2008, *NAR* 36:7100): a weakly-pairing seed is largely
off-target free even where its matches exist. HUSH computes the seed-duplex
melting temperature for both strands and rejects a candidate above the threshold
(21.5 C by default).

Two things about this filter are worth knowing.

It is **deliberately strict**, and on a very GC-rich CDS such as HLA it can
reject every candidate. That is the documented, expected behaviour, not a bug;
relax it with `--no-seed`.

Its **physics is contested**. The seed duplex does not form in solution: Ago
reshapes its guide's binding properties, and seed pairing proceeds 50-200x faster
than annealing of the equivalent naked RNAs (Salomon et al. 2015, *Cell* 162:84;
Wee et al. 2012, *Cell* 151:1055). HUSH keeps the filter because it is the
industry standard and reproduces siDirect, but treats it as a calibrated
empirical proxy rather than a mechanism.

Consequently the seed filter's internals are **not tunable**: the Freier
parameter set, the initiation constant, the concentration and the salt are
pinned, and `--params` never reaches them. Only the threshold moves. siDirect's
21.5 C cut-off is meaningful only against exactly those numbers.

## Transcriptome off-target

**On by default.** `core.offtarget`, via `bowtie1` in `-v` mode over a prebuilt
index.

This is the safety floor and the one gate that consults reality rather than a
sequence-intrinsic proxy: it asks whether the guide has a near-match anywhere
else in the human transcriptome. A candidate is rejected when its nearest
off-target carries fewer than `--offtarget-min-mismatch` mismatches (2 by
default, which matches siDirect on the PDCD1 benchmark).

### What the design owns

"Anywhere else" has to be defined, and for a gene family it cannot be defined by
one accession. A guide is meant to match its target; on HLA it is also meant to
match the family, and the reference carries the locus more than once. HLA-A is
five RefSeq transcripts (`NM_002116`, `NM_001242758`, `NR_001434`,
`XM_041680767`, `XM_041680768`) before any alternate-haplotype assembly copy is
counted. A gate that excludes one of them measures on-target homology, not
off-target risk: with nothing excluded, every candidate on HLA-A\*01, A\*02 and
A\*03 is rejected, and naming a single accession by hand lifts seven production
scenarios from 3 to 5 surviving candidates out of 107. That is not a quality
signal, it is the wrong question.

So the owned set is a set (`core.offtarget.SelfTranscripts`), filled from three
sources, unioned:

* `--offtarget-self-accession`, a comma-separated list of accessions, matched
  version-insensitively (`NM_002116` covers `NM_002116.8`).
* the RefSeq accession in the target's own FASTA header, when it has one. IMGT
  allele headers such as `HLA-A02:01:01:01` do not, which is why this route
  alone left HLA runs with nothing excluded.
* homology resolution (on by default, `--offtarget-resolve-self`): the target,
  every panel member and any `--offtarget-self-fasta` entry is tiled into its
  k-mers, and a reference transcript is owned when at least
  `--offtarget-self-coverage` percent of one input's distinct k-mers align to it
  exactly. Coverage is judged per input, so a short panel member is not diluted
  by a long target.

Panel members seed resolution because they are intended targets: a pan-allele
design means to hit them. This does not weaken an allele-specific design -
discrimination between alleles is the [allele stage](#allele-targeting)'s
job, enforced at a resolution the transcriptome does not have, since every
HLA-A allele maps onto the same RefSeq transcripts.

The default coverage of 10% comes from the measured distribution on the four
bundled HLA-A alleles. Exact 19-mer coverage of the reference by one allele is
sharply clustered: 5 transcripts above 25% (HLA-A itself, the same five on all
four alleles), roughly 100-140 between 10% and 25% (alternate-haplotype copies
of the locus and the closest class I paralogs), then a tail of transcripts
sharing one or two k-mers, which are genuine off-targets and stay in. Raising
the threshold above 25% narrows the set to the target gene; lowering it below 5%
starts absorbing incidental matches. Sequence is the only signal available, so
an alternate-haplotype copy of the target gene and a true paralog are separated
only by how much of the target they carry - a bowtie index stores accessions,
not gene symbols.

On the seven original wet-lab scenarios this is the difference between a gate
that rejects everything and a gate that reports real homology: the A\*02
target-specific design goes from 0 surviving candidates to 36, and the wet-lab
pick at region 497-515 passes with its nearest off-target 2 mismatches away and
no hit closer than that.

**An empty owned set is never silent.** A run with the filter on and nothing
excluded prints a `NO SELF EXCLUDED` warning naming the consequence, because a
silent "nothing is self" is what made this invisible in the first place.

The reference is built once with `hush reference build`, which fetches the
RefSeq human mRNA set, runs `bowtie-build`, and writes a manifest recording the
RefSeq release, source URLs, FASTA sha256, bowtie version and build date, so a
build is auditable. Designs then find it automatically.

bowtie1 was chosen over BLAST because `-v` mode is an end-to-end,
mismatch-bounded search, which is exactly the question being asked; BLAST's
local alignment is not.

**Running unreferenced is possible but never silent.** HUSH supports it because
you may want to supply a different reference genome, but the CLI prints an
`UNREFERENCED RUN` warning naming the consequence and the design TSV carries an
`offtarget_checked` column, so an empty hit count can never be read as "no
off-targets found".

## Homopolymer runs

**On by default.** Rejects a core containing a single-nucleotide run longer than
`--max-run` (3). Poly-runs hurt siRNA stability and cause polymerase
read-through. The gate is base-agnostic and inspects the sense core only; the
U-specific check that matters for expressed hairpins happens after grafting, in
[stage 4](shrna-conversion.md).

## Opt-in gates

**GC window** (`--gc`, 30-52%) - a standalone counterpart to Reynolds' internal
GC criterion, which keeps its published form and stays pinned.

**Internal hairpin** (`--hairpin`, Tm <= 20 C) - rejects a sense strand that
folds into a stable self-structure. The standalone counterpart to Reynolds'
internal-repeat criterion. `core.hairpin` treats intramolecular melting as
concentration-independent (`Tm = dH/dS`) after Groebe & Uhlenbeck 1988, *NAR*
16:11725, whose loop destabilisation values it interpolates.

**Duplex-end asymmetry** (`--asymmetry`) - requires the duplex's looser end to
favour loading the guide. RISC keeps the strand whose 5' end sits at the less
stable end (Khvorova et al. 2003, *Cell* 115:209; Schwarz et al. 2003, *Cell*
115:199). The filter uses the **sign** of the differential only - there is
deliberately no magnitude threshold, because turning kcal/mol into a cut-off
needs calibration data HUSH does not have.

Thermodynamics throughout use the Xia et al. 1998 (*Biochemistry* 37:14719)
nearest-neighbor parameters by default, with Freier et al. 1986 (*PNAS* 83:9373)
available via `--params`; the Tm formulation follows SantaLucia & Hicks 2004.

## Allele targeting

Supplying `--panel` engages stage 3. Each candidate-length window is evaluated
for whether it **hits** the alleles you want silenced and **discriminates**
against those you want spared, and that verdict is ANDed with every gate above.

Four intents: `--pan-allele` (conserved sites across the family, the default),
`--allele-specific ID` (single one member out), `--target-specific` (spare the
whole panel, for when the target is loaded apart from its reference panel), and
explicit `--hit`/`--spare` sets.

**Coverage is required.** A window where a targeted or spared allele is uncovered
cannot satisfy the request: missing data is never read as evidence.

**A window is scored by its weakest spared allele.** An siRNA is only as specific
as the allele it distinguishes least well, so `minDiff` and `disc` describe that
allele; `--explain` shows the rest.

### Discrimination ranking

When a request spares something, survivors are ordered by how well they
discriminate. The key is a **lexicographic tuple of booleans, not a weighted
score** - mismatch count is deliberately excluded from the per-mismatch key. It
rests on exactly three papers, applied in this order:

1. **Ohnishi et al. 2008** (*PLoS ONE* 3:e2248) - a second mismatch, even inside
   the seed, markedly improves discrimination. This is why mismatch count sits
   at the top of the ordering.
2. **Huang et al. 2009** (*NAR* 37:7560) - Table 1, a lookup keyed by the exact
   base change at the exact position. The most specific evidence available, and
   the first tier of the per-mismatch tuple (`H+` / `H-` / `H?`).
3. **Schwarz et al. 2006** (*PLoS Genetics* 2:e140) - the general principles that
   order what Huang endorses: purine:purine mismatches discriminate best, and
   mismatches 3' of the seed discriminate robustly while seed mismatches are
   overwhelmed at high siRNA concentration. Also the source of the consistently
   discriminating positions, including the P5-in-the-seed anomaly, which HUSH
   reproduces intact rather than smoothing over.

The `disc` readout encodes this as e.g. `P16 post H+ S pp`: guide position,
region (`post` / `seed` / `anch`), Huang status, position tier (`S` strongest,
`C` consistent, `-` neither), and mismatch class (`pp` purine:purine, `pu`/`up`/
`uu`, `?` ambiguous). `--explain` adds each allele's difference count and flags
`agree=no` where Huang and the regional finding disagree.

## Validation

HUSH is benchmarked against **siDirect** (Naito et al. 2004, 2009), which it
positions itself as a successor to: on matched settings the rule pre-screen
reproduces its picks candidate-for-candidate. `tools/compare_sidirect.py` diffs
a design TSV against a siDirect export.

When comparing in `natural` overhang mode, compare on the **guide** or the target
position, not the passenger: the guide equals siDirect's listed guide, but the
passenger is offset because its overhang comes from the downstream flank.

## Not implemented

A **sequence-motif constraint** on `design` (e.g. forcing a base at a core
position) remains unimplemented. Its original motivation was Pol III initiation
for shRNA, which stage 4 now addresses as the `starts_with_purine` QC flag
rather than a hard gate - see [shrna-conversion.md](shrna-conversion.md).

A **machine-learning efficacy layer** is out of scope. The prior-art position is
that modern ML tools confirm efficacy rather than generate candidates; HUSH
generates.
