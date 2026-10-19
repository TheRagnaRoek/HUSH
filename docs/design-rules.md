# Design rules

What each gate tests, why it exists, and what it rests on. Implementation lives
in `core/` and `msa/`; the flags are in [cli.md](cli.md).

Selection is a **filter, not a score**. Every criterion below is an independent,
optional gate, and a candidate is selected only if it clears every enabled one.
Nothing is combined into a single number. A composite score hides which property
disqualified a candidate, and the three rule sets below were each published with
their own thresholds, which a weighted blend would invalidate.

## The duplex

HUSH's canonical duplex has a **19-nt paired core and symmetric 2-nt 3'
overhangs**, so the synthesised strands are 21 nt (Elbashir et al. 2001, *Nature*
411:494; *EMBO J* 20:6877).

`--length` sets the length of the **core**, not of the strand. All three rule sets
are defined on the core, so scoring follows the source papers regardless of the
overhang policy. Overhangs are stored separately from the core, and no gate
scores them.

With the default `natural` policy, each overhang is taken from the flanks of the
target. The passenger overhang is the 2 nt 3' of the core. The guide overhang is
the reverse complement of the 2 nt 5' of the core, which gives the offset duplex
of Elbashir et al. If a flank extends past the end of the transcript, `dtdt` is
used for that overhang. `dtdt` is 2'-deoxythymidine, which is tolerated at the 3'
end and is the common synthetic default; `uu` is ribo-uridine.

The two overhangs come from opposite flanks, so **the guide overhang is not the
reverse complement of the passenger overhang**. Code that reconstructs a duplex
must read both overhangs.

## The three rule sets

All three are from 2004 and are reimplemented in `core.rules` as published. By
default a candidate passes if any enabled rule set accepts it (`--mode any`).

Positions are 1-based from the 5' end of the named strand. The sense strand is
the passenger (the mRNA target site), the antisense strand is the guide, both
5'->3'. Each rule runs on the strand it is defined for.

**Ui-Tei** (Ui-Tei et al. 2004, *NAR* 32:936) is a strict four-rule pass/fail
filter. All four rules must hold:

1. A or U at the guide 5' end.
2. G or C at the sense 5' end.
3. At least 4 A/U among the first 7 nt of the guide.
4. No G/C run longer than 9 nt.

Rules 1-3 express one principle: RISC keeps the strand whose 5' end lies at the
less stable end of the duplex (Khvorova et al. 2003; Schwarz et al. 2003). The
constants are fixed for a 19-mer and are not rescaled for longer cores.
Ambiguity codes count as neither A/U nor G/C. Two of the rules pull in opposite
directions on seed stability: rule 2 keeps the passenger seed GC-rich, while
rule 3 lowers the guide seed Tm. HUSH documents this and does not resolve it.

**Reynolds** (Reynolds et al. 2004, *Nat Biotechnol* 22:326) is a score out of 10,
implemented as in the original algorithm and not as the common simplification of
eight equal binary criteria. On the sense strand:

| Criterion | Points |
|---|---|
| GC content 30-52% | +1 |
| No internal repeat | +1 |
| A at 3 | +1 |
| U at 10 | +1 |
| A at 19 | +1 |
| A/U base in positions 15-19 | +1 per base (0 to +5) |
| G at 13 | -1 |
| G or C at 19 | -1 |

The range is -2 to +10 and the default cut-off is 6. Position 19 contributes
three terms: A scores +2, U scores +1, and G or C scores -1. "No internal repeat"
means that the most stable hairpin of the sense strand melts at or below 20 C
(see [Opt-in gates](#opt-in-gates)); the verdict of an external folding tool can
replace it. Reynolds published no per-siRNA score table, so there is no
reproduction fixture and the tests check the arithmetic. On cores longer than
19 nt, the 3' criteria (positions 15-19 and 19) follow the actual 3' end, and
positions 3, 10 and 13 remain absolute.

**Amarzguioui** (Amarzguioui & Prydz 2004, *BBRC* 316:1050) is an additive score
on the sense strand: `score = motif points + DT3`.

- Motifs: +1 each for G/C at 1, A at 6 and A/U at 19; -1 each for U at 1 and G
  at 19.
- DT3 is the A/U count in the last 3 nt minus the A/U count in the first 3 nt
  (3' end minus 5' end, range -3 to +3). It enters the score as its own value,
  not through a lookup.

The sense 3' end is the same physical end as the guide 5' end, so a positive DT3
favours loading of the guide, and functional siRNAs score positive. The sign is
easy to invert in an implementation. The default cut-off is 3, inclusive; the
authors describe it as conservative.

Table 3 of the paper (46 siRNAs) is the exact-reproduction fixture.
`tests/test_rules.py` reproduces every DT3, motif set and total, the reported
class means (2.6 functional, 0 non-functional) and the selection result (14
selected, 11 functional) from `tests/data/amarzguioui_table3.tsv`.

Known limitation: on cores longer than 19 nt, the A/U at 19 and G at 19 motifs
read absolute position 19, while DT3 follows the actual 3' end, so the two use
different reference positions. A 19-mer is not affected.

## Seed-duplex off-target Tm

**On by default.** Implemented in `core.seed` after Naito et al. 2009 (*BMC
Bioinformatics* 10:392).

Seed-mediated off-target silencing correlates with the thermodynamic stability of
the seed duplex (Ui-Tei et al. 2008, *NAR* 36:7100): a weakly pairing seed causes
few off-target effects even where matching sites exist. HUSH computes the
seed-duplex melting temperature for both strands and rejects a candidate unless
both seeds are below the threshold (21.5 C by default).

**Seed and formula.** The seed is positions 2-8 of a strand (1-based from its 5'
end). Position 1 lies in the Argonaute MID pocket and does not pair with the
target. Both strands are filtered, because the passenger is loaded often enough
to cause its own seed-dependent off-target silencing. The melting temperature
follows Naito et al. 2009:

```
Tm = 1000 * dH / (S0 + dS + R ln(Ct/4)) - 273.15 + 16.6 log10[Na+]
```

dH and dS are the Freier 1986 stacking sums of the seed, without initiation,
terminal A-U or symmetry terms. S0 = -10.8 cal/(mol K) is the helix initiation
entropy, R = 1.987 cal/(mol K), Ct = 100 uM and [Na+] = 100 mM. The factor 4
always applies. A seed passes when its Tm is strictly below the threshold, and a
candidate passes when both seeds pass.

The equation differs from the general [model](#nearest-neighbor-thermodynamics)
in three ways: it uses the Freier parameter set, it has no duplex terms, and it
uses the rounded gas constant of the paper (`core.thermo` uses 1.9872). The
21.5 C cut-off is empirical and is only valid with these numbers. With Xia
stacking values and the same cut-off, 13% of the 16,384 possible heptamers change
verdict; if the cut-off is re-matched to the same pass rate, 2.3% change.

**Validation.** The Tm of all 16,384 heptamers reproduces the published figures:
4,488 (27.4%) are below 21.5 C, and the values range from -12.0 to 59.6 C (the
paper reports -12 to 60 C). Any wrong constant changes the count
(`tests/test_seed.py`).

**Strictness.** The filter is strict. On a very GC-rich CDS such as HLA it can
reject every candidate. This is expected behaviour. Relax it with `--no-seed`.

**Physical basis.** The seed duplex does not form in solution: Argonaute changes
the binding properties of its guide, and seed pairing is 50-200 times faster than
annealing of the corresponding free RNAs (Salomon et al. 2015, *Cell* 162:84; Wee
et al. 2012, *Cell* 151:1055). HUSH keeps the filter because it is the standard
published filter (Naito et al. 2009) and HUSH reproduces its figures. HUSH treats
it as a calibrated empirical proxy, not as a mechanism.

**Fixed parameters.** The internals of the seed filter cannot be changed: the
Freier parameter set, the initiation constant, the concentration and the salt
concentration are fixed, and `--params` does not affect them. Only the threshold
can be changed. The 21.5 C cut-off is only valid with exactly these values.

## Transcriptome off-target

**On by default.** Implemented in `core.offtarget` with `bowtie1` in `-v` mode
over a prebuilt index.

This is the only gate that searches a reference transcriptome instead of using a
sequence-based proxy. It checks whether the guide has a near match elsewhere in
the human transcriptome. A candidate is rejected when its nearest off-target has
fewer than `--offtarget-min-mismatch` mismatches (default 2, which matches
siDirect on the PDCD1 benchmark).

**Both strands are searched.** Either strand of the duplex can be loaded into
RISC, so each candidate contributes two reads. Both are the 19-nt core only; the
overhangs are not searched.

- The **guide side** searches the sense core. The guide silences transcripts
  complementary to it, which are the transcripts that carry the sense sequence.
- The **passenger side** searches the antisense core, for the transcripts that the
  passenger would silence if it were loaded instead of the guide.

Only the forward strand of each transcript is searched (`--norc`), because
transcripts are single-stranded sense mRNA. The design's own transcripts (see
below) are removed from the hits of both sides. A candidate fails when the
nearest remaining hit on either strand has fewer mismatches than the threshold.

The gate uses the mismatch count only, not the positions of the mismatches,
because the seed filter already accounts for the seed region. The search depth
(default 3, the `-v` maximum) exceeds what the gate needs (threshold minus 1), so
hit counts are reported down to 3 mismatches. The TSV and API report the nearer of
the two strands' nearest hits, and the hit counts summed over both strands.

### What the design owns

The gate must exclude the design's own transcripts, and for a gene family one
accession is not enough. A guide is meant to match its target. On HLA it is also
meant to match the family, and the reference contains the locus more than once:
HLA-A has five RefSeq transcripts (`NM_002116`, `NM_001242758`, `NR_001434`,
`XM_041680767`, `XM_041680768`) before alternate-haplotype assembly copies are
counted. If any of them counts as an off-target, the gate measures homology to the
target instead of off-target risk. With nothing excluded, every candidate on
HLA-A\*01, A\*02 and A\*03 is rejected. Naming a single accession by hand
increases the surviving candidates across seven production scenarios from 3 to 5
out of 107.

The owned transcripts (`core.offtarget.SelfTranscripts`) are the union of three
sources:

* `--offtarget-self-accession`: a comma-separated list of accessions, matched
  without the version (`NM_002116` matches `NM_002116.8`).
* The RefSeq accession in the FASTA header of the target, if it has one. IMGT
  allele headers such as `HLA-A02:01:01:01` have none, so this source alone left
  HLA runs with nothing excluded.
* Homology resolution (on by default, `--offtarget-resolve-self`): the target,
  every panel member and every `--offtarget-self-fasta` entry is split into
  k-mers, and a reference transcript is owned when at least
  `--offtarget-self-coverage` percent of the distinct k-mers of one input match it
  exactly. Coverage is computed per input, so a short panel member is not diluted
  by a long target.

Panel members are used for resolution because they are intended targets: a
pan-allele design is meant to hit them. This does not weaken an allele-specific
design, because discrimination between alleles is done by the
[allele stage](#allele-targeting), at a resolution that the transcriptome lacks
(every HLA-A allele maps to the same RefSeq transcripts).

The default coverage of 10% comes from the measured distribution on the four
bundled HLA-A alleles. Exact 19-mer coverage of the reference by one allele is
clustered: 5 transcripts above 25% (HLA-A itself, the same five for all four
alleles); roughly 100-140 between 10% and 25% (alternate-haplotype copies of the
locus and the closest class I paralogs); and a tail of transcripts that share one
or two k-mers, which are real off-targets and stay in the gate. A threshold above
25% restricts the set to the target gene; below 5% it includes incidental
matches. Sequence is the only available signal, because a bowtie index stores
accessions and not gene symbols. An alternate-haplotype copy of the target gene
and a true paralog are therefore distinguished only by how much of the target
they contain.

On the seven original wet-lab scenarios, this changes the gate from rejecting
everything to reporting real homology: the A\*02 target-specific design goes from
0 surviving candidates to 36, and the wet-lab pick at region 497-515 passes, with
its nearest off-target 2 mismatches away and no closer hit.

**Empty owned set.** If the filter is on and nothing is excluded, the run prints a
`NO SELF EXCLUDED` warning that states the consequence. Without the warning this
condition is not visible in the output.

The reference is built once with `hush reference build`. The command downloads
the RefSeq human mRNA set, runs `bowtie-build`, and writes a manifest with the
RefSeq release, source URLs, FASTA sha256, bowtie version and build date, so the
build can be audited. Designs then find the index automatically.

bowtie1 was chosen over BLAST because its `-v` mode performs an end-to-end search
bounded by a mismatch count, which matches the question the gate asks. BLAST uses
local alignment.

**Unreferenced runs.** A run without the built reference is possible, for example
to use a different reference genome. The CLI then prints an `UNREFERENCED RUN`
warning, and the design TSV has an `offtarget_checked` column, so an empty hit
count cannot be read as "no off-targets found".

## Homopolymer runs

**On by default.** A core with a single-nucleotide run longer than `--max-run`
(default 3) is rejected. Homopolymer runs reduce siRNA stability and cause
polymerase read-through. The gate does not depend on the base and checks the sense
core only. The check for U runs in expressed hairpins is done after grafting, in
[stage 4](shrna-conversion.md).

## Opt-in gates

**GC window** (`--gc`, 30-52%) is a tunable gate alongside the GC criterion of
Reynolds, which keeps its published form and cannot be changed.

**Internal hairpin** (`--hairpin`, Tm <= 20 C) rejects a sense strand that folds
into a stable self-structure. It is the tunable counterpart to the internal-repeat
criterion of Reynolds. The model is described under
[Nearest-neighbor thermodynamics](#nearest-neighbor-thermodynamics).

**Duplex-end asymmetry** (`--asymmetry`) requires the less stable end of the
duplex to favour loading of the guide. RISC keeps the strand whose 5' end lies at
the less stable end (Khvorova et al. 2003, *Cell* 115:209; Schwarz et al. 2003,
*Cell* 115:199). The filter uses only the **sign** of the differential. It has no
magnitude threshold, because converting kcal/mol into a cut-off requires
calibration data that HUSH does not have.

## Nearest-neighbor thermodynamics

`core.thermo` converts an RNA sequence into (dH, dS); free energy and Tm follow. A
duplex is the sum of its stacking steps (each overlapping dinucleotide and its
partner on the opposite strand) plus duplex terms. Only complementary nucleotide
pairs are modelled, in the 1 M NaCl reference state of the published parameters.
There is no salt correction, and mismatches, G-U pairs and dangling ends are not
included.

**Parameter sets.** Xia et al. 1998 (*Biochemistry* 37:14719) is the default and
Freier et al. 1986 (*PNAS* 83:9373) is available via `--params`, which applies to
the asymmetry and hairpin read-outs only. A published threshold is only valid
against the parameter set it was calibrated on (the cut-off of the seed filter
uses Freier), so the caller always names the set explicitly and it is never
replaced implicitly. Only 10 of the 16 steps are independent: read from the other
strand, step XY becomes its reverse complement with the same energy. The code
stores the 10 published values and derives the other 6. Units are dH in kcal/mol
and dS in cal/(mol K).

**Two entry points.**

- A complete duplex adds the duplex terms of Xia to the stacking sum: initiation
  (dH +3.61, dS -1.5), +3.72 and +10.5 per helix end closed by an A-U pair, and
  dS -1.4 for a self-complementary strand. These terms do not transfer to
  another parameter set, so this entry point uses Xia only.
- Stacking only has no duplex terms. It is used for windows inside a longer helix,
  whose edges are not helix ends: the asymmetry windows, hairpin stems and the
  seed Tm.

Using the wrong entry point shifts every value by the duplex terms.

**Formulas.** dG = dH - T*dS, reported at 37 C as dG37. Tm = dH / (dS + R
ln(Ct/x)) with R = 1.9872 cal/(mol K), Ct the total strand concentration, x = 4
for two distinct strands and x = 1 for a self-complementary strand (SantaLucia &
Hicks 2004). The default Ct is 1 uM, which is a convention. A published score uses
its own conditions: the seed filter uses its own equation at 100 uM, and hairpin
melting has no concentration term.

**Hairpins.** The internal-hairpin gate and the "no internal repeat" criterion of
Reynolds use `core.hairpin`. It is an approximation, not a reproduction of the
folding method of the Reynolds paper, which the paper does not specify. The
internal structure of a strand is its most stable simple hairpin: a perfectly
paired stem of at least 3 bp closing one loop of at least 3 nt, anywhere in the
strand (the minimums are conventions, not values from the sources). The cut-off
(20 C, inclusive) is compared with the Tm of the fold with the lowest dG37, not
the highest Tm over all folds. A strand that folds into nothing passes.

- **Stem:** stacking only, Xia by default (`--params`), without initiation or
  terminal A-U terms. G-U pairs, mismatches and bulges are not modelled, so
  structure may be slightly under-estimated.
- **Loop:** the destabilisation dG37 of Groebe & Uhlenbeck 1988 (*NAR* 16:11725,
  Table 3), given for homopolymer loops of size 3, 4, 5, 7 and 9. It is averaged
  over the bases of the loop and interpolated linearly in size, clamped at 3 and
  9. G and ambiguity codes take the mean of A, C and U. The loop penalty replaces
  the initiation term, as in the derivation of the table.
- **Tm:** the table gives dG37 only, so the loop is treated as purely entropic
  (dH_loop = 0, dS_loop = -dG37_loop / T37). Tm = dH / dS, without a concentration
  term, because intramolecular melting does not depend on concentration
  (Groebe & Uhlenbeck, 1 to 10 uM).

**Validation.** Tm agrees with the independent implementation of the same Xia
table in Biopython (`RNA_NN2`) within 0.01 C. Biopython rounds the gas constant
to 1.987.

## Allele targeting

Supplying `--panel` enables stage 3. Each window of candidate length is checked
for whether it **hits** the alleles to be silenced and **discriminates** against
the alleles to be spared. The result is combined (AND) with every gate above.

Four intents are available: `--pan-allele` (conserved sites across the family,
the default), `--allele-specific ID` (hit one member and spare the rest),
`--target-specific` (spare the whole panel, for a target that is loaded
separately from its reference panel), and explicit `--hit` / `--spare` sets.

**Coverage is required.** A window in which a targeted or spared allele is not
covered cannot satisfy the request: missing data is never treated as evidence.

**A window is scored by its weakest spared allele.** An siRNA is only as specific
as the allele it distinguishes least well, so `minDiff` and `disc` describe that
allele; `--explain` shows the other alleles.

The allele gate counts differences (at least `--min-spare-mismatches` per spared
allele). It does not weigh them by type or position; placement enters only
through the [ranking](#discrimination-ranking).

### Alignment

Each panel member is aligned pairwise to the target (a star alignment with the
target as the fixed centre). The alignment is then projected onto the coordinates
of the target: a row as long as the target, where position i holds the residue of
the member at target position i. A deletion in the member becomes a gap in the
row (uncovered). Residues of the member opposite a target gap have no target
position, so they are kept as insertions attached to the next position. Window
and candidate positions therefore coincide, and no conversion between alignment
columns and target positions is needed.

The centre is fixed because candidates are enumerated from the target. The
classical center-star method (Gusfield) chooses the centre that minimises the
summed pairwise cost; that does not apply here, so no optimality guarantee holds.
The star alignment provides a common coordinate frame and does not guarantee
alignment quality.

The scores are project choices, not published constants: match +2, mismatch -1,
gap open -10, gap extension -0.5, terminal gaps free. An indel costs about as much
as four mismatches. A 1-nt deletion is placed as an internal gap only when roughly
four or more positions follow it (measured on random sequence). Nearer to a
sequence end, the aligner shifts the alignment instead and leaves that end
uncovered. Free terminal gaps let a member that covers only part of the target
extend past the ends without penalty, and its uncovered positions appear as gaps.
Members are converted to the molecule type of the target first, so T and U never
mismatch. Among equal-scoring alignments the first is taken: results are
reproducible, but the placement of an indel within a repeat is arbitrary.

**Limits.** The method suits near-identical members. The bundled HLA-A alleles
are 95-98% identical to A*01:01 (coding sequence, no indels). Accuracy decreases
as members diverge, because a single reference sequence cannot represent the
relationships between panel members. Gaps ignore the reading frame, so this
aligner can open frame-shifting gaps that a codon-aware aligner would avoid. The
high gap cost makes this rare among closely related members. Divergent families
would need a progressive aligner behind the same `AlignmentBackend` interface.

### Discrimination ranking

When a request spares an allele, the windows that pass the allele gate are ordered
by how reliably their difference from the spared allele is expected to
discriminate. Every such window already differs from each spared allele; the
windows differ in the position and type of that difference. The ranking never
removes a window. A high rank does not guarantee that the spared allele escapes
silencing: Du et al. 2005 (*NAR* 33:1671) found that 73% of singly mismatched
targets were still significantly silenced. A difference is necessary for
discrimination but not sufficient.

The key is a **lexicographic tuple of five booleans, not a weighted score**. The
criteria are checked in order. The first criterion that separates two windows
decides, and a deficit in an earlier criterion cannot be offset by a later one.
Nothing is weighted, so nothing needs calibration. Positions are guide positions,
1-based from the 5' end.

| # | Criterion | True when | Source |
|---|---|---|---|
| 1 | Listed for this base change | Table 1 lists the substitution at that position | Huang et al. 2009 (*NAR* 37:7560) |
| 2 | Inside the low-tolerance window | P9-P16 | Du 2005, Schwarz, Huang |
| 3 | At a replicated position | P10, P13 or P16 | Schwarz et al. 2006 (*PLoS Genet* 2:e140) |
| 4 | 3' of the seed | position 9 or later | Schwarz |
| 5 | Purine:purine mismatch | guide and spared base both purines | Schwarz, Du |

**Why this order.**

1. Criterion 1 is the only one that depends on the actual base substitution and
   not only on the position.
2. P9-P15 is the low-tolerance region of Du (target-site positions 5-11; for a
   19-nt site, guide = 20 - target, checked against the printed sequences of the
   paper). P16 is added although Du place it in medium tolerance: it is the
   best-replicated position of Schwarz, and Huang list it in 7 of 12 rows. 54 of
   the 55 position entries of Huang lie in P9-P16.
3. Schwarz report positions in five figures. P13 and P16 appear in all five and
   P10 in four (table below). P16 also passed their generality test: ten
   independent siRNA-mRNA pairs against two targets, a higher IC50 for the
   mismatched combination in every pair and at least 20-fold in seven.
4. A seed mismatch blocks binding and is overcome by a high siRNA concentration;
   a mismatch 3' of the seed blocks catalysis and is not. Criterion 4 ranks below
   the window because Du place P16-P19 in medium and high tolerance.
5. Criterion 5 is last and is a judgement: the design advice of both papers is
   positional, and the comparison across criteria (purine:purine at P12 against a
   milder class at P16) has not been measured. It is the criterion most open to
   reordering; changing its position requires changing one line in `_single_key`.

| Schwarz figure | System | Positions |
|---|---|---|
| 1B | Drosophila lysate, G:G mismatch | P5 P9 P10 P13 P14 P15 P16 |
| 3 | HEK293, discriminating in both systems | P5 P10 P12 P13 P16 |
| 3C | G:U wobble, 3' of the seed | P8 P11 P13 P14 P16 |
| 6A | purine:purine scan, P10 scaffold | P4 P7 P9 P10 P11 P13 P16 |
| 6B | purine:purine scan, P4 scaffold | P3 P4 P5 P9 P10 P11 P12 P13 P16 |

A window is ranked by its best-placed difference from its weakest spared allele.
Ties are broken by target position.

**The number of differences is not ranked.** All sources measure single
mismatches, and none states whether several poorly placed differences outweigh one
well-placed difference. A window with one difference listed by Huang therefore
ranks above a window with three differences that are only post-seed. This does not
mean that the first window is the better design. The count and the remaining
positions are reported next to the rank (`minDiff`, the `+P11,P4` cell,
`--explain`, TSV `disc_all_positions`, API `all_guide_positions`) so that the user
can weigh them, and two or more differences are flagged for reference only.
Ohnishi et al. 2008 (*PLoS ONE* 3:e2248) show that a second mismatch engineered
into the siRNA improves discrimination. This is a design intervention and not a
natural difference between alleles, so it is not used here.

**Limits of the ranking.**

- The order of the five criteria is only partly supported by evidence. Criterion 1
  first is well founded (specific information before general information). The
  window before single positions, and position before class, are judgements.
- Table 1 of Huang and the regional finding disagree in places: A>C lists P8,
  which lies inside the seed. They are compared and not reconciled. `agree` is
  yes when Huang list a post-seed position, or do not list a seed or anchor
  position; otherwise it is no. Most `agree=no` cases are therefore post-seed
  positions that Table 1 does not list.
- Schwarz found strong discrimination at P3-P5 with one siRNA whose seed is the
  most stable of their set and which loaded its sense strand. These seed positions
  are not used.
- Not modelled: structure within purine:purine (G:G and A:G were more selective
  than A:A and G:A at P10, and only A:G held up at 20 nM), and the finding of Du
  that A:C, C:A and the G:U wobble are best tolerated and therefore the weakest
  discriminators. They rank only as "not purine:purine".
- Guide position 1 does not pair with the target (region `anch`). It ranks with
  the seed and, like any difference, counts toward the allele gate: a spared allele
  that differs only at P1 passes the gate with the lowest rank.
- Ambiguity codes have no Huang row and no purine or pyrimidine identity. They give
  `H?` and `?`.
- The twelve rows of Huang are asserted verbatim in `tests/test_discrimination.py`,
  but they were transcribed from a secondary source and have not been checked
  against the paper (as of 2026-10-08). Schwarz and Du were read in full.

The `disc` readout encodes one mismatch, for example `P16 post H+ S pp`: guide
position, region (`post` / `seed` / `anch`), Huang status (`H+` listed, `H-` not
listed, `H?` ambiguous), position tier (`S` = P10, P13 or P16; `W` = inside
P9-P16; `-` outside), and class (`pp` purine:purine, `pu` / `up` / `uu`, `?`
ambiguous). `--explain` adds the difference count of each allele and marks
`agree=no` where Huang and the regional finding disagree.

## Validation

HUSH was compared against siDirect (Naito et al. 2004, 2009). On matched settings
the rule pre-screen reproduces the siDirect picks candidate for candidate. The
siDirect export, the PDCD1 CDS and the matching HUSH design are in `data/testing/`.

In `natural` overhang mode, compare on the guide or on the target position, not on
the passenger. The guide equals the guide listed by siDirect; the passenger is
offset because its overhang comes from the downstream flank.

## Not implemented

A **sequence-motif constraint** on `design` (for example, requiring a base at a
core position) is not implemented. Its original motivation was Pol III initiation
for shRNA, which stage 4 now handles as the `starts_with_purine` QC flag instead
of a hard gate (see [shrna-conversion.md](shrna-conversion.md)).

An **insertion-aware allele gate** is not implemented. Residues that a member
carries between two target positions are recorded and attached to the next
position, but they do not enter either gate. A targeted allele with an insertion
inside the window still counts as identical, and a spared allele that differs only
by an insertion counts as not discriminating, although in both cases the guide
would form a bulge. An insertion counts as inside a window only when it lies
strictly between two positions of the window. Insertions are not shown in the CLI,
TSV or API output. Scoring a bulge would require mismatch and bulge thermodynamics
that `core.thermo` does not contain.

**Ambiguity codes** in a panel member are compared as plain characters: N against A
is a difference, and N against N is a match. A panel from IPD-IMGT/HLA contains no
ambiguity codes; a hand-assembled panel may.

A **machine-learning efficacy layer** is out of scope. Current machine-learning
tools score the efficacy of existing candidates and do not generate candidates,
whereas HUSH generates candidates.
