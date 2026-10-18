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
runs past a transcript end falls back to `dtdt` there. `dtdt` is 2'-deoxythymidine,
tolerated at the 3' end and the common synthetic default; `uu` is ribo-uridine.

Because the two overhangs come from opposite flanks, **the guide's overhang is
not the reverse-complement of the passenger's**. Anything reconstructing a duplex
must read both, not derive one from the other.

## The three rule sets

All three are from 2004 and are reimplemented faithfully in `core.rules`. They
are combined by **union** (`--mode any`) by default, each weighted equally.

Positions are 1-based from the 5' end of the named strand. The sense strand is
the passenger (the mRNA target site), the antisense strand is the guide, both
5'->3'. Each rule runs on the strand it is defined for.

**Ui-Tei** (Ui-Tei et al. 2004, *NAR* 32:936) - a strict four-rule pass/fail
filter. All four must hold:

1. A or U at the guide 5' end.
2. G or C at the sense 5' end.
3. At least 4 A/U among the first 7 nt of the guide.
4. No G/C run longer than 9 nt.

Rules 1-3 are three statements of one principle: RISC keeps the strand whose 5'
end sits at the looser end of the duplex (Khvorova et al. 2003; Schwarz et al.
2003). The constants are fixed for a 19-mer and are not rescaled for longer
cores. Ambiguity codes count as neither A/U nor G/C. The rule set is also the
source of a tension HUSH documents rather than resolves: rule 2 keeps the
passenger seed GC-rich while rule 3 pulls the guide seed Tm down.

**Reynolds** (Reynolds et al. 2004, *Nat Biotechnol* 22:326) - a score out of 10,
implemented following the original algorithm rather than the widespread "eight
equal binary criteria" simplification. On the sense strand:

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

The range is -2 to +10 and the default cut-off is 6. Position 19 carries three
terms at once: A scores +2, U +1, G or C -1. "No internal repeat" means the most
stable hairpin of the sense strand melts at 20 C or below (see
[Opt-in gates](#opt-in-gates)); a verdict from an external folder can replace
it. Reynolds published no per-siRNA score table, so there is no reproduction
fixture and the tests check the arithmetic. On cores longer than 19 nt the 3'
parts (positions 15-19 and 19) follow the real 3' end, while positions 3, 10 and
13 stay absolute.

**Amarzguioui** (Amarzguioui & Prydz 2004, *BBRC* 316:1050) - an additive score,
sense strand: `score = motif points + DT3`.

- Motifs: +1 each for G/C at 1, A at 6 and A/U at 19; -1 each for U at 1 and G
  at 19.
- DT3 is the A/U count in the last 3 nt minus the count in the first 3 nt (3'
  end minus 5' end, range -3 to +3). It enters as its own value, not through a
  lookup.

The sense 3' end is the same physical end as the guide 5' end, so a positive DT3
favours the guide and functional siRNAs score positive. The sign is the easiest
thing to invert. The default cut-off is 3, inclusive; the authors call it
conservative.

Its Table 3 (46 siRNAs) is HUSH's exact-reproduction validation fixture:
`tests/test_rules.py` reproduces every DT3, motif set and total, the reported
class means (2.6 functional, 0 non-functional) and the selection result (14
selected, 11 functional) from `tests/data/amarzguioui_table3.tsv`.

Open item: on cores longer than 19 nt the A/U at 19 and G at 19 motifs read
absolute position 19, while DT3 follows the real 3' end, so the two disagree. A
19-mer is unaffected.

## Seed-duplex off-target Tm

**On by default.** `core.seed`, after Naito et al. 2009 (*BMC Bioinformatics*
10:392).

Seed-mediated off-target silencing tracks the thermodynamic stability of the
seed duplex (Ui-Tei et al. 2008, *NAR* 36:7100): a weakly-pairing seed is largely
off-target free even where its matches exist. HUSH computes the seed-duplex
melting temperature for both strands and rejects a candidate unless both seeds
are below the threshold (21.5 C by default).

**Seed and formula.** The seed is positions 2-8 of a strand (1-based from its 5'
end). Position 1 sits in the Argonaute MID pocket and does not pair with the
target. Both strands are filtered, because the passenger is loaded often enough
to drive its own seed-dependent off-target silencing. The melting temperature is
the equation of Naito et al. 2009:

```
Tm = 1000 * dH / (S0 + dS + R ln(Ct/4)) - 273.15 + 16.6 log10[Na+]
```

dH and dS are the Freier 1986 stacking sums of the seed, with no initiation,
terminal A-U or symmetry terms. S0 = -10.8 cal/(mol K) is the helix initiation
entropy, R = 1.987 cal/(mol K), Ct = 100 uM and [Na+] = 100 mM. The factor 4
applies unconditionally. A seed passes when its Tm is strictly below the
threshold, and a candidate passes when both seeds do.

This differs from the general [model](#nearest-neighbor-thermodynamics) in three
deliberate ways: the Freier set, no duplex terms, and the published rounded gas
constant (`core.thermo` uses 1.9872). The 21.5 C cut-off is empirical and only
meaningful with these numbers. With Xia stacking values and the same cut-off, 13%
of the 16,384 possible heptamers change verdict, and 2.3% if the cut-off is
re-matched to the same pass rate.

**Validation.** The Tm of all 16,384 heptamers reproduces the published figures:
4,488 (27.4%) fall below 21.5 C, with values from -12.0 to 59.6 C (the paper
reports -12 to 60 C). Any wrong constant moves the count
(`tests/test_seed.py`).

Two things about this filter are worth knowing.

It is **deliberately strict**, and on a very GC-rich CDS such as HLA it can
reject every candidate. That is the documented, expected behaviour, not a bug;
relax it with `--no-seed`.

Its **physics is contested**. The seed duplex does not form in solution: Ago
reshapes its guide's binding properties, and seed pairing proceeds 50-200x faster
than annealing of the equivalent naked RNAs (Salomon et al. 2015, *Cell* 162:84;
Wee et al. 2012, *Cell* 151:1055). HUSH keeps the filter because it is the
standard published filter (Naito et al. 2009) and HUSH reproduces its figures,
but treats it as a calibrated empirical proxy rather than a mechanism.

Consequently the seed filter's internals are **not tunable**: the Freier
parameter set, the initiation constant, the concentration and the salt are
pinned, and `--params` never reaches them. Only the threshold moves. The 21.5 C
cut-off is meaningful only against exactly those numbers.

## Transcriptome off-target

**On by default.** `core.offtarget`, via `bowtie1` in `-v` mode over a prebuilt
index.

This is the safety floor and the one gate that consults reality rather than a
sequence-intrinsic proxy: it asks whether the guide has a near-match anywhere
else in the human transcriptome. A candidate is rejected when its nearest
off-target carries fewer than `--offtarget-min-mismatch` mismatches (2 by
default, which matches siDirect on the PDCD1 benchmark).

**Both strands are searched.** Either strand of the duplex can be loaded into
RISC, so each candidate contributes two reads. Both are the 19-nt core only, and
the overhangs are not searched.

- The **guide side** searches the sense core. The guide silences transcripts
  complementary to it, which are the transcripts that carry the sense sequence.
- The **passenger side** searches the antisense core, for the transcripts the
  passenger would silence if it were loaded by mistake.

Only the forward strand of each transcript is searched (`--norc`), because
transcripts are single-stranded sense mRNA. The design's own transcripts (see
below) are removed from the hits of both sides. A candidate fails when the
nearest remaining hit on either strand has fewer mismatches than the threshold.

The gate is flat: only the mismatch count matters, not where the mismatches fall,
because the seed filter already weighs the seed region. The search depth
(default 3, the `-v` maximum) exceeds what the gate needs (threshold minus 1), so
hit counts are reported down to 3 mismatches. The TSV and API report the nearer
of the two strands' nearest hits, and hit counts summed over both strands.

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
internal-repeat criterion. The model is described under
[Nearest-neighbor thermodynamics](#nearest-neighbor-thermodynamics).

**Duplex-end asymmetry** (`--asymmetry`) - requires the duplex's looser end to
favour loading the guide. RISC keeps the strand whose 5' end sits at the less
stable end (Khvorova et al. 2003, *Cell* 115:209; Schwarz et al. 2003, *Cell*
115:199). The filter uses the **sign** of the differential only - there is
deliberately no magnitude threshold, because turning kcal/mol into a cut-off
needs calibration data HUSH does not have.

## Nearest-neighbor thermodynamics

`core.thermo` turns an RNA sequence into (dH, dS); free energy and Tm follow. A
duplex is the sum of its stacking steps (each overlapping dinucleotide and its
partner on the opposite strand) plus duplex terms. Only complementary nucleotide
pairs are modelled, in the 1 M NaCl reference state of the published parameters.
There is no salt correction, and mismatches, G-U pairs and dangling ends are not
included.

**Parameter sets.** Xia et al. 1998 (*Biochemistry* 37:14719) is the default and
Freier et al. 1986 (*PNAS* 83:9373) is available via `--params`, which applies to
the asymmetry and hairpin read-outs only. A published threshold is only valid
against the set it was calibrated on (the seed filter's cut-off uses Freier), so
the set is always named by the caller and never swapped silently. Only 10 of the
16 steps are independent: read from the other strand, step XY becomes its reverse
complement with the same energy. The code stores the 10 published values and
derives the other 6. Units are dH in kcal/mol and dS in cal/(mol K).

**Two entry points.**

- A complete duplex adds Xia's duplex terms to the stacking sum: initiation
  (dH +3.61, dS -1.5), +3.72 and +10.5 per helix end closed by an A-U pair, and
  dS -1.4 for a self-complementary strand. These terms do not transfer to
  another parameter set, so this entry point is Xia only.
- Stacking only has no duplex terms. It is for windows inside a longer helix,
  whose edges are not helix ends: the asymmetry windows, hairpin stems and the
  seed Tm.

Using the wrong one shifts every number by the duplex terms.

**Formulas.** dG = dH - T*dS, reported at 37 C as dG37. Tm = dH / (dS + R
ln(Ct/x)) with R = 1.9872 cal/(mol K), Ct the total strand concentration, x = 4
for two distinct strands and 1 for a self-complementary strand (SantaLucia &
Hicks 2004). The default Ct is 1 uM, a convention. A published score supplies its
own: the seed filter uses its own equation at 100 uM, and hairpin melting has no
concentration term.

**Hairpins.** The internal-hairpin gate and Reynolds' "no internal repeat"
criterion use `core.hairpin`. It is a transparent approximation, not a
reproduction of the folder that the Reynolds paper leaves unspecified. A strand's
internal structure is its most stable simple hairpin: a perfectly paired stem of
at least 3 bp closing one loop of at least 3 nt, anywhere in the strand (the
minimums are conventions, not from the sources). The cut-off (20 C, inclusive) is
compared with the Tm of the fold with the lowest dG37, not the highest Tm over
all folds. A strand that folds into nothing passes.

- **Stem:** stacking only, Xia by default (`--params`), with no initiation or
  terminal A-U terms. G-U pairs, mismatches and bulges are not modelled, so
  structure may be slightly under-estimated.
- **Loop:** the destabilisation dG37 of Groebe & Uhlenbeck 1988 (*NAR* 16:11725,
  Table 3), given for homopolymer loops of size 3, 4, 5, 7 and 9. It is averaged
  over the loop's bases and interpolated linearly in size, clamped at 3 and 9. G
  and ambiguity codes take the mean of A, C and U. The loop penalty takes the
  role of initiation, as in how the table was derived.
- **Tm:** the table gives dG37 only, so the loop is treated as purely entropic
  (dH_loop = 0, dS_loop = -dG37_loop / T37). Tm = dH / dS with no concentration
  term, because intramolecular melting does not depend on concentration
  (Groebe & Uhlenbeck, 1 to 10 uM).

**Validation.** Tm agrees with Biopython's independent implementation of the same
Xia table (`RNA_NN2`) within 0.01 C. Its gas constant is rounded to 1.987.

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

The allele gate counts differences (at least `--min-spare-mismatches` per spared
allele). It does not weigh them by type or position; placement enters only
through the [ranking](#discrimination-ranking).

### Alignment

Each panel member is aligned to the target on its own, a star with the target as
its fixed centre. It is then projected onto the target's coordinates: a row as
long as the target, where position i holds the member's residue at target
position i. A deletion in the member becomes a gap in the row (uncovered).
Residues the member has opposite a target gap have no target position, so they
are kept as insertions attached to the next position. Window and candidate
positions therefore coincide, with no conversion between alignment columns and
target positions.

The centre is fixed because candidates are enumerated from the target. The
classical center-star method (Gusfield) picks the centre that minimises the
summed pairwise cost; that does not apply here, so no optimality is claimed. The
star gives a consistent frame, not a quality guarantee.

The scores are project choices, not published constants: match +2, mismatch -1,
gap open -10, gap extension -0.5, terminal gaps free. An indel costs about as
much as four mismatches. A 1-nt deletion is placed as an internal gap only when
roughly four or more positions follow it (measured on random sequence). Nearer a
sequence end the aligner shifts instead and leaves that end uncovered. Free
terminal gaps let a partial member hang off the ends, and its uncovered
positions show as gaps. Members are converted to the target's molecule type
first, so T and U never mismatch. Among equal-scoring alignments the first is
taken: results are reproducible, but the placement of an indel within a repeat is
arbitrary.

**Limits.** The method suits near-identical members. The bundled HLA-A alleles
are 95-98% identical to A*01:01 (coding sequence, no indels). Accuracy drops as
members diverge, because one centre cannot represent the relationships within
the panel. Gaps ignore the reading frame, so this aligner can open
frame-shifting gaps that a codon-aware one would avoid. The high gap cost makes
that rare among closely related members. Divergent families would need a
progressive aligner behind the same `AlignmentBackend` interface.

### Discrimination ranking

When a request spares something, the windows that pass the allele gate are
ordered by how reliably their difference to the spared allele should work. Every
such window already differs from each spared allele; what differs is where. The
ranking never removes a window, and a high rank does not predict silence on the
spared allele: Du et al. 2005 (*NAR* 33:1671) found that 73% of singly
mismatched targets were still significantly silenced. A difference is necessary
for discrimination, not sufficient.

The key is a **lexicographic tuple of five booleans, not a weighted score**. The
criteria are consulted in order, the first that separates two windows decides,
and a deficit in an earlier one cannot be offset by a later one. Nothing is
weighted, so nothing needs calibrating. Positions are guide positions, 1-based
from the 5' end.

| # | Criterion | True when | Source |
|---|---|---|---|
| 1 | Listed for this base change | Table 1 lists the substitution at that position | Huang et al. 2009 (*NAR* 37:7560) |
| 2 | Inside the low-tolerance window | P9-P16 | Du 2005, Schwarz, Huang |
| 3 | At a replicated position | P10, P13 or P16 | Schwarz et al. 2006 (*PLoS Genet* 2:e140) |
| 4 | 3' of the seed | position 9 or later | Schwarz |
| 5 | Purine:purine mismatch | guide and spared base both purines | Schwarz, Du |

**Why this order.**

1. It is the only criterion conditioned on the actual substitution, not on
   position alone.
2. P9-P15 is Du's low-tolerance region (target-site positions 5-11; for a 19-nt
   site guide = 20 - target, checked against the paper's printed sequences). P16
   is added although Du place it in medium tolerance: it is Schwarz's
   best-replicated position and Huang list it in 7 of 12 rows. 54 of Huang's 55
   position entries lie in P9-P16.
3. Schwarz report positions in five figures; P13 and P16 appear in all five and
   P10 in four (table below). P16 also passed their generality test: ten
   independent siRNA-mRNA pairs against two targets, a higher IC50 for the
   mismatched combination in every pair and at least 20-fold in seven.
4. A seed mismatch blocks binding and is overwhelmed by high siRNA concentration;
   one 3' of the seed blocks catalysis and is not. It ranks below the window
   because Du place P16-P19 in medium and high tolerance.
5. Last, as a judgement: both papers' design advice is positional, and the
   comparison across axes (purine:purine at P12 against a milder class at P16)
   is not measured. It is the entry most open to reordering; moving it means
   moving one line in `_single_key`.

| Schwarz figure | System | Positions |
|---|---|---|
| 1B | Drosophila lysate, G:G mismatch | P5 P9 P10 P13 P14 P15 P16 |
| 3 | HEK293, discriminating in both systems | P5 P10 P12 P13 P16 |
| 3C | G:U wobble, 3' of the seed | P8 P11 P13 P14 P16 |
| 6A | purine:purine scan, P10 scaffold | P4 P7 P9 P10 P11 P13 P16 |
| 6B | purine:purine scan, P4 scaffold | P3 P4 P5 P9 P10 P11 P12 P13 P16 |

A window is ranked by its best-placed difference to its weakest spared allele;
ties break by target position.

**The number of differences is not ranked.** All sources measure single
mismatches, and none says whether several poorly placed differences outweigh one
well-placed one. A window with one Huang-listed difference therefore ranks above
one with three merely post-seed differences; that is not a claim that the first
is the better design. The count and the remaining positions are reported next to
the rank (`minDiff`, the `+P11,P4` cell, `--explain`, TSV `disc_all_positions`,
API `all_guide_positions`) for you to weigh, and two or more differences are
flagged for reference only. Ohnishi et al. 2008 (*PLoS ONE* 3:e2248) show that a
second mismatch engineered into the siRNA improves discrimination; that is a
design intervention, not a natural difference between alleles, so it is not used
here.

**Limits of the ranking.**

- The order of the five entries is only partly evidenced. Huang first is well
  founded (specific beats general). The window before single positions, and
  position before class, are judgements.
- Huang's Table 1 and the regional finding disagree in places: A>C lists P8,
  inside the seed. They are compared, not reconciled. `agree` is yes when Huang
  lists a post-seed position, or does not list a seed or anchor position;
  otherwise it is no. Most `agree=no` cases are therefore post-seed positions
  that Table 1 omits.
- Schwarz's strong discrimination at P3-P5 comes from one siRNA whose seed is the
  most stable of their set and which loaded its sense strand, so those seed
  positions are not used.
- Not modelled: structure within purine:purine (G:G and A:G were more selective
  than A:A and G:A at P10, and only A:G held up at 20 nM), and Du's finding that
  A:C, C:A and the G:U wobble are best tolerated, which makes them the weakest
  discriminators. They rank only as "not purine:purine".
- Guide position 1 does not pair with the target (region `anch`). It ranks with
  the seed and, like any difference, counts toward the allele gate: a spared
  allele that differs only at P1 passes the gate with the lowest rank.
- Ambiguity codes have no Huang row and no purine or pyrimidine identity. They
  give `H?` and `?` rather than a guess.
- Huang's twelve rows are asserted verbatim in `tests/test_discrimination.py`,
  but they were transcribed from a secondary source and not checked against the
  paper (as of 2026-10-08). Schwarz and Du were read in full.

The `disc` readout encodes one mismatch, e.g. `P16 post H+ S pp`: guide position,
region (`post` / `seed` / `anch`), Huang status (`H+` listed, `H-` not listed,
`H?` ambiguous), position tier (`S` = P10, P13 or P16; `W` = inside P9-P16; `-`
outside), and class (`pp` purine:purine, `pu` / `up` / `uu`, `?` ambiguous).
`--explain` adds each allele's difference count and flags `agree=no` where Huang
and the regional finding disagree.

## Validation

HUSH is benchmarked against **siDirect** (Naito et al. 2004, 2009), which it
positions itself as a successor to: on matched settings the rule pre-screen
reproduces its picks candidate-for-candidate. The siDirect export and the
matching HUSH design are kept in `data/testing/`.

When comparing in `natural` overhang mode, compare on the **guide** or the target
position, not the passenger: the guide equals siDirect's listed guide, but the
passenger is offset because its overhang comes from the downstream flank.

## Not implemented

A **sequence-motif constraint** on `design` (e.g. forcing a base at a core
position) remains unimplemented. Its original motivation was Pol III initiation
for shRNA, which stage 4 now addresses as the `starts_with_purine` QC flag
rather than a hard gate - see [shrna-conversion.md](shrna-conversion.md).

An **insertion-aware allele gate** is not implemented. Residues a member carries
between two target positions are recorded, attached to the next position, but do
not enter either gate. A targeted allele with an insertion inside the window
still counts as identical, and a spared allele that differs only by an insertion
counts as not discriminating, although in both cases the guide would bulge. An
insertion counts as inside a window only when it lies strictly between two of its
positions. Insertions are not shown in the CLI, TSV or API output. Scoring a
bulge would need mismatch and bulge thermodynamics that `core.thermo` does not
carry.

**Ambiguity codes** in a panel member are compared as plain characters: N
against A is a difference, N against N a match. A panel from IPD-IMGT/HLA has
none; a hand-assembled panel may.

A **machine-learning efficacy layer** is out of scope. The prior-art position is
that modern ML tools confirm efficacy rather than generate candidates; HUSH
generates.
