# shRNA conversion

Stage 4 carries a validated siRNA candidate into a cloneable stem-loop shRNA.
Implementation: `shrna/convert.py` and `shrna/loops.py`.

## What comes out

The bare hairpin insert as top and bottom DNA oligos, in the
sense-loop-antisense orientation, with a Pol III poly-T terminator and **no**
vector-specific restriction tails: the wet lab adds the sticky ends for whichever
vector it clones into.

The design follows the canonical pSUPER system (Brummelkamp et al. 2002,
*Science* 296:550): a Pol III promoter, an inverted repeat, a loop, the
complementary repeat, and a poly-T terminator.

## Geometry

Both strands are encoded **in full** - paired core plus its 3' overhang - so the
hairpin expresses exactly the guide the design stage validated:

```
sense core | sense overhang | loop | guide core | guide overhang | terminator
|___________________ body (folded) ______________|______ 3' tail ______|
```

Three regions, and the distinction between them is what the QC read-outs rest on:

- The **sense overhang** is unpaired where it sits, between stem and loop, so it
  reads as part of the loop.
- The **guide overhang** is a single-stranded 3' end, so it joins the terminator
  in the tail and stays out of the fold search.
- Only the **core** is stem.

This is a correction. Earlier versions built the stem from the 19-nt core alone
and dropped both overhangs, so the hairpin did not encode the 21-nt guide the
upstream gates had validated. If you are comparing against output from before
that change, every construct differs.

**dTdT overhangs are encoded as ribo-U.** An expressed hairpin cannot carry
2'-deoxythymidine; the expressed equivalent of that position is a ribo-U, and the
ordered DNA oligo carries T either way. This also matters in `natural` mode,
which falls back to dTdT at transcript ends - without the normalisation those
candidates could not convert at all.

## Loops

Selectable with `--loop`. All three are documented choices, not arbitrary:

| Name | Sequence | Source |
|---|---|---|
| `jensen` (default) | `UGUGCUU` | Jensen et al. 2012, *PLoS ONE* 7:e43095 |
| `psuper` | `UUCAAGAGA` | Brummelkamp et al. 2002, *Science* 296:550 |
| `trc` | `CUCGAG` | Moffat et al. 2006, *Cell* 124:1283 (pLKO.1 / TRC) |

`jensen` is the default because it was selected from randomised retroviral
libraries for robust knockdown with the **least stem-sequence dependence** -
the safest property for an automated tool that does not know the stem in
advance. `trc` embeds an XhoI site (`CTCGAG` as DNA).

Loop length matters: 8-nt loops perform best overall and a 4-nt loop drops below
50% knockdown (Wang et al. 2013, *PLoS ONE* 8:e56110).

## QC read-outs

Each construct is **scored, not merely stamped** - the difference from tools that
template a pick into a cassette without checking what it folds into.

### folds_as_designed

The body is folded (`core.hairpin`) and the most stable structure compared
against the intended one: stem equal to the core, loop equal to the sense
overhang plus the loop sequence. `False` means a competing structure is more
stable.

One interaction is worth knowing. The `jensen` loop ends in `UU`, so an **`AA`
sense overhang base-pairs with it** and pulls 2 nt of loop into the stem,
extending it to 21 bp and shifting the Dicer processing register. `AA` is a
common passenger overhang under the `natural` policy, so this fires regularly.
It is reported, not corrected: choosing a different loop is the fix.

Stem length sits on a processing boundary in any case - stems of 19 bp and below
behave distinctly as "sshRNA" and can be extremely potent (Ge et al. 2010, *RNA*
16:106).

### The terminator

The deliberate terminator is a run of T (U in the transcript) after the guide
overhang: 6 by default, at least 4 when included (the length is not checked with
`--no-terminator`). T4 is the minimal Pol III termination signal and full
efficiency needs six or more. Termination within the run is heterogeneous,
leaving a U-tail of 1 to 6 nt (Gao et al. 2018). The terminator is excluded from
the premature-terminator check by construction.

### premature_terminator_free

**This one is a filter, not just a flag.** Pol III terminates within a run of
four or more T residues (Bogenhagen & Brown 1981, *Cell* 24:261), with full
efficiency from six; Gao et al. 2018 (*Mol Ther Nucleic Acids* 10:36) confirmed
this in vivo on the human U6, 7SK and H1 promoters. So a `UUUU` anywhere in the
hairpin body truncates the transcript before the hairpin completes.

Grafting can build such a run where no upstream gate can see it: the core's last
base, the overhang, and the loop's first base are three different objects that
only meet after the graft. The upstream homopolymer gate runs before grafting,
is base-agnostic, and inspects the sense core only.

Constructs carrying one are **withheld by default**, counted in the CLI summary,
the API response and the TSV, and recoverable with
`--keep-premature-terminators`.

The check covers the **body only**. A U run where the guide overhang meets the
poly-U terminator is the terminator doing its job at the 3' end, not a premature
one, so the tail is excluded by construction.

Note that retaining overhangs makes this fire more often: with `uu` or `dtdt`
overhangs, a core ending in `U` gives `U + UU + U` across the stem/loop junction.
That is correct behaviour, which is why the withheld count is always reported
rather than the list silently shrinking.

### starts_with_purine

A human U6 or 7SK promoter initiates efficiently and with a precise start site on
A or G at +1 (Gao et al. 2017, *Transcription* 8:275), so a C or U start is
flagged. The widely repeated rule "U6 requires a G" is too narrow, so HUSH checks
for a purine instead. This is a **QC flag, not a gate**: it reports rather than drops,
because the right response depends on the promoter the wet lab uses.

### max_homopolymer_run

Longest single-nucleotide run in the body. The 3' tail's own poly-U is excluded
by construction.

## Scope

**Dicer-processed shRNA only.** AgoshRNA - a very short stem with a small loop,
processed by Argonaute-2 rather than Dicer and yielding only the guide strand
(Gao et al. 2019, *RNA Biology* 16:890) - is out of scope and not exposed as a
mode.

The **miR-30 / miR-E pri-miRNA route** (Fellmann et al. 2013, *Cell Reports*
5:1704) is also out of scope: it needs research effort and wet-lab
infrastructure beyond this project, despite being the route being commercialised
for this application area.

Orientation is **sense-loop-antisense only**, modelled as an enum so
antisense-loop-sense can be added without changing the construct's shape.
