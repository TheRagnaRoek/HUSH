# shRNA conversion

Stage 4 converts a validated siRNA candidate into a cloneable stem-loop shRNA.
Implementation: `shrna/convert.py` and `shrna/loops.py`.

## Output

The output is the hairpin insert as top and bottom DNA oligos, in
sense-loop-antisense orientation, with a Pol III poly-T terminator. It has **no**
vector-specific restriction tails; the sticky ends for the target vector are added
in the lab.

The design follows the canonical pSUPER system (Brummelkamp et al. 2002,
*Science* 296:550): a Pol III promoter, an inverted repeat, a loop, the
complementary repeat, and a poly-T terminator.

## Geometry

Both strands are encoded **in full**, each as paired core plus 3' overhang, so the
hairpin expresses exactly the guide that the design stage validated:

```
sense core | sense overhang | loop | guide core | guide overhang | terminator
|___________________ body (folded) ______________|______ 3' tail ______|
```

The QC read-outs distinguish three regions:

- The **sense overhang** lies between stem and loop and is unpaired, so it is
  treated as part of the loop.
- The **guide overhang** is a single-stranded 3' end, so it is part of the 3' tail
  with the terminator and is excluded from the fold search.
- Only the **core** forms the stem.

Earlier versions built the stem from the 19-nt core alone and omitted both
overhangs, so the hairpin did not encode the 21-nt guide that the upstream gates
had validated. Output from before this change differs for every construct.

**dTdT overhangs are encoded as ribo-U.** An expressed hairpin cannot contain
2'-deoxythymidine, so the position is encoded as ribo-U; the ordered DNA oligo
contains T in both cases. This also applies in `natural` mode, which falls back to
dTdT at transcript ends. Without the conversion, those candidates could not be
converted.

## Loops

Three loops are available with `--loop`, each from a published source:

| Name | Sequence | Source |
|---|---|---|
| `jensen` (default) | `UGUGCUU` | Jensen et al. 2012, *PLoS ONE* 7:e43095 |
| `psuper` | `UUCAAGAGA` | Brummelkamp et al. 2002, *Science* 296:550 |
| `trc` | `CUCGAG` | Moffat et al. 2006, *Cell* 124:1283 (pLKO.1 / TRC) |

`jensen` is the default because it was selected from randomised retroviral
libraries for robust knockdown with the **least dependence on the stem sequence**.
This suits an automated tool that does not know the stem in advance. `trc`
contains an XhoI site (`CTCGAG` as DNA).

Loop length affects knockdown: 8-nt loops perform best overall, and a 4-nt loop
gives less than 50% knockdown (Wang et al. 2013, *PLoS ONE* 8:e56110).

## QC read-outs

Each construct is folded and checked. It is not inserted into a template without
checking its structure.

### folds_as_designed

The body is folded (`core.hairpin`) and the most stable structure is compared with
the intended one: stem equal to the core, loop equal to the sense overhang plus the
loop sequence. `False` means that a competing structure is more stable.

**Known interaction.** The `jensen` loop ends in `UU`, so an `AA` sense overhang
pairs with it. This pulls 2 nt of the loop into the stem, extends the stem to
21 bp and shifts the Dicer processing site. `AA` is a common passenger overhang
under the `natural` policy, so this occurs often. HUSH reports it and does not
correct it; choosing a different loop avoids it.

Stem length also matters for processing: stems of 19 bp or shorter behave
differently ("sshRNA") and can be very potent (Ge et al. 2010, *RNA* 16:106).

### The terminator

The terminator is a run of T (U in the transcript) after the guide overhang: 6 by
default, at least 4 when included (the length is not checked with
`--no-terminator`). T4 is the minimal Pol III termination signal, and full
efficiency needs six or more. Termination within the run is heterogeneous and
leaves a U-tail of 1 to 6 nt (Gao et al. 2018). The terminator is not part of the
premature-terminator check.

### premature_terminator_free

This read-out is also a filter. Pol III terminates within a run of four or more T
residues (Bogenhagen & Brown 1981, *Cell* 24:261), with full efficiency from six;
Gao et al. 2018 (*Mol Ther Nucleic Acids* 10:36) confirmed this in vivo on the
human U6, 7SK and H1 promoters. A `UUUU` anywhere in the hairpin body therefore
truncates the transcript before the hairpin is complete.

Grafting can create such a run at a position that no upstream gate checks: the
last base of the core, the overhang and the first base of the loop are only
adjacent after grafting. The upstream homopolymer gate runs before grafting, does
not depend on the base, and checks the sense core only.

Constructs with such a run are **withheld by default**. The count appears in the
CLI summary, the API response and the TSV, and `--keep-premature-terminators`
retains the constructs.

The check covers the **body only**. A U run where the guide overhang meets the
poly-U terminator belongs to the intended terminator and is not premature, so the
3' tail is excluded.

Because overhangs are now retained, this check triggers more often: with `uu` or
`dtdt` overhangs, a core ending in `U` gives `U + UU + U` across the stem/loop
junction. This is correct behaviour. The withheld count is always reported so that
the list does not shrink without notice.

### starts_with_purine

A human U6 or 7SK promoter initiates efficiently and with a precise start site on
A or G at +1 (Gao et al. 2017, *Transcription* 8:275), so a C or U start is
flagged. The widely repeated rule "U6 requires a G" is too narrow, so HUSH checks
for a purine. This read-out is a **flag, not a gate**: constructs are reported and
not removed, because the appropriate response depends on the promoter used in the
lab.

### max_homopolymer_run

The longest single-nucleotide run in the body. The poly-U run of the 3' tail is
excluded.

## Scope

**Only Dicer-processed shRNA is supported.** AgoshRNA has a very short stem with a
small loop, is processed by Argonaute-2 instead of Dicer and yields only the guide
strand (Gao et al. 2019, *RNA Biology* 16:890). It is out of scope and is not
available as a mode.

The **miR-30 / miR-E pri-miRNA route** (Fellmann et al. 2013, *Cell Reports*
5:1704) is also out of scope: it would require research and wet-lab infrastructure
beyond this project, although it is the route being commercialised for this
application area.

Only the sense-loop-antisense orientation is implemented. It is an enum, so
antisense-loop-sense can be added without changing the structure of the construct.
