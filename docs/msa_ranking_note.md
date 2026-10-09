# Design note: ranking allele-specific candidates by mismatch placement

Status: **implemented.** `msa/discrimination.py` ranks candidate windows by where
their difference to a spared allele sits on the guide strand. This note records
the decisions behind it, and in particular the one thing the ranking refuses to
do, so neither has to be re-derived.

Date of decision: 2026-10-08. Supersedes the first version of the ranking, which
is discussed in §6.

---

## 1. What problem this solves, and what it does not

In protein families like HLA you often find only a single difference between two
alleles inside a region suitable for RNAi. Established tools do not help here:
siDirect requires several mismatches before it calls an off-target ruled out. So
every window that reaches this stage already *formally* discriminates — it
matches each allele to be silenced and differs from each allele to be spared.

What differs between those windows is how reliably the difference is expected to
work. The published literature gives an ordering over **mismatch placements**. It
does not give an ordering over whole candidates, and this module does not pretend
otherwise:

- **Ranking never removes a candidate.** Filtering happens earlier, in
  `msa/windows.py`. A badly ranked window is still a window that satisfies the
  request.
- **A high rank is not a prediction of silence on the spared allele.** Du et al.
  found 73% of their singly mismatched targets still significantly silenced. A
  difference is necessary for discrimination, not sufficient.

## 2. The three sources

| | design | what it gives us |
|---|---|---|
| Schwarz et al. 2006 | tiled siRNAs against a SOD1 point mutation, five assays | the seed/post-seed mechanism; the positions that replicate; purine:purine as a class |
| Huang et al. 2009 | every mismatch at every position | Table 1: a lookup keyed by the exact base change |
| Du et al. 2005 | all 57 single-nucleotide permutations of the *target site* | low / medium / high tolerance regions |

**Du's coordinates.** Du reports *target-site* positions, so they have to be
converted. The paper prints both strands of siCD46:

```
target site  5'-CTTATTGGAGAGAGCACGA-3'   (19 nt)
guide        5'-TCGTGCTCTCTCCAATAAG-3'   = exact reverse complement
```

so `guide = 20 - target`, and Du's low-tolerance region (target 5–11) is guide
**P9–P15**. `tests/test_discrimination.py::TestLowToleranceWindow` asserts the
reverse-complement relationship rather than trusting the arithmetic.

**Where they converge.** P9–P15 is low tolerance in Du, holds P9, P10, P12, P13
and P14 of Schwarz's named positions, and contains the bulk of Huang's table (54
of its 55 position entries fall in P9–P16).

**Where they disagree, and how we resolve it.**

- *P16.* Schwarz's best-replicated position — it discriminated in all five of
  their assays and carries their own generality test: ten independent
  siRNA–mRNA pairs against two targets, a higher IC₅₀ for the mismatched
  combination in every pair and at least 20-fold in seven. Huang list it in 7 of
  12 rows. Du place it one tier lower, in medium tolerance. We keep it inside the
  window, on the strength of two sources against one, and record the dissent here
  and in the module docstring.
- *P3, P4, P5.* Schwarz saw strong discrimination at these seed positions, but
  from a single siRNA whose seed pairing they describe as the most
  thermodynamically stable of their set — and their own strand-loading analysis
  (Table 1) shows that siRNA loaded predominantly its **sense** strand, so its
  apparent discrimination is confounded. Not used as a criterion. An earlier
  version of this module did credit P5; see §6.

## 3. The ranking

Ordinal, non-compensatory, lexicographic. Criteria are consulted in a fixed
order; the first that separates two candidates decides; a deficit in a
higher-priority criterion cannot be offset by a lower-priority one. No score is
assigned, so no weights have to be calibrated.

| # | criterion | condition | source |
|---|---|---|---|
| 1 | listed for this exact base change | the substitution appears at that position in Table 1 | Huang 2009 |
| 2 | inside the low-tolerance window | P9–P16 | all three |
| 3 | at a position replicated in every assay | P10, P13, P16 | Schwarz 2006 |
| 4 | elsewhere 3′ of the seed | P9–P19 | Schwarz 2006 |
| 5 | mismatch class | purine:purine > other > A·C, C·A, G·U wobble | Schwarz 2006; Du 2005 |

Criterion 1 is first because it is the only one conditioned on the actual
substitution rather than on position alone. Criterion 4 ranks *below* the window
because Du put P16–P19 in the medium- and high-tolerance regions: "anywhere 3′ of
the seed" carries less than the window does, even though the underlying mechanism
(seed mismatches block binding and are overwhelmed by concentration; mismatches
3′ of it block catalysis and are not) is what makes the window plausible in the
first place. Criterion 5 is last and is the one entry a reader may reasonably
want moved; `DISCRIMINATION_ORDER_NOTE` says what that costs.

The ranking is applied to the best-placed difference to the **weakest spared
allele** — a maximin, matching `min_spare_mismatches`. A design is only as
allele-specific as its least-well-separated off-allele, so an average would hide
exactly the case that matters.

Positions are guide-strand, 1-based from the 5′ end; a mismatch is written
guide·spared-allele.

## 4. What is reported but never ranked

**The number of differences is not a criterion.** All three sources measure
single mismatches. None of them says whether three poorly placed differences
outweigh one well-placed difference, so any threshold would assert something the
literature does not support.

The consequence is real and worth stating plainly: a window with one
Huang-listed difference will rank above a window with three differences that are
merely post-seed. That is not a claim that the first is the better design. It is
the ranking saying only what it can say — this one's best-placed difference is
better placed — and handing the rest to the user.

Which means the rest has to actually reach the user. The rank rests on one
mismatch, so every renderer carries the whole set alongside it:

- `WindowEvaluation.deciding_mismatches` — every difference to the weakest
  spared allele, in ranking order
- CLI table — `minDiff` plus the `placement` cell, which appends `+P11,P4`
- CLI TSV — `disc_all_positions`
- CLI `--explain` / shell `assess` — every spared allele, every position
- API — `DiscriminationResult.all_guide_positions`

`MULTI_MISMATCH = 2` survives as a reporting flag only, for comparison with the
siDirect convention. It must not re-enter `discrimination_key` without a source.

## 5. Scope and deferrals

- Within purine:purine there is further structure — G:G and A:G were more
  selective than A:A and G:A at P10, and A:G alone held up at 20 nM — which is
  not modelled.
- Du's finding that A·C, C·A and the G·U wobble are the *best* tolerated, i.e.
  the weakest discriminators, is not modelled either; they currently rank only as
  "not purine:purine". Demoting them explicitly is the most obvious next
  refinement, and it is well sourced.
- Huang's table and the regional finding disagree in places (A>C lists guide
  position 8, inside the seed). Reported via `MismatchAssessment.sources_agree`,
  not resolved.
- Guide position 1 sits in Argonaute's MID pocket and does not pair, so a
  "mismatch" there is not one. Reported as `GuideRegion.ANCHOR`; it ranks with
  the seed because no measurement separates them.
- IUPAC ambiguity codes yield `None` rather than a guess.
- **Huang's twelve table rows have not been checked against the PDF.** They are
  reproduced in `HUANG_RECOMMENDED_POSITIONS` and asserted verbatim in the tests,
  but the source for them is secondary. Schwarz and Du have been read in full.

## 6. What changed, and why

The first version of this module ranked differently. The changes, and what forced
each one:

| was | is | why |
|---|---|---|
| count ≥ 4 as the top tier, count again as the final tiebreaker | count not ranked at all | No source supports 4, or 3. The threshold had been derived from `BOWTIE_MAX_V = 3` in `core/offtarget.py`, which is a **limitation of bowtie1's `-v` mode**, not a biological finding. Ohnishi 2008 had been cited for "two mismatches beat one", but there the second mismatch is *engineered into the siRNA* and mismatches both alleles, and the natural SNP is central — it does not support the claim. |
| `STRONGEST_POSITIONS = {10, 16}` | `{10, 13, 16}` | Counting Schwarz's five reported sets: P13 and P16 appear in all five, P10 in four. P13 had been in the lower tier. |
| `CONSISTENT_POSITIONS` included P5 | P5 dropped; set is `{9, 12, 14}`, reported but not a ranking tier | Schwarz's seed result comes from one siRNA with an unusually stable seed, and Table 1 shows that siRNA loaded its sense strand. |
| region (3′ of seed) ranked second | window P9–P16 second, region fourth | Du's regions, once converted, are far more specific than "3′ of the seed" and agree with the other two. |
| purine:purine attributed to the ordering convenience, then briefly to Du | attributed to Schwarz, corroborated by Du | Schwarz test it as a class and it is their abstract's second finding. Du never use the term but report the same pattern per base pair. |
| only the best-placed mismatch surfaced | every difference surfaced | Once count left the ranking, collapsing the output to the best mismatch meant the user could see neither the count's detail nor the other positions, and so could not make the judgement the module had just declined to make. |

Public names (`discrimination_key`, `--sort discrimination`, the API field names)
were **deliberately not renamed**, although "placement" would describe them
better. They are part of the service and HTTP contract the web client speaks, and
renaming them would have meant touching the GUI. The CLI column header is the one
place the better word appears.

## 7. For the web client

Nothing under `web/` or `service/` was touched, and **the app keeps working
as-is**: the response shape is backwards compatible and `sort=discrimination` is
still the string the service expects. But a few labels now say something that is
no longer true, and there is one new field worth showing.

### What changed on the wire

`DiscriminationResult` (`api/schemas.py` → `web/src/api/types.ts:132`) gains one
field:

```ts
all_guide_positions: number[];   // NEW, defaults to [] - every difference to
                                 // this allele, in ranking order
```

Two existing fields changed **meaning** without changing shape, which is the part
that needs a human eye:

| field | was | is |
|---|---|---|
| `strongest_position` | P10 or P16 | P10, P13 or P16 |
| `consistent_position` | Schwarz's second tier (P5, P9, P12, P13, P14) | **inside the P9–P16 window** |
| `multi_mismatch` | a ranking criterion | reporting only — never ranked |

`consistent_position` kept its name purely so existing clients keep parsing. It
is not load-bearing anywhere else; renaming it on the wire is a two-line change
in both repos if that is preferred.

Unchanged: `/design` request params, `AlleleResult`, `PerAlleleDiscrimination`,
`spare_min_mismatch`, `diff_count`.

### Roughly what to do in the client

Rough, not prescriptive.

**a) Fix the stale labels.** `web/src/lib/disc.ts`:

```ts
export const TIER = {
  S: "a strongest position (P10 / P16)",   // -> P10 / P13 / P16
  C: "a merely consistent position",       // -> inside the P9-P16 window
  "-": "neither tier",                     // -> outside the window
};
```

`tierOf()` needs no logic change, only the strings. The CLI now writes `W` rather
than `C` for the same thing, if the two should read alike. The file header
comment points at `docs/design_options.md` and should point here instead.

**b) Show every difference, not just the deciding one.** This is the real work,
and the reason the new field exists. `DiscriminationBreakdown.tsx` currently
prints `{diff_count} differences in window` and then grades only the best one.
With `all_guide_positions`:

```
HLA-A*03:01    3 differences    P11 · P15 · P4     [deciding allele]
               └ graded chips for P11, the one the rank rests on
```

Ranking order, not numeric order — the first entry is always the deciding one. A
muted style for entries 2..n would say "these exist but did not decide", which is
exactly right.

**c) Say somewhere that count is not ranked.** One line of hint text under the
candidate list, e.g.:

> Ranked by where the difference sits on the guide. The number of differences is
> shown but not ranked — no study measures one against the other.

Without it, a user who sees a 1-difference candidate above a 3-difference one
will read it as a bug. `CandidateList.tsx:31` already shows the `minDiff` chip,
so the two numbers sit next to each other and that sentence closes the gap. §4
above is the argument behind it.

**d) Optional, low effort.** `ScoreReport.tsx:115` renders the one-line token
(`P16 post H+ S pp`); appending `+P11,P4` matches what the CLI now prints.

**e) Mock data.** `web/src/lib/sample.ts:37` builds a `DiscriminationResult` by
hand and needs `all_guide_positions` added, or the sample view shows an empty
list.

### What does not need touching

The sort control (`"discrimination"` is unchanged), `AlignmentTrack.tsx`
(`discProse()` keeps working; only its wording changes via `disc.ts`), and
anything in `Funnel`, `GateRow` or `TargetingStep` — filtering is unchanged.

---

## References

- Schwarz, D.S. et al. (2006). Designing siRNA that distinguish between genes
  that differ by a single nucleotide. *PLoS Genetics* 2(9), e140.
  doi:10.1371/journal.pgen.0020140
- Huang, H. et al. (2009). Profiling of mismatch discrimination in RNAi enabled
  rational design of allele-specific siRNAs. *Nucleic Acids Research* 37(22),
  7560–7569. doi:10.1093/nar/gkp835
- Du, Q. et al. (2005). A systematic analysis of the silencing effects of an
  active siRNA at all single-nucleotide mismatched target sites. *Nucleic Acids
  Research* 33(5), 1671–1677. doi:10.1093/nar/gki312
