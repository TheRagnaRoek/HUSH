# From siRNA to shRNA: the stem-loop conversion

HUSH pipeline stage 4. A validated siRNA candidate is carried into a cloneable
stem-loop shRNA, and the resulting hairpin is scored rather than merely stamped
into a template. Citations are keyed to the groupings in `references.md`
(section 8, shRNA scaffold engineering). Implementation: `shrna.convert`,
`shrna.loops`. Rationale and the decision record: `docs/shrna_design_note.md`.

---

## 1. Why convert at all

- A synthetic siRNA duplex is transfected and silences transiently; it is
  consumed and diluted as cells divide.
- A stem-loop shRNA is **transcribed from a vector**, so silencing is
  continuously replenished and can be made stable and selectable.
- For an allele-specific HLA knockdown intended for cell-therapy work, that
  distinction matters: the reagent must keep working in a population expanded
  over weeks, not for a two-day assay.
- The conversion is therefore not a formatting step. It changes how the guide is
  produced, which introduces an entirely new set of constraints: those imposed
  by RNA Polymerase III and by Dicer.

---

## 2. What a stem-loop shRNA is

- The canonical form is the pSUPER design: an **RNA Pol III promoter** (U6 or
  H1) transcribes an **inverted repeat - loop - inverted repeat**, closed by a
  **poly-T terminator**.
- The transcript folds back on itself into a hairpin: the two inverted repeats
  pair to form the stem, the intervening sequence stays unpaired as the loop.
- After nuclear export, **Dicer removes the loop**, regenerating an siRNA-like
  duplex that loads into RISC. From that point the mechanism is the one the
  siRNA design rules already address.
- Pol III is used rather than Pol II because it transcribes short, precisely
  ended, non-polyadenylated RNAs - exactly the requirement for a defined
  hairpin.

Cite: Brummelkamp et al. 2002 (pSUPER); Anthony & Cantlon 2007 (review context).

---

## 3. Scope: the established stem-loop scaffold only

- HUSH commits to the stem-loop scaffold and **excludes the miR-30 / miR-E
  pri-miRNA route**, descoped on 2026-09-12.
- Reason: miR-E demands research effort plus wet-lab infrastructure the team does
  not have, namely single-copy expression and validation that the backbone is
  processed correctly.
- miR-E is kept in the reference set as context, not as a build target. It is
  the stronger scaffold for single-copy work and is being commercialised for
  this exact application area, so it belongs in future work rather than nowhere.

Cite: Fellmann et al. 2013 (miR-E); US patent 11384335 (B2M shRNAmiR, miR-E
backbone).

---

## 4. The stem is the validated core, unchanged

- The 19-bp duplex core from stage 3 becomes the stem **verbatim**.
- Consequence, and the main design claim of this stage: the guide that was
  validated upstream for seed-Tm, asymmetry and off-target is **exactly** the
  guide the hairpin expresses. No re-shift, no re-validation.
- The alternative was to extend the stem to 21-29 bp so the hairpin acts as a
  **Dicer substrate**. Dicer then cuts at its own measured register, shifting
  the mature guide by a few nucleotides, and every upstream result would have to
  be recomputed for a guide that is no longer the one that was scored.
- Holding the stem at 19 bp keeps a single guide definition across the whole
  pipeline. The cost is that the construct sits on a processing boundary
  (section 10).

Cite: Elbashir, Harborth et al. 2001 (the 19-bp paired core the stem reproduces);
Ge et al. 2010 (stem length structure-activity).

Implementation: `shrna.convert.convert`.

---

## 5. Loop choice

- The loop measurably affects knockdown, so HUSH offers a **small registry of
  established loops** rather than one hard-coded choice:

| Name | Sequence (RNA) | Length | Source |
|------|----------------|--------|--------|
| `jensen` (default) | UGUGCUU | 7 nt | Jensen et al. 2012 |
| `psuper` | UUCAAGAGA | 9 nt | Brummelkamp et al. 2002 (pSUPER / pSHAG) |
| `trc` | CUCGAG | 6 nt | Moffat et al. 2006 (pLKO.1 / TRC; embeds XhoI) |

- **Why Jensen is the default:** it was selected from randomized libraries for
  robust knockdown with the **least stem-sequence dependence**. That is the
  decisive property for an automated tool, which must graft a loop onto a stem
  it has not seen and cannot test.
- Length evidence: 8-nt loops performed best overall (TTCTAGAA, TTGGCCAA), a
  4-nt loop fell below 50% knockdown, and TTCAAGAGA outperformed CTCGAG. Useful
  range 4-11 nt, with 4-9 preferred.
- Loops must not contain **runs of 4 or more T/U**, which Pol III reads as
  termination (section 7).
- `trc` is retained because its XhoI site is convenient for pLKO.1 cloning, not
  because it is the best performer. Note that the loop and its flanks are
  documented in the TRC manual rather than in Moffat et al. itself.

Cite: Jensen et al. 2012 (default loop); Wang et al. 2013 (loop length and
comparison; note this paper was long miscited as "Liu et al."); Moffat et al.
2006 and Stewart et al. 2003 (the CTCGAG loop).

Implementation: `shrna.loops` (`LOOPS`, `DEFAULT_LOOP`, `get_loop`).

---

## 6. Orientation

- Two arrangements are in use: **sense-loop-antisense (SA)**, the pSUPER
  convention, and antisense-loop-sense (AS).
- The literature gives **no universal winner**; the optimum is target-dependent.
- HUSH emits **SA only**, but models orientation as an enum so AS can be added
  without changing the shape of the construct.
- Being explicit about this is the point: an arbitrary choice recorded as
  arbitrary is honest, whereas emitting both without evidence would double the
  cloning burden for no predicted gain.

Implementation: `shrna.convert.Orientation`.

---

## 7. Pol III constraints

Two constraints come from the polymerase rather than from RNAi, and both are
checked per construct.

- **Purine at +1.** Human U6 and 7SK promoters initiate efficiently and with a
  precise start site on **A or G**. Reported as `starts_with_purine`; a C or U
  start is flagged.
  - Worth stating plainly: the widely repeated rule is "U6 requires a G". The
    primary literature shows A works too, so HUSH deliberately checks for a
    purine rather than inheriting the narrower folklore rule.
- **Premature termination.** Pol III terminates within a run of **4 or more T**,
  with full efficiency at 6 or more. The hairpin body (stem plus loop) must
  therefore contain no such run, **including at the two loop/stem junctions**,
  where a run can be created by concatenation and is invisible to the upstream
  per-strand homopolymer gate.
  - Reported as `premature_terminator_free`.
- The **deliberate** terminator is appended at the 3' end: 6 T by default,
  minimum 4, and it is excluded from the junction check by construction.

Cite: Gao et al. 2017 (+1 start site, A or G); Bogenhagen & Brown 1981 (the
4-T termination signal); Gao et al. 2018 (in-vivo confirmation on U6, 7SK and
H1; T4 minimal, 6 or more for full efficiency).

Implementation: `shrna.convert` (`_POL3_TERMINATOR_RUN`,
`MIN_TERMINATOR_LENGTH`, `DEFAULT_TERMINATOR_LENGTH`).

---

## 8. Scoring the hairpin

This is what separates HUSH from the naive templating of earlier tools.

- Prior art (SSD, InvivoGen siRNA Wizard) **stamps** siRNA picks into a hairpin
  template. SSD publishes neither its loop sequences nor any folding rule.
- A grafted hairpin is a new molecule, and nothing guarantees it folds into the
  structure that was drawn. HUSH therefore **folds the construct and checks**.
- `folds_as_designed` is true only when the most stable fold of the stem plus
  loop is the **intended** one: stem length equal to the full core, closing the
  intended loop. Anything else is a misfold risk, since a construct that prefers
  a competing structure will not be processed as planned.
- The poly-T terminator is excluded from the fold search: it is a
  single-stranded tail, not part of the hairpin, and including it would let the
  search find structures that do not exist in the folded transcript.
- Reported alongside: `max_homopolymer_run` over the stem and loop, and the two
  Pol III flags from section 7.
- These are **legible flags, not a black-box score**. A team can see which
  constraint a construct failed and decide whether it matters for their vector.

Cite: Silva et al. (SSD, the templating precedent); InvivoGen siRNA Wizard.

Implementation: `shrna.convert` (`folds_as_designed`), `core.hairpin`
(`most_stable_hairpin`).

---

## 9. What the converter emits

- The transcribed **hairpin RNA**, 5'->3', terminator included as U's.
- The orderable **DNA oligos**: `insert_top` and `insert_bottom`, where the
  bottom strand is the reverse complement of the top. Annealed, they form the
  double-stranded insert.
- The source candidate's coordinates (`region`, `start`, `end`), so every
  construct is traceable back to a position on the target.
- **No vector-specific restriction overhangs.** The wet lab adds sticky-end
  tails for whichever vector it clones into.
  - Rationale: overhangs are a property of the chosen backbone, not of the
    design. Emitting them would silently commit the team to one vector and
    make the output wrong for every other.
- Consequence: restriction-site avoidance stays with the wet lab, while the
  purine-at-+1 check lands here as a QC flag.

Implementation: `shrna.convert.ShrnaConstruct`.

---

## 10. The processing boundary, and limitations

- A 19-bp stem sits on a **processing boundary**. A very short stem (around
  18-19 bp) combined with a small loop (3-5 nt) can be cleaved by
  **Argonaute-2 instead of Dicer** - the AgoshRNA regime.
- AgoshRNA is attractive in principle: it yields **only the guide strand**,
  eliminating passenger-mediated off-target entirely. But it follows different
  design rules, including its own 5'-terminal nucleotide preferences.
- HUSH commits to the **Dicer-processed** regime and does not expose AgoshRNA as
  a mode. The boundary is documented so the choice is deliberate rather than
  accidental, which is the whole reason this section exists.
- Honest limitations of the current stage:
  - Stems <= 19 bp ("sshRNA") are known to behave distinctly and can be
    extremely potent; HUSH does not exploit this.
  - No orientation sweep (SA only), no cassette assembly (no promoter or
    vector context), no multimeric constructs.
  - The fold check confirms the intended **structure**, not knockdown. It
    removes constructs that cannot work as drawn; it does not predict which of
    the survivors works best.

Cite: Gao et al. 2019 and Berkhout & Herrera-Carrillo 2021 (AgoshRNA);
Ge et al. 2010 (sshRNA potency).

---

## 11. References

All entries for this section, with identifiers and a note on what HUSH takes
from each, are in `references.md` section 8 (shRNA scaffold engineering).
Section 5 of that file covers the siRNA duplex geometry the stem reproduces.
