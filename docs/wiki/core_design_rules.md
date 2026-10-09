# Core siRNA design rules

Skeleton only. Each section carries a one-line note on what belongs there and the
references to cite in it; replace the notes with prose as it is written. Citations
are keyed to the groupings in `references.md`.

---

## 1. Why rules at all

*To write:* the 2004 problem - a randomly chosen 21-mer silences its target only
a fraction of the time, and these rule sets were the first reproducible way to
improve that hit rate. State up front that HUSH uses them as a pre-screen filter,
not as a predictor of knockdown.

---

## 2. The mechanism the rules are proxies for

*To write:* asymmetric RISC loading. The strand whose 5' end sits at the
thermodynamically looser end of the duplex is the one Argonaute retains, so every
rule below is one way of biasing that choice toward the guide.

Cite: Khvorova et al. 2003 (internal stability profiles, 5-nt windows);
Schwarz et al. 2003 (complex assembly, terminal base pairs).

---

## 3. The duplex HUSH designs

*To write:* short and factual - a 19-bp paired core with symmetric 2-nt 3'
overhangs, i.e. 21-nt synthesized strands, and the strand naming convention
(sense/passenger, antisense/guide, both stored 5'->3'). This section exists so the
positional numbering used in sections 4-6 is unambiguous.

Cite: Elbashir, Harborth et al. 2001 (canonical duplex);
Elbashir, Martinez et al. 2001 (2-nt 3' overhang optimal, dTdT tolerated).

Implementation: `core.sirna`.

---

## 4. Ui-Tei 2004

*To write:* the four rules and their all-or-nothing semantics.

- (i) A/U at the guide 5' terminus
- (ii) G/C at the sense 5' terminus
- (iii) at least 4 A/U within the guide's 5'-terminal 7 nt
- (iv) no G/C stretch longer than 9 nt

Note worth making: rules (i)-(iii) are three statements of the same asymmetry
principle. Also flag the tension with the seed filter - rule (ii) keeps the
passenger seed GC-rich while rule (iii) pulls the guide seed Tm down - and link
to the seed section rather than arguing it here.

Cite: Ui-Tei et al. 2004.

Implementation: `core.rules` (`UI_TEI_AU_WINDOW`, `UI_TEI_MIN_AU_IN_WINDOW`,
`UI_TEI_MAX_GC_RUN`).

---

## 5. Reynolds 2004

*To write:* the eight criteria, derived from 180 siRNAs against two genes - GC
30-52%, A/U richness at positions 15-19, absence of internal repeat, A at 3,
A at 19, U at 10, no G/C at 19, no G at 13.

Then a clearly marked aside on the scoring, which most secondary sources get
wrong: the original is a score out of 10 with an additive A/U count over positions
15-19 and two negative penalties, not eight equal binary points. Position 19
carries three terms at once, giving an A(+2) / U(+1) / G,C(-1) gradient. Range
-2 to +10, threshold >= 6.

Then the internal-repeat criterion and how HUSH supplies it (hairpin Tm <= 20 C,
overridable with an external folder's verdict).

Cite: Reynolds et al. 2004; Xia et al. 1998 and Groebe & Uhlenbeck 1988 if the
hairpin calculation is explained here rather than linked.

Implementation: `core.rules`, `core.hairpin`.

---

## 6. Amarzguioui & Prydz 2004

*To write:* the additive motif score - S1 (G/C at 1), A6 (A at 6), W19 (A/U at
19) positive; U1 (U at 1), G19 (G at 19) negative - plus the raw duplex-end A/U
differential, threshold >= 3.

Spell out the differential explicitly: it enters as its own value (-3 to +3), not
as a lookup, counted over the three terminal base pairs at each end with the sense
strand as reference, measured 3' end minus 5' end. The sign is the easiest thing
in this whole section to invert, so state why functional siRNAs score positive.

Validation to claim here: all 46 siRNAs of the paper's Table 3 reproduce exactly -
every differential, motif and total - along with the reported class means and
selection statistic.

Cite: Amarzguioui & Prydz 2004.

Implementation: `core.rules`; fixture `tests/data/amarzguioui_table3.tsv`.

---

## 7. What the three rule sets agree on

*To write:* a short comparison table, plus a paragraph on how the same insight is
encoded three ways (binary gate, weighted score, continuous differential).

| Principle | Ui-Tei | Reynolds | Amarzguioui |
| --- | --- | --- | --- |
| Guide 5' end looseness | | | |
| GC content | | | |
| 3'-terminal A/U | | | |
| Internal structure | | | |
| Positional motifs | | | |
| Output type | pass/fail | score out of 10 | additive score |

Decide while writing: the longer-core behaviour (20-24 nt, where 3'-anchored
criteria follow the real 3' terminus while positions 3, 10 and 13 stay at their
absolute 1-based positions) applies to both Reynolds and Amarzguioui. Either
repeat it in both sections or put it in a shared note here.

---

## 8. How HUSH combines them

*To write:* a filter, not a ranking - why the rule sets narrow the candidate field
rather than order it. The default paradigm, and the distinction between
user-configurable options and fixed rule constants (rule thresholds from the
papers are not editable; selection thresholds are).

Implementation: `core.selection`.

---

## 9. Limitations and the modern view

*To write:* the rules are correlative, calibrated on small gene panels, and
superseded for ranking by ML efficacy models. Argonaute reshapes its guide's
binding properties, so naked-RNA thermodynamics is not what happens inside the
complex, and guide position 1 is anchored in the MID pocket rather than pairing at
all. None of this invalidates the rules as empirical filters - it explains why
HUSH keeps them as an interpretable, cheap pre-screen and leaves ranking to a
later stage.

Cite: Salomon et al. 2015; Wee et al. 2012. Forward-link to the ML scoring
section.

---

## 10. References

Full entries with identifiers and "used for" notes live in `references.md`;
keep this list in sync with sections 1 of that file (classical rules), 2
(mechanism), 4 (thermodynamics) and 5 (duplex geometry).

Classical rule sets:

- Ui-Tei, K., Naito, Y., Takahashi, F., Haraguchi, T., Ohki-Hamazaki, H., Juni, A.,
  Ueda, R. & Saigo, K. (2004). Guidelines for the selection of highly effective
  siRNA sequences for mammalian and chick RNA interference. *Nucleic Acids
  Research* 32(3), 936-948. PMID 14769950
- Reynolds, A., Leake, D., Boese, Q., Scaringe, S., Marshall, W.S. & Khvorova, A.
  (2004). Rational siRNA design for RNA interference. *Nature Biotechnology* 22,
  326-330. doi:10.1038/nbt936
- Amarzguioui, M. & Prydz, H. (2004). An algorithm for selection of functional
  siRNA sequences. *Biochemical and Biophysical Research Communications* 316(4),
  1050-1058. doi:10.1016/j.bbrc.2004.02.157

Mechanism:

- Khvorova, A., Reynolds, A. & Jayasena, S.D. (2003). Functional siRNAs and miRNAs
  exhibit strand bias. *Cell* 115(2), 209-216. doi:10.1016/s0092-8674(03)00801-8
- Schwarz, D.S., Hutvagner, G., Du, T., Xu, Z., Aronin, N. & Zamore, P.D. (2003).
  Asymmetry in the assembly of the RNAi enzyme complex. *Cell* 115(2), 199-208.
  doi:10.1016/s0092-8674(03)00759-1
- Salomon, W.E., Jolly, S.M., Moore, M.J., Zamore, P.D. & Serebrov, V. (2015).
  Single-molecule imaging reveals that Argonaute reshapes the binding properties
  of its nucleic acid guides. *Cell* 162(1), 84-95. doi:10.1016/j.cell.2015.06.029
- Wee, L.M., Flores-Jasso, C.F., Salomon, W.E. & Zamore, P.D. (2012). Argonaute
  divides its RNA guide into domains with distinct functions and RNA-binding
  properties. *Cell* 151(5), 1055-1067. doi:10.1016/j.cell.2012.10.036

Duplex geometry:

- Elbashir, S.M., Harborth, J., Lendeckel, W., Yalcin, A., Weber, K. & Tuschl, T.
  (2001). Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured
  mammalian cells. *Nature* 411, 494-498. doi:10.1038/35078107
- Elbashir, S.M., Martinez, J., Patkaniowska, A., Lendeckel, W. & Tuschl, T.
  (2001). Functional anatomy of siRNAs for mediating efficient RNAi in Drosophila
  melanogaster embryo lysate. *EMBO Journal* 20, 6877-6888.
  doi:10.1093/emboj/20.23.6877

Thermodynamics (cite only if restated here rather than linked):

- Xia, T., SantaLucia, J., Burkard, M.E., Kierzek, R., Schroeder, S.J., Jiao, X.,
  Cox, C. & Turner, D.H. (1998). Thermodynamic parameters for an expanded
  nearest-neighbor model for formation of RNA duplexes with Watson-Crick base
  pairs. *Biochemistry* 37(42), 14719-14735. doi:10.1021/bi9809425
- SantaLucia, J. & Hicks, D. (2004). The thermodynamics of DNA structural motifs.
  *Annual Review of Biophysics and Biomolecular Structure* 33, 415-440.
  doi:10.1146/annurev.biophys.32.110601.141800
- Groebe, D.R. & Uhlenbeck, O.C. (1988). Characterization of RNA hairpin loop
  stability. *Nucleic Acids Research* 16(24), 11725-11735.
  doi:10.1093/nar/16.24.11725
