# Design note: siRNA to shRNA conversion via stem-loop scaffolding

Status: **scope fixed, references gathered, design decisions pending.** This note
records the literature basis for HUSH pipeline stage 4 (`shrna/`, currently a
placeholder package), the design parameters the literature fixes, and the open
decisions that must be resolved before the module is specified. It exists so the
conversion step can be built without re-deriving the background.

Date: 2026-09-12.

---

## 1. Scope decision

HUSH commits to the **established stem-loop scaffold only**. The miR-30 / miR-E
pri-miRNA context named as pillar 2 of the concept dossier
([`HUSH_concept_dossier.md`](HUSH_concept_dossier.md), line 38) is **explicitly
out of scope** for the current timeline: it demands both research effort and
wet-lab infrastructure (single-copy expression, backbone processing validation)
the team does not have available. The Fellmann (2013) miR-E work is kept in the
reference set for context only, not as a build target.

This narrows stage 4 to: take a validated siRNA candidate from `selection`
(guide/passenger strands, both stored 5'->3', per `core.sirna` conventions) and
emit a scored stem-loop shRNA construct, rather than naive template-stamping (the
gap left by SSD / siRNA Wizard).

## 2. What the literature fixes about the stem-loop scaffold

A vector-expressed stem-loop shRNA is a well-defined object. The canonical form
is the pSUPER design (Brummelkamp et al. 2002): an RNA Pol III promoter (U6 or
H1) transcribes an inverted repeat - loop - inverted repeat, closed by a poly-T
terminator (5-6 T). After nuclear export Dicer removes the loop, regenerating an
siRNA-like duplex that loads into RISC.

Parameters the literature pins down:

- **Stem length = the siRNA core.** 19-29 bp. A 19-bp stem directly mimics a
  classical siRNA (pSUPER). Stems of 21-29 bp act as Dicer substrates and get
  processed, shifting the exact guide by a few nt. Ge et al. (2009) is the key
  structure-activity study: stems <= 19 bp ("sshRNA") have distinct behaviour and
  can be extremely potent.
- **Loop.** 4-11 nt, short (4-9) preferred. Runs of T/U in the loop must be
  avoided - Pol III reads them as termination. Canonical loops:
    - `TTCAAGAGA` (9 nt) - pSUPER / pSHAG.
    - `CTCGAG` (6 nt, embeds an XhoI site) - pLKO.1 / the TRC library.
    - `TGTGCTT` (7 nt) - Jensen et al. (2012), selected from randomized libraries
      as giving robust knockdown with the least stem-sequence dependence; the
      most design-tool-friendly default for that reason.
    - Liu et al. (2013) found 8-nt loops best overall (`TTCTAGAA`, `TTGGCCAA`), a
      4-nt loop dropping below 50% knockdown, and TTCAAGAGA beating CTCGAG.
- **Orientation.** Both sense-loop-antisense (SA) and antisense-loop-sense (AS)
  are used; the optimum is target-dependent, with no universal winner. pSUPER
  convention is sense-loop-antisense.
- **Pol III initiation.** The U6 promoter strongly prefers a G at the +1
  position. This is the real basis for the "force a base at the graft junction"
  requirement deferred into this module in
  [`design_options.md`](design_options.md) (Planned / deferred section).

## 3. One boundary condition: the AgoshRNA regime

Our siRNA cores are 19-mers, which sits at a processing boundary. A very short
stem (~18-19 bp) combined with a small loop (3-5 nt) can be processed by
Argonaute-2 instead of Dicer ("AgoshRNA", Gao et al. 2019; Berkhout &
Herrera-Carrillo 2021). This yields only the guide strand - eliminating passenger
off-target entirely - but follows different design rules. Because we sit on this
boundary, the module must decide deliberately whether it targets Dicer- or
Ago-processing rather than landing there by accident.

## 4. What we improve on

SSD (Silva et al.) is the dossier's siRNA->shRNA precedent. The paper publishes
neither loop sequences nor folding rules; it stamps its siRNA picks into
mono/bi/trivalent cassettes. HUSH's differentiator (dossier pillar 2) is to
*score* the resulting hairpin. The machinery already exists in the codebase:
`core.hairpin` (self-fold Tm), `core.asymmetry`, and the guide/passenger strands
emitted by `selection`.

## 5. Open design decisions (to resolve before spec)

1. **Loop set** - single evidence-backed default (`TGTGCTT`) vs. a small
   selectable set (pSUPER / TRC / Jensen).
2. **Stem** - use the 19-mer core as-is, or extend toward a Dicer-substrate stem
   (determines whether the guide stays fixed or shifts on processing).
3. **Orientation** - fix sense-loop-antisense, or emit both.
4. **Processing regime** - commit to Dicer-processed shRNA and treat AgoshRNA as
   out of scope, or expose it as a mode.
5. **Cloning context** - output the hairpin RNA only, or a full cloneable cassette
   (promoter + insert + terminator + restriction overhangs). Decides whether the
   Pol III 5'-G constraint and restriction-site avoidance land here.

## 6. References

### Founding scaffold and reviews
- Brummelkamp, T.R., Bernards, R. & Agami, R. (2002). A system for stable
  expression of short interfering RNAs in mammalian cells. *Science* 296, 550-553.
  DOI: 10.1126/science.1068999
- Anthony, R.V. & Cantlon, J.D. (2006). Ribonucleic acid interference: a new
  approach to the in vivo study of gene function. *J Anim Sci* 85(13 Suppl),
  E18-19. DOI: 10.2527/jas.2006-486

### Stem length and loop sequence structure-activity
- Ge, Q. et al. (2009). Minimal-length short hairpin RNAs: the relationship of
  structure and RNAi activity. *RNA* 16(1), 106-117.
- Jensen, S.M.R. et al. (2012). Functional selection of shRNA loops from
  randomized retroviral libraries. *PLoS ONE* 7(8), e43095.
  DOI: 10.1371/journal.pone.0043095
- Liu, Y.P. et al. (2013). A simple and robust vector-based shRNA expression
  system used for RNA interference. *PLoS ONE* 8(2), e56110.
  DOI: 10.1371/journal.pone.0056110

### AgoshRNA / Pol III initiation boundary
- Gao, Z., Berkhout, B. & Herrera-Carrillo, E. (2019). Boosting AgoshRNA activity
  by optimized 5'-terminal nucleotide selection. *RNA Biol* 16(7), 890-898.
  DOI: 10.1080/15476286.2019.1599259
- Berkhout, B. & Herrera-Carrillo, E. (2021). Design and evaluation of AgoshRNAs
  with 3'-terminal HDV ribozymes. *Methods Mol Biol* 2167, 225-252.
  DOI: 10.1007/978-1-0716-0716-9_13

### Templating precedent and scaffold engineering (context)
- Silva et al. SSD - a free software for designing multimeric mono-, bi- and
  trivalent shRNAs. PMC7197978.
- Fellmann, C. et al. (2013). An optimized microRNA backbone for effective
  single-copy RNAi (miR-E). *Cell Reports* 5(6), 1704-1713. Out of scope; kept
  for context.
