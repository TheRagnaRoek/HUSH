# HUSH design notes

Design rationale, literature provenance, formulae, and validation records for
the stage-2 siRNA design engine. This is the long-form companion to the module
docstrings: the modules say *what* they do and the critical caller contract; the
*why*, the derivations, and the "checked against the paper's own figures" records
live here so the code stays readable.

Related documents:

- [`HUSH_concept_dossier.md`](HUSH_concept_dossier.md) - the project vision and
  prior-art landscape.
- [`seed_offtarget_design_note.md`](seed_offtarget_design_note.md) - the specific
  critique of the seed-Tm filter and the proposed replacement.
- [`design_options.md`](design_options.md) - a one-page user reference to every
  `hush design` / `SelectionCriteria` option, its default and effect.
- The `README.md` - user-facing usage.

Contents:

- [Nearest-neighbor duplex thermodynamics (`core.thermo`)](#nearest-neighbor-duplex-thermodynamics-corethermo)
- [Seed-duplex Tm off-target filter (`core.seed`)](#seed-duplex-tm-off-target-filter-coreseed)
- [Duplex-end asymmetry (`core.asymmetry`)](#duplex-end-asymmetry-coreasymmetry)
- [Intramolecular hairpin Tm (`core.hairpin`)](#intramolecular-hairpin-tm-corehairpin)
- [Classical rule sets (`core.rules`)](#classical-rule-sets-corerules)
- [Candidate selection (`core.selection`)](#candidate-selection-coreselection)

---

## Nearest-neighbor duplex thermodynamics (`core.thermo`)

### The primitive

Every thermodynamic idea the design stage relies on - RISC loading asymmetry,
seed-duplex off-target filtering, allele-discriminating mismatch penalties -
reduces to one primitive: the free energy of a short RNA helix. `core.thermo` is
that primitive and nothing more. It converts a sequence into `(dH, dS)`, from
which dG at any temperature and the melting temperature follow.

The model is the nearest-neighbor (NN) approximation: the stability of a duplex
is the sum of the stacking contributions of its overlapping dinucleotide steps,
plus a per-duplex initiation term and end corrections.

### Two entry points, deliberately separate

- `duplex_thermo()` - a complete, isolated duplex. Includes initiation, the
  terminal A-U penalties, and the self-complementarity correction. This is what a
  free-standing oligo duplex (e.g. an isolated seed duplex) needs.
- `stacking_terms()` - the summed stacking contributions alone, with no
  initiation or end corrections. This is what a window *inside* a longer helix
  needs, where the window's boundaries are not real helix ends (e.g. comparing
  the two ends of an siRNA duplex for loading asymmetry, or scoring a hairpin
  stem).

Using the wrong one silently shifts every number by the initiation term, so the
choice is made explicit at the call site rather than hidden behind a flag.

### Two parameter sets, named explicitly

Two NN parameter sets are provided: Xia et al. (1998), the default and current
standard, and Freier et al. (1986). This is not redundancy. A published threshold
is only valid against the parameter set it was calibrated on - siDirect's 21.5 C
seed-Tm cut-off was fit to Freier numbers (see the seed section) - so the caller
names the set explicitly rather than getting whichever is currently considered
best. Scoring a Freier-calibrated threshold with Xia values, or vice versa,
silently shifts every number relative to the cut-off.

### How the tables are stored

A nearest-neighbor step is written 5'-XY-3' paired with 3'-X'Y'-5'. Because the
partner strand is implied by Watson-Crick pairing, the top-strand dinucleotide
alone identifies the step. Only 10 of the 16 steps are independent: reading a
duplex from the other strand turns step XY into its reverse-complement, and the
two must carry the same energy (5'-CA-3'/3'-GU-5' *is* 5'-UG-3'/3'-AC-5' seen from
the other side). The code stores the 10 published values and derives the other 6
by that symmetry (`_expand_by_symmetry`), so no value is typed twice and the two
spellings of a step can never disagree. Values are `(dH kcal/mol, dS cal/mol/K)`,
measured in 1 M NaCl.

### Scope and deferrals

- Watson-Crick pairs only. Mismatch, G-U wobble, and dangling-end parameters are
  not included; allele-discrimination scoring will need them.
- No salt correction. Both parameter sets are reported in their 1 M NaCl
  reference state. A published score that specifies its own salt condition applies
  its own correction (see the seed section).
- No siRNA-specific thresholds. Those live with the score that defines them.

### Validation

The Tm agrees with Biopython's independent implementation of the same Xia table
(`RNA_NN2`) to within 0.007 C across a range of sequences; the residual is
entirely Biopython's gas constant rounded to 1.987 against the accurate 1.9872.
Enforced in `tests/test_thermo.py`.

### References

- Primary: Xia, T. et al. (1998). Thermodynamic parameters for an expanded
  nearest-neighbor model for formation of RNA duplexes with Watson-Crick base
  pairs. *Biochemistry* 37(42), 14719-14735. doi:10.1021/bi9809425
- Freier, S.M. et al. (1986). Improved free-energy parameters for predictions of
  RNA duplex stability. *PNAS* 83(24), 9373-9377.
- SantaLucia, J. & Hicks, D. (2004). The thermodynamics of DNA structural motifs.
  *Annu Rev Biophys Biomol Struct* 33, 415-440. (Tm formulation.)

---

## Seed-duplex Tm off-target filter (`core.seed`)

### What the filter is for

An siRNA silences its intended target through near-perfect complementarity, but
it also behaves like a microRNA: the seven nucleotides at positions 2-8 from the
5' end of the RISC-loaded strand - the *seed* - are enough to knock down any
transcript carrying the complementary site in its 3' UTR. There are only 16384
possible 7-mers against tens of thousands of 3' UTRs, so seed matches cannot be
avoided by sequence choice alone.

Ui-Tei et al. (2008) showed the way out: the strength of the off-target effect
tracks the thermodynamic stability of the seed-target duplex, so a seed that
pairs *weakly* is largely off-target free even where its matches exist. siDirect
2.0 (Naito et al. 2009) turned that into a filter at Tm < 21.5 C.

Both strands are filtered. The passenger strand is loaded into RISC often enough
to drive its own seed-dependent off-target silencing, so siDirect requires both
strands' seed duplexes to fall below the threshold.

### The exact formula (Naito et al. 2009)

    Tm = 1000 * dH / (A + dS + R * ln(Ct / 4)) - 273.15 + 16.6 * log10([Na+])

with `A` the helix initiation constant (-10.8), `R = 1.987` cal/(mol*K), `Ct` the
total strand concentration (100 uM), and `[Na+]` fixed at 100 mM.

### Why it is not built on `core.thermo`

Three things differ from HUSH's general nearest-neighbor engine, and all three
are deliberate - the 21.5 C threshold is only meaningful in combination with
exactly these choices:

1. **Freier et al. (1986) parameters, not Xia et al. (1998).** The initiation
   constant of -10.8 is Freier's. The 21.5 C threshold is an empirical cut-off
   calibrated on Freier numbers; recomputing the Tm with the newer, better Xia
   parameters and comparing against the same cut-off would silently change which
   candidates pass. (See `seed_offtarget_design_note.md`: the two parameter sets
   re-rank only ~2.3% of heptamers, so the swap is a near-no-op that would
   nonetheless invalidate the threshold.)
2. **No terminal A-U penalty, no symmetry correction, and a division by 4
   applied unconditionally.** So the seed Tm is built from
   `core.thermo.stacking_terms()` with the Freier table, not from
   `duplex_thermo()`.
3. **The paper's rounded gas constant `1.987`,** rather than `core.thermo`'s more
   precise `1.9872`. Reproducing a published number requires the published
   constant.

### Validation

Computing this Tm for all 4^7 = 16384 possible seed heptamers reproduces the
paper's own figures exactly: 4488 of them (27.4%) fall below 21.5 C, spanning
-12.0 C to 59.6 C (the paper reports "ranged from -12 C to 60 C"). This single
check pins down every constant at once - the parameter set, the -10.8 initiation
term, the 100 uM concentration, the 100 mM salt correction, and the absence of
the terminal-A-U and symmetry corrections. Enforced in `tests/test_seed.py`.

### Status and known limitations

The filter is kept siDirect-faithful on purpose, as the benchmark to beat, not as
the intended final answer. It is believed to be over-conservative, and on a
GC-rich polymorphic target it can reject every candidate (it returns nothing for
the HLA-A*01:01 CDS in `data/examples`). The full critique, the supporting
measurements, the proposed expected-off-target-load replacement, and the
experiment that would justify building it are recorded in
[`seed_offtarget_design_note.md`](seed_offtarget_design_note.md). **Read that
before changing any constant in `core.seed`.**

Out of scope: which strand is preferentially loaded (that is `core.asymmetry`),
and the near-perfect-match homology search that forms siDirect's third selection
step.

### References

- Primary: Naito, Y. et al. (2009). siDirect 2.0. *BMC Bioinformatics* 10, 392.
  doi:10.1186/1471-2105-10-392
- Ui-Tei, K. et al. (2008). Thermodynamic stability and Watson-Crick base pairing
  in the seed duplex... *Nucleic Acids Research* 36(22), 7100-7109.
  doi:10.1093/nar/gkn902
- Freier, S.M. et al. (1986). *PNAS* 83(24), 9373-9377.

---

## Duplex-end asymmetry (`core.asymmetry`)

### What it measures

An siRNA duplex has two ends and RISC only keeps one strand. Which one it keeps
is decided thermodynamically: the helicase unwinds from whichever end is easier
to open, and the strand whose 5' end sits at that looser end is retained as the
guide (Khvorova et al. 2003; Schwarz et al. 2003). The 2004 rule sets all encode
this indirectly - Ui-Tei through A/U counts at the guide 5' terminus, Amarzguioui
through the terminal-triplet A/U differential. This module measures it directly,
in kcal/mol:

    differential = dG(guide 5' end window) - dG(sense 5' end window)

A **positive** differential means the guide's 5' end is the less stable end (dG
is less negative there), the configuration that favours loading the guide - the
orientation you want. A negative differential warns that the passenger strand is
the thermodynamically preferred guide, which both wastes the duplex and puts the
passenger's seed in charge of off-target behaviour.

### Why `stacking_terms()`, not `duplex_thermo()`

These windows are the ends of a longer helix, not free duplexes. Their outer edge
is a real helix end but their inner edge is not, and adding a duplex initiation
term plus a terminal A-U penalty to both windows would not cancel cleanly.
`core.thermo` draws that distinction explicitly; this module is the reason it
does.

### Scope and deferrals

- The window length is a parameter (default 4 bp), not a published constant.
  Khvorova computed internal-stability profiles over 5-nt windows; Schwarz
  emphasised the few terminal base pairs; later tools differ again. No single
  canonical window exists, so none is asserted.
- **No threshold is defined.** The module returns a measured quantity and its
  sign, not a pass/fail verdict. Turning kcal/mol into a cut-off needs
  calibration against efficacy data HUSH does not have, and the sign alone
  already carries the mechanistic claim.
- Parameters are Xia 1998, 1 M NaCl, uncorrected.

### Validation

On the real HLA CDS the sign of the thermodynamic differential agrees with
Amarzguioui's A/U-count proxy in 95.2% of the candidates where the proxy is
decisive - the expected agreement between a quantity and its published stand-in,
and a check that neither sign convention is inverted. Enforced in
`tests/test_asymmetry.py`.

### References

- Primary: Khvorova, A., Reynolds, A., Jayasena, S.D. (2003). Functional siRNAs
  and miRNAs exhibit strand bias. *Cell* 115(2), 209-216.
- Schwarz, D.S. et al. (2003). Asymmetry in the assembly of the RNAi enzyme
  complex. *Cell* 115(2), 199-208.

---

## Intramolecular hairpin Tm (`core.hairpin`)

### What it supplies

Reynolds et al. (2004) award a point when an siRNA has no stable internal
structure - operationally, when the most stable hairpin it can fold into melts at
Tm <= 20 C (see the rules section, "absence of internal repeat"). This module
supplies that Tm.

### A hairpin is not a free duplex

Its melting is not the duplex melting `core.thermo` and `core.seed` compute:

- It melts **intramolecularly**. Groebe & Uhlenbeck showed the transition is
  concentration-independent (identical profiles from 1 to 10 uM), so the Tm has
  no concentration term: `Tm = dH / dS`.
- It carries a **loop penalty**. The closing loop destabilises the fold by an
  amount that depends mostly on loop size (Groebe & Uhlenbeck, Table 3), which the
  stem's base-pair stacking does not capture.

### The model (a transparent approximation)

This is deliberately a transparent approximation, not a reproduction of whatever
folding tool Reynolds used (which is not specified in their paper and cannot be
reproduced exactly). It is built from parts HUSH already trusts, with every
assumption named:

- **Stem:** perfect Watson-Crick stacking only, from `stacking_terms()` with the
  Xia (1998) parameters. No initiation and no terminal-A-U term - the loop penalty
  plays initiation's role, matching how Groebe & Uhlenbeck derived the loop data
  (they subtracted the bare stem stacking). G-U wobble and mismatched/bulged stems
  are not modelled, so internal structure may be slightly *under*-estimated - the
  conservative direction for a filter that rewards its absence.
- **Loop:** destabilising dG37 from Groebe & Uhlenbeck Table 3, weighted by the
  loop's actual A/C/U composition and linearly interpolated across loop size.
- **Tm:** because Table 3 gives the loop as dG37 *only* (not dH and dS
  separately), the loop destabilisation is taken to be purely entropic
  (`dH_loop = 0`, so `dS_loop = -dG37_loop / T37`). This is the one unavoidable
  modelling assumption, and it is the standard one: loops raise the conformational
  entropy cost of folding without themselves forming or breaking base pairs.

The "internal structure" is taken to be the single most stable (minimum-dG37)
simple hairpin; its Tm is what the criterion compares against 20 C. Loop and stem
minimums (3 nt, 3 bp) are steric conventions, not from the source, and are
tunable.

### References

- Primary: Groebe, D.R. & Uhlenbeck, O.C. (1988). Characterization of RNA hairpin
  loop stability. *Nucleic Acids Research* 16(24), 11725-11735.
- Reynolds, A. et al. (2004). Rational siRNA design for RNA interference.
  *Nature Biotechnology* 22, 326-330.

---

## Classical rule sets (`core.rules`)

Three empirical rule sets from 2004 form the backbone of rational siRNA design.
They all encode thermodynamic asymmetry - the duplex should unwind more easily at
the guide (antisense) 5' end so RISC loads the guide rather than the passenger -
but express it differently. Strand conventions match `core.sirna`: sense/passenger
and antisense/guide are both stored 5'->3'; each rule is evaluated on whichever
strand it is naturally defined for.

### Ui-Tei (2004)

A strict 4-rule pass/fail filter (A/U at the guide 5' terminus; G/C at the sense
5' terminus; >= 4 A/U in the guide's 5'-terminal 7 nt; no G/C stretch longer than
9 nt). "Highly effective" requires all four.

### Reynolds (2004)

Scored **out of 10**, following the original algorithm rather than the common
"8 equal binary criteria" simplification. All parts are evaluated on the sense
strand, 5'->3':

    +1 each (binary):  GC content 30-52%; absence of internal repeat;
                       'A' at position 19; 'A' at position 3; 'U' at position 10.
    +1 per A/U base at positions 15-19  (additive, 0 to +5 - NOT a threshold).
    -1 each (penalty): 'G' at position 13; 'G' or 'C' at position 19.

so the score runs from -2 to +10, and >= 6 predicts a functional siRNA.

Position 19 (the 3' terminus) carries three of these at once - the additive 15-19
count, the 'A' at 19 bonus, and the G/C at 19 penalty - giving a clean
A(+2) / U(+1) / G,C(-1) gradient. The additive 3'-terminal A/U term and the two
penalties are the two places the common simplification diverges from the original;
the out-of-10 form is validated by construction and internal consistency (there is
no published per-siRNA Reynolds score table to reproduce).

The "absence of internal repeat" criterion is a hairpin-Tm calculation supplied by
`core.hairpin` (Tm <= 20 C); it can be overridden with an external folder's
verdict via `reynolds(candidate, no_internal_repeat=...)`.

Longer cores (20-24 nt): the original criteria are defined for a 19-mer. The
3'-anchored parts (the 15-19 window, position 19) follow the actual 3' terminus,
while the internal positional parts (positions 3, 10, 13) stay at their absolute
1-based positions from the 5' end.

### Amarzguioui & Prydz (2004)

A simple additive, criterion-based score (like Reynolds), not a weighted
regression. Positions are 1-based on the sense strand, 5'->3':

    Positive motifs (+1): S1 = G/C at 1;  A6 = A at 6;  W19 = A/U at 19.
    Negative motifs (-1): U1 = U at 1;    G19 = G at 19.

The duplex-end asymmetry enters as the **raw differential itself**, not a lookup:
the paper incorporates it "in the form of the actual value of the A/U differential
(ranging from -3 to +3)", so `score = motif points + differential`. The
differential (the paper's "DT3") is counted over the three terminal base pairs at
each end of the duplex, sense strand as reference, and measured **3' end minus 5'
end**. The sign is easy to get backwards: functional siRNAs score positive,
because an A/U-rich sense 3' end is the same physical end as the guide's 5' end, so
unwinding starts there and RISC retains the guide. A score of >= 3 predicts a
functional siRNA.

**Validation:** all 46 siRNAs of the paper's own Table 3 reproduce exactly - every
DT3, every motif, and every total score - along with the reported class means (2.6
functional / 0 non-functional) and selection statistic (14 selected, 11
functional). The fixture is `tests/data/amarzguioui_table3.tsv`; the check is in
`tests/test_rules.py`.

### References

- Primary (per rule): Ui-Tei, K. et al. (2004). *Nucleic Acids Research* 32(3),
  936-948. | Reynolds, A. et al. (2004). *Nature Biotechnology* 22, 326-330. |
  Amarzguioui, M. & Prydz, H. (2004). *Biochem Biophys Res Commun* 316(4),
  1050-1058. doi:10.1016/j.bbrc.2004.02.157

---

## Candidate selection (`core.selection`)

### Why a filter, not a ranking

This is where a user's design choices are applied to the enumerated candidates.
It is deliberately **not** a ranking: there is no single combined score, because
the signals (three functionality rule sets, a GC window, a seed-duplex off-target
Tm, a loading-asymmetry sign) live on incommensurable scales and no defensible way
exists to weight one against another without efficacy data HUSH does not have.
Inventing weights would be a black box that needs tuning and validation; a
transparent, option-rich filter - in the spirit of siDirect - does not. Every
candidate passes or does not, and keeps its full per-criterion breakdown so the
reason for any verdict is visible.

### The two groups of criteria

- **The three rule sets** (Ui-Tei, Reynolds, Amarzguioui). Each toggles
  independently; the enabled ones combine by a mode: `ALL` (must satisfy every
  enabled set) or `ANY` (kept if any enabled set passes it - the union of what the
  selected algorithms would each select). With no rule set enabled, the rule step
  imposes no constraint.
- **Independent filters**, each optional, ANDed on top: a GC-content window, the
  seed-duplex off-target Tm, a standalone internal-structure (self-hairpin) gate,
  and the loading asymmetry.

### Customisable values are options, not rule edits

The GC range, seed-Tm threshold, hairpin-Tm cut-off, and rule score thresholds are
options here, not edits to the rules themselves. Two of them shadow a criterion
that lives fixed inside a published rule, and the split is deliberate:

- the standalone GC filter is separate from Reynolds' internal GC criterion, which
  stays fixed at the published 30-52%;
- the standalone hairpin gate is separate from Reynolds' internal-repeat
  sub-criterion, which stays pinned to Xia 1998 at Tm <= 20 C.

In both cases the rule keeps its published form (fidelity to the source), and the
standalone filter is the tunable one. `core.hairpin` and `core.asymmetry` both now
take a `params` argument so the standalone gates can be scored with either
nearest-neighbor set, while `reynolds()` keeps calling them at their Xia defaults.

### Advanced thermodynamic knobs, and where the parameter set may not go

Three options are "advanced": the asymmetry window (default 4 bp; an empirical
range, not a single canonical value) and temperature, and the nearest-neighbor
parameter set (`ThermoParams`, Xia 1998 / Freier 1986). The parameter set is
scoped to the asymmetry and hairpin readouts **only**. It is deliberately not
wired into the seed filter: siDirect's 21.5 C cut-off is calibrated on Freier
numbers and is meaningful only in combination with the seed module's own pinned
internals (see the seed section), so exposing a parameter-set switch there would
let a caller silently invalidate the threshold. Same reasoning as the pinned rule
criteria above - a published threshold travels with the parameter set it was fit
on.

### Ordering

Nothing is ranked or reordered: selected candidates are returned in target-position
order, exactly as enumerated. Ordering by any single criterion is left to the
caller.
