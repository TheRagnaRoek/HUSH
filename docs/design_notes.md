# HUSH design notes

Design rationale, literature provenance, formulae, and validation records for
the stage-2 siRNA design engine and the stage-3 MSA/allele stage. This is the
long-form companion to the module docstrings: the modules say *what* they do and
the critical caller contract; the *why*, the derivations, and the "checked
against the paper's own figures" records live here so the code stays readable.

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
- [Transcriptome homology off-target filter (`core.offtarget`)](#transcriptome-homology-off-target-filter-coreofftarget)
- [Duplex-end asymmetry (`core.asymmetry`)](#duplex-end-asymmetry-coreasymmetry)
- [Intramolecular hairpin Tm (`core.hairpin`)](#intramolecular-hairpin-tm-corehairpin)
- [Classical rule sets (`core.rules`)](#classical-rule-sets-corerules)
- [Candidate selection (`core.selection`)](#candidate-selection-coreselection)
- [Duplex geometry: core length and 3' overhangs (`core.sirna`)](#duplex-geometry-core-length-and-3-overhangs-coresirna)
- [MSA / allele stage (`msa`)](#msa--allele-stage-msa)

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

- Primary: Xia, T., SantaLucia, J., Burkard, M.E., Kierzek, R., Schroeder, S.J., Jiao, X., Cox, C., Turner, D.H. (1998). Thermodynamic parameters for an expanded
  nearest-neighbor model for formation of RNA duplexes with Watson-Crick base
  pairs. *Biochemistry* 37(42), 14719-14735. doi:10.1021/bi9809425
- Freier, S.M., Kierzek, R., Jaeger, J.A., Sugimoto, N., Caruthers, M.H., Neilson, T., Turner, D.H. (1986). Improved free-energy parameters for predictions of
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
- Ui-Tei, K., Naito, Y., Nishi, K., Juni, A., Saigo, K. (2008). Thermodynamic stability and Watson-Crick base pairing
  in the seed duplex... *Nucleic Acids Research* 36(22), 7100-7109.
  doi:10.1093/nar/gkn902
- Freier, S.M., Kierzek, R., Jaeger, J.A., Sugimoto, N., Caruthers, M.H., Neilson, T., Turner, D.H. (1986). *PNAS* 83(24), 9373-9377.

---

## Transcriptome homology off-target filter (`core.offtarget`)

### What it is, and how it differs from the seed filter

siDirect's third selection step is a *homology* search: it looks each guide up
against the human transcriptome and drops candidates that have a near-perfect
match somewhere other than the intended target. This is a different question from
`core.seed`. The seed filter asks a sequence-intrinsic question ("does this seed
pair too stably?") and never looks at any transcript; this module asks a database
question ("where else in the transcriptome could this guide bind, and how well?").
The two are orthogonal and compose: a candidate can be seed-clean yet have a
near-perfect full-length off-target, or vice versa.

The motivating gap: benchmarked on the PDCD1 CDS with matched functional rules,
HUSH kept a **superset** of siDirect's picks. The extras were exactly the
candidates siDirect had removed for near-perfect off-targets - decisively, a
candidate at CDS 650-668 that clears every functional rule and composition gate,
so no rule difference could explain the disagreement. Only a transcriptome search
closes it.

### The search: bowtie1 `-v`, the open equivalent of GGGenome

siDirect relies on GGGenome, a fast k-mismatch full-text search that is not open.
The open, offline counterpart is **bowtie1 in `-v` mode**: end-to-end alignment
allowing at most N mismatches, reporting every hit. It is *not* BLAST - BLAST is a
local heuristic that can miss short near-matches and reports HSPs rather than a
clean full-guide mismatch count, so it cannot reproduce siDirect's "minimum number
of mismatches" column reliably. bowtie1 (Artistic License) is also lighter than
bowtie2 (GPLv3) and, unlike bowtie2, keeps the fixed-mismatch `-v` mode (capped at
3) this needs. The reference is **RefSeq human mRNA**, pinned per build; its NM_/NR_
accessions match HUSH's benchmark fixtures, so on-target exclusion by accession is
exact. Build the index with `tools/build_offtarget_db.py`, which records a manifest
(release, source URLs, checksums, bowtie version) for reproducibility.

### Two backends behind one protocol

An `OffTargetSearcher` protocol has two real implementations. `BruteForceSearcher`
is a pure-Python Hamming scan - correct by construction, no external tool - for
small references (a curated allele panel, and the tests); it is the oracle the
bowtie backend is checked against. `BowtieSearcher` shells out to bowtie (via
stdlib `subprocess`, so importing the module never needs the binary) and, through
`precompute()`, aligns every strand of every candidate in a **single** bowtie run,
returning a lookup-backed `PrecomputedSearcher` so the per-candidate selection loop
pays no process cost. Both engines are ungapped, Hamming-distance and
forward-orientation, so their results agree exactly.

### Strand semantics and on-target exclusion

A guide silences transcripts complementary to it, i.e. those carrying its target
site (the **sense** core); the passenger, if mis-loaded, silences those carrying
the **antisense** core. So each candidate contributes two reads - the guide side
searches the sense core, the passenger side the antisense core - forward strand
only (`--norc`), because transcripts are single-stranded sense mRNA. The
candidate's own source transcript is excluded by accession (version aside),
matching siDirect; it is auto-detected from the target FASTA header and
overridable.

### The gate, and its calibration

The gate rejects a candidate whose nearest non-self off-target (on either strand)
has **fewer than** `offtarget_min_mismatch` mismatches. The default is **2**:
on the PDCD1 benchmark siDirect keeps candidates whose nearest off-target is >= 2
mismatches away and drops only the near-perfect (0/1) ones, so a threshold of 2
reproduces its decision - a threshold of 3 would wrongly drop siDirect's own
2-mismatch picks. The search depth defaults to 3, deeper than the gate needs, so
the reported nearest-mismatch and per-level counts line up with siDirect's own
columns (which reach 3). The gate is **flat**, not position-weighted: seed-region
weighting is already `core.seed`'s job, so the two off-target angles stay separate.

### Relevance to allele-aware design

For a polymorphic family like HLA the generic framing is only partly right:
cross-reactivity *within* the intended allele family may be the design goal
(pan-allele knockdown), not a hazard, so a single global gate that treats every
non-source hit as off-target cannot express it. v1 excludes only the source
accession; a later `msa/` integration is where "intended vs unintended
cross-reactivity" belongs (the same distinction the seed design note raises). The
`BruteForceSearcher` over a small curated panel is the natural tool for that
targeted question.

### Validation

Unit tests plant off-targets at known mismatch levels in a tiny synthetic
reference and check the mismatch counts, on-target exclusion (version-insensitive),
both-strand detection, the threshold boundary, and the bowtie output parser; a
skipped-if-absent parity test asserts `BowtieSearcher` and `BruteForceSearcher`
agree on a freshly built index. Enforced in `tests/test_offtarget.py`. The
transcriptome-scale acceptance check (HUSH-only candidates drop to zero against
siDirect on PDCD1) is run with a built index via `tools/compare_sidirect.py`.

### Scope and deferrals

Genomic (non-transcript) off-targets are out of scope, as is any expression-weighted
off-target model (the proposal in `seed_offtarget_design_note.md` is a separate,
later research task). The full decision record - bowtie1 vs BLAST, RefSeq pinning,
strand semantics, the HLA caveat - is in
[`offtarget_design_note.md`](offtarget_design_note.md).

### References

- Primary: Naito, Y. et al. (2009). siDirect 2.0. *BMC Bioinformatics* 10, 392.
  doi:10.1186/1471-2105-10-392 (the transcriptome off-target selection step).
- Langmead, B., Trapnell, C., Pop, M., Salzberg, S.L. (2009). Ultrafast and
  memory-efficient alignment of short DNA sequences to the human genome. *Genome
  Biology* 10, R25. doi:10.1186/gb-2009-10-3-r25 (bowtie).
- RefSeq: O'Leary, N.A. et al. (2016). Reference sequence (RefSeq) database at
  NCBI. *Nucleic Acids Research* 44(D1), D733-D745. doi:10.1093/nar/gkv1189

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
- Schwarz, D.S., Hutvagner, G., Du, T., Xu, Z., Aronin, N., Zamore, P.D. (2003). Asymmetry in the assembly of the RNAi enzyme
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
- Reynolds, A., Leake, D., Boese, Q., Scaringe, S., Marshall, W.S., Khvorova, A. (2004). Rational siRNA design for RNA interference.
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

- Primary (per rule): Ui-Tei, K., Naito, Y., Takahashi, F., Haraguchi, T., Ohki-Hamazaki, H., Juni, A., Ueda, R., Saigo, K. (2004). *Nucleic Acids Research* 32(3),
  936-948. | Reynolds, A., Leake, D., Boese, Q., Scaringe, S., Marshall, W.S., Khvorova, A. (2004). *Nature Biotechnology* 22, 326-330. |
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
  seed-duplex off-target Tm, the transcriptome homology off-target filter
  (`core.offtarget`, which needs a reference database handle), a standalone
  internal-structure (self-hairpin) gate, a homopolymer-run gate (reject a
  single-nucleotide run over `max_homopolymer_run`, default 3, e.g. `GGGG` - judged
  on the core only), and the loading asymmetry.

### The default paradigm

`SelectionCriteria`'s defaults encode a best-practice screen rather than the most
permissive or the most conservative one. All three rule sets are enabled and
combined by **`ANY`** (each classical algorithm is treated as an equal vote, so a
candidate any one of them endorses survives - `ALL` would demand three
independent 2004 heuristics agree, which erases almost everything). The **seed**
off-target Tm gate is on because it is the siDirect industry standard, and the
**homopolymer** gate is on because poly-runs are disruptive to both siRNA
stability and polymerase read-through - both are cheap, widely-agreed quality
gates. The GC, hairpin and asymmetry gates default off: GC is counter-productively
strict on the GC-rich polymorphic families HUSH targets, and asymmetry is a niche,
conservative signal not found in comparable tools. The transcriptome off-target
filter defaults **on at the CLI** (`hush design`) because off-targets can be
harmful, but it needs a reference index; with none configured the CLI skips it
with a warning rather than erroring, and `SelectionCriteria` itself keeps it off
because the filter has a mandatory searcher argument with no sensible default. The
one deliberate cost of this paradigm: the seed filter is strict enough to reject
every candidate on a very GC-rich CDS (HLA), where the user is expected to relax it
per-target with `--no-seed` (the GUI will surface this hint).

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

`core.selection` ranks nothing: it is a pure filter, and the candidates it selects
are returned in target-position order, exactly as enumerated. The one ordering the
pipeline does impose lives a stage later, in the allele module: a request that
spares alleles is ranked by discrimination (see the MSA/allele section and
`msa.discrimination`), which reorders survivors but never admits or rejects one.
Ordering by any single siRNA criterion is still left to the caller.

---

## Duplex geometry: core length and 3' overhangs (`core.sirna`)

### The standard: a 19-nt core with symmetric 2-nt 3' overhangs

HUSH's canonical duplex is a **19-nt fully paired core carrying a 2-nt 3' overhang on
each strand** - 21-nt synthesized strands, the Elbashir/Tuschl design. Two numbers get
conflated in the literature and it is worth separating them: **19** is the paired core
(the RISC-relevant, thermodynamically active region), **21** is the strand a vendor
synthesizes (core + overhang), and **23** is a different measurement again (an mRNA
target-site footprint, or a Dicer-substrate design). `--length` sets the **core** and
defaults to 19.

Keeping the core at 19 is not just convention: the three rule sets (Ui-Tei, Reynolds,
Amarzguioui) are all *defined and validated on the 19-nt duplex*, with position-specific
criteria indexed against it. Enlarging the core would feed those scorers inputs their
authors never specified. So the overhang is modelled as a separate 2-nt feature, tracked
apart from the core, and **no filter (rules, GC, seed, hairpin, homopolymer) ever sees
it** - the seed keys off guide positions 2-8 from the 5' end, which the 3' overhang does
not touch.

### Overhang policy

`OverhangPolicy` offers three choices (`--overhang`):

- **`natural`** (default): the overhang is the target's own flanking sequence. The
  passenger's 3' overhang is the 2 nt immediately 3' of the core on the mRNA; the guide's
  3' overhang is the reverse-complement of the 2 nt immediately 5' of the core. This is
  the true *offset* Tuschl duplex - the two strands share the 19-nt core but their
  overhangs come from opposite flanks, so they are not a blunt reverse-complement pair.
- **`dtdt`**: two 2'-deoxythymidines on both strands (the common synthetic overhang),
  rendered lowercase `tt` to flag the DNA residues in an otherwise-RNA strand.
- **`uu`**: two ribo-uridines on both strands.

A candidate whose natural flank runs past a transcript end (the first/last couple of
sliding-window positions) has no flank to draw from and falls back to `dtdt` there.
Downstream, the MSA/allele module restricts candidates to conserved interior windows, so
this edge fallback is a rare safety net rather than a common path.

### Consequence for siDirect benchmarking

In `natural` mode the guide strand as ordered (`guide_full` = core + upstream-derived
overhang) equals the reverse-complement of the full core-plus-upstream window - i.e.
exactly siDirect's listed guide. So adding the correct offset geometry leaves the
functionally decisive strand reproducing siDirect. The passenger, by contrast, is offset
(its overhang is the downstream flank) and no longer matches siDirect's blunt-listed
passenger, which is the intended physical duplex, not a defect. Benchmarks should compare
on the guide or the target position.

### References

- Elbashir SM, Harborth J, Lendeckel W, Yalcin A, Weber K, Tuschl T. Duplexes of
  21-nucleotide RNAs mediate RNA interference in cultured mammalian cells. *Nature* 411,
  494-498 (2001). doi:10.1038/35078107

---

## MSA / allele stage (`msa`)

The stage that makes HUSH allele-aware. Given the target and a panel of related
sequences (the motivating case is an HLA allele family), it decides, per target
position and per candidate window, whether an siRNA built from the target would
hit the alleles the user wants silenced and discriminate against the ones they
want spared. The per-unit rationale lives in the module docstrings
(`msa.alignment`, `msa.align_star`, `msa.variation`, `msa.windows`,
`msa.targeting`); this section records the design decisions that span the stage
and the choice made when wiring it into `hush design`.

### The coordinate contract, and why it holds the stage together

Every panel member is projected onto **ungapped target coordinates**: it becomes
a string exactly as long as the target, index *i* carrying whatever residue the
member has at target position *i*, with columns the alignment inserted into the
target removed (kept as insertion payloads, not dropped). This is the single
decision the rest of the stage rests on. Because the frame is the target's own,
window position *p* on the target is the same *p* that `core.sirna` enumerates a
candidate at - the two lists join on position without either side knowing about
the other. No unit ever converts between alignment columns and target positions,
because after projection there are no alignment columns left to convert.

`PanelAlignment` enforces the invariant (every row exactly the target's width)
rather than trusting it, so a backend that violates it fails loudly at the
hand-off instead of silently shifting which positions look polymorphic.

### Three states, kept distinct: match, mismatch, uncovered

`msa.variation` reduces each member to one of three states per position, and the
distinction between *mismatch* and *uncovered* is load-bearing, not cosmetic.
"Uncovered" (a deletion, or an alignment that does not reach) is not a weak form
of "different": a window in which an allele the user wants to **hit** has no
residues says nothing about that allele, and one in which an allele to be
**spared** is missing cannot be claimed to discriminate against it. Folding the
two together would turn missing data into evidence - a silent error, not a design
choice - so coverage is required by both gates rather than treated as harmless.

### Two gates, and why a filter and not a ranking

`msa.windows` applies an `AlleleTargeting` request under two gates, in the order a
rejection is reported: the targeted alleles must be covered and identical across
the window (the guide is built from the target, so a difference there is a
mismatch against the very allele it should silence), and every spared allele must
be covered and carry at least `min_spare_mismatches` differences. The four named
intents - pan-allele, allele-specific, target-specific, subset - are not four
types but one structure (a *hit* set and a disjoint *spare* set) with different
values, so a single code path serves all of them and `spare` is stored explicitly
rather than inferred as "whatever else is in the panel" (a panel is data that can
change; an implicit spare set would silently mean something different after a
member is added). Target-specific is the empty-*hit* case: candidates are built
from the loaded target, so it is silenced by construction, and hitting no panel
member while sparing the whole panel is a coherent, meaningful request (the mirror
of pan-allele, which hits everything and spares nothing).

The **gate** ranks nothing - every window passes or does not, and keeps its full
per-allele breakdown. The **ordering** of survivors is a separate layer:
`msa.discrimination` places each difference against published measurements of
which mismatches, at which guide positions, block silencing of the spared allele
(Huang et al. 2009 for the exact base change, Schwarz et al. 2006 for region and
position, Ohnishi et al. 2008 for a second mismatch). It is a lexicographic tuple
of published facts, not a weighted score - no calibration data is invented - and a
window is ranked by its **weakest** spared allele, the "weakest link" rule
`min_spare_mismatches` already follows. A guide that sits *3' of the seed*
discriminates robustly while a seed mismatch is overwhelmed by concentration, so
where a difference sits still matters enormously; the positions are recorded in
target-window coordinates (`guide_position()` converts to the guide strand, which
runs the opposite way). Insertions are recorded but not yet scored: an insertion
would make the guide bulge on that member, a strong discriminating feature, but
scoring it needs parameters not yet carried.

### The star backend, and its scope

`msa.align_star` aligns every member pairwise against the target, the target as
the fixed centre of the star (Gusfield's center-star construction). The classical
method picks the centre to minimise summed pairwise cost; here the centre is
fixed by the pipeline - candidates are enumerated from the target - so the
construction's factor-2 bound does not apply and no optimality is claimed. What it
buys is a consistent frame, not a quality guarantee. It is adequate while members
are near-identical (an HLA locus runs ~98-99% identity) and degrades as they
diverge, because one centre cannot represent the relationships a guide tree would.
The backend is named through `AlignmentBackend`, so a progressive aligner can be
slotted in behind it without touching `variation`, `windows`, or the CLI. Gap
costs default to making an indel buy back at least four mismatches before it pays
for itself, which keeps the aligner from scattering frame-shifting gaps to chase a
few mismatches - the failure mode that would move which positions look
polymorphic; it is an approximation (no codon awareness), not a constraint.

### Wiring into `hush design`: the allele gate lives at the app layer

The allele stage is engaged from `hush design --panel`, not as a separate command,
so one run does allele-aware design end to end. The integration decision worth
recording is **where the gate lives**: the dependency direction is `msa -> core`,
so `core.selection` must never import `msa`. The CLI (`app.cli`) is the composition
root that may import both, and it is there that the allele verdict is ANDed with
the siRNA selection verdict - a candidate is selected only if `CandidateEvaluation.
selected` **and** the window at its position `satisfies` the targeting. The window
length is bound to the siRNA core length (`--length`) so the two lists stay
one-for-one; the candidate is looked up by its start position in a dict built once
from `evaluate_all`. The criteria block and pass-each-gate funnel gain an
`allele:` line; a discriminating request adds two per-candidate columns describing
its weakest spared allele - `minDiff` (that allele's difference count) and `disc`
(its best mismatch's full criteria: guide position, region, Huang listing,
position tier, class) - and is ranked by discrimination by default (`--sort auto`;
`--sort position` forces target order). `--explain`, or `assess <region>` in the
shell, prints the full per-allele breakdown for every spared allele, not just the
weakest; the TSV gains `allele_hit_max_mismatch`, `allele_spare_min_mismatch`,
`allele_spare_count`, and the deciding mismatch's `disc_*` criteria. Keeping the
gate at the app layer means the allele stage and the siRNA engine each stay
testable in isolation, and the funnel reports the allele gate as just one more
independent gate rather than a special case.

### Scope and deferrals

- Residues are compared as characters; IUPAC ambiguity codes are not resolved (an
  IPD-IMGT/HLA panel is unambiguous, a hand-assembled one may not be).
- Differences and insertions are counted, not weighted (see above) - the position
  data needed for weighting is recorded against the day `core.thermo` carries
  mismatch thermodynamics.
- Only a star backend exists; divergent families want a progressive one.
- Elbashir SM, Martinez J, Patkaniowska A, Lendeckel W, Tuschl T. Functional anatomy of
  siRNAs for mediating efficient RNAi in Drosophila melanogaster embryo lysate. *EMBO J*
  20, 6877-6888 (2001). doi:10.1093/emboj/20.23.6877 - establishes the 2-nt 3' overhang
  as most efficient and that 2'-deoxy (dTdT) is tolerated at the 3' end.
