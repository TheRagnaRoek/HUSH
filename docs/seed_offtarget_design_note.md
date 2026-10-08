# Design note: the seed-duplex Tm filter, and what should eventually replace it

Status: **decision recorded, alternative deferred.** `core/seed.py` implements
the siDirect 2.0 filter faithfully and stays that way for now. This note records
why the current approach is believed to be flawed, what the alternative is, and
what would have to happen before the alternative could be built. It exists so
that work can be picked up later without re-deriving the argument.

Date of decision: 2026-07-24.

---

## 1. What was decided, and why

`core/seed.py` reproduces the siDirect 2.0 seed-duplex Tm filter exactly:
Freier et al. (1986) nearest-neighbor parameters, helix initiation constant
-10.8, total strand concentration 100 uM, `[Na+]` 100 mM via a
`16.6*log10[Na]` term, gas constant 1.987, no terminal A-U penalty and no
symmetry correction, threshold Tm < 21.5 C on **both** strands.

It is verified against the figures published in the siDirect 2.0 paper: of all
4^7 = 16384 possible seed heptamers, 4488 (27.4%) fall below 21.5 C, spanning
-12.0 C to 59.6 C ("ranged from -12 C to 60 C"). See `tests/test_seed.py`.

This is deliberately kept unchanged, for three reasons:

1. **It is the baseline to beat.** Any claim that HUSH improves on siDirect
   requires running siDirect's algorithm exactly, not approximately.
2. **It labels existing experimental data.** siRNAs already designed and tested
   with siDirect carry siDirect's verdict as metadata; that verdict is only
   usable as a label if it can be reproduced.
3. **Reproducibility.** A published comparison will be asked to show the
   incumbent was implemented correctly.

## 2. A rejected alternative: swapping in newer parameters

The obvious "modernisation" is to replace Freier 1986 with Xia 1998 (the set
`core/thermo.py` defaults to). This was evaluated and **rejected as a no-op.**

Measured across all 16384 heptamers:

| Quantity | Freier 1986 | Xia 1998 |
|---|---|---|
| Range | -12.0 to 59.6 C | -31.3 to 61.1 C |
| Threshold selecting the same 27.4% | 21.5 C | 22.31 C |

- Pearson r = 0.988, Spearman rho = 0.990 between the two Tm scales.
- Under quantile-matched thresholds only **2.3%** of heptamers change verdict.

Conclusions:

- Swapping the parameter set does **not** address the over-conservatism
  described in section 3. It re-ranks almost nothing.
- Swapping the table while keeping the 21.5 C threshold is strictly worse than
  either consistent choice, because 21.5 C is an empirical cut-off calibrated on
  Freier numbers. This is why `core/thermo.py` exposes both tables by explicit
  name and refuses to have a single "current" default.

The reason the two sets agree so closely is given in section 3.1: they share the
same systematic error.

## 3. Why the current filter is believed to be flawed

### 3.1 A free-duplex model is the wrong physics, and both parameter sets share the error

Freier and Xia both measure **free RNA duplexes in solution**. The seed duplex
does not form in solution - it forms inside Argonaute, which displays guide
nucleotides g2-g8 in a pre-organized near-helical conformation that removes much
of the entropic cost of pairing. Seed pairing proceeds 50-200x faster than
annealing of the equivalent naked RNAs, and guide-loaded Ago behaves more like
an RNA-binding protein than like free RNA.

Both tables are therefore proxies carrying the same directional error, which is
exactly why substituting one for the other moves only 2.3% of heptamers.

### 3.2 The filter controls affinity per site but ignores site count, and the two are anticorrelated

The siDirect 2.0 paper reports (its Figure 4) that seeds with **lower** Tm are
**more** frequently present in human 3' UTRs than seeds with higher Tm. The
filter therefore selects for seeds that bind each site weakly while
systematically enriching for seeds that have more sites to bind. The paper
surfaces this observation but does not resolve it.

The quantity that plausibly determines the actual off-target burden is:

```
expected off-target load ~ (seed-matched sites in EXPRESSED transcripts)
                           x (repression per site)
```

siDirect controls only the second factor. This predicts a specific failure mode:
**a high-Tm seed that is rare in the expressed transcriptome can be genuinely
clean**, and would be wrongly rejected. Anecdotal experimental evidence for
exactly this failure mode motivated this note (see section 5).

### 3.3 A hard conjunctive gate can return the empty set

The filter is a binary gate ANDed with the Ui-Tei rule filter. On a GC-rich,
polymorphic target this can eliminate every candidate. Measured on
`data/examples/hla_a01.fasta` (HLA-A*01:01 CDS, 1098 nt, 63.8% GC, 1080
candidates):

| Step | Surviving candidates |
|---|---|
| Ui-Tei functional (siDirect step 1) | 85 |
| ... of which seed-clean on both strands | **0** |
| Seed-clean on both strands (any Ui-Tei verdict) | 2 |

The mechanism of the interaction is specific and worth recording:

| | median guide seed Tm | median passenger seed Tm |
|---|---|---|
| All candidates | 37.9 C | 37.9 C |
| Ui-Tei passers | 28.4 C | 38.2 C |

Ui-Tei rule (iii) requires A/U richness in the guide's 5'-terminal 7 nt, which
overlaps the seed, so it does pull guide seed Tm down (16 of 85 clear 21.5 C).
But Ui-Tei rule (ii) requires G/C at the sense 5' terminus, which on a 64% GC
gene keeps the passenger seed GC-rich: **0 of 85** clear it. The two filters are
in direct tension on GC-rich targets.

This is consistent with siDirect's own caveat that it cannot design for 5.3% of
RefSeq mRNAs under default parameters, with multigene families named as the
typical example. HLA is such a family. For HUSH's core use case a gate that can
return nothing is the wrong construct; a ranked score is needed so that a best
available candidate always exists, with the accepted risk stated explicitly.

## 4. The proposed alternative

Replace the binary gate with a **continuous, calibrated estimate of expected
off-target load**, combining the two factors siDirect separates:

1. **Per-site repression strength** - a seed-duplex stability term. Use Xia 1998
   here, not because it is newer but because a new model requires its own
   calibration anyway, so there is no reason to inherit the older parameters.
   Consider whether an Ago-context correction (section 3.1) can be fitted rather
   than assumed.
2. **Site abundance** - the number of seed-complementary sites in the 3' UTRs of
   transcripts **expressed in the relevant cell type**. Expression-weighting is
   the part existing tools omit; DSIR and Dharmacon siDESIGN compute raw seed
   frequency, siDirect computes Tm, and nothing combines them.

Output a score and a rank, not a pass/fail. Retain the siDirect gate as a
selectable mode for benchmarking.

Novelty position: the combination and its calibration, not a new mechanism.
This should be stated honestly - it is a reframing and recalibration advance.

### Relevance to allele-aware design

For HLA the generic off-target framing is partly inapplicable: cross-reactivity
within the allele family may be the design goal (pan-allele targeting) rather
than a hazard. A single global gate cannot express that. The off-target model
and the allele-awareness module (`msa/`) should share a notion of "intended vs.
unintended cross-reactivity" rather than treating all non-target hits alike.

## 5. What gates implementation: the validation experiment

Do not build the model before running this test. Required inputs:

- Guide and passenger sequences of siRNAs already designed and experimentally
  tested, **including those that behaved as siDirect predicted**, not only the
  exceptions.
- siDirect's verdict for each (pass/fail, and the seed Tm value).
- The off-target readout and how it was measured. Transcriptome-wide
  (microarray/RNA-seq) carries far more weight than a small qPCR panel.
- The cell type, which defines the expressed transcriptome and therefore the
  abundance term.

The test: among siRNAs that siDirect **rejected** on seed Tm but that measured
clean, are their seeds **rare** in the expressed transcriptome relative to
seeds that siDirect accepted? If yes, section 3.2 is supported and the model is
worth building. If no, the over-conservatism has some other cause and this
direction should be dropped cheaply.

### Known risks

- **Small, biased n.** Sequences retained by a lab are usually the ones that
  worked, so the sample is not random with respect to the outcome.
- **Weak readout.** "Low off-target" from a targeted qPCR panel does not
  establish transcriptome-wide cleanliness.
- **Real infrastructure cost.** The abundance term needs a 3' UTR database plus
  cell-type expression data - not a weekend's work, and a maintenance burden.
- **Calibration data.** Moving from a gate to a score requires labelled data to
  fit against; the published siRNA efficacy datasets may not carry usable
  off-target labels.

## 6. References

- Naito, Y., Yoshimura, J., Morishita, S., Ui-Tei, K. (2009). siDirect 2.0:
  updated software for designing functional siRNA with reduced seed-dependent
  off-target effect. BMC Bioinformatics 10, 392. doi:10.1186/1471-2105-10-392
- Ui-Tei, K., Naito, Y., Nishi, K., Juni, A., Saigo, K. (2008). Thermodynamic stability and Watson-Crick base
  pairing in the seed duplex are major determinants of the efficiency of the
  siRNA-based off-target effect. Nucleic Acids Research 36(22), 7100-7109.
  doi:10.1093/nar/gkn902
- Freier, S.M., Kierzek, R., Jaeger, J.A., Sugimoto, N., Caruthers, M.H., Neilson, T., Turner, D.H. (1986). Improved free-energy parameters for predictions of
  RNA duplex stability. PNAS 83(24), 9373-9377.
- Xia, T., SantaLucia, J., Burkard, M.E., Kierzek, R., Schroeder, S.J., Jiao, X., Cox, C., Turner, D.H. (1998). Thermodynamic parameters for an expanded
  nearest-neighbor model for formation of RNA duplexes with Watson-Crick base
  pairs. Biochemistry 37(42), 14719-14735. doi:10.1021/bi9809425
- Salomon, W.E., Jolly, S.M., Moore, M.J., Zamore, P.D., Serebrov, V. (2015).
  Single-molecule imaging reveals that Argonaute reshapes the binding properties
  of its nucleic acid guides. Cell 162(1), 84-95.
- Wee, L.M., Flores-Jasso, C.F., Salomon, W.E., Zamore, P.D. (2012). Argonaute
  divides its RNA guide into domains with distinct functions and RNA-binding
  properties. Cell 151(5), 1055-1067.
