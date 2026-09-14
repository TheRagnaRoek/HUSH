"""HUSH MSA / allele module (pipeline stage 3).

Multiple-sequence alignment of the target against a related sequence panel
(e.g. an HLA allele family), with per-position variation and window-level
selection logic for conserved-region (pan-allele) vs. polymorphic
(allele-specific) targeting.

Every panel member is projected onto ungapped target coordinates, so a window on
the target lines up one-for-one with a siRNA candidate at the same position. The
stage answers, per window, whether a candidate built from the target would hit
the alleles the user wants hit and discriminate against the ones they want
spared. It is a filter, not a ranking - differences are counted, never weighted.

Submodules:

* alignment  - the alignment result type (PanelAlignment / AlignedRow) every
               backend produces and every consumer reads.
* align_star - the star-alignment backend (each member aligned pairwise to the
               target; Gusfield center-star construction).
* variation  - per-position match / mismatch / uncovered comparison of the
               aligned panel against the target.
* windows    - allele evaluation of every candidate window (the hit / spare
               gates), plus region merging.
* targeting  - the request model (AlleleTargeting): which alleles to hit and
               which to spare, with pan-allele / allele-specific / subset
               constructors.
* msa_errors - the stage's exception vocabulary, all deriving from HushError.

The stage is wired into the CLI through ``hush design --panel`` (see app.cli):
the allele gate is ANDed with the siRNA selection filters at the app layer, so
core.selection never depends on this module.
"""
