"""Intramolecular RNA hairpin stability - Reynolds' internal-structure criterion.

Supplies the Tm of the most stable hairpin an siRNA can fold into; Reynolds
(2004) awards a point when that is <= 20 C (no stable internal structure).

A hairpin melts INTRAMOLECULARLY, so unlike the duplex Tm in core.thermo/core.seed
it has no concentration term: Tm = dH / dS. The model is a transparent
approximation (not a reproduction of Reynolds' unspecified folder), with named
assumptions: perfect Watson-Crick stems via stacking_terms() (Xia 1998); loop
destabilisation from Groebe & Uhlenbeck (1988) Table 3, composition-weighted and
size-interpolated; and - because Table 3 gives loops as dG37 only - loop
destabilisation treated as purely entropic (dH_loop = 0). WC stems only (no G-U
or bulges), so internal structure is if anything under-estimated - the
conservative direction for a filter that rewards its absence.

Primary source: Groebe & Uhlenbeck 1988 (loop stability); Reynolds et al. 2004
(the criterion). Full model and assumptions in docs/design_notes.md (hairpin
section).
"""

from __future__ import annotations

from dataclasses import dataclass

from core.sequence import BioSequence
from core.thermo import (
    REFERENCE_TEMPERATURE_C,
    ZERO_CELSIUS_K,
    gibbs_free_energy,
    stacking_terms,
)

# --------------------------------------------------------------------------- #
# Geometry and threshold
# --------------------------------------------------------------------------- #

# Sterically an RNA hairpin loop needs at least 3 unpaired nucleotides, and a
# stem needs at least 3 base pairs to be a credible fold (2 gives a single
# stacking step and would flag noise). Both are conventions - Groebe & Uhlenbeck
# used a fixed stem, and Reynolds did not specify these - so they are tunable.
MIN_LOOP = 3
MIN_STEM = 3

# The Reynolds cut-off: internal structure melting at or below this is treated
# as absent (Reynolds et al. 2004). The comparison is inclusive ("<= 20 C").
DEFAULT_INTERNAL_TM_THRESHOLD_C = 20.0

_T37_K = REFERENCE_TEMPERATURE_C + ZERO_CELSIUS_K

_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


# --------------------------------------------------------------------------- #
# Groebe & Uhlenbeck (1988) loop stability
# --------------------------------------------------------------------------- #
#
# Table 3, destabilising free energy of a homopolymer loop at 37 C, 1 M NaCl.
# The paper tabulates -dG37; the values below are dG37 itself (positive =
# destabilising), for loop sizes 3, 4, 5, 7 and 9.

_LOOP_DG37: dict[str, dict[int, float]] = {
    "A": {3: 4.9, 4: 3.6, 5: 3.5, 7: 4.0, 9: 5.6},
    "C": {3: 6.6, 4: 6.9, 5: 7.1, 7: 8.3, 9: 9.1},
    "U": {3: 4.2, 4: 3.0, 5: 3.1, 7: 3.8, 9: 3.2},
}
_LOOP_SIZES = (3, 4, 5, 7, 9)


def _interp_row(row: dict[int, float], size: int) -> float:
    """Loop dG37 for one composition at a given size, linearly interpolated.

    Sizes below the smallest / above the largest tabulated value are clamped to
    the nearest measured size (the data does not extend past 9).
    """
    if size <= _LOOP_SIZES[0]:
        return row[_LOOP_SIZES[0]]
    if size >= _LOOP_SIZES[-1]:
        return row[_LOOP_SIZES[-1]]
    if size in row:
        return row[size]
    lo = max(s for s in _LOOP_SIZES if s < size)
    hi = min(s for s in _LOOP_SIZES if s > size)
    frac = (size - lo) / (hi - lo)
    return row[lo] + frac * (row[hi] - row[lo])


def loop_penalty_dg37(loop: str) -> float:
    """Destabilising dG37 (kcal/mol, positive) of a hairpin loop.

    Weighted by the loop's actual A/C/U composition and interpolated across its
    size. A base outside {A, C, U} - a G (which Groebe & Uhlenbeck excluded to
    avoid aggregation) or an ambiguity code - is scored as the mean of the three
    measured compositions at that size.
    """
    size = len(loop)
    if size == 0:
        return 0.0
    total = 0.0
    for base in loop:
        if base in _LOOP_DG37:
            total += _interp_row(_LOOP_DG37[base], size)
        else:
            total += sum(_interp_row(_LOOP_DG37[b], size) for b in "ACU") / 3.0
    return total / size


# --------------------------------------------------------------------------- #
# Hairpin fold
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class HairpinFold:
    """One intramolecular hairpin: a Watson-Crick stem closing a single loop.

    Fields:
        stem_start:  0-based index of the 5' stem arm's first base.
        stem_length: base pairs in the stem.
        loop:        the unpaired loop, 5'->3'.
        dh: fold enthalpy in kcal/mol (stem stacking; loop taken as dH = 0).
        ds: fold entropy in cal/(mol*K) (stem stacking plus the loop's entropic
            penalty).
        dg37: fold free energy at 37 C in kcal/mol (more negative = more stable).
        tm: intramolecular melting temperature in degrees C (dH / dS).
    """

    stem_start: int
    stem_length: int
    loop: str
    dh: float
    ds: float
    dg37: float
    tm: float

    def __repr__(self) -> str:
        return (
            f"HairpinFold(stem={self.stem_length}bp@{self.stem_start}, "
            f"loop={self.loop!r}, dg37={self.dg37:.2f}, tm={self.tm:.1f}C)"
        )


def _normalize(seq: BioSequence | str) -> str:
    if isinstance(seq, BioSequence):
        return seq.as_rna().seq
    return "".join(str(seq).split()).upper().replace("T", "U")


def _arms_pair(arm5: str, arm3: str) -> bool:
    """True when the two stem arms are a perfect antiparallel WC match."""
    n = len(arm3)
    for k, base in enumerate(arm5):
        partner = arm3[n - 1 - k]
        if _COMPLEMENT.get(base) != partner:
            return False
    return True


def _fold(rna: str, stem_start: int, stem_length: int, loop: str) -> HairpinFold:
    dh, ds_stem = stacking_terms(rna[stem_start : stem_start + stem_length])
    loop_dg37 = loop_penalty_dg37(loop)
    # Loop entropic penalty: dG37 = -T * dS  (dH_loop = 0), dS in cal/(mol*K).
    ds_loop = -loop_dg37 * 1000.0 / _T37_K
    ds = ds_stem + ds_loop
    dg37 = gibbs_free_energy(dh, ds_stem, REFERENCE_TEMPERATURE_C) + loop_dg37
    tm = (dh * 1000.0) / ds - ZERO_CELSIUS_K
    return HairpinFold(
        stem_start=stem_start,
        stem_length=stem_length,
        loop=loop,
        dh=dh,
        ds=ds,
        dg37=dg37,
        tm=tm,
    )


def most_stable_hairpin(
    seq: BioSequence | str,
    min_stem: int = MIN_STEM,
    min_loop: int = MIN_LOOP,
) -> HairpinFold | None:
    """Return the minimum-dG37 simple hairpin foldable within 'seq', or None.

    Searches every perfect Watson-Crick stem (>= min_stem bp) closing a single
    loop (>= min_loop nt) and returns the most stable one by free energy at 37 C.
    None means no such hairpin exists - the sequence has no internal structure of
    the kind the criterion looks for.
    """
    rna = _normalize(seq)
    n = len(rna)
    best: HairpinFold | None = None
    for stem_start in range(n):
        max_stem = (n - stem_start - min_loop) // 2
        for stem_length in range(min_stem, max_stem + 1):
            arm5 = rna[stem_start : stem_start + stem_length]
            loop_start = stem_start + stem_length
            for loop_length in range(min_loop, n - loop_start - stem_length + 1):
                arm3_start = loop_start + loop_length
                arm3 = rna[arm3_start : arm3_start + stem_length]
                if not _arms_pair(arm5, arm3):
                    continue
                fold = _fold(rna, stem_start, stem_length, rna[loop_start:arm3_start])
                if best is None or fold.dg37 < best.dg37:
                    best = fold
    return best


def internal_structure_tm(
    seq: BioSequence | str,
    min_stem: int = MIN_STEM,
    min_loop: int = MIN_LOOP,
) -> float | None:
    """Melting temperature (degrees C) of the most stable internal hairpin.

    None when no hairpin can form, i.e. there is no internal structure to melt.
    """
    fold = most_stable_hairpin(seq, min_stem=min_stem, min_loop=min_loop)
    return None if fold is None else fold.tm


def has_no_internal_repeat(
    seq: BioSequence | str,
    threshold_c: float = DEFAULT_INTERNAL_TM_THRESHOLD_C,
    min_stem: int = MIN_STEM,
    min_loop: int = MIN_LOOP,
) -> bool:
    """True when 'seq' has no stable internal structure (Reynolds criterion).

    The most stable internal hairpin must melt at or below 'threshold_c' (20 C by
    default). A sequence that folds into nothing trivially satisfies this.
    """
    tm = internal_structure_tm(seq, min_stem=min_stem, min_loop=min_loop)
    return tm is None or tm <= threshold_c
