"""Nearest-neighbor duplex thermodynamics - the engine every other unit builds on.

Converts a short RNA sequence into (dH, dS), from which dG at any temperature and
the melting temperature follow, via the nearest-neighbor stacking model.

Two entry points, deliberately separate - using the wrong one silently shifts
every number by the initiation term, so the choice is explicit at the call site:

* duplex_thermo()   - a complete, isolated duplex (initiation, terminal A-U
  penalties and the self-complementarity correction included).
* stacking_terms()  - the stacking contributions alone, for a window INSIDE a
  longer helix whose boundaries are not real helix ends (loading asymmetry, a
  hairpin stem).

Two parameter sets are exposed by name, Xia 1998 (default) and Freier 1986,
because a published threshold is only valid against the set it was calibrated on;
the caller picks. Watson-Crick pairs only, 1 M NaCl, no salt correction.

Primary source: Xia et al. 1998. Model, the two-parameter-set rationale, the
symmetry-derivation of the tables, and the Biopython cross-check are in
docs/design_notes.md (thermo section).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from core.errors import SequenceValidationError
from core.sequence import BioSequence

# --------------------------------------------------------------------------- #
# Physical constants and conventions
# --------------------------------------------------------------------------- #

# Gas constant in cal/(mol*K), matching the cal-based dS units of the tables
# (dH is in kcal/mol throughout, so it is scaled by 1000 where the two meet).
R_CAL = 1.9872

ZERO_CELSIUS_K = 273.15

# dG is conventionally quoted at 37 C in this literature; dg37 is the shorthand.
REFERENCE_TEMPERATURE_C = 37.0

# Default total strand concentration (mol/L) for Tm - a convention, not a
# property of the sequence. A published score passes its own (core.seed uses
# siDirect's 100 uM).
DEFAULT_TOTAL_STRAND_CONC_M = 1e-6

# The four unambiguous RNA bases; ambiguity codes have no NN parameters.
RNA_WC_BASES = frozenset("ACGU")

_RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


# --------------------------------------------------------------------------- #
# Nearest-neighbor parameters
# --------------------------------------------------------------------------- #
#
# Keyed by the top-strand 5'->3' dinucleotide (the partner is implied by WC
# pairing). Only these 10 steps are independent; _expand_by_symmetry derives the
# other 6 by reverse-complement, so no value is typed twice. Values are
# (dH kcal/mol, dS cal/mol/K), 1 M NaCl. See docs/design_notes.md (thermo).

_XIA_1998_PUBLISHED: dict[str, tuple[float, float]] = {
    "AA": (-6.82, -19.0),   # AA/UU
    "AU": (-9.38, -26.7),   # AU/AU
    "UA": (-7.69, -20.5),   # UA/UA
    "CU": (-10.48, -27.1),  # CU/AG
    "CA": (-10.44, -26.9),  # CA/UG
    "GU": (-11.40, -29.5),  # GU/AC
    "GA": (-12.44, -32.5),  # GA/UC
    "CG": (-10.64, -26.7),  # CG/CG
    "GG": (-13.39, -32.7),  # GG/CC
    "GC": (-14.88, -36.9),  # GC/GC
}

# The older Freier et al. (1986) RNA parameters, same 10 independent steps. Kept
# because siDirect's seed filter is calibrated on them (see core.seed); not the
# default. Mixing a Freier-calibrated threshold with Xia values (or vice versa)
# shifts every number relative to the cut-off.
_FREIER_1986_PUBLISHED: dict[str, tuple[float, float]] = {
    "AA": (-6.6, -18.4),   # AA/UU
    "AU": (-5.7, -15.5),   # AU/AU
    "UA": (-8.1, -22.6),   # UA/UA
    "CU": (-7.6, -19.2),   # CU/AG
    "CA": (-10.5, -27.8),  # CA/UG
    "GU": (-10.2, -26.2),  # GU/AC
    "GA": (-13.3, -35.5),  # GA/UC
    "CG": (-8.0, -19.4),   # CG/CG
    "GG": (-12.2, -29.7),  # GG/CC
    "GC": (-14.2, -34.9),  # GC/GC
}

# Duplex initiation, applied once per duplex.
INIT_DH, INIT_DS = 3.61, -1.5

# Applied once per helix end that is closed by an A-U pair (destabilising: an
# A-U end frays more readily than a G-C end).
TERMINAL_AU_DH, TERMINAL_AU_DS = 3.72, 10.5

# Entropic correction applied once when the strand is its own reverse-
# complement, which halves the number of distinguishable duplex orientations.
SYMMETRY_DH, SYMMETRY_DS = 0.0, -1.4


def _reverse_complement(step: str) -> str:
    return "".join(_RNA_COMPLEMENT[base] for base in reversed(step))


def _expand_by_symmetry(
    published: dict[str, tuple[float, float]],
) -> dict[str, tuple[float, float]]:
    """Complete the 16-step table from the 10 published, independent steps."""
    table = dict(published)
    for step, values in published.items():
        table.setdefault(_reverse_complement(step), values)
    return table


#: dH/dS for all 16 nearest-neighbor steps, keyed by the top-strand 5'->3' pair.
NN_PARAMS_XIA_1998 = _expand_by_symmetry(_XIA_1998_PUBLISHED)
NN_PARAMS_FREIER_1986 = _expand_by_symmetry(_FREIER_1986_PUBLISHED)

#: The default parameter set for general-purpose stability work.
NN_PARAMS = NN_PARAMS_XIA_1998


class ThermoParams(str, Enum):
    """A nearest-neighbor parameter set, selectable by name.

    A hashable stand-in for the raw tables above, so callers (e.g. a frozen
    SelectionCriteria) can carry a choice of parameter set without embedding an
    unhashable dict. Resolve to the actual table with `.table`.

    XIA_1998 is the modern default. FREIER_1986 exists because published
    thresholds calibrated on it (siDirect's seed Tm) are only valid against it -
    swapping tables silently shifts every number relative to such a cut-off, so
    the set is named explicitly rather than chosen for the caller.
    """

    XIA_1998 = "xia"
    FREIER_1986 = "freier"

    @property
    def table(self) -> dict[str, tuple[float, float]]:
        """The nearest-neighbor table this name selects."""
        if self is ThermoParams.FREIER_1986:
            return NN_PARAMS_FREIER_1986
        return NN_PARAMS_XIA_1998


# --------------------------------------------------------------------------- #
# Input normalisation
# --------------------------------------------------------------------------- #


def as_wc_rna(seq: BioSequence | str) -> str:
    """Return 'seq' as an unambiguous 5'->3' RNA string.

    Accepts a BioSequence (DNA is transcribed; a protein raises MolTypeError) or
    a bare string (whitespace stripped, uppercased, T read as U). Any residue
    outside ACGU - including N and the IUPAC ambiguity codes - has no NN
    parameter and is rejected.
    """
    if isinstance(seq, BioSequence):
        text = seq.as_rna().seq
    else:
        text = "".join(str(seq).split()).upper().replace("T", "U")

    bad = sorted(set(text) - RNA_WC_BASES)
    if bad:
        raise SequenceValidationError(
            f"Nearest-neighbor thermodynamics need unambiguous RNA (ACGU); "
            f"found {bad} in {text!r:.40}"
        )
    return text


def is_self_complementary(seq: BioSequence | str) -> bool:
    """True when the strand is its own reverse-complement (a palindrome duplex).

    Such a strand pairs with a copy of itself, which changes both the entropy
    (symmetry correction) and the concentration term of the Tm equation.
    """
    rna = as_wc_rna(seq)
    return bool(rna) and rna == _reverse_complement(rna)


# --------------------------------------------------------------------------- #
# Free energy and melting temperature
# --------------------------------------------------------------------------- #


def gibbs_free_energy(
    dh: float, ds: float, temperature_c: float = REFERENCE_TEMPERATURE_C
) -> float:
    """dG = dH - T*dS in kcal/mol, for dH in kcal/mol and dS in cal/(mol*K)."""
    return dh - (temperature_c + ZERO_CELSIUS_K) * ds / 1000.0


def melting_temperature(
    dh: float,
    ds: float,
    total_strand_conc: float = DEFAULT_TOTAL_STRAND_CONC_M,
    self_complementary: bool = False,
) -> float:
    """Melting temperature in degrees C from duplex dH and dS.

        Tm = dH / (dS + R * ln(Ct / x))

    with Ct the total strand concentration and x = 4 for two different strands
    present in equal amounts, or x = 1 for a self-complementary strand (which
    has no distinct partner to be diluted against).

    Args:
        dh: duplex enthalpy in kcal/mol (initiation and end terms included).
        ds: duplex entropy in cal/(mol*K).
        total_strand_conc: total strand concentration in mol/L (both strands).
        self_complementary: whether the strand pairs with a copy of itself.

    Raises:
        ValueError: if the concentration is not positive, or if the denominator
            vanishes (no physical melting temperature for those inputs).
    """
    if total_strand_conc <= 0:
        raise ValueError(
            f"Total strand concentration must be positive, got {total_strand_conc}."
        )
    effective = total_strand_conc / (1.0 if self_complementary else 4.0)
    denominator = ds + R_CAL * math.log(effective)
    if denominator == 0:
        raise ValueError(
            "Tm is undefined: dS plus the concentration term is exactly zero."
        )
    return (dh * 1000.0) / denominator - ZERO_CELSIUS_K


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class DuplexThermo:
    """Thermodynamics of one complete, perfectly Watson-Crick paired RNA duplex.

    dh and ds describe the whole duplex: stacking, initiation, terminal A-U
    penalties and the symmetry correction are already summed in. They are in the
    Xia reference state (1 M NaCl), uncorrected for any other salt condition.

    Fields:
        sequence: the top strand, 5'->3', as normalised RNA. The partner strand
            is its reverse-complement and is not stored.
        dh: duplex enthalpy in kcal/mol.
        ds: duplex entropy in cal/(mol*K).
        terminal_au: how many of the two helix ends are closed by an A-U pair.
        self_complementary: whether the strand is its own reverse-complement.
    """

    sequence: str
    dh: float
    ds: float
    terminal_au: int
    self_complementary: bool

    def dg(self, temperature_c: float = REFERENCE_TEMPERATURE_C) -> float:
        """Duplex free energy in kcal/mol at 'temperature_c' (more negative = more stable)."""
        return gibbs_free_energy(self.dh, self.ds, temperature_c)

    @property
    def dg37(self) -> float:
        """Duplex free energy in kcal/mol at 37 C, the conventional reference."""
        return self.dg(REFERENCE_TEMPERATURE_C)

    def tm(
        self, total_strand_conc: float = DEFAULT_TOTAL_STRAND_CONC_M
    ) -> float:
        """Melting temperature in degrees C at the given total strand concentration."""
        return melting_temperature(
            self.dh,
            self.ds,
            total_strand_conc=total_strand_conc,
            self_complementary=self.self_complementary,
        )

    def __repr__(self) -> str:
        return (
            f"DuplexThermo(seq={self.sequence!r}, dg37={self.dg37:.2f} kcal/mol, "
            f"dH={self.dh:.2f}, dS={self.ds:.1f}, terminal_au={self.terminal_au})"
        )


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def stacking_terms(
    seq: BioSequence | str,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> tuple[float, float]:
    """Sum the nearest-neighbor stacking terms of a strand: (dH, dS).

    This is the helix interior only - no initiation, no terminal A-U penalty, no
    symmetry correction. Use it to compare windows WITHIN a longer duplex, where
    the window edges are not real helix ends; use duplex_thermo() for a duplex
    that actually exists on its own.

    A strand of n bases has n-1 steps; a single base has none and sums to zero.

    Args:
        seq: the strand as a BioSequence (DNA is transcribed to RNA) or a string.
        params: which nearest-neighbor table to sum with. Defaults to Xia 1998;
            pass NN_PARAMS_FREIER_1986 when reproducing a Freier-calibrated
            published score.

    Raises:
        SequenceValidationError: if the sequence contains anything but ACGU.
    """
    rna = as_wc_rna(seq)
    dh = 0.0
    ds = 0.0
    for i in range(len(rna) - 1):
        step_dh, step_ds = params[rna[i : i + 2]]
        dh += step_dh
        ds += step_ds
    return dh, ds


def duplex_thermo(seq: BioSequence | str) -> DuplexThermo:
    """Full nearest-neighbor thermodynamics of the duplex formed by 'seq'.

    'seq' is the top strand 5'->3'; its partner is taken to be the perfect
    reverse-complement. Stacking terms are summed, then duplex initiation, one
    terminal A-U penalty per A-U-closed end, and the symmetry correction if the
    strand is self-complementary.

    Xia 1998 throughout - unlike stacking_terms() this takes no parameter-set
    argument, because the initiation and terminal-A-U terms applied here are
    Xia's own and do not carry over to another set (Freier 1986, for one,
    defines no terminal A-U penalty at all).

    Args:
        seq: the strand as a BioSequence (DNA is transcribed to RNA) or a string.

    Raises:
        SequenceValidationError: if the sequence contains anything but ACGU.
        ValueError: if the sequence is shorter than 2 nt - a duplex needs at
            least one nearest-neighbor step to be defined.
    """
    rna = as_wc_rna(seq)
    if len(rna) < 2:
        raise ValueError(
            f"A nearest-neighbor duplex needs at least 2 nt; got {len(rna)}."
        )

    dh, ds = stacking_terms(rna)

    dh += INIT_DH
    ds += INIT_DS

    terminal_au = sum(1 for end in (rna[0], rna[-1]) if end in ("A", "U"))
    dh += terminal_au * TERMINAL_AU_DH
    ds += terminal_au * TERMINAL_AU_DS

    self_complementary = rna == _reverse_complement(rna)
    if self_complementary:
        dh += SYMMETRY_DH
        ds += SYMMETRY_DS

    return DuplexThermo(
        sequence=rna,
        dh=dh,
        ds=ds,
        terminal_au=terminal_au,
        self_complementary=self_complementary,
    )
