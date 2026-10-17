"""Nearest-neighbor RNA duplex thermodynamics (Xia 1998, Freier 1986).

duplex_thermo() handles a complete duplex; stacking_terms() sums the stacking
terms only, for windows inside a longer helix. Model, parameters and formulas:
docs/design-rules.md (Nearest-neighbor thermodynamics).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from core.errors import SequenceValidationError
from core.sequence import BioSequence

# cal/(mol*K); dH is in kcal/mol, so it is scaled by 1000 where the two meet.
R_CAL = 1.9872

ZERO_CELSIUS_K = 273.15

REFERENCE_TEMPERATURE_C = 37.0

# Default total strand concentration for Tm, in mol/L.
DEFAULT_TOTAL_STRAND_CONC_M = 1e-6

# Ambiguity codes have no nearest-neighbor parameters.
RNA_WC_BASES = frozenset("ACGU")

_RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}

# Keyed by the top-strand 5'->3' dinucleotide. Only these 10 steps are
# independent; _expand_by_symmetry derives the other 6.
# Values are (dH kcal/mol, dS cal/mol/K) at 1 M NaCl.
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

# Freier 1986 is kept because the seed filter's cut-off is calibrated on it.
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

# Xia 1998 duplex terms: initiation (once per duplex), penalty per helix end
# closed by A-U, and the entropy correction for a self-complementary strand.
INIT_DH, INIT_DS = 3.61, -1.5
TERMINAL_AU_DH, TERMINAL_AU_DS = 3.72, 10.5
SYMMETRY_DH, SYMMETRY_DS = 0.0, -1.4


def _reverse_complement(step: str) -> str:
    return "".join(_RNA_COMPLEMENT[base] for base in reversed(step))


def _expand_by_symmetry(
    published: dict[str, tuple[float, float]],
) -> dict[str, tuple[float, float]]:
    """Complete the 16-step table from the 10 independent published steps."""
    table = dict(published)
    for step, values in published.items():
        table.setdefault(_reverse_complement(step), values)
    return table


NN_PARAMS_XIA_1998 = _expand_by_symmetry(_XIA_1998_PUBLISHED)
NN_PARAMS_FREIER_1986 = _expand_by_symmetry(_FREIER_1986_PUBLISHED)

# Default parameter set.
NN_PARAMS = NN_PARAMS_XIA_1998


class ThermoParams(str, Enum):
    """Named nearest-neighbor parameter set; hashable, unlike the raw tables."""

    XIA_1998 = "xia"
    FREIER_1986 = "freier"

    @property
    def table(self) -> dict[str, tuple[float, float]]:
        """The nearest-neighbor table for this set."""
        if self is ThermoParams.FREIER_1986:
            return NN_PARAMS_FREIER_1986
        return NN_PARAMS_XIA_1998


def as_wc_rna(seq: BioSequence | str) -> str:
    """The sequence as a 5'->3' ACGU string.

    Strings are uppercased with T read as U. Anything outside ACGU, ambiguity
    codes included, raises SequenceValidationError.
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
    """True when the strand equals its own reverse complement."""
    rna = as_wc_rna(seq)
    return bool(rna) and rna == _reverse_complement(rna)


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
    """Melting temperature in degrees C from duplex dH (kcal/mol) and dS (cal/mol/K).

        Tm = dH / (dS + R * ln(Ct / x))

    Ct is the total strand concentration (mol/L); x = 4 for two distinct
    strands, 1 for a self-complementary strand. Raises ValueError for a
    non-positive concentration or a zero denominator.
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


@dataclass(frozen=True)
class DuplexThermo:
    """Thermodynamics of one complete, perfectly paired RNA duplex (1 M NaCl).

    sequence is the top strand; dh (kcal/mol) and ds (cal/mol/K) already include
    initiation, terminal A-U and symmetry terms. terminal_au counts the helix
    ends closed by an A-U pair.
    """

    sequence: str
    dh: float
    ds: float
    terminal_au: int
    self_complementary: bool

    def dg(self, temperature_c: float = REFERENCE_TEMPERATURE_C) -> float:
        """Free energy in kcal/mol at temperature_c."""
        return gibbs_free_energy(self.dh, self.ds, temperature_c)

    @property
    def dg37(self) -> float:
        """Free energy in kcal/mol at 37 C."""
        return self.dg(REFERENCE_TEMPERATURE_C)

    def tm(
        self, total_strand_conc: float = DEFAULT_TOTAL_STRAND_CONC_M
    ) -> float:
        """Melting temperature in degrees C."""
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


def stacking_terms(
    seq: BioSequence | str,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> tuple[float, float]:
    """Summed nearest-neighbor stacking terms (dH, dS) of a strand.

    Helix interior only: no initiation, terminal A-U or symmetry terms. Use it
    for windows inside a longer helix. A single base sums to zero.
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
    """Thermodynamics of the duplex of seq with its perfect complement.

    Uses Xia 1998 throughout: its initiation and terminal A-U terms do not
    transfer to another parameter set. Raises ValueError below 2 nt.
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
