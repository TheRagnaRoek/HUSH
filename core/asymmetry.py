"""Duplex-end thermodynamic asymmetry (Khvorova 2003, Schwarz 2003).

    differential = dG(guide 5' end window) - dG(sense 5' end window)

Positive means the guide's 5' end is the less stable end. The windows are helix
ends, so only stacking terms are used (core.thermo.stacking_terms). No
magnitude threshold is defined; see docs/design-rules.md.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.sirna import SiRNACandidate
from core.thermo import (
    NN_PARAMS,
    REFERENCE_TEMPERATURE_C,
    gibbs_free_energy,
    stacking_terms,
)

# Terminal base pairs compared at each end; a convention, not a published value.
DEFAULT_WINDOW_BP = 4
MIN_WINDOW_BP = 2


def _validate_window(window: int, core_length: int) -> None:
    if not isinstance(window, bool) and not isinstance(window, int):
        raise TypeError(f"window must be an int, got {type(window).__name__}")
    if window < MIN_WINDOW_BP:
        raise ValueError(
            f"A stacking window needs at least {MIN_WINDOW_BP} base pairs to "
            f"contain one nearest-neighbor step; got {window}."
        )
    if 2 * window > core_length:
        raise ValueError(
            f"Two {window}-bp windows overlap in a {core_length}-bp core; the "
            f"two duplex ends would share base pairs and the differential would "
            f"be meaningless. Use a window of at most {core_length // 2}."
        )


@dataclass(frozen=True)
class EndAsymmetry:
    """Stacking dG (kcal/mol) at both duplex ends and their differential."""

    window: int
    temperature_c: float
    guide_5p_dg: float
    sense_5p_dg: float
    differential: float

    @property
    def favours_guide(self) -> bool:
        """True when the guide's 5' end is the less stable end; zero is False."""
        return self.differential > 0

    @property
    def preferred_strand(self) -> str:
        """Predicted loaded strand: "guide", "passenger" or "none"."""
        if self.differential > 0:
            return "guide"
        if self.differential < 0:
            return "passenger"
        return "none"

    def __repr__(self) -> str:
        return (
            f"EndAsymmetry(differential={self.differential:+.2f} kcal/mol, "
            f"preferred={self.preferred_strand}, "
            f"guide_5p={self.guide_5p_dg:.2f}, sense_5p={self.sense_5p_dg:.2f}, "
            f"window={self.window}bp)"
        )


def terminal_stacking_dg(
    strand: str,
    window: int = DEFAULT_WINDOW_BP,
    temperature_c: float = REFERENCE_TEMPERATURE_C,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> float:
    """Stacking dG (kcal/mol) of the first 'window' bases of a strand."""
    dh, ds = stacking_terms(strand[:window], params)
    return gibbs_free_energy(dh, ds, temperature_c)


def end_asymmetry(
    candidate: SiRNACandidate,
    window: int = DEFAULT_WINDOW_BP,
    temperature_c: float = REFERENCE_TEMPERATURE_C,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> EndAsymmetry:
    """Asymmetry of a candidate, from the 5' ends of its guide and sense strands.

    Raises ValueError for a window under MIN_WINDOW_BP or one whose two windows
    would overlap in the core.
    """
    _validate_window(window, candidate.length)

    guide_5p_dg = terminal_stacking_dg(
        candidate.antisense.seq, window, temperature_c, params
    )
    sense_5p_dg = terminal_stacking_dg(
        candidate.sense.seq, window, temperature_c, params
    )

    return EndAsymmetry(
        window=window,
        temperature_c=temperature_c,
        guide_5p_dg=guide_5p_dg,
        sense_5p_dg=sense_5p_dg,
        differential=guide_5p_dg - sense_5p_dg,
    )
