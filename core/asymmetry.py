"""Duplex-end thermodynamic asymmetry - which strand RISC loads.

RISC keeps the strand whose 5' end sits at the duplex's less stable end. This
module measures that directly, in kcal/mol, rather than through the A/U proxies
the 2004 rule sets use:

    differential = dG(guide 5' end window) - dG(sense 5' end window)

A POSITIVE differential means the guide's 5' end is the looser end - the
orientation that loads the intended strand. Negative warns the passenger is the
thermodynamically preferred guide.

Built from stacking_terms() (not duplex_thermo()): these windows are the ends of
a longer helix, so initiation and terminal-A-U terms do not belong. The window
length is a parameter (default 4 bp; no canonical value exists), and NO THRESHOLD
IS DEFINED - the module returns a measured quantity and its sign, since turning
kcal/mol into a cut-off needs efficacy data HUSH does not have.

Primary source: Khvorova et al. 2003; Schwarz et al. 2003. Rationale, the
stacking_terms choice, and the proxy-agreement check are in docs/design_notes.md
(asymmetry section).
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

# --------------------------------------------------------------------------- #
# Window policy
# --------------------------------------------------------------------------- #

# Number of terminal base pairs compared at each end. A convention, not a
# published constant - see the module docstring. A window of n base pairs spans
# n-1 nearest-neighbor steps, so 2 is the smallest meaningful value.
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


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class EndAsymmetry:
    """Thermodynamic asymmetry between the two ends of one siRNA duplex.

    Both free energies are stacking-only sums over the terminal window, in
    kcal/mol at the stated temperature, and both describe the SAME duplex read
    from its two ends - so they are directly comparable.

    Fields:
        window: number of terminal base pairs compared at each end.
        temperature_c: temperature the free energies were evaluated at.
        guide_5p_dg: stacking dG of the window at the guide's 5' end.
        sense_5p_dg: stacking dG of the window at the sense strand's 5' end
            (the same physical end as the guide's 3' end).
        differential: guide_5p_dg - sense_5p_dg. Positive favours the guide.
    """

    window: int
    temperature_c: float
    guide_5p_dg: float
    sense_5p_dg: float
    differential: float

    @property
    def favours_guide(self) -> bool:
        """True when the guide's 5' end is the less stable end of the duplex.

        This is the configuration that loads the intended strand. Exactly zero
        counts as False: a symmetric duplex expresses no preference, so there is
        nothing favouring the guide.
        """
        return self.differential > 0

    @property
    def preferred_strand(self) -> str:
        """Which strand this duplex is predicted to load: guide/passenger/none."""
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


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def terminal_stacking_dg(
    strand: str,
    window: int = DEFAULT_WINDOW_BP,
    temperature_c: float = REFERENCE_TEMPERATURE_C,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> float:
    """Stacking free energy (kcal/mol) of the 5'-terminal 'window' bases.

    Stacking terms only - no initiation, no terminal A-U penalty - because this
    window is the end of a longer helix, not a duplex in its own right. 'params'
    is the nearest-neighbor table (default Xia 1998).
    """
    dh, ds = stacking_terms(strand[:window], params)
    return gibbs_free_energy(dh, ds, temperature_c)


def end_asymmetry(
    candidate: SiRNACandidate,
    window: int = DEFAULT_WINDOW_BP,
    temperature_c: float = REFERENCE_TEMPERATURE_C,
    params: dict[str, tuple[float, float]] = NN_PARAMS,
) -> EndAsymmetry:
    """Measure the duplex-end thermodynamic asymmetry of a candidate.

    Compares the terminal stacking free energy at each end of the duplex, taking
    each end from the 5' side of the strand that starts there: the guide's 5'
    end from the antisense strand, and the opposite end from the sense strand's
    5' end.

    Args:
        candidate: the siRNA candidate to measure.
        window: terminal base pairs to compare at each end (default 4).
        temperature_c: temperature for dG (default 37 C).
        params: nearest-neighbor table for the stacking dG (default Xia 1998).

    Returns:
        An EndAsymmetry whose 'differential' is positive when the guide's 5' end
        is the less stable end - the orientation that loads the guide.

    Raises:
        ValueError: if the window is smaller than MIN_WINDOW_BP, or large enough
            that the two windows would overlap in the core.
        SequenceValidationError: if either strand contains anything but ACGU.
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
