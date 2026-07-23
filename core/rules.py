"""Classical siRNA design rules (pipeline stage 2, rule pre-screen).

Three empirical rule sets from 2004 form the backbone of rational siRNA design.
They all encode thermodynamic asymmetry - the duplex should unwind more easily
at the guide (antisense) 5' end so RISC loads the guide rather than the
passenger - but express it differently:

* Ui-Tei    - a strict 4-rule pass/fail filter (implemented here).
* Reynolds  - an 8-point additive checklist (later unit).
* Amarzguioui - a compact weighted regression score (later unit).

Strand conventions match core.sirna: sense/passenger and antisense/guide are
both stored 5'->3'. Rules are evaluated on whichever strand they are naturally
defined for and documented per rule.

Reference:
    Ui-Tei, K. et al. (2004). Guidelines for the selection of highly effective
    siRNA sequences for mammalian and chick RNA interference.
    Nucleic Acids Research 32(3), 936-948.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.sirna import SiRNACandidate

# --------------------------------------------------------------------------- #
# Ui-Tei constants
# --------------------------------------------------------------------------- #

# The guide 5'-terminal window Ui-Tei rule (iii) inspects for A/U richness, and
# the minimum A/U count required within it. These are fixed values from the
# original guidelines, defined for the standard 19-mer; they are not rescaled
# for longer cores.
UI_TEI_AU_WINDOW = 7
UI_TEI_MIN_AU_IN_WINDOW = 4

# Ui-Tei rule (iv): a GC stretch must not be LONGER than this many nt, so a run
# of length UI_TEI_MAX_GC_RUN + 1 or more is a violation.
UI_TEI_MAX_GC_RUN = 9

_AU = frozenset("AU")
_GC = frozenset("GC")


def _longest_run(seq: str, chars: frozenset[str]) -> int:
    """Length of the longest maximal run of characters drawn from 'chars'."""
    longest = 0
    current = 0
    for base in seq:
        if base in chars:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class UiTeiResult:
    """Outcome of the four Ui-Tei rules for one candidate.

    Each field is the pass/fail of one rule. A candidate is "highly effective"
    under Ui-Tei only when all four hold (see 'satisfied').

    Fields:
        guide_5p_au:    (i)   A/U at the guide (antisense) 5' terminus.
        sense_5p_gc:    (ii)  G/C at the sense (passenger) 5' terminus.
        guide_5p_au_rich: (iii) >= 4 A/U within the guide's 5'-terminal 7 nt.
        no_long_gc_run: (iv)  no GC stretch longer than 9 nt.
    """

    guide_5p_au: bool
    sense_5p_gc: bool
    guide_5p_au_rich: bool
    no_long_gc_run: bool

    @property
    def satisfied(self) -> bool:
        """True when all four rules pass (Ui-Tei 'highly effective')."""
        return (
            self.guide_5p_au
            and self.sense_5p_gc
            and self.guide_5p_au_rich
            and self.no_long_gc_run
        )

    @property
    def score(self) -> int:
        """Number of rules passed (0-4); a convenience for ranking/inspection."""
        return sum(
            (
                self.guide_5p_au,
                self.sense_5p_gc,
                self.guide_5p_au_rich,
                self.no_long_gc_run,
            )
        )

    def __repr__(self) -> str:
        return (
            f"UiTeiResult(satisfied={self.satisfied}, score={self.score}/4, "
            f"guide_5p_au={self.guide_5p_au}, sense_5p_gc={self.sense_5p_gc}, "
            f"guide_5p_au_rich={self.guide_5p_au_rich}, "
            f"no_long_gc_run={self.no_long_gc_run})"
        )


# --------------------------------------------------------------------------- #
# Ui-Tei evaluation
# --------------------------------------------------------------------------- #


def ui_tei(candidate: SiRNACandidate) -> UiTeiResult:
    """Evaluate the four Ui-Tei rules for a candidate.

    Rules (guide = antisense, sense = passenger; both 5'->3'):
        (i)   guide[0] is A or U.
        (ii)  sense[0] is G or C.
        (iii) at least 4 of guide[0:7] are A or U.
        (iv)  the longest G/C run is not longer than 9 nt.

    Rule (iv) is strand-symmetric (a GC run in the sense strand is a CG run of
    equal length in the antisense), so it is computed once on the sense strand.
    Ambiguous residues (N and IUPAC codes) belong to neither the A/U nor the G/C
    set, so they simply fail the base-identity rules rather than being guessed.
    """
    guide = candidate.antisense.seq
    sense = candidate.sense.seq

    guide_5p_au = bool(guide) and guide[0] in _AU
    sense_5p_gc = bool(sense) and sense[0] in _GC
    au_in_window = sum(base in _AU for base in guide[:UI_TEI_AU_WINDOW])
    guide_5p_au_rich = au_in_window >= UI_TEI_MIN_AU_IN_WINDOW
    no_long_gc_run = _longest_run(sense, _GC) <= UI_TEI_MAX_GC_RUN

    return UiTeiResult(
        guide_5p_au=guide_5p_au,
        sense_5p_gc=sense_5p_gc,
        guide_5p_au_rich=guide_5p_au_rich,
        no_long_gc_run=no_long_gc_run,
    )
