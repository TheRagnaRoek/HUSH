"""Classical siRNA rule sets: Ui-Tei, Reynolds and Amarzguioui (all 2004).

Strands are 5'->3'. Each rule runs on the strand it is defined for. Rationale
and combination policy: docs/design-rules.md.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.hairpin import has_no_internal_repeat
from core.sirna import SiRNACandidate

# Ui-Tei rule definitions: docs/design-rules.md (The three rule sets).
UI_TEI_AU_WINDOW = 7
UI_TEI_MIN_AU_IN_WINDOW = 4
UI_TEI_MAX_GC_RUN = 9

_AU = frozenset("AU")
_GC = frozenset("GC")


def _longest_run(seq: str, chars: frozenset[str]) -> int:
    """Length of the longest run of characters from chars."""
    longest = 0
    current = 0
    for base in seq:
        if base in chars:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


@dataclass(frozen=True)
class UiTeiResult:
    """Pass/fail of each Ui-Tei rule; a candidate qualifies when all four pass."""

    guide_5p_au: bool
    sense_5p_gc: bool
    guide_5p_au_rich: bool
    no_long_gc_run: bool

    @property
    def satisfied(self) -> bool:
        """True when all four rules pass."""
        return (
            self.guide_5p_au
            and self.sense_5p_gc
            and self.guide_5p_au_rich
            and self.no_long_gc_run
        )

    @property
    def score(self) -> int:
        """Number of rules passed (0-4)."""
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


def ui_tei(candidate: SiRNACandidate) -> UiTeiResult:
    """Evaluate the four Ui-Tei rules (docs/design-rules.md, The three rule sets)."""
    guide = candidate.antisense.seq
    sense = candidate.sense.seq

    guide_5p_au = bool(guide) and guide[0] in _AU
    sense_5p_gc = bool(sense) and sense[0] in _GC
    au_in_window = sum(base in _AU for base in guide[:UI_TEI_AU_WINDOW])
    guide_5p_au_rich = au_in_window >= UI_TEI_MIN_AU_IN_WINDOW
    # A G/C run has the same length on both strands, so one strand suffices.
    no_long_gc_run = _longest_run(sense, _GC) <= UI_TEI_MAX_GC_RUN

    return UiTeiResult(
        guide_5p_au=guide_5p_au,
        sense_5p_gc=sense_5p_gc,
        guide_5p_au_rich=guide_5p_au_rich,
        no_long_gc_run=no_long_gc_run,
    )


# Reynolds scoring: docs/design-rules.md (The three rule sets).
REYNOLDS_MIN_LENGTH = 19
REYNOLDS_GC_MIN = 0.30
REYNOLDS_GC_MAX = 0.52
REYNOLDS_AU_3P_WINDOW = 5
REYNOLDS_POS_A_INTERNAL = 3
REYNOLDS_POS_U = 10
REYNOLDS_POS_G_PENALTY = 13
REYNOLDS_DEFAULT_THRESHOLD = 6


def _base_at(seq: str, pos_1based: int) -> str:
    """Residue at a 1-based position from the 5' end."""
    return seq[pos_1based - 1]


@dataclass(frozen=True)
class ReynoldsResult:
    """Components of the Reynolds score; the two g_/gc_ fields are penalties."""

    gc_in_range: bool
    no_internal_repeat: bool
    a_at_pos19: bool
    a_at_pos3: bool
    u_at_pos10: bool
    au_count_3p: int
    g_at_pos13: bool
    gc_at_pos19: bool

    @property
    def positive_score(self) -> int:
        """Points from the positive criteria."""
        return (
            self.au_count_3p
            + (1 if self.gc_in_range else 0)
            + (1 if self.a_at_pos19 else 0)
            + (1 if self.a_at_pos3 else 0)
            + (1 if self.u_at_pos10 else 0)
            + (1 if self.no_internal_repeat else 0)
        )

    @property
    def penalty(self) -> int:
        """Points deducted (0, 1 or 2)."""
        return (1 if self.g_at_pos13 else 0) + (1 if self.gc_at_pos19 else 0)

    @property
    def score(self) -> int:
        """Positive points minus penalties."""
        return self.positive_score - self.penalty

    @property
    def max_possible(self) -> int:
        """Highest attainable score: the 3' A/U window plus five binary criteria."""
        return REYNOLDS_AU_3P_WINDOW + 5

    def meets_threshold(self, threshold: int = REYNOLDS_DEFAULT_THRESHOLD) -> bool:
        """True when the score reaches threshold."""
        return self.score >= threshold

    def __repr__(self) -> str:
        return (
            f"ReynoldsResult(score={self.score}/{self.max_possible}, "
            f"meets_threshold={self.meets_threshold()})"
        )


def reynolds(
    candidate: SiRNACandidate,
    no_internal_repeat: bool | None = None,
) -> ReynoldsResult:
    """Score a candidate on its sense strand.

    The internal-repeat criterion comes from core.hairpin unless
    no_internal_repeat is given. Raises ValueError for cores under 19 nt.
    """
    sense = candidate.sense.seq
    if len(sense) < REYNOLDS_MIN_LENGTH:
        raise ValueError(
            f"Reynolds criteria need a core of at least {REYNOLDS_MIN_LENGTH} nt; "
            f"got {len(sense)} nt."
        )

    if no_internal_repeat is None:
        no_internal_repeat = has_no_internal_repeat(sense)

    gc = candidate.gc_content
    gc_in_range = REYNOLDS_GC_MIN <= gc <= REYNOLDS_GC_MAX

    tail = sense[-REYNOLDS_AU_3P_WINDOW:]
    au_count_3p = sum(base in _AU for base in tail)

    three_prime = sense[-1]
    a_at_pos19 = three_prime == "A"
    gc_at_pos19 = three_prime in _GC

    a_at_pos3 = _base_at(sense, REYNOLDS_POS_A_INTERNAL) == "A"
    u_at_pos10 = _base_at(sense, REYNOLDS_POS_U) == "U"
    g_at_pos13 = _base_at(sense, REYNOLDS_POS_G_PENALTY) == "G"

    return ReynoldsResult(
        gc_in_range=gc_in_range,
        no_internal_repeat=no_internal_repeat,
        a_at_pos19=a_at_pos19,
        a_at_pos3=a_at_pos3,
        u_at_pos10=u_at_pos10,
        au_count_3p=au_count_3p,
        g_at_pos13=g_at_pos13,
        gc_at_pos19=gc_at_pos19,
    )


# Amarzguioui scoring: docs/design-rules.md (The three rule sets). Motifs at 1,
# 6 and 19 use absolute positions; DT3 follows the 3' end (see the open item).
AMARZGUIOUI_MIN_LENGTH = 19
AMARZGUIOUI_TERMINAL_BP = 3

AMARZGUIOUI_POS_S1 = 1
AMARZGUIOUI_POS_A6 = 6
AMARZGUIOUI_POS_W19 = 19

AMARZGUIOUI_DIFFERENTIAL_MIN = -AMARZGUIOUI_TERMINAL_BP
AMARZGUIOUI_DIFFERENTIAL_MAX = AMARZGUIOUI_TERMINAL_BP

# Inclusive cut-off; the authors describe it as conservative.
AMARZGUIOUI_DEFAULT_THRESHOLD = 3


@dataclass(frozen=True)
class AmarzguiouiResult:
    """Components of the Amarzguioui score; u1 and g19 are negative motifs."""

    end_au_differential: int
    s1_gc: bool
    a6_a: bool
    w19_au: bool
    u1: bool
    g19: bool

    @property
    def motif_score(self) -> int:
        """Sum of the five motif contributions."""
        return (
            (1 if self.s1_gc else 0)
            + (1 if self.a6_a else 0)
            + (1 if self.w19_au else 0)
            + (-1 if self.u1 else 0)
            + (-1 if self.g19 else 0)
        )

    @property
    def score(self) -> int:
        """Motif points plus the end differential."""
        return self.motif_score + self.end_au_differential

    def meets_threshold(
        self, threshold: int = AMARZGUIOUI_DEFAULT_THRESHOLD
    ) -> bool:
        """True when the score reaches threshold."""
        return self.score >= threshold

    def __repr__(self) -> str:
        return (
            f"AmarzguiouiResult(score={self.score}, "
            f"motif_score={self.motif_score}, "
            f"end_au_differential={self.end_au_differential:+d})"
        )


def amarzguioui(candidate: SiRNACandidate) -> AmarzguiouiResult:
    """Score a candidate on its sense strand. Raises ValueError for cores under 19 nt."""
    sense = candidate.sense.seq
    if len(sense) < AMARZGUIOUI_MIN_LENGTH:
        raise ValueError(
            f"Amarzguioui criteria need a core of at least "
            f"{AMARZGUIOUI_MIN_LENGTH} nt; got {len(sense)} nt."
        )

    pos1 = _base_at(sense, AMARZGUIOUI_POS_S1)
    pos6 = _base_at(sense, AMARZGUIOUI_POS_A6)
    pos19 = _base_at(sense, AMARZGUIOUI_POS_W19)

    s1_gc = pos1 in _GC
    a6_a = pos6 == "A"
    w19_au = pos19 in _AU
    u1 = pos1 == "U"
    g19 = pos19 == "G"

    au_5p = sum(base in _AU for base in sense[:AMARZGUIOUI_TERMINAL_BP])
    au_3p = sum(base in _AU for base in sense[-AMARZGUIOUI_TERMINAL_BP:])
    end_au_differential = au_3p - au_5p

    return AmarzguiouiResult(
        end_au_differential=end_au_differential,
        s1_gc=s1_gc,
        a6_a=a6_a,
        w19_au=w19_au,
        u1=u1,
        g19=g19,
    )
