"""Classical siRNA design rules (pipeline stage 2, rule pre-screen).

Three empirical rule sets from 2004 form the backbone of rational siRNA design.
They all encode thermodynamic asymmetry - the duplex should unwind more easily
at the guide (antisense) 5' end so RISC loads the guide rather than the
passenger - but express it differently:

* Ui-Tei    - a strict 4-rule pass/fail filter.
* Reynolds  - an 8-point additive checklist.
* Amarzguioui - a compact weighted regression score (later unit).

Strand conventions match core.sirna: sense/passenger and antisense/guide are
both stored 5'->3'. Rules are evaluated on whichever strand they are naturally
defined for and documented per rule.

References:
    Ui-Tei, K. et al. (2004). Guidelines for the selection of highly effective
    siRNA sequences for mammalian and chick RNA interference.
    Nucleic Acids Research 32(3), 936-948.
    Reynolds, A. et al. (2004). Rational siRNA design for RNA interference.
    Nature Biotechnology 22, 326-330.
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


# --------------------------------------------------------------------------- #
# Reynolds constants
# --------------------------------------------------------------------------- #

# All Reynolds criteria are evaluated on the sense strand, 5'->3'. The original
# criteria are defined for a 19-mer; for a longer core (20-24 nt) we anchor the
# 3'-end rules (II, IV, VII) to the actual 3' terminus and keep the internal
# rules (V, VI, VIII) at their absolute 1-based positions from the 5' end.
REYNOLDS_MIN_LENGTH = 19

# Criterion I: GC content window (fraction, inclusive bounds).
REYNOLDS_GC_MIN = 0.30
REYNOLDS_GC_MAX = 0.52

# Criterion II: at least this many A/U in the 3'-terminal window.
REYNOLDS_AU_3P_WINDOW = 5
REYNOLDS_MIN_AU_3P = 3

# Internal positional criteria (1-based positions from the 5' end).
REYNOLDS_POS_A_INTERNAL = 3   # V:    'A' at position 3
REYNOLDS_POS_U = 10           # VI:   'U' at position 10
REYNOLDS_POS_NOT_G = 13       # VIII: not 'G' at position 13

# Score at or above which Reynolds predicts a functional siRNA. The classic
# cutoff is 6 (out of 8). While criterion III is deferred the effective maximum
# is 7 (see ReynoldsResult.max_possible); the cutoff itself is unchanged.
REYNOLDS_DEFAULT_THRESHOLD = 6


def _base_at(seq: str, pos_1based: int) -> str:
    """Return the residue at a 1-based position from the 5' end."""
    return seq[pos_1based - 1]


# --------------------------------------------------------------------------- #
# Reynolds result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ReynoldsResult:
    """Outcome of the eight Reynolds criteria for one candidate.

    Unlike Ui-Tei this is a graded checklist, not a gate: each satisfied
    criterion contributes one point and a total of at least
    REYNOLDS_DEFAULT_THRESHOLD predicts a functional siRNA.

    Criterion III (low internal secondary structure) needs a folding/Tm
    calculation that lives in the later thermodynamics unit. Until it is
    supplied it is left as None ("not assessed") rather than assumed pass or
    fail: it is excluded from 'score' and 'max_possible' and named in 'pending',
    so the checklist never silently under- or over-counts.

    Fields (Roman numerals are the original Reynolds criteria):
        gc_in_range:      I    GC content within [30%, 52%].
        au_rich_3p:       II   >= 3 A/U in the sense 3'-terminal 5 nt.
        low_internal_structure: III  low internal secondary structure (deferred).
        a_at_3p_end:      IV   'A' at the sense 3' terminus.
        a_at_pos3:        V    'A' at position 3.
        u_at_pos10:       VI   'U' at position 10.
        not_gc_at_3p_end: VII  sense 3' terminus is not G/C.
        not_g_at_pos13:   VIII position 13 is not 'G'.
    """

    gc_in_range: bool
    au_rich_3p: bool
    low_internal_structure: bool | None
    a_at_3p_end: bool
    a_at_pos3: bool
    u_at_pos10: bool
    not_gc_at_3p_end: bool
    not_g_at_pos13: bool

    def _criteria(self) -> list[tuple[str, bool | None]]:
        return [
            ("gc_in_range", self.gc_in_range),
            ("au_rich_3p", self.au_rich_3p),
            ("low_internal_structure", self.low_internal_structure),
            ("a_at_3p_end", self.a_at_3p_end),
            ("a_at_pos3", self.a_at_pos3),
            ("u_at_pos10", self.u_at_pos10),
            ("not_gc_at_3p_end", self.not_gc_at_3p_end),
            ("not_g_at_pos13", self.not_g_at_pos13),
        ]

    @property
    def pending(self) -> list[str]:
        """Names of criteria not yet assessed (value is None)."""
        return [name for name, value in self._criteria() if value is None]

    @property
    def score(self) -> int:
        """Points earned among the assessed criteria."""
        return sum(1 for _, value in self._criteria() if value is True)

    @property
    def max_possible(self) -> int:
        """Number of criteria actually assessed (excludes pending ones)."""
        return sum(1 for _, value in self._criteria() if value is not None)

    def meets_threshold(self, threshold: int = REYNOLDS_DEFAULT_THRESHOLD) -> bool:
        """True when 'score' reaches 'threshold' (default 6)."""
        return self.score >= threshold

    def __repr__(self) -> str:
        pending = ",".join(self.pending) or "none"
        return (
            f"ReynoldsResult(score={self.score}/{self.max_possible}, "
            f"meets_threshold={self.meets_threshold()}, pending={pending})"
        )


# --------------------------------------------------------------------------- #
# Reynolds evaluation
# --------------------------------------------------------------------------- #


def reynolds(
    candidate: SiRNACandidate,
    low_internal_structure: bool | None = None,
) -> ReynoldsResult:
    """Evaluate the Reynolds criteria for a candidate (sense strand, 5'->3').

    Criterion III (low internal secondary structure) is not computed here; pass
    'low_internal_structure' once the thermodynamics unit can assess it. Left as
    None it is reported as pending and excluded from the score.

    Raises:
        ValueError: if the sense strand is shorter than REYNOLDS_MIN_LENGTH; the
        criteria are only defined for a standard-length (>= 19 nt) core.
    """
    sense = candidate.sense.seq
    if len(sense) < REYNOLDS_MIN_LENGTH:
        raise ValueError(
            f"Reynolds criteria need a core of at least {REYNOLDS_MIN_LENGTH} nt; "
            f"got {len(sense)} nt."
        )

    gc = candidate.gc_content
    gc_in_range = REYNOLDS_GC_MIN <= gc <= REYNOLDS_GC_MAX

    tail = sense[-REYNOLDS_AU_3P_WINDOW:]
    au_rich_3p = sum(base in _AU for base in tail) >= REYNOLDS_MIN_AU_3P

    three_prime = sense[-1]
    a_at_3p_end = three_prime == "A"
    not_gc_at_3p_end = three_prime not in _GC

    a_at_pos3 = _base_at(sense, REYNOLDS_POS_A_INTERNAL) == "A"
    u_at_pos10 = _base_at(sense, REYNOLDS_POS_U) == "U"
    not_g_at_pos13 = _base_at(sense, REYNOLDS_POS_NOT_G) != "G"

    return ReynoldsResult(
        gc_in_range=gc_in_range,
        au_rich_3p=au_rich_3p,
        low_internal_structure=low_internal_structure,
        a_at_3p_end=a_at_3p_end,
        a_at_pos3=a_at_pos3,
        u_at_pos10=u_at_pos10,
        not_gc_at_3p_end=not_gc_at_3p_end,
        not_g_at_pos13=not_g_at_pos13,
    )
