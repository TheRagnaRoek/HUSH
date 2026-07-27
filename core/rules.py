"""Classical siRNA design rules (stage 2 rule pre-screen).

Three empirical 2004 rule sets, each encoding the same thermodynamic-asymmetry
insight (the duplex should unwind more easily at the guide 5' end) differently:

* Ui-Tei      - a strict 4-rule pass/fail filter.
* Reynolds    - an additive score out of 10 (binary criteria, an additive
  3'-terminal A/U count, and two penalties).
* Amarzguioui - an additive criterion score: five motifs plus the duplex-end A/U
  differential.

Strand conventions match core.sirna: sense/passenger and antisense/guide are both
stored 5'->3'; each rule is evaluated on whichever strand it is naturally defined
for. The exact scoring of each set - and where Reynolds and Amarzguioui differ
from the common simplifications, plus the Table-3 validation - is in
docs/design_notes.md (rules section).

Primary sources: Ui-Tei et al. 2004; Reynolds et al. 2004; Amarzguioui & Prydz
2004.
"""

from __future__ import annotations

from dataclasses import dataclass

from core.hairpin import has_no_internal_repeat
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

# Reynolds is scored OUT OF 10 (not eight equal binary criteria), on the sense
# strand 5'->3': five +1 binary criteria, +1 per A/U base at positions 15-19
# (additive), and -1 for 'G' at 13 and for 'G'/'C' at 19; score runs -2..+10,
# and >= 6 predicts functional. Longer cores anchor the 3' parts to the actual
# terminus and keep the internal ones at absolute positions. Full scoring and the
# position-19 gradient: docs/design_notes.md (rules section).
REYNOLDS_MIN_LENGTH = 19

# GC content window (fraction, inclusive bounds).
REYNOLDS_GC_MIN = 0.30
REYNOLDS_GC_MAX = 0.52

# 3'-terminal A/U: +1 per A/U base across the 3'-terminal 5 nt (positions 15-19
# on a 19-mer). Additive, not a threshold.
REYNOLDS_AU_3P_WINDOW = 5

# Positional parts, 1-based. The internal ones are counted from the 5' end.
REYNOLDS_POS_A_INTERNAL = 3   # 'A' at position 3   (+1)
REYNOLDS_POS_U = 10           # 'U' at position 10  (+1)
REYNOLDS_POS_G_PENALTY = 13   # 'G' at position 13  (-1)
# Position 19 is handled as the 3' terminus: 'A' -> +1, 'G'/'C' -> -1.

# Score at or above which Reynolds predicts a functional siRNA: 6 of a possible
# 10. The cutoff is unchanged whether the internal-repeat criterion is computed
# by HUSH (the default) or supplied from an external folder.
REYNOLDS_DEFAULT_THRESHOLD = 6


def _base_at(seq: str, pos_1based: int) -> str:
    """Return the residue at a 1-based position from the 5' end."""
    return seq[pos_1based - 1]


# --------------------------------------------------------------------------- #
# Reynolds result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class ReynoldsResult:
    """Outcome of the Reynolds (2004) scoring for one candidate.

    Unlike Ui-Tei this is a graded score, not a gate, and it runs from -2 to
    +10: five +1 binary criteria, an additive 3'-terminal A/U count worth up to
    +5, and two -1 penalties. A total of at least REYNOLDS_DEFAULT_THRESHOLD
    predicts a functional siRNA.

    The 'absence of internal repeat' criterion is a hairpin-Tm calculation
    (core.hairpin); reynolds() computes it by default but it can be supplied from
    an external folder instead.

    Fields:
        gc_in_range:       GC content within [30%, 52%]                    (+1).
        no_internal_repeat: no stable internal hairpin (Tm <= 20 C)        (+1).
        a_at_pos19:        'A' at the sense 3' terminus (position 19)      (+1).
        a_at_pos3:         'A' at position 3                               (+1).
        u_at_pos10:        'U' at position 10                              (+1).
        au_count_3p:       number of A/U bases in the 3'-terminal 5 nt  (0..+5).
        g_at_pos13:        'G' at position 13 (penalty, -1 when True).
        gc_at_pos19:       'G'/'C' at the 3' terminus (penalty, -1 when True).
    """

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
        """Points earned from the positive criteria (excludes penalties)."""
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
        """Points deducted by the two negative criteria (0, 1 or 2)."""
        return (1 if self.g_at_pos13 else 0) + (1 if self.gc_at_pos19 else 0)

    @property
    def score(self) -> int:
        """Net Reynolds score: positive points minus penalties."""
        return self.positive_score - self.penalty

    @property
    def max_possible(self) -> int:
        """Best attainable score: the additive 3' window (5) plus five binary
        positives; penalties never raise the ceiling."""
        return REYNOLDS_AU_3P_WINDOW + 5

    def meets_threshold(self, threshold: int = REYNOLDS_DEFAULT_THRESHOLD) -> bool:
        """True when 'score' reaches 'threshold' (default 6)."""
        return self.score >= threshold

    def __repr__(self) -> str:
        return (
            f"ReynoldsResult(score={self.score}/{self.max_possible}, "
            f"meets_threshold={self.meets_threshold()})"
        )


# --------------------------------------------------------------------------- #
# Reynolds evaluation
# --------------------------------------------------------------------------- #


def reynolds(
    candidate: SiRNACandidate,
    no_internal_repeat: bool | None = None,
) -> ReynoldsResult:
    """Evaluate the Reynolds score for a candidate (sense strand, 5'->3').

    The 'absence of internal repeat' criterion is computed from the sense strand
    via core.hairpin (most stable internal hairpin melting at Tm <= 20 C). Pass
    'no_internal_repeat' explicitly to override that with an external folder's
    verdict.

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


# --------------------------------------------------------------------------- #
# Amarzguioui constants
# --------------------------------------------------------------------------- #
#
# Additive criterion score (sense strand, 1-based): positive motifs S1 (G/C at 1),
# A6 (A at 6), W19 (A/U at 19) each +1; negative motifs U1 (U at 1), G19 (G at 19)
# each -1; plus the raw duplex-end A/U differential ("DT3"), so
# score = motif points + differential. DT3 is measured 3' END MINUS 5' END over
# the three terminal base pairs - the sign is easy to invert; functional siRNAs
# score POSITIVE. All 46 siRNAs of the paper's Table 3 reproduce exactly (DT3,
# motifs, total). Derivation, sign rationale and validation: docs/design_notes.md
# (rules section).

AMARZGUIOUI_MIN_LENGTH = 19
AMARZGUIOUI_TERMINAL_BP = 3

# Motif positions, 1-based on the sense strand.
AMARZGUIOUI_POS_S1 = 1
AMARZGUIOUI_POS_A6 = 6
AMARZGUIOUI_POS_W19 = 19

# Bounds of the differential: each end contributes 0-3 A/U, so the difference
# spans -3 to +3 inclusive.
AMARZGUIOUI_DIFFERENTIAL_MIN = -AMARZGUIOUI_TERMINAL_BP
AMARZGUIOUI_DIFFERENTIAL_MAX = AMARZGUIOUI_TERMINAL_BP

# Score at or above which the paper calls an siRNA functional (inclusive). It
# reads as a minimally positive differential (+1) plus two of the three positive
# motifs. The authors call this cut-off deliberately conservative.
AMARZGUIOUI_DEFAULT_THRESHOLD = 3


# --------------------------------------------------------------------------- #
# Amarzguioui result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AmarzguiouiResult:
    """Outcome of the Amarzguioui & Prydz (2004) criterion-based score.

    Complete: the five motif criteria and the duplex-end differential are all
    evaluated here, and 'score' is their sum.

    Fields:
        end_au_differential: the paper's DT3 - A/U count among the three
            terminal base pairs at the sense 3' end MINUS that at the 5' end.
            Spans -3 to +3; positive favours functionality, and the value enters
            'score' directly rather than through a lookup.
        s1_gc:  G/C at sense position 1   (positive motif S1, +1).
        a6_a:   A at sense position 6     (positive motif A6, +1).
        w19_au: A/U at sense position 19  (positive motif W19, +1).
        u1:     U at sense position 1     (negative motif U1, -1 when True).
        g19:    G at sense position 19    (negative motif G19, -1 when True).
    """

    end_au_differential: int
    s1_gc: bool
    a6_a: bool
    w19_au: bool
    u1: bool
    g19: bool

    @property
    def motif_score(self) -> int:
        """Sum of the five motif contributions (+1 positive, -1 negative)."""
        return (
            (1 if self.s1_gc else 0)
            + (1 if self.a6_a else 0)
            + (1 if self.w19_au else 0)
            + (-1 if self.u1 else 0)
            + (-1 if self.g19 else 0)
        )

    @property
    def score(self) -> int:
        """Total score: motif points plus the duplex-end differential."""
        return self.motif_score + self.end_au_differential

    def meets_threshold(
        self, threshold: int = AMARZGUIOUI_DEFAULT_THRESHOLD
    ) -> bool:
        """True when 'score' reaches 'threshold' (default 3, inclusive)."""
        return self.score >= threshold

    def __repr__(self) -> str:
        return (
            f"AmarzguiouiResult(score={self.score}, "
            f"motif_score={self.motif_score}, "
            f"end_au_differential={self.end_au_differential:+d})"
        )


# --------------------------------------------------------------------------- #
# Amarzguioui evaluation
# --------------------------------------------------------------------------- #


def amarzguioui(candidate: SiRNACandidate) -> AmarzguiouiResult:
    """Evaluate the Amarzguioui & Prydz criteria for a candidate (sense strand).

    Computes the five motif criteria and the duplex-end A/U differential, and
    sums them into 'score'. Nothing is deferred.

    Args:
        candidate: the siRNA candidate to score.

    Raises:
        ValueError: if the sense strand is shorter than AMARZGUIOUI_MIN_LENGTH.
    """
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

    # DT3: 3' end minus 5' end, three terminal base pairs at each end.
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
