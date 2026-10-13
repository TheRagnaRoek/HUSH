"""Allele evaluation of candidate windows - the MSA stage's filter step.

Applies a user's AlleleTargeting to every window of the target and reports, per
window, whether it can hit what the user wants hit and discriminate against what
they want spared. Deliberately **not** a ranking, for the same reason
core.selection is not one: the signals here are counts of differences, and
turning them into a single weighted score would need efficacy data HUSH does not
have. Every window passes or does not, and keeps its full per-allele breakdown so
any verdict is legible.

Two gates, in the order a rejection is reported:

* the targeted alleles must be covered across the whole window and carry no
  difference from the target - the guide is built from the target, so a
  difference is a mismatch against the very allele it is meant to silence,
* every spared allele must be covered and carry at least
  ``min_spare_mismatches`` differences inside the window.

Coverage is required rather than optional. A window in which a targeted allele
has no residues says nothing about that allele, and one in which a spared allele
is missing cannot be claimed to discriminate against it; treating missing data as
either evidence or harmless would be a silent error, not a design choice.

Windows are indexed by their 0-based start on the target mRNA and their length -
the same coordinates core.sirna gives its candidates - so the two lists can be
joined on position without either side knowing about the other.

Scope and deferrals:

* Differences are counted, not weighted. Where a difference sits inside the
  duplex matters a great deal - a mismatch in the guide seed discriminates far
  more strongly than one at the 3' end - but weighting them needs the mismatch
  thermodynamics core.thermo does not yet carry. The positions are recorded so
  the weighting can be added without re-deriving anything.
* Insertions are recorded when their anchor falls inside a window, but they do
  not enter the verdict. An insertion would make the guide bulge on that allele,
  which is a strong discriminating feature; scoring it needs the same missing
  parameters.
* Positions are reported in **target-window** coordinates. Guide coordinates are
  available through guide_position(); see its docstring for why the two run in
  opposite directions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Iterator

from core.sirna import DEFAULT_LENGTH
from msa.discrimination import MismatchAssessment, assess_mismatch
from msa.discrimination import discrimination_key as _allele_discrimination_key
from msa.discrimination import ordered_by_placement as _ordered_by_placement
from msa.msa_errors import MsaError
from msa.targeting import AlleleTargeting
from msa.variation import PanelVariation, SiteState

# How many differences a spared allele must carry inside a window before the
# window is credited with discriminating against it. One is the minimum that
# means anything; a higher setting buys robustness against a single sequencing
# error or a single ambiguous alignment column.
DEFAULT_MIN_SPARE_MISMATCHES = 1


# --------------------------------------------------------------------------- #
# Criteria
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class WindowCriteria:
    """What the user asks of a window.

    Parameters:
        targeting: Which alleles must be hit and which spared.
        length: Window length in nucleotides. Defaults to the siRNA core length
            so windows line up with core.sirna's candidates one for one.
        min_spare_mismatches: Differences required per spared allele.
    """

    targeting: AlleleTargeting
    length: int = DEFAULT_LENGTH
    min_spare_mismatches: int = DEFAULT_MIN_SPARE_MISMATCHES

    def __post_init__(self) -> None:
        if self.length < 1:
            raise MsaError(f"Window length must be at least 1, got {self.length}.")
        if self.min_spare_mismatches < 1:
            raise MsaError(
                "min_spare_mismatches must be at least 1; requiring zero "
                "differences would credit every window with discriminating."
            )


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class WindowEvaluation:
    """The outcome of applying one WindowCriteria to one window.

    The dicts cover only the alleles the request constrains; an allele the
    request says nothing about is not examined and does not appear. Positions in
    them are 0-based **within the window**, not on the target.

    'satisfies' is the overall verdict; 'reasons' lists, in plain language, what
    a rejected window failed (empty when it satisfies).

    Parameters:
        position: 0-based start of the window on the target mRNA.
        length: Window length in nucleotides.
        mismatches: Allele id -> window positions where it differs.
        uncovered: Allele id -> window positions it has no residue for.
        insertions: Allele id -> window positions an insertion is anchored at.
            Recorded only; see the module's deferrals.
        max_hit_mismatches: Worst difference count among the targeted alleles.
            Zero when every targeted allele matches throughout.
        min_spare_mismatches: Weakest discrimination among the spared alleles -
            the count for whichever spared allele differs least. None when the
            request spares nothing.
        hit_pass: Whether every targeted allele is covered and identical.
        spare_pass: Whether every spared allele is covered and differs enough.
            None when the request spares nothing.
    """

    position: int
    length: int
    mismatches: dict[str, tuple[int, ...]]
    uncovered: dict[str, tuple[int, ...]]
    insertions: dict[str, tuple[int, ...]]
    max_hit_mismatches: int
    min_spare_mismatches: int | None
    hit_pass: bool
    spare_pass: bool | None
    satisfies: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)
    mismatch_assessments: dict[str, tuple[MismatchAssessment, ...]] = field(
        default_factory=dict
    )

    @property
    def end(self) -> int:
        """0-based, exclusive end of the window on the target."""
        return self.position + self.length

    @property
    def discrimination_key(self) -> tuple | None:
        """How well this window discriminates, for ranking. Higher is better.

        The weakest spared allele decides, the same "weakest link" rule
        min_spare_mismatches follows: a window is only as discriminating as the
        allele it distinguishes least well.

        None when the request spares nothing - pan-allele targeting asks for no
        discrimination, so there is none to rank. Callers that sort must decide
        what to do with that rather than being handed a misleading zero.

        See msa.discrimination for what the entries rest on.
        """
        if not self.mismatch_assessments:
            return None
        return min(
            _allele_discrimination_key(assessments)
            for assessments in self.mismatch_assessments.values()
        )

    @property
    def spared_alleles_by_discrimination(
        self,
    ) -> list[tuple[str, tuple[MismatchAssessment, ...]]]:
        """Every spared allele with its mismatches, weakest discrimination first.

        The first entry is the allele that sets the window's rank (the weakest
        link discrimination_key reports); the rest follow in improving order. A
        renderer that wants the whole picture - not just the deciding allele -
        walks this. Empty when the request spares nothing.
        """
        return sorted(
            self.mismatch_assessments.items(),
            key=lambda item: _allele_discrimination_key(item[1]),
        )

    @property
    def deciding_mismatches(self) -> tuple[MismatchAssessment, ...]:
        """Every difference to the weakest spared allele, best-placed first.

        The rank rests on the first element alone (deciding_mismatch), so
        without this the remaining differences would never reach the user. They
        matter: a window with three mediocre differences may well be the safer
        design than one with a single ideal difference, and msa.discrimination
        deliberately does not decide that - no source measures it. Showing the
        whole set is how that decision is handed back.

        Empty when nothing is spared.
        """
        if not self.mismatch_assessments:
            return ()
        weakest = min(
            self.mismatch_assessments.values(), key=_allele_discrimination_key
        )
        return _ordered_by_placement(weakest)

    @property
    def deciding_mismatch(self) -> MismatchAssessment | None:
        """The single mismatch this window's rank rests on.

        The best-placed difference from the weakest spared allele - the one that
        sets discrimination_key. Reporting it lets a user see *why* a window
        ranks where it does without unpacking every assessment; see
        deciding_mismatches for the rest of them.

        None when nothing is spared, or when the weakest allele has no
        difference at all (in which case the window does not discriminate and
        will not have been selected).
        """
        ordered = self.deciding_mismatches
        return ordered[0] if ordered else None

    @property
    def region(self) -> str:
        """1-based inclusive coordinate label, matching SiRNACandidate.region."""
        return f"{self.position + 1}-{self.end}"

    def __repr__(self) -> str:
        verdict = "satisfies" if self.satisfies else "rejected"
        return (
            f"WindowEvaluation(region={self.region}, {verdict}"
            + ("" if self.satisfies else f", reasons={list(self.reasons)}")
            + ")"
        )


# --------------------------------------------------------------------------- #
# Coordinate conversion
# --------------------------------------------------------------------------- #


def guide_position(window_position: int, length: int) -> int:
    """Convert a target-window position to a position on the guide strand.

    The guide is the reverse complement of the target window, so the two run in
    opposite directions: window position 0 pairs with the guide's 3' end, and
    the guide's 5' end pairs with the last window position. Both are 0-based.

    This matters because the seed - guide positions 1..7 counting from 0, the
    published 2-8 counting from 1 - therefore sits at the **3' end of the target
    window**, which is the opposite of where one expects it. The conversion is
    kept as its own function so the single place it can be got backwards is
    covered by its own tests.

    Args:
        window_position: 0-based position inside the target window.
        length: Window length in nucleotides.

    Returns:
        The 0-based position on the guide strand, 5'->3'.

    Raises:
        IndexError: When window_position lies outside the window.
    """
    if not 0 <= window_position < length:
        raise IndexError(
            f"Window position {window_position} outside 0..{length - 1}"
        )
    return length - 1 - window_position


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def evaluate_window(
    variation: PanelVariation, criteria: WindowCriteria, position: int
) -> WindowEvaluation:
    """Apply 'criteria' to the window starting at 'position'.

    Only the alleles the request constrains are examined; on a panel of
    thousands, a request over a handful of alleles costs a handful of slices.
    """
    _check_window(variation, criteria, position)
    targeting = criteria.targeting
    end = position + criteria.length

    mismatches: dict[str, tuple[int, ...]] = {}
    uncovered: dict[str, tuple[int, ...]] = {}
    insertions: dict[str, tuple[int, ...]] = {}

    for allele in sorted(targeting.constrained):
        states = variation.states_for(allele)[position:end]
        mismatches[allele] = _positions_of(states, SiteState.MISMATCH)
        uncovered[allele] = _positions_of(states, SiteState.UNCOVERED)
        anchors = variation.alignment[allele].insertions
        insertions[allele] = tuple(
            anchor - position for anchor in sorted(anchors) if position < anchor < end
        )

    hit_pass = all(
        not mismatches[allele] and not uncovered[allele] for allele in targeting.hit
    )
    max_hit_mismatches = max(
        (len(mismatches[allele]) for allele in targeting.hit), default=0
    )

    if targeting.spare:
        spare_counts = {a: len(mismatches[a]) for a in targeting.spare}
        min_spare_mismatches: int | None = min(spare_counts.values())
        spare_pass: bool | None = all(
            not uncovered[allele]
            and spare_counts[allele] >= criteria.min_spare_mismatches
            for allele in targeting.spare
        )
    else:
        min_spare_mismatches = None
        spare_pass = None

    satisfies = hit_pass and (spare_pass is not False)
    reasons = (
        ()
        if satisfies
        else _reasons(criteria, mismatches, uncovered, hit_pass, spare_pass)
    )

    return WindowEvaluation(
        position=position,
        length=criteria.length,
        mismatches=mismatches,
        uncovered=uncovered,
        insertions=insertions,
        max_hit_mismatches=max_hit_mismatches,
        min_spare_mismatches=min_spare_mismatches,
        hit_pass=hit_pass,
        spare_pass=spare_pass,
        satisfies=satisfies,
        reasons=reasons,
        mismatch_assessments=_assess(variation, criteria, position, mismatches),
    )


def _assess(
    variation: PanelVariation,
    criteria: WindowCriteria,
    position: int,
    mismatches: dict[str, tuple[int, ...]],
) -> dict[str, tuple[MismatchAssessment, ...]]:
    """Place every difference from a spared allele against the published evidence.

    Only spared alleles are assessed: discrimination is a statement about what an
    siRNA fails to silence, so a difference from a targeted allele is a defect,
    not a feature, and is already reported through hit_pass.
    """
    target_seq = variation.target.seq
    assessments: dict[str, tuple[MismatchAssessment, ...]] = {}
    for allele in sorted(criteria.targeting.spare):
        row = variation.alignment[allele].bases
        assessments[allele] = tuple(
            assess_mismatch(
                window_position=offset,
                window_length=criteria.length,
                target_base=target_seq[position + offset],
                spared_base=row[position + offset],
            )
            for offset in mismatches[allele]
        )
    return assessments


def evaluate_windows(
    variation: PanelVariation, criteria: WindowCriteria
) -> Iterator[WindowEvaluation]:
    """Lazily evaluate every window of the target, satisfying or not.

    Inputs are validated eagerly - before iteration begins - so a target shorter
    than the window, or a request naming an allele the alignment does not carry,
    raises at the call rather than partway through consumption. Same split as
    core.sirna.enumerate_candidates.
    """
    last_start = variation.width - criteria.length
    if last_start < 0:
        raise MsaError(
            f"Target is {variation.width} nt, shorter than the window length "
            f"{criteria.length}; no window can be evaluated."
        )
    _require_alleles_present(variation, criteria)
    return _window_stream(variation, criteria, last_start)


def _window_stream(
    variation: PanelVariation, criteria: WindowCriteria, last_start: int
) -> Iterator[WindowEvaluation]:
    """Generator over the windows; assumes inputs are validated."""
    for position in range(last_start + 1):
        yield evaluate_window(variation, criteria, position)


def select(
    variation: PanelVariation, criteria: WindowCriteria
) -> list[WindowEvaluation]:
    """Every window that satisfies the request, in target-position order.

    For the full set of evaluations, satisfying or not, map evaluate_windows
    over the variation directly.
    """
    return [
        evaluation
        for evaluation in evaluate_windows(variation, criteria)
        if evaluation.satisfies
    ]


def regions(evaluations: Iterable[WindowEvaluation]) -> Iterator[tuple[int, int]]:
    """Merge satisfying windows into contiguous target spans.

    Takes evaluations in target order and yields (start, end) pairs, 0-based and
    end-exclusive, covering every stretch of the target that satisfying windows
    reach. Overlapping windows are merged; a gap between them ends a region.
    """
    current: list[int] | None = None
    for evaluation in evaluations:
        if not evaluation.satisfies:
            continue
        if current is not None and evaluation.position <= current[1]:
            current[1] = max(current[1], evaluation.end)
            continue
        if current is not None:
            yield current[0], current[1]
        current = [evaluation.position, evaluation.end]
    if current is not None:
        yield current[0], current[1]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _check_window(
    variation: PanelVariation, criteria: WindowCriteria, position: int
) -> None:
    """Reject a window that does not fit on the target, or an unknown allele."""
    if position < 0 or position + criteria.length > variation.width:
        raise IndexError(
            f"Window {position}..{position + criteria.length} does not fit a "
            f"target of {variation.width} nt."
        )
    _require_alleles_present(variation, criteria)


def _require_alleles_present(
    variation: PanelVariation, criteria: WindowCriteria
) -> None:
    """Reject a request naming an allele the alignment does not carry."""
    missing = sorted(criteria.targeting.constrained - set(variation.ids))
    if missing:
        raise MsaError(
            f"Allele(s) {missing} are named in the request but absent from the "
            "alignment."
        )


def _positions_of(states: str, state: SiteState) -> tuple[int, ...]:
    """Window positions carrying a given state."""
    return tuple(index for index, char in enumerate(states) if char == state.value)


def _reasons(
    criteria: WindowCriteria,
    mismatches: dict[str, tuple[int, ...]],
    uncovered: dict[str, tuple[int, ...]],
    hit_pass: bool,
    spare_pass: bool | None,
) -> tuple[str, ...]:
    """Human-readable list of what a rejected window failed."""
    targeting = criteria.targeting
    out: list[str] = []

    if not hit_pass:
        differing = sorted(a for a in targeting.hit if mismatches[a])
        missing = sorted(a for a in targeting.hit if uncovered[a])
        if differing:
            out.append(f"targeted allele(s) differ here: {differing}")
        if missing:
            out.append(f"targeted allele(s) not covered here: {missing}")

    if spare_pass is False:
        weak = sorted(
            a
            for a in targeting.spare
            if not uncovered[a] and len(mismatches[a]) < criteria.min_spare_mismatches
        )
        missing = sorted(a for a in targeting.spare if uncovered[a])
        if weak:
            out.append(
                f"fewer than {criteria.min_spare_mismatches} difference(s) from "
                f"spared allele(s): {weak}"
            )
        if missing:
            out.append(f"spared allele(s) not covered here: {missing}")

    return tuple(out)
