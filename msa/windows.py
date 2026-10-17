"""Allele evaluation of candidate windows.

For each window of the target, two gates: targeted alleles must be covered and
identical to the target across the window, and every spared allele must be
covered with at least min_spare_mismatches differences. Missing coverage never
passes. Windows use the same 0-based coordinates as core.sirna candidates.
See docs/design-rules.md (Allele targeting). Differences are counted, not
weighted, and insertions are recorded but do not enter the verdict (Not
implemented).
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

# Differences a spared allele needs inside a window to count as discriminated.
DEFAULT_MIN_SPARE_MISMATCHES = 1


@dataclass(frozen=True)
class WindowCriteria:
    """Targeting request, window length (default: siRNA core length) and the
    differences required per spared allele."""

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


@dataclass(frozen=True)
class WindowEvaluation:
    """Verdict and per-allele breakdown for one window.

    The dicts cover only the alleles the request constrains, with positions
    0-based within the window. 'max_hit_mismatches' is the worst count among
    targeted alleles; 'min_spare_mismatches' the weakest among spared ones
    (None when nothing is spared). 'hit_pass' and 'spare_pass' are the two
    gates (spare_pass is None when nothing is spared). 'reasons' names what a
    rejected window failed.
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
        """0-based exclusive end of the window."""
        return self.position + self.length

    @property
    def discrimination_key(self) -> tuple | None:
        """Ranking key of the weakest spared allele (higher is better).

        None when nothing is spared. See msa.discrimination.
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
        """Spared alleles with their mismatches, weakest discrimination first."""
        return sorted(
            self.mismatch_assessments.items(),
            key=lambda item: _allele_discrimination_key(item[1]),
        )

    @property
    def deciding_mismatches(self) -> tuple[MismatchAssessment, ...]:
        """All differences to the weakest spared allele, best-placed first.

        The rank rests on the first one; the rest are shown because no source
        weighs several poor mismatches against one good one. Empty when nothing
        is spared.
        """
        if not self.mismatch_assessments:
            return ()
        weakest = min(
            self.mismatch_assessments.values(), key=_allele_discrimination_key
        )
        return _ordered_by_placement(weakest)

    @property
    def deciding_mismatch(self) -> MismatchAssessment | None:
        """The best-placed difference to the weakest spared allele, which sets the rank.

        None when nothing is spared or that allele has no difference.
        """
        ordered = self.deciding_mismatches
        return ordered[0] if ordered else None

    @property
    def region(self) -> str:
        """1-based inclusive label, as SiRNACandidate.region."""
        return f"{self.position + 1}-{self.end}"

    def __repr__(self) -> str:
        verdict = "satisfies" if self.satisfies else "rejected"
        return (
            f"WindowEvaluation(region={self.region}, {verdict}"
            + ("" if self.satisfies else f", reasons={list(self.reasons)}")
            + ")"
        )


def guide_position(window_position: int, length: int) -> int:
    """0-based guide position paired with a 0-based window position.

    The guide is the reverse complement of the window, so the directions are
    opposite: the seed lies at the 3' end of the target window. Raises
    IndexError for a position outside the window.
    """
    if not 0 <= window_position < length:
        raise IndexError(
            f"Window position {window_position} outside 0..{length - 1}"
        )
    return length - 1 - window_position


def evaluate_window(
    variation: PanelVariation, criteria: WindowCriteria, position: int
) -> WindowEvaluation:
    """Evaluate the window starting at position; only constrained alleles are read."""
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
    """Assess the differences to spared alleles only.

    A difference to a targeted allele is a defect, reported through hit_pass.
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
    """Lazily evaluate every window; inputs are validated before the first yield."""
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
    for position in range(last_start + 1):
        yield evaluate_window(variation, criteria, position)


def select(
    variation: PanelVariation, criteria: WindowCriteria
) -> list[WindowEvaluation]:
    """Windows that satisfy the request, in target order."""
    return [
        evaluation
        for evaluation in evaluate_windows(variation, criteria)
        if evaluation.satisfies
    ]


def regions(evaluations: Iterable[WindowEvaluation]) -> Iterator[tuple[int, int]]:
    """Merge satisfying windows, given in target order, into (start, end) spans.

    0-based, end-exclusive. Overlapping windows merge; a gap ends a region.
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


def _check_window(
    variation: PanelVariation, criteria: WindowCriteria, position: int
) -> None:
    """Reject a window outside the target or an allele not in the alignment."""
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
    """Plain-language list of what a rejected window failed."""
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
