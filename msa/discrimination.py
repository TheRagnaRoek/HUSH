"""Placement of a mismatch against published discrimination data.

Each difference between the target and a spared allele is assessed with Huang
2009 (Table 1, NAR 37:7560), Schwarz 2006 (PLoS Genet 2:e140) and Du 2005 (low
tolerance region). The key is a tuple of booleans compared lexicographically,
not a weighted score; the number of differences is not in it. Ranking never
removes a candidate. See docs/design-rules.md (Discrimination ranking, including
its limits).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence

from core.seed import SEED_END, SEED_START
from msa.msa_errors import MsaError

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CTU")

_COMPLEMENT = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}


def _to_dna(base: str) -> str:
    """Uppercase, with U read as T."""
    base = base.upper()
    return "T" if base == "U" else base


def _is_unambiguous(base: str) -> bool:
    return _to_dna(base) in {"A", "C", "G", "T"}


# Huang 2009, Table 1, verbatim. A row is (spared base, targeted base); the guide
# carries the complement of the targeted base and so mismatches the spared one.
# Positions are 1-based guide positions where that mismatch was least tolerated.
HUANG_RECOMMENDED_POSITIONS: dict[tuple[str, str], frozenset[int]] = {
    ("A", "C"): frozenset({8, 9, 10, 11, 16}),
    ("A", "G"): frozenset({10, 11, 12, 16}),
    ("A", "T"): frozenset({9, 10, 14, 15, 16}),
    ("C", "A"): frozenset({11, 12, 13, 14, 16}),
    ("C", "G"): frozenset({9, 10, 12, 13, 16}),
    ("C", "T"): frozenset({9, 10, 11, 12, 14}),
    ("G", "A"): frozenset({10, 13, 15, 16}),
    ("G", "C"): frozenset({10, 11, 12, 13}),
    ("G", "T"): frozenset({9, 10, 11, 12, 14}),
    ("T", "A"): frozenset({9, 11, 12, 15}),
    ("T", "C"): frozenset({10, 11, 12, 13}),
    ("T", "G"): frozenset({10, 11, 12, 13, 16}),
}

# Schwarz 2006 positions that held up across their five figures: P13 and P16 in
# all five, P10 in four. P16 also passed their ten-pair generality test.
STRONGEST_POSITIONS = frozenset({10, 13, 16})

# Other positions Schwarz names, excluding the seed. Reported, not a ranking
# tier, since the window below already covers them.
CONSISTENT_POSITIONS = frozenset({9, 12, 14})

# Du 2005 low-tolerance region: target-site positions 5-11, which is guide
# P9-P15 (guide = 20 - target, checked against the paper's printed 19-mer pair).
DU_LOW_TOLERANCE = frozenset(range(9, 16))

# Ranking window: Du's region plus P16 (rationale in the docs).
LOW_TOLERANCE_WINDOW = DU_LOW_TOLERANCE | frozenset({16})

# Differences at which a window is flagged against the two-mismatch convention
# of Naito et al. 2009. Reported only; the count is never a ranking criterion.
MULTI_MISMATCH = 2


class MismatchClass(str, Enum):
    """Guide residue : spared residue, as the papers name it (siRNA nucleotide first).

    Purine:purine disrupts the helix most and discriminates best.
    """

    PURINE_PURINE = "purine:purine"
    PURINE_PYRIMIDINE = "purine:pyrimidine"
    PYRIMIDINE_PURINE = "pyrimidine:purine"
    PYRIMIDINE_PYRIMIDINE = "pyrimidine:pyrimidine"


class GuideRegion(str, Enum):
    """Region of the guide: ANCHOR (position 1), SEED or POST_SEED (3' of the seed)."""

    ANCHOR = "anchor"
    SEED = "seed"
    POST_SEED = "post_seed"


def guide_region(guide_position: int) -> GuideRegion:
    """Region of a 1-based guide position."""
    if guide_position < 1:
        raise MsaError(f"Guide positions are 1-based; got {guide_position}.")
    if guide_position == 1:
        return GuideRegion.ANCHOR
    if SEED_START <= guide_position <= SEED_END:
        return GuideRegion.SEED
    return GuideRegion.POST_SEED


def classify_mismatch(guide_base: str, other_base: str) -> MismatchClass | None:
    """Mismatch class of a guide/target base pair, or None if either is ambiguous."""
    if not (_is_unambiguous(guide_base) and _is_unambiguous(other_base)):
        return None
    guide_purine = _to_dna(guide_base) in PURINES
    other_purine = _to_dna(other_base) in PURINES
    if guide_purine and other_purine:
        return MismatchClass.PURINE_PURINE
    if guide_purine:
        return MismatchClass.PURINE_PYRIMIDINE
    if other_purine:
        return MismatchClass.PYRIMIDINE_PURINE
    return MismatchClass.PYRIMIDINE_PYRIMIDINE


def huang_recommends(
    spared_base: str, target_base: str, guide_position: int
) -> bool | None:
    """Whether Huang's Table 1 lists this position for this base change.

    None when either base is ambiguous or the two bases are equal.
    """
    if not (_is_unambiguous(spared_base) and _is_unambiguous(target_base)):
        return None
    key = (_to_dna(spared_base), _to_dna(target_base))
    positions = HUANG_RECOMMENDED_POSITIONS.get(key)
    if positions is None:  # identical bases: not a mismatch, not a table row
        return None
    return guide_position in positions


@dataclass(frozen=True)
class MismatchAssessment:
    """One difference between the target and a spared allele, placed in the evidence.

    window_position is 0-based in the window; guide_position is 1-based on the
    guide, 5'->3', and runs the opposite way. guide_base is the complement of
    target_base. The ranking is derived by discrimination_key, not stored.
    """

    window_position: int
    guide_position: int
    target_base: str
    spared_base: str
    guide_base: str
    huang_recommended: bool | None
    mismatch_class: MismatchClass | None
    region: GuideRegion

    @property
    def mutation_type(self) -> str:
        """The Huang Table 1 row, e.g. 'A>C'."""
        return f"{_to_dna(self.spared_base)}>{_to_dna(self.target_base)}"

    @property
    def sources_agree(self) -> bool:
        """Whether Huang's table and the Schwarz regional finding agree.

        False for a Huang position inside the seed, or an omitted one 3' of it.
        """
        if self.huang_recommended is None:
            return False
        return self.huang_recommended == (self.region is GuideRegion.POST_SEED)

    def __repr__(self) -> str:
        cls = "?" if self.mismatch_class is None else self.mismatch_class.value
        return (
            f"MismatchAssessment(P{self.guide_position}, {self.mutation_type}, "
            f"{cls}, {self.region.value}, huang={self.huang_recommended})"
        )


def guide_position_1based(window_position: int, window_length: int) -> int:
    """1-based guide position paired with a 0-based window position.

    The 0-based counterpart is msa.windows.guide_position. Raises MsaError for a
    position outside the window.
    """
    if not 0 <= window_position < window_length:
        raise MsaError(
            f"Window position {window_position} outside 0..{window_length - 1}"
        )
    return window_length - window_position


def assess_mismatch(
    window_position: int, window_length: int, target_base: str, spared_base: str
) -> MismatchAssessment:
    """Assess one window-relative difference; MsaError for a position outside the window."""
    position = guide_position_1based(window_position, window_length)
    guide_base = _COMPLEMENT.get(_to_dna(target_base), target_base.upper())
    return MismatchAssessment(
        window_position=window_position,
        guide_position=position,
        target_base=target_base.upper(),
        spared_base=spared_base.upper(),
        guide_base=guide_base,
        huang_recommended=huang_recommends(spared_base, target_base, position),
        mismatch_class=classify_mismatch(guide_base, spared_base),
        region=guide_region(position),
    )


def _single_key(assessment: MismatchAssessment) -> tuple[bool, bool, bool, bool, bool]:
    """The five ranking booleans of one mismatch, in the order of the docs table."""
    return (
        assessment.huang_recommended is True,
        assessment.guide_position in LOW_TOLERANCE_WINDOW,
        assessment.guide_position in STRONGEST_POSITIONS,
        assessment.region is GuideRegion.POST_SEED,
        assessment.mismatch_class is MismatchClass.PURINE_PURINE,
    )


def best_mismatch(
    assessments: Sequence[MismatchAssessment],
) -> MismatchAssessment | None:
    """The most discriminating mismatch against one allele, or None."""
    if not assessments:
        return None
    return max(assessments, key=_single_key)


def ordered_by_placement(
    assessments: Sequence[MismatchAssessment],
) -> tuple[MismatchAssessment, ...]:
    """All mismatches against one allele, best-placed first."""
    return tuple(sorted(assessments, key=_single_key, reverse=True))


def discrimination_key(
    assessments: Sequence[MismatchAssessment],
) -> tuple[bool, bool, bool, bool, bool]:
    """Ranking key of one spared allele: its best-placed mismatch; higher is better.

    No mismatch at all gives the lowest key.
    """
    if not assessments:
        return (False, False, False, False, False)
    return _single_key(best_mismatch(assessments))
