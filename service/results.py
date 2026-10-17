"""Results returned by the service layer.

They carry the full computed state and no formatting. See
docs/architecture.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from core.selection import CandidateEvaluation, SelectionCriteria
from core.sequence import BioSequence, DesignInput
from core.sirna import OverhangPolicy
from msa.targeting import AlleleTargeting
from msa.windows import DEFAULT_MIN_SPARE_MISMATCHES, WindowEvaluation
from shrna.convert import ShrnaConstruct, ShrnaOptions


@dataclass
class InspectResult:
    """A loaded target and its optional panel."""

    design_input: DesignInput

    @property
    def target(self) -> BioSequence:
        return self.design_input.target

    @property
    def has_panel(self) -> bool:
        return self.design_input.has_panel


@dataclass
class DesignResult:
    """Outcome of one design run.

    'evaluations' covers every candidate in position order; 'selected' is the
    subset that cleared every gate, in the requested sort order. 'targeting' and
    'window_by_pos' are set only for a run with a panel.
    """

    target: BioSequence
    length: int
    overhang: OverhangPolicy
    criteria: SelectionCriteria
    evaluations: list[CandidateEvaluation]
    selected: list[CandidateEvaluation]
    targeting: AlleleTargeting | None = None
    window_by_pos: dict[int, WindowEvaluation] | None = None
    min_spare_mismatches: int = DEFAULT_MIN_SPARE_MISMATCHES
    offtarget_skipped: bool = False


@dataclass
class ConvertResult:
    """shRNA constructs built from candidates.

    'dropped' holds the constructs withheld by a conversion filter, so a
    renderer can report them.
    """

    constructs: list[ShrnaConstruct]
    options: ShrnaOptions
    dropped: list[ShrnaConstruct] = field(default_factory=list)

    @property
    def total_converted(self) -> int:
        """Constructs built, before any conversion filter."""
        return len(self.constructs) + len(self.dropped)
