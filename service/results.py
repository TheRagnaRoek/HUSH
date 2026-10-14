"""Structured results returned by the HUSH service layer.

Lossless carriers: each result holds the full computed state a renderer needs,
whether that renderer is the CLI (terminal tables and TSV) or the future HTTP
API (JSON). No formatting decisions live here - a result never knows about
terminals, colours, or wire formats.
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
    """A preview of a design input: the loaded target and its optional panel."""

    design_input: DesignInput

    @property
    def target(self) -> BioSequence:
        return self.design_input.target

    @property
    def has_panel(self) -> bool:
        return self.design_input.has_panel


@dataclass
class DesignResult:
    """The full outcome of one design run over the loaded target (+ panel).

    'evaluations' is every candidate with its per-criterion breakdown, in
    target-position order; 'selected' is the subset that cleared every enabled
    gate, already ordered per the request's sort. 'window_by_pos' and
    'targeting' are populated only for an allele-aware (panel) run.
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
    """The stem-loop shRNA constructs built from a set of siRNA candidates.

    'constructs' is what survived conversion; 'dropped' is every construct
    withheld by a conversion filter, kept so a renderer can say what went
    and why rather than silently returning a shorter list.
    """

    constructs: list[ShrnaConstruct]
    options: ShrnaOptions
    dropped: list[ShrnaConstruct] = field(default_factory=list)

    @property
    def total_converted(self) -> int:
        """Constructs built, before any conversion filter."""
        return len(self.constructs) + len(self.dropped)
