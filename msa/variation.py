"""Per-position comparison of an aligned panel against the target.

Each member gets a string of SiteState characters, one per target position.
UNCOVERED is kept apart from MISMATCH: missing data is never evidence
(docs/design-rules.md, Allele targeting). Nothing is scored here. Residues are
compared as plain characters and insertions are not examined (docs/design-rules.md,
Not implemented).
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator

from core.sequence import BioSequence
from msa.alignment import GAP, PanelAlignment
from msa.msa_errors import AlignmentError


class SiteState(str, Enum):
    """A member's state at one target position; values are the string characters.

    UNCOVERED: no residue (deletion, or the alignment does not reach).
    """

    MATCH = "="
    MISMATCH = "X"
    UNCOVERED = GAP


@dataclass(frozen=True)
class PanelVariation:
    """Every member compared to the target.

    states maps a member id to its state string, in panel order and target
    coordinates; windows are taken by slicing.
    """

    alignment: PanelAlignment
    states: dict[str, str]

    def __post_init__(self) -> None:
        object.__setattr__(self, "states", dict(self.states))
        self._validate()

    @classmethod
    def from_alignment(cls, alignment: PanelAlignment) -> "PanelVariation":
        """Compare every aligned row against the target."""
        target_seq = alignment.target.seq
        states = {
            row.id: _compare(target_seq, row.bases) for row in alignment
        }
        return cls(alignment=alignment, states=states)

    def _validate(self) -> None:
        expected = set(self.alignment.ids)
        got = set(self.states)
        if got != expected:
            missing = sorted(expected - got)
            extra = sorted(got - expected)
            raise AlignmentError(
                f"State rows do not match the alignment: missing {missing}, "
                f"unexpected {extra}."
            )
        allowed = {state.value for state in SiteState}
        for row_id, row in self.states.items():
            if len(row) != self.width:
                raise AlignmentError(
                    f"State row {row_id!r} has {len(row)} positions but the "
                    f"target has {self.width}."
                )
            bad = sorted(set(row) - allowed)
            if bad:
                raise AlignmentError(
                    f"State row {row_id!r} contains unknown state character(s) "
                    f"{bad}; expected one of {sorted(allowed)}."
                )

    @property
    def target(self) -> BioSequence:
        return self.alignment.target

    @property
    def width(self) -> int:
        """Number of target positions."""
        return self.alignment.width

    @property
    def ids(self) -> list[str]:
        """Member ids in panel order."""
        return self.alignment.ids

    def __len__(self) -> int:
        return len(self.states)

    def __contains__(self, member_id: str) -> bool:
        return member_id in self.states

    def states_for(self, member_id: str) -> str:
        """The full state row of one member."""
        try:
            return self.states[member_id]
        except KeyError:
            raise AlignmentError(
                f"No state row for {member_id!r}. Present: {sorted(self.states)}"
            ) from None

    def state_at(self, member_id: str, position: int) -> SiteState:
        """One member's state at one target position."""
        self._check_position(position)
        return SiteState(self.states_for(member_id)[position])

    def column(self, position: int) -> dict[str, SiteState]:
        """Every member's state at one target position."""
        self._check_position(position)
        return {
            member_id: SiteState(row[position])
            for member_id, row in self.states.items()
        }

    def _check_position(self, position: int) -> None:
        if not 0 <= position < self.width:
            raise IndexError(
                f"Target position {position} outside 0..{self.width - 1}"
            )

    def polymorphic_positions(self) -> Iterator[int]:
        """Positions where at least one member mismatches; uncovered does not count."""
        rows = list(self.states.values())
        for position in range(self.width):
            if any(row[position] == SiteState.MISMATCH for row in rows):
                yield position

    def conserved_positions(self) -> Iterator[int]:
        """Positions every member covers and matches."""
        rows = list(self.states.values())
        for position in range(self.width):
            if all(row[position] == SiteState.MATCH for row in rows):
                yield position

    def __repr__(self) -> str:
        return (
            f"PanelVariation(target={self.target.id!r}, width={self.width}, "
            f"members={len(self.states)})"
        )


def _compare(target_seq: str, member_bases: str) -> str:
    """One state character per position; both strings are in target coordinates."""
    states: list[str] = []
    for target_char, member_char in zip(target_seq, member_bases):
        if member_char == GAP:
            states.append(SiteState.UNCOVERED.value)
        elif member_char == target_char:
            states.append(SiteState.MATCH.value)
        else:
            states.append(SiteState.MISMATCH.value)
    return "".join(states)
