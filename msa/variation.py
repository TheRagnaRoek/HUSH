"""Per-position comparison of an aligned panel against the target (stage 3).

Turns a PanelAlignment into the question every later step actually asks: at each
target position, does this member carry the target's residue, a different one,
or none at all?

The three states are kept apart on purpose. "Uncovered" is not a weak form of
"different": a candidate window in which a member the user wants to *hit* is not
covered says nothing about that member, and a window in which a member to be
*spared* is not covered cannot be claimed to discriminate against it. Folding
the two together would silently turn missing data into evidence.

Layout mirrors PanelAlignment: one string per member, exactly as long as the
target, one character per position. Storing states as a string rather than a
list of objects keeps a whole panel cheap - an IPD-IMGT/HLA locus can hold
thousands of members over a coding sequence of a few thousand bases - and lets a
window be taken with a slice.

Nothing is scored or judged here. Whether a set of differences satisfies a
user's request is msa.windows' question; this module only reports what is there.

Scope and deferrals:

* Residues are compared as characters. IUPAC ambiguity codes are not resolved,
  so 'N' against 'A' counts as a difference and 'N' against 'N' as a match. A
  panel drawn from IPD-IMGT/HLA is unambiguous; a hand-assembled one may not be.
* Insertions carried by AlignedRow are not examined. They have no target
  position, so they cannot appear in a per-position view; a window that spans
  their anchor has to consult the row itself.
* Case is not normalised - BioSequence already uppercases on construction.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator

from core.sequence import BioSequence
from msa.alignment import GAP, PanelAlignment
from msa.msa_errors import AlignmentError

# --------------------------------------------------------------------------- #
# The three states
# --------------------------------------------------------------------------- #


class SiteState(str, Enum):
    """What one panel member does at one target position.

    The values double as the characters used in the packed state strings, so a
    state can be read straight out of a row with SiteState(row[position]).

    MATCH:     the member carries the target's residue.
    MISMATCH:  the member carries a different residue.
    UNCOVERED: the member has no residue here - a deletion, or the alignment
               does not reach. Not the same as a difference.
    """

    MATCH = "="
    MISMATCH = "X"
    UNCOVERED = GAP


# --------------------------------------------------------------------------- #
# The comparison
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PanelVariation:
    """Every panel member compared to the target, position by position.

    Parameters:
        alignment: The alignment this was derived from. Kept so a consumer can
            reach the residues themselves, the insertions, and the target.
        states: Sequence id -> one SiteState character per target position, in
            panel order and target coordinates.
    """

    alignment: PanelAlignment
    states: dict[str, str]

    def __post_init__(self) -> None:
        object.__setattr__(self, "states", dict(self.states))
        self._validate()

    # -- construction ------------------------------------------------------ #

    @classmethod
    def from_alignment(cls, alignment: PanelAlignment) -> "PanelVariation":
        """Compare every aligned row against the target."""
        target_seq = alignment.target.seq
        states = {
            row.id: _compare(target_seq, row.bases) for row in alignment
        }
        return cls(alignment=alignment, states=states)

    # -- validation -------------------------------------------------------- #

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

    # -- accessors --------------------------------------------------------- #

    @property
    def target(self) -> BioSequence:
        """The target every member was compared against."""
        return self.alignment.target

    @property
    def width(self) -> int:
        """Number of target positions."""
        return self.alignment.width

    @property
    def ids(self) -> list[str]:
        """Member ids, in panel order."""
        return self.alignment.ids

    def __len__(self) -> int:
        return len(self.states)

    def __contains__(self, member_id: str) -> bool:
        return member_id in self.states

    def states_for(self, member_id: str) -> str:
        """The full state row of one member, one character per target position.

        Slice it to get a window: ``variation.states_for(id)[start:start + n]``.
        """
        try:
            return self.states[member_id]
        except KeyError:
            raise AlignmentError(
                f"No state row for {member_id!r}. Present: {sorted(self.states)}"
            ) from None

    def state_at(self, member_id: str, position: int) -> SiteState:
        """What one member does at one target position."""
        self._check_position(position)
        return SiteState(self.states_for(member_id)[position])

    def column(self, position: int) -> dict[str, SiteState]:
        """What every member does at one target position, in panel order."""
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

    # -- queries ----------------------------------------------------------- #

    def polymorphic_positions(self) -> Iterator[int]:
        """Target positions where at least one member differs from the target.

        Uncovered positions do not count: missing data is not a difference.
        """
        rows = list(self.states.values())
        for position in range(self.width):
            if any(row[position] == SiteState.MISMATCH for row in rows):
                yield position

    def conserved_positions(self) -> Iterator[int]:
        """Target positions every member covers and matches.

        The complement of polymorphic_positions() only for fully covered panels;
        a position no member covers is neither conserved nor polymorphic.
        """
        rows = list(self.states.values())
        for position in range(self.width):
            if all(row[position] == SiteState.MATCH for row in rows):
                yield position

    def __repr__(self) -> str:
        return (
            f"PanelVariation(target={self.target.id!r}, width={self.width}, "
            f"members={len(self.states)})"
        )


# --------------------------------------------------------------------------- #
# The comparison itself
# --------------------------------------------------------------------------- #


def _compare(target_seq: str, member_bases: str) -> str:
    """One state character per target position, comparing residue by residue.

    Both strings are already in target coordinates and of equal length - the
    invariant PanelAlignment enforces - so this is a straight positional walk.
    """
    states: list[str] = []
    for target_char, member_char in zip(target_seq, member_bases):
        if member_char == GAP:
            states.append(SiteState.UNCOVERED.value)
        elif member_char == target_char:
            states.append(SiteState.MATCH.value)
        else:
            states.append(SiteState.MISMATCH.value)
    return "".join(states)
