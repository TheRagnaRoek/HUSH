"""Alignment result types for the allele stage.

Each panel member is projected onto ungapped target coordinates: a string as
long as the target, where index i is the member's residue at target position i.
A position holds a residue, GAP (member has none there; not the same as a
difference) or, between positions, an insertion kept in AlignedRow.insertions.
Insertions are recorded, not scored (docs/design-rules.md, Not implemented).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Iterator

from core.sequence import BioSequence
from msa.msa_errors import AlignmentError

# Marks an uncovered position; accepted by core.sequence as a gap character.
GAP = "-"


class AlignmentBackend(str, Enum):
    """Method that produced an alignment (hashable, so it can sit in frozen objects).

    STAR: every member aligned pairwise to the target (msa.align_star).
    """

    STAR = "star"


@dataclass(frozen=True)
class AlignedRow:
    """One panel member in target coordinates.

    bases has one character per target position (GAP where uncovered).
    insertions maps a target position i to the residues inserted between
    positions i-1 and i; key 0 is before the target, len(target) after it.
    """

    id: str
    bases: str
    insertions: dict[int, str] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.bases)

    def covers(self, position: int) -> bool:
        """True when the member has a residue at this target position."""
        return self.bases[position] != GAP

    @property
    def covered(self) -> int:
        """Number of covered target positions."""
        return len(self.bases) - self.bases.count(GAP)

    @property
    def coverage(self) -> float:
        """Fraction of the target covered; 0.0 for an empty target."""
        return self.covered / len(self.bases) if self.bases else 0.0

    def __repr__(self) -> str:
        extra = f", insertions={len(self.insertions)}" if self.insertions else ""
        return f"AlignedRow(id={self.id!r}, coverage={self.coverage:.0%}{extra})"


@dataclass(frozen=True)
class PanelAlignment:
    """Every panel member projected onto one target, keyed by id in panel order.

    Every row must be exactly as long as the target; this is validated.
    backend is None for an alignment assembled by hand.
    """

    target: BioSequence
    rows: dict[str, AlignedRow]
    backend: AlignmentBackend | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", dict(self.rows))
        self._validate()

    @classmethod
    def from_rows(
        cls,
        target: BioSequence,
        rows: Iterable[AlignedRow],
        backend: AlignmentBackend | None = None,
    ) -> "PanelAlignment":
        """Build from rows in panel order; duplicate ids raise AlignmentError."""
        mapping: dict[str, AlignedRow] = {}
        for row in rows:
            if row.id in mapping:
                raise AlignmentError(
                    f"Duplicate row id {row.id!r} in alignment; each panel member "
                    "must appear exactly once."
                )
            mapping[row.id] = row
        return cls(target=target, rows=mapping, backend=backend)

    def _validate(self) -> None:
        width = len(self.target)
        for row_id, row in self.rows.items():
            if row.id != row_id:
                raise AlignmentError(
                    f"Row filed under {row_id!r} carries id {row.id!r}."
                )
            if len(row) != width:
                raise AlignmentError(
                    f"Row {row_id!r} has {len(row)} positions but the target has "
                    f"{width}. Every row must be projected onto ungapped target "
                    "coordinates."
                )
            bad = sorted(p for p in row.insertions if not 0 <= p <= width)
            if bad:
                raise AlignmentError(
                    f"Row {row_id!r} anchors insertion(s) at target position(s) "
                    f"{bad}, outside 0..{width}."
                )

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[AlignedRow]:
        return iter(self.rows.values())

    def __contains__(self, row_id: str) -> bool:
        return row_id in self.rows

    def __getitem__(self, row_id: str) -> AlignedRow:
        try:
            return self.rows[row_id]
        except KeyError:
            raise AlignmentError(
                f"No aligned row for {row_id!r}. Present: {sorted(self.rows)}"
            ) from None

    def get(self, row_id: str) -> AlignedRow | None:
        """The row for this id, or None."""
        return self.rows.get(row_id)

    @property
    def width(self) -> int:
        """Number of target positions."""
        return len(self.target)

    @property
    def ids(self) -> list[str]:
        """Row ids in panel order."""
        return list(self.rows)

    def column(self, position: int) -> dict[str, str]:
        """Each member's residue at one target position (GAP if uncovered).

        The target's own base is not included.
        """
        if not 0 <= position < self.width:
            raise IndexError(
                f"Target position {position} outside 0..{self.width - 1}"
            )
        return {row_id: row.bases[position] for row_id, row in self.rows.items()}

    def __repr__(self) -> str:
        via = f", backend={self.backend.value!r}" if self.backend else ""
        return (
            f"PanelAlignment(target={self.target.id!r}, width={self.width}, "
            f"rows={len(self.rows)}{via})"
        )
