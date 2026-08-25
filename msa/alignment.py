"""The alignment result the MSA / allele stage speaks in (pipeline stage 3).

Every panel member is projected onto **ungapped target coordinates**: each member
becomes a string exactly as long as the target, where index *i* holds the residue
that member carries at target position *i*. Columns the alignment inserted into
the target are removed during projection, so no downstream unit ever converts
between alignment columns and target positions.

Three states are kept distinct at each target position:

* a residue - matching the target or not (this module does not judge which; that
  is msa.variation's job),
* GAP - the member has no residue there. Not the same as "different": a window
  where a targeted allele is uncovered is unusable, not good or bad,
* an insertion - extra residues with no target coordinate at all, recorded in
  AlignedRow.insertions rather than dropped.

Nothing is aligned here. The backends live in their own modules (msa.align_star,
and whatever follows) and are named by AlignmentBackend, so one can be swapped
for another without touching this module or the layers above it.

Scope and deferrals:

* No alignment algorithm, no scoring, no verdict.
* Insertions are recorded but not interpreted. An insertion would make an siRNA
  bulge on that member, which is a discriminating feature; weighting it needs
  the mismatch thermodynamics core.thermo does not yet carry.
* The projection is lossy by construction: residues a member carries outside the
  target's span are kept only as insertion payloads, without alignment detail.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Iterator

from core.sequence import BioSequence
from msa.msa_errors import AlignmentError

# The character marking a target position the member does not cover. Matches the
# gap alphabet core.sequence already accepts, so a projected row can be wrapped
# in a BioSequence without further translation.
GAP = "-"


# --------------------------------------------------------------------------- #
# Backend identity
# --------------------------------------------------------------------------- #


class AlignmentBackend(str, Enum):
    """Which method produced an alignment, selectable and recordable by name.

    A hashable stand-in for the backend itself, so a caller can carry a choice
    of method through a frozen criteria object and into a report without
    embedding a function. Each name corresponds to one module under msa/.

    STAR: pairwise alignment of every member against the target, the target
        acting as the common reference (msa.align_star). Adequate while panel
        members are near-identical; it degrades as they diverge.
    """

    STAR = "star"


# --------------------------------------------------------------------------- #
# One panel member, projected onto the target
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AlignedRow:
    """One panel member expressed in target coordinates.

    Parameters:
        id:    The member's sequence id, as it appears in the panel.
        bases: One character per target position, in target order. A residue
               where the member covers that position, GAP where it does not.
               Always the same length as the target.
        insertions: Extra residues the member carries that the target does not.
               Keyed by the target position they sit *before*: key ``i`` means
               "between target positions i-1 and i". ``0`` is before the target
               starts, ``len(target)`` after it ends.
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
        """Number of target positions this member has a residue for."""
        return len(self.bases) - self.bases.count(GAP)

    @property
    def coverage(self) -> float:
        """Fraction of the target this member covers (0.0 for an empty target)."""
        return self.covered / len(self.bases) if self.bases else 0.0

    def __repr__(self) -> str:
        extra = f", insertions={len(self.insertions)}" if self.insertions else ""
        return f"AlignedRow(id={self.id!r}, coverage={self.coverage:.0%}{extra})"


# --------------------------------------------------------------------------- #
# The whole panel, projected onto the target
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class PanelAlignment:
    """Every panel member projected onto one ungapped target sequence.

    The hand-off from an alignment backend to msa.variation. The invariant every
    backend must uphold - enforced here rather than trusted - is that each row is
    exactly as long as the target, so index *i* means target position *i* in
    every row without exception.

    Rows are keyed by sequence id because that is how every consumer addresses
    them; panel order is preserved by the mapping's insertion order.

    Parameters:
        target:  The target the panel was aligned against, ungapped and in its
                 original coordinates.
        rows:    Sequence id -> AlignedRow, in panel order.
        backend: Which method produced this alignment, or None for an alignment
                 assembled by hand (tests, imported files).
    """

    target: BioSequence
    rows: dict[str, AlignedRow]
    backend: AlignmentBackend | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", dict(self.rows))
        self._validate()

    # -- construction ------------------------------------------------------ #

    @classmethod
    def from_rows(
        cls,
        target: BioSequence,
        rows: Iterable[AlignedRow],
        backend: AlignmentBackend | None = None,
    ) -> "PanelAlignment":
        """Build an alignment from rows in panel order.

        The convenient entry point for backends, which naturally produce rows one
        after another rather than a mapping. Duplicate ids are rejected here
        rather than silently overwriting each other.
        """
        mapping: dict[str, AlignedRow] = {}
        for row in rows:
            if row.id in mapping:
                raise AlignmentError(
                    f"Duplicate row id {row.id!r} in alignment; each panel member "
                    "must appear exactly once."
                )
            mapping[row.id] = row
        return cls(target=target, rows=mapping, backend=backend)

    # -- validation -------------------------------------------------------- #

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

    # -- accessors --------------------------------------------------------- #

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self) -> Iterator[AlignedRow]:
        """Iterate over the rows, in panel order."""
        return iter(self.rows.values())

    def __contains__(self, row_id: str) -> bool:
        return row_id in self.rows

    def __getitem__(self, row_id: str) -> AlignedRow:
        """The row for this sequence id, raising a legible error if absent."""
        try:
            return self.rows[row_id]
        except KeyError:
            raise AlignmentError(
                f"No aligned row for {row_id!r}. Present: {sorted(self.rows)}"
            ) from None

    def get(self, row_id: str) -> AlignedRow | None:
        """The row for this sequence id, or None when it is absent."""
        return self.rows.get(row_id)

    @property
    def width(self) -> int:
        """Number of target positions - the length of every row."""
        return len(self.target)

    @property
    def ids(self) -> list[str]:
        """Row ids, in panel order."""
        return list(self.rows)

    def column(self, position: int) -> dict[str, str]:
        """The residue every member carries at one target position.

        GAP for members that do not cover it. The target's own base is not
        included - read it from ``alignment.target.seq[position]``.
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
