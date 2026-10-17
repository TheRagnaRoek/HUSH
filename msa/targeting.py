"""Allele targeting request: which panel members a candidate must hit and spare.

All intents (pan-allele, allele-specific, subset, target-specific) are one
structure of two disjoint id sets, built by named constructors. The mode is
derived from the sets, never stored. See docs/design-rules.md (Allele targeting).
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Iterable

from core.sequence import SequencePanel
from msa.msa_errors import PanelError, TargetingError


class TargetingMode(str, Enum):
    """Intent label derived from a request's sets, used for reporting."""

    ALLELE_SPECIFIC = "allele_specific"
    PAN_ALLELE = "pan_allele"
    SUBSET = "subset"
    TARGET_SPECIFIC = "target_specific"


@dataclass(frozen=True)
class AlleleTargeting:
    """Panel member ids a candidate must match (hit) and discriminate against (spare).

    The sets must be disjoint and not both empty. hit may be empty only for
    target-specific targeting. spare is stored explicitly rather than taken as
    the rest of the panel, so a request keeps its meaning when the panel
    changes; the constructors fill it from the panel at request time.
    """

    hit: frozenset[str]
    spare: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        # Frozen dataclass: normalise any iterable to frozensets via object.__setattr__.
        object.__setattr__(self, "hit", frozenset(self.hit))
        object.__setattr__(self, "spare", frozenset(self.spare))
        self._validate()

    def _validate(self) -> None:
        if not self.hit and not self.spare:
            raise TargetingError(
                "AlleleTargeting requires at least one allele in 'hit' or 'spare'; "
                "a request that constrains nothing has no meaning. To target the "
                "loaded target alone, spare the whole panel (target-specific)."
            )
        overlap = self.hit & self.spare
        if overlap:
            raise TargetingError(
                f"Allele(s) {sorted(overlap)} appear in both 'hit' and 'spare'. "
                "A candidate cannot be required to match and to discriminate "
                "against the same sequence."
            )

    @property
    def mode(self) -> TargetingMode:
        """Intent label: no hit is target-specific, no spare is pan-allele."""
        if not self.hit:
            return TargetingMode.TARGET_SPECIFIC
        if not self.spare:
            return TargetingMode.PAN_ALLELE
        if len(self.hit) == 1:
            return TargetingMode.ALLELE_SPECIFIC
        return TargetingMode.SUBSET

    @property
    def constrained(self) -> frozenset[str]:
        """Every allele in hit or spare."""
        return self.hit | self.spare

    @property
    def discriminating(self) -> bool:
        """True when any allele is spared."""
        return bool(self.spare)

    def unconstrained(self, panel: SequencePanel) -> frozenset[str]:
        """Panel members in neither set; empty for the named constructors."""
        return frozenset(panel.ids) - self.constrained

    @classmethod
    def allele_specific(cls, allele: str, panel: SequencePanel) -> "AlleleTargeting":
        """Hit one allele; spare every other member."""
        ids = _panel_ids(panel)
        _require_present({allele}, ids)
        return cls(hit={allele}, spare=ids - {allele})

    @classmethod
    def pan_allele(cls, panel: SequencePanel) -> "AlleleTargeting":
        """Hit every member; spare none."""
        return cls(hit=_panel_ids(panel), spare=frozenset())

    @classmethod
    def subset(
        cls, alleles: Iterable[str], panel: SequencePanel
    ) -> "AlleleTargeting":
        """Hit the given alleles; spare the remaining members."""
        wanted = frozenset(alleles)
        ids = _panel_ids(panel)
        _require_present(wanted, ids)
        return cls(hit=wanted, spare=ids - wanted)

    @classmethod
    def target_specific(cls, panel: SequencePanel) -> "AlleleTargeting":
        """Hit no member; spare the whole panel (the target is loaded apart from it)."""
        return cls(hit=frozenset(), spare=_panel_ids(panel))

    def __repr__(self) -> str:
        return (
            f"AlleleTargeting(mode={self.mode.value}, "
            f"hit={sorted(self.hit)}, spare={sorted(self.spare)})"
        )


# Panels can hold thousands of alleles, so error messages list at most this many ids.
_MAX_REPORTED_IDS = 10


def _panel_ids(panel: SequencePanel) -> frozenset[str]:
    """The panel's ids; blank or duplicate ids raise PanelError."""
    ids = panel.ids
    blank = sum(1 for seq_id in ids if not seq_id)
    if blank:
        raise PanelError(
            f"Panel {panel.name!r} has {blank} sequence(s) without an id; "
            "allele targeting addresses panel members by id."
        )
    duplicates = sorted(
        seq_id for seq_id, count in Counter(ids).items() if count > 1
    )
    if duplicates:
        raise PanelError(
            f"Panel {panel.name!r} contains duplicate sequence id(s) "
            f"{_truncate(duplicates)}; targeting by id would be ambiguous."
        )
    return frozenset(ids)


def _require_present(wanted: frozenset[str], available: frozenset[str]) -> None:
    """Raise TargetingError for ids that are not in the panel."""
    missing = wanted - available
    if missing:
        raise TargetingError(
            f"Allele(s) {_truncate(sorted(missing))} are not in the panel "
            f"({len(available)} member(s))."
        )


def _truncate(ids: list[str]) -> str:
    """Render ids for an error message, capped at _MAX_REPORTED_IDS."""
    if len(ids) <= _MAX_REPORTED_IDS:
        return str(ids)
    shown = ids[:_MAX_REPORTED_IDS]
    return f"{shown} (+{len(ids) - _MAX_REPORTED_IDS} more)"
