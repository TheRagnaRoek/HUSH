"""Allele targeting requests - what the user wants a candidate to hit and spare.

The MSA/allele stage answers one question: which sites on the target mRNA let an
siRNA discriminate between the members of a panel the way the user asked for.
This module holds the *request* side of that question - the domain model the rest
of the stage is parameterised by. It performs no alignment and computes no
verdicts; those live in the alignment and evaluation units.

Three user intents are commonly named, but they are not three different things -
they are one structure with different values:

* allele-specific - hit exactly one panel member, spare the rest.
* pan-allele      - hit every panel member, spare none.
* subset          - hit a chosen group, spare the rest.

All three reduce to two disjoint sets of allele ids: the ones a candidate must
match ('hit') and the ones it must be able to discriminate against ('spare').
Modelling them as one type rather than three keeps a single code path through the
stage, and admits combinations the three named cases do not cover (e.g. sparing
only part of the panel and leaving the remainder unconstrained). The named cases
are provided as constructors, not as separate types.

TargetingMode is *derived* from the two sets rather than stored, so a request's
reported mode can never disagree with its contents.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Iterable

from core.sequence import SequencePanel
from msa.msa_errors import PanelError, TargetingError

# --------------------------------------------------------------------------- #
# Mode
# --------------------------------------------------------------------------- #


class TargetingMode(str, Enum):
    """The named intent a request corresponds to, derived from its sets.

    ALLELE_SPECIFIC: exactly one allele is targeted.
    PAN_ALLELE:      nothing is spared - the candidate should hit everything asked
                     for and discriminate against nothing.
    SUBSET:          several alleles are targeted while others are spared.
    """

    ALLELE_SPECIFIC = "allele_specific"
    PAN_ALLELE = "pan_allele"
    SUBSET = "subset"


# --------------------------------------------------------------------------- #
# Request
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AlleleTargeting:
    """Which panel members a candidate site must hit, and which it must spare.

    Immutable, mirroring the rest of the HUSH domain model. Both fields are
    normalised to frozensets on construction, so any iterable of ids may be
    passed. The two sets must be disjoint: an allele cannot be both a target and
    a bystander to protect.

    'spare' is stored explicitly rather than being taken as "whatever else is in
    the panel". A panel is data that can change; a request that silently meant
    something different after a sequence was added to the panel would be a trap.
    The constructors below fill it in from a panel at the moment the request is
    made, which is when the user's intent is actually known.

    Parameters:
        hit:   Ids of the panel members the candidate must match. Never empty.
        spare: Ids the candidate must be able to discriminate against. May be
               empty (pan-allele targeting).
    """

    hit: frozenset[str]
    spare: frozenset[str] = frozenset()

    def __post_init__(self) -> None:
        # Accept any iterable of ids; store frozensets. The dataclass is frozen,
        # so normalisation goes through object.__setattr__ - the same pattern
        # core.sequence.BioSequence uses to clean its sequence string.
        object.__setattr__(self, "hit", frozenset(self.hit))
        object.__setattr__(self, "spare", frozenset(self.spare))
        self._validate()

    # -- validation -------------------------------------------------------- #

    def _validate(self) -> None:
        if not self.hit:
            raise TargetingError(
                "AlleleTargeting requires at least one allele in 'hit'; a request "
                "that targets nothing has no meaning."
            )
        overlap = self.hit & self.spare
        if overlap:
            raise TargetingError(
                f"Allele(s) {sorted(overlap)} appear in both 'hit' and 'spare'. "
                "A candidate cannot be required to match and to discriminate "
                "against the same sequence."
            )

    # -- derived properties ------------------------------------------------ #

    @property
    def mode(self) -> TargetingMode:
        """The named intent this request corresponds to.

        Derived, never stored: a request is defined by its sets, and the mode is
        only a label for reporting.
        """
        if not self.spare:
            return TargetingMode.PAN_ALLELE
        if len(self.hit) == 1:
            return TargetingMode.ALLELE_SPECIFIC
        return TargetingMode.SUBSET

    @property
    def constrained(self) -> frozenset[str]:
        """Every allele the request says something about (hit or spare)."""
        return self.hit | self.spare

    @property
    def discriminating(self) -> bool:
        """True when the request asks a candidate to tell alleles apart."""
        return bool(self.spare)

    def unconstrained(self, panel: SequencePanel) -> frozenset[str]:
        """Panel members this request says nothing about.

        Non-empty only for requests built by hand; the constructors below assign
        every panel member to one of the two sets.
        """
        return frozenset(panel.ids) - self.constrained

    # -- named constructors ------------------------------------------------ #

    @classmethod
    def allele_specific(cls, allele: str, panel: SequencePanel) -> "AlleleTargeting":
        """Hit exactly 'allele'; spare every other member of 'panel'."""
        ids = _panel_ids(panel)
        _require_present({allele}, ids)
        return cls(hit={allele}, spare=ids - {allele})

    @classmethod
    def pan_allele(cls, panel: SequencePanel) -> "AlleleTargeting":
        """Hit every member of 'panel'; spare nothing.

        Selects conserved sites: a candidate satisfying this request matches all
        panel members and is not asked to discriminate against anything.
        """
        return cls(hit=_panel_ids(panel), spare=frozenset())

    @classmethod
    def subset(
        cls, alleles: Iterable[str], panel: SequencePanel
    ) -> "AlleleTargeting":
        """Hit the given 'alleles'; spare the remaining members of 'panel'."""
        wanted = frozenset(alleles)
        ids = _panel_ids(panel)
        _require_present(wanted, ids)
        return cls(hit=wanted, spare=ids - wanted)

    def __repr__(self) -> str:
        return (
            f"AlleleTargeting(mode={self.mode.value}, "
            f"hit={sorted(self.hit)}, spare={sorted(self.spare)})"
        )


# --------------------------------------------------------------------------- #
# Panel helpers
# --------------------------------------------------------------------------- #


# How many missing ids to name before truncating the message. A panel pulled
# from IPD-IMGT/HLA can hold thousands of alleles, so 'available' is never listed
# in full.
_MAX_REPORTED_IDS = 10


def _panel_ids(panel: SequencePanel) -> frozenset[str]:
    """The panel's sequence ids, rejecting blank or duplicated ones.

    Ids are the only handle a request has on a panel member, so ambiguity here
    would silently make a request mean something other than what was asked. The
    checks are linear in panel size: an IPD-IMGT/HLA locus can hold thousands of
    alleles, and this runs on every request.
    """
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
    """Raise if any requested id is absent from the panel."""
    missing = wanted - available
    if missing:
        raise TargetingError(
            f"Allele(s) {_truncate(sorted(missing))} are not in the panel "
            f"({len(available)} member(s))."
        )


def _truncate(ids: list[str]) -> str:
    """Render an id list for an error message, capped at _MAX_REPORTED_IDS."""
    if len(ids) <= _MAX_REPORTED_IDS:
        return str(ids)
    shown = ids[:_MAX_REPORTED_IDS]
    return f"{shown} (+{len(ids) - _MAX_REPORTED_IDS} more)"
