"""siRNA candidate model and enumeration (stage 2).

One candidate per window along the target. Nothing is scored here. The duplex
geometry and overhang policies are described in docs/design-rules.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator

from core.errors import MolTypeError
from core.sequence import BioSequence

# Core length excludes the overhang; the range is enforced, not clamped.
DEFAULT_LENGTH = 19
CANONICAL_MIN_LENGTH = 19
CANONICAL_MAX_LENGTH = 24


def _validate_length(length: int) -> None:
    if not isinstance(length, int):
        raise TypeError(f"siRNA length must be an int, got {type(length).__name__}")
    if not CANONICAL_MIN_LENGTH <= length <= CANONICAL_MAX_LENGTH:
        raise ValueError(
            f"siRNA core length {length} is outside the canonical range "
            f"{CANONICAL_MIN_LENGTH}-{CANONICAL_MAX_LENGTH} nt."
        )


OVERHANG_LENGTH = 2

# Lowercase marks the DNA residues of a dTdT overhang inside an RNA strand.
_DTDT_OVERHANG = "tt"
_UU_OVERHANG = "UU"


class OverhangPolicy(str, Enum):
    """3' overhang source: target flanks (NATURAL), dTdT or UU."""

    NATURAL = "natural"
    DTDT = "dtdt"
    UU = "uu"


def _overhangs(
    rna: BioSequence, position: int, length: int, policy: OverhangPolicy
) -> tuple[str, str]:
    """(passenger_overhang, guide_overhang), 5'->3', for one site."""
    if policy is OverhangPolicy.DTDT:
        return _DTDT_OVERHANG, _DTDT_OVERHANG
    if policy is OverhangPolicy.UU:
        return _UU_OVERHANG, _UU_OVERHANG

    core_end = position + length
    tail = rna.seq[core_end : core_end + OVERHANG_LENGTH]
    passenger = tail if len(tail) == OVERHANG_LENGTH else _DTDT_OVERHANG
    if position >= OVERHANG_LENGTH:
        guide = rna.subseq(position - OVERHANG_LENGTH, position).reverse_complement().seq
    else:
        guide = _DTDT_OVERHANG
    return passenger, guide


@dataclass(frozen=True)
class SiRNACandidate:
    """One siRNA duplex against a single site of the target mRNA.

    position is the 0-based start on the mRNA. sense is the passenger core
    (the mRNA site), antisense the guide core (its reverse complement).
    The overhangs lie outside the paired core.
    """

    position: int
    length: int
    sense: BioSequence
    antisense: BioSequence
    sense_overhang: str = ""
    guide_overhang: str = ""

    @property
    def end(self) -> int:
        """0-based exclusive end of the site."""
        return self.position + self.length

    @property
    def region(self) -> str:
        """1-based inclusive label, e.g. "497-515"."""
        return f"{self.position + 1}-{self.end}"

    @property
    def guide(self) -> BioSequence:
        """Alias for antisense."""
        return self.antisense

    @property
    def passenger(self) -> BioSequence:
        """Alias for sense."""
        return self.sense

    @property
    def guide_full(self) -> str:
        """Guide core plus its overhang, 5'->3'."""
        return self.antisense.seq + self.guide_overhang

    @property
    def passenger_full(self) -> str:
        """Passenger core plus its overhang, 5'->3'."""
        return self.sense.seq + self.sense_overhang

    @property
    def gc_content(self) -> float:
        """GC fraction of the core."""
        return self.sense.gc_content

    def __repr__(self) -> str:
        return (
            f"SiRNACandidate(region={self.region}, len={self.length}, "
            f"sense={self.sense.seq!r}, antisense={self.antisense.seq!r})"
        )


def _prepare_target(target: BioSequence, length: int) -> BioSequence:
    """Validate the inputs and return the target as RNA."""
    _validate_length(length)
    if not target.is_nucleotide:
        raise MolTypeError(
            f"Target {target.label!r} is {target.moltype.value}; siRNA design "
            "requires a nucleotide (mRNA/CDS) target."
        )
    rna = target.as_rna()
    if len(rna) < length:
        raise ValueError(
            f"Target {rna.label!r} is {len(rna)} nt, shorter than the requested "
            f"siRNA length {length} nt; no candidates can be enumerated."
        )
    return rna


def enumerate_candidates(
    target: BioSequence,
    length: int = DEFAULT_LENGTH,
    overhang: OverhangPolicy = OverhangPolicy.NATURAL,
) -> Iterator[SiRNACandidate]:
    """Yield a candidate for every window of 'length' nt, in position order.

    Inputs are validated before the first candidate is yielded.
    """
    rna = _prepare_target(target, length)
    return _candidate_windows(rna, length, overhang)


def _candidate_windows(
    rna: BioSequence, length: int, overhang: OverhangPolicy
) -> Iterator[SiRNACandidate]:
    last_start = len(rna) - length
    for position in range(last_start + 1):
        sense = rna.subseq(position, position + length)
        antisense = sense.reverse_complement()
        sense_oh, guide_oh = _overhangs(rna, position, length, overhang)
        yield SiRNACandidate(
            position=position,
            length=length,
            sense=sense,
            antisense=antisense,
            sense_overhang=sense_oh,
            guide_overhang=guide_oh,
        )
