"""siRNA design core (pipeline stage 2): candidate model and enumeration.

This is the substrate the rest of the design stage consumes. Nothing is scored
here - this module only turns a target mRNA into the full set of possible siRNA
duplex candidates, one per sliding window. The classical rule pre-screen
(Ui-Tei / Reynolds / Amarzguioui), thermodynamic asymmetry, and ranking are
layered on top in later units.

Conventions used throughout:

* The duplex 'core' is the fully complementary region (default 19 nt), the
  substrate every downstream filter (rules, GC, seed, ...) operates on.
* Each strand also carries a 2-nt 3' overhang, tracked separately from the core
  so the filters never see it. Default overhangs are mRNA-derived (the real
  offset Tuschl duplex); dTdT / UU are selectable (see OverhangPolicy).
* sense / passenger strand  - matches the mRNA target site, 5'->3'.
* antisense / guide strand  - reverse-complement of the target site, 5'->3'
  (the strand loaded into RISC to recognise the mRNA).

Both strands are stored as validated RNA BioSequence objects so downstream code
can call the existing nucleotide operations (gc_content, seed extraction, ...)
without re-parsing.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterator

from core.errors import MolTypeError
from core.sequence import BioSequence

# --------------------------------------------------------------------------- #
# Length policy
# --------------------------------------------------------------------------- #

# The complementary core length, excluding any 3' overhang. siRNA duplexes
# canonically run 19-24 nt; 19 is the siDirect-style default. This is a
# user-tunable parameter - the range is enforced (not silently accepted) so a
# typo or an out-of-range value surfaces immediately, and can be relaxed later.
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


# --------------------------------------------------------------------------- #
# Overhang policy
# --------------------------------------------------------------------------- #

# Canonical synthetic siRNAs carry symmetric 2-nt 3' overhangs (Elbashir et al.
# 2001, EMBO J 20:6877). They sit outside the paired core, so no filter sees
# them - they only define the strand a user would order.
OVERHANG_LENGTH = 2

# dTdT is 2'-deoxythymidine; lowercase marks the DNA residues in an RNA strand.
_DTDT_OVERHANG = "tt"
_UU_OVERHANG = "UU"


class OverhangPolicy(str, Enum):
    """Which 3' overhang to append to each strand's paired core.

    NATURAL: the target-derived flanks. The passenger overhang is the 2 nt just
        3' of the core on the mRNA; the guide overhang is the reverse-complement
        of the 2 nt just 5' of the core (the true offset duplex). A candidate
        whose flank runs past a transcript end falls back to dTdT there.
    DTDT: two 2'-deoxythymidines on both strands - the common synthetic default.
    UU: two ribo-uridines on both strands.
    """

    NATURAL = "natural"
    DTDT = "dtdt"
    UU = "uu"


def _overhangs(
    rna: BioSequence, position: int, length: int, policy: OverhangPolicy
) -> tuple[str, str]:
    """Return (passenger_overhang, guide_overhang), 5'->3', for one site."""
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


# --------------------------------------------------------------------------- #
# Candidate model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SiRNACandidate:
    """One siRNA duplex candidate against a single target mRNA site.

    Immutable, mirroring BioSequence. Coordinates are relative to the target
    mRNA the candidate was enumerated from.

    Attributes:
        position:  0-based start of the target site on the mRNA.
        length:    Length of the complementary core (nt).
        sense:     Passenger strand core (RNA), equal to the mRNA target site.
        antisense: Guide strand core (RNA), the reverse-complement of the site.
        sense_overhang / guide_overhang: the 2-nt 3' overhangs, 5'->3', outside
            the paired core (empty only if a candidate is built without them).
    """

    position: int
    length: int
    sense: BioSequence
    antisense: BioSequence
    sense_overhang: str = ""
    guide_overhang: str = ""

    @property
    def end(self) -> int:
        """0-based, exclusive end of the target site on the mRNA."""
        return self.position + self.length

    @property
    def region(self) -> str:
        """1-based inclusive coordinate label, matching subseq's convention."""
        return f"{self.position + 1}-{self.end}"

    @property
    def guide(self) -> BioSequence:
        """Alias for the antisense strand (the RISC-loaded guide)."""
        return self.antisense

    @property
    def passenger(self) -> BioSequence:
        """Alias for the sense strand."""
        return self.sense

    @property
    def guide_full(self) -> str:
        """Guide strand as ordered, 5'->3': paired core plus its 3' overhang."""
        return self.antisense.seq + self.guide_overhang

    @property
    def passenger_full(self) -> str:
        """Passenger strand as ordered, 5'->3': paired core plus its 3' overhang."""
        return self.sense.seq + self.sense_overhang

    @property
    def gc_content(self) -> float:
        """GC fraction of the duplex core (identical for both strands)."""
        return self.sense.gc_content

    def __repr__(self) -> str:
        return (
            f"SiRNACandidate(region={self.region}, len={self.length}, "
            f"sense={self.sense.seq!r}, antisense={self.antisense.seq!r})"
        )


# --------------------------------------------------------------------------- #
# Enumeration
# --------------------------------------------------------------------------- #


def _prepare_target(target: BioSequence, length: int) -> BioSequence:
    """Validate inputs and return the target as RNA, ready for windowing."""
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
    """Yield every siRNA candidate site by sliding a window along the target.

    The target is validated as a nucleotide sequence and converted to RNA. A
    window of 'length' nt is then slid base by base from the 5' end; each window
    position gives one candidate whose sense strand is the mRNA site and whose
    antisense strand is that site's reverse-complement. Each candidate also gets
    2-nt 3' overhangs per 'overhang' (see OverhangPolicy), tracked apart from the
    core so no downstream filter is affected by them.

    Input is validated eagerly (before iteration begins) so an invalid target or
    length raises immediately rather than partway through consumption.

    Args:
        target: The target mRNA/CDS as a BioSequence (DNA is transcribed to RNA).
        length: Complementary core length in nt (default 19, range 19-24).
        overhang: Which 3' overhang to append (default mRNA-derived NATURAL).

    Yields:
        SiRNACandidate objects, 5'->3', in order of increasing position.
    """
    rna = _prepare_target(target, length)
    return _candidate_windows(rna, length, overhang)


def _candidate_windows(
    rna: BioSequence, length: int, overhang: OverhangPolicy
) -> Iterator[SiRNACandidate]:
    """Generator over the prepared RNA target; assumes inputs are validated."""
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
