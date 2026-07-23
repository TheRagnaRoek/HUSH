"""siRNA design core (pipeline stage 2): candidate model and enumeration.

This is the substrate the rest of the design stage consumes. Nothing is scored
here - this module only turns a target mRNA into the full set of possible siRNA
duplex candidates, one per sliding window. The classical rule pre-screen
(Ui-Tei / Reynolds / Amarzguioui), thermodynamic asymmetry, and ranking are
layered on top in later units.

Conventions used throughout:

* The duplex 'core' is the fully complementary region (default 19 nt). The 2-nt
  3' overhangs are a design choice handled later, not baked into candidates.
* sense / passenger strand  - matches the mRNA target site, 5'->3'.
* antisense / guide strand  - reverse-complement of the target site, 5'->3'
  (the strand loaded into RISC to recognise the mRNA).

Both strands are stored as validated RNA BioSequence objects so downstream code
can call the existing nucleotide operations (gc_content, seed extraction, ...)
without re-parsing.
"""

from __future__ import annotations

from dataclasses import dataclass
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
        sense:     Passenger strand (RNA), equal to the mRNA target site.
        antisense: Guide strand (RNA), the reverse-complement of the site.
    """

    position: int
    length: int
    sense: BioSequence
    antisense: BioSequence

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
    target: BioSequence, length: int = DEFAULT_LENGTH
) -> Iterator[SiRNACandidate]:
    """Yield every siRNA candidate site by sliding a window along the target.

    The target is validated as a nucleotide sequence and converted to RNA. A
    window of 'length' nt is then slid base by base from the 5' end; each window
    position gives one candidate whose sense strand is the mRNA site and whose
    antisense strand is that site's reverse-complement.

    Input is validated eagerly (before iteration begins) so an invalid target or
    length raises immediately rather than partway through consumption.

    Args:
        target: The target mRNA/CDS as a BioSequence (DNA is transcribed to RNA).
        length: Complementary core length in nt (default 19, range 19-24).

    Yields:
        SiRNACandidate objects, 5'->3', in order of increasing position.
    """
    rna = _prepare_target(target, length)
    return _candidate_windows(rna, length)


def _candidate_windows(rna: BioSequence, length: int) -> Iterator[SiRNACandidate]:
    """Generator over the prepared RNA target; assumes inputs are validated."""
    last_start = len(rna) - length
    for position in range(last_start + 1):
        sense = rna.subseq(position, position + length)
        antisense = sense.reverse_complement()
        yield SiRNACandidate(
            position=position,
            length=length,
            sense=sense,
            antisense=antisense,
        )
