"""siRNA to shRNA conversion (stage 4).

Grafts both strands of a candidate, overhangs included, onto a stem-loop and
emits the bare hairpin insert as DNA oligos. Each construct is folded and
checked (core.hairpin). See docs/shrna-conversion.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from core.hairpin import HairpinFold, most_stable_hairpin
from core.selection import longest_homopolymer_run
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate
from shrna.loops import DEFAULT_LOOP, LOOPS, LoopScaffold

# Pol III stops at a run of four or more U (Bogenhagen & Brown 1981).
_POL3_TERMINATOR_RUN = "UUUU"

# T4 is the minimal terminator; full efficiency from T6 (Gao et al. 2018).
MIN_TERMINATOR_LENGTH = 4
DEFAULT_TERMINATOR_LENGTH = 6


def _as_encoded_overhang(overhang: str) -> str:
    """Overhang as encoded in the hairpin: dT becomes ribo-U."""
    return overhang.upper().replace("T", "U")


class Orientation(str, Enum):
    """Strand order of the hairpin; only SA (sense-loop-antisense) is built."""

    SA = "sa"


class StemOverhang(str, Enum):
    """Fate of the sense overhang between stem and loop.

    KEEP encodes it, so the hairpin carries the full 21-nt sense strand. TRIM
    drops it, so the loop joins the core directly. The guide overhang is always
    kept as the 3' end.
    """

    KEEP = "keep"
    TRIM = "trim"


def get_stem_overhang(spec: str) -> StemOverhang:
    """The StemOverhang for a name; ValueError if unknown."""
    try:
        return StemOverhang(spec.strip().lower())
    except ValueError:
        choices = ", ".join(m.value for m in StemOverhang)
        raise ValueError(f"unknown stem overhang {spec!r}; choose from {choices}.")


@dataclass(frozen=True)
class ShrnaOptions:
    """Loop, orientation, stem overhang, terminator and premature-terminator policy."""

    loop: LoopScaffold = LOOPS[DEFAULT_LOOP]
    orientation: Orientation = Orientation.SA
    stem_overhang: StemOverhang = StemOverhang.KEEP
    include_terminator: bool = True
    terminator_length: int = DEFAULT_TERMINATOR_LENGTH
    drop_premature_terminators: bool = True

    def __post_init__(self) -> None:
        if not isinstance(self.loop, LoopScaffold):
            raise TypeError(
                f"loop must be a LoopScaffold, got {type(self.loop).__name__}"
            )
        if not isinstance(self.stem_overhang, StemOverhang):
            raise TypeError(
                "stem_overhang must be a StemOverhang, got "
                f"{type(self.stem_overhang).__name__}"
            )
        if self.orientation is not Orientation.SA:
            raise ValueError(
                f"only the SA orientation is implemented, got {self.orientation}."
            )
        if self.include_terminator and self.terminator_length < MIN_TERMINATOR_LENGTH:
            raise ValueError(
                f"terminator_length must be >= {MIN_TERMINATOR_LENGTH} when the "
                f"terminator is included, got {self.terminator_length}."
            )


@dataclass(frozen=True)
class ShrnaConstruct:
    """One construct with its QC read-outs (described in docs/shrna-conversion.md).

    start and end are the 0-based core coordinates on the target. hairpin_rna
    includes the terminator as U. insert_top and insert_bottom are the DNA
    oligos, bottom being the reverse complement of top. fold is the most stable
    fold of the body, or None.
    """

    region: str
    start: int
    end: int
    loop_name: str
    orientation: Orientation
    hairpin_rna: BioSequence
    insert_top: BioSequence
    insert_bottom: BioSequence
    fold: HairpinFold | None
    folds_as_designed: bool
    starts_with_purine: bool
    premature_terminator_free: bool
    max_homopolymer_run: int
    sense_overhang: str = ""
    guide_overhang: str = ""

    def __repr__(self) -> str:
        ok = "folds-as-designed" if self.folds_as_designed else "MISFOLD-RISK"
        return (
            f"ShrnaConstruct(region={self.region}, loop={self.loop_name}, "
            f"{ok}, 5'={'A/G' if self.starts_with_purine else 'C/U'})"
        )


def convert(
    sense: BioSequence,
    antisense: BioSequence,
    *,
    region: str,
    start: int,
    end: int,
    sense_overhang: str = "",
    guide_overhang: str = "",
    options: ShrnaOptions | None = None,
) -> ShrnaConstruct:
    """Build a construct from the core strands (RNA, 5'->3') and their overhangs.

    The overhangs are separate because the sense overhang lands between stem
    and loop, and the guide overhang in the 3' tail, so the stem is the core.
    With stem_overhang=TRIM the sense overhang is left out of the hairpin.
    """
    opts = options or ShrnaOptions()
    sense_rna = sense.as_rna()
    antisense_rna = antisense.as_rna()
    trim = opts.stem_overhang is StemOverhang.TRIM
    sense_oh = "" if trim else _as_encoded_overhang(sense_overhang)
    guide_oh = _as_encoded_overhang(guide_overhang)

    # The sense overhang is unpaired, so it counts as loop; the guide overhang
    # is a single-stranded 3' end and joins the terminator.
    intended_loop = sense_oh + opts.loop.seq
    body = sense_rna.seq + intended_loop + antisense_rna.seq
    terminator = "U" * opts.terminator_length if opts.include_terminator else ""
    hairpin_str = body + guide_oh + terminator

    hairpin_rna = BioSequence(
        hairpin_str, id=f"shRNA:{region}", moltype=MolType.RNA
    )
    insert_top = hairpin_rna.reverse_transcribe()
    insert_bottom = insert_top.reverse_complement()

    # Fold the body only; the 3' tail is single-stranded.
    fold = most_stable_hairpin(body)
    folds_as_designed = (
        fold is not None
        and fold.stem_length == len(sense_rna)
        and fold.loop == intended_loop
    )

    run, _ = longest_homopolymer_run(body)

    return ShrnaConstruct(
        region=region,
        start=start,
        end=end,
        loop_name=opts.loop.name,
        orientation=opts.orientation,
        hairpin_rna=hairpin_rna,
        insert_top=insert_top,
        insert_bottom=insert_bottom,
        fold=fold,
        folds_as_designed=folds_as_designed,
        starts_with_purine=hairpin_str[:1] in ("A", "G"),
        premature_terminator_free=_POL3_TERMINATOR_RUN not in body,
        max_homopolymer_run=run,
        sense_overhang=sense_oh,
        guide_overhang=guide_oh,
    )


def convert_candidate(
    candidate: SiRNACandidate, options: ShrnaOptions | None = None
) -> ShrnaConstruct:
    """Convert a SiRNACandidate."""
    return convert(
        candidate.sense,
        candidate.antisense,
        region=candidate.region,
        start=candidate.position,
        end=candidate.end,
        sense_overhang=candidate.sense_overhang,
        guide_overhang=candidate.guide_overhang,
        options=options,
    )
