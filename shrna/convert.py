"""siRNA to shRNA conversion (pipeline stage 4): stem-loop scaffold.

Carries a validated siRNA candidate into a cloneable stem-loop shRNA. The
19-bp duplex core becomes the stem unchanged, so the guide the design stage
validated (seed-Tm, asymmetry, off-target) is exactly the guide the hairpin
expresses - no re-shift, no re-validation. The construct is emitted in the
sense-loop-antisense (SA) orientation as the bare hairpin insert (top and
bottom DNA oligos), with a Pol III poly-T terminator but without any
vector-specific restriction overhangs: the wet-lab adds sticky-end tails for
whichever vector they clone into.

Unlike the naive templating of earlier tools, each construct is scored: the
hairpin is folded (core.hairpin) to confirm it forms the intended stem-loop
rather than a competing structure, and a few legible QC flags (Pol III 5'-G
initiation, loop terminator-safety, homopolymer runs) are attached.

Rationale and references in docs/shrna_design_note.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from core.hairpin import HairpinFold, most_stable_hairpin
from core.selection import longest_homopolymer_run
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate
from shrna.loops import DEFAULT_LOOP, LOOPS, LoopScaffold

# A run of >= 4 U/T in the transcribed strand is read by Pol III as a
# termination signal. The hairpin body (stem + loop, before the deliberate 3'
# terminator) must not contain one - including at the loop/stem junctions, which
# the upstream homopolymer gate cannot see.
_POL3_TERMINATOR_RUN = "UUUU"

# Minimum poly-T terminator length: Pol III needs a run of at least ~4-5 T.
MIN_TERMINATOR_LENGTH = 4
DEFAULT_TERMINATOR_LENGTH = 6


class Orientation(str, Enum):
    """Strand order of the hairpin.

    SA: sense - loop - antisense (the pSUPER convention). The only orientation
    built for now; modelled as an enum so antisense-loop-sense can be added
    later without changing the construct's shape.
    """

    SA = "sa"


@dataclass(frozen=True)
class ShrnaOptions:
    """How to build the hairpin from a duplex core.

    Attributes:
        loop: the loop scaffold to graft (default: Jensen UGUGCUU).
        orientation: strand order (SA only for now).
        include_terminator: append the Pol III poly-T terminator.
        terminator_length: length of that poly-T run (>= 4 when included).
    """

    loop: LoopScaffold = LOOPS[DEFAULT_LOOP]
    orientation: Orientation = Orientation.SA
    include_terminator: bool = True
    terminator_length: int = DEFAULT_TERMINATOR_LENGTH

    def __post_init__(self) -> None:
        if not isinstance(self.loop, LoopScaffold):
            raise TypeError(
                f"loop must be a LoopScaffold, got {type(self.loop).__name__}"
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
    """One stem-loop shRNA construct, with its fold-based QC read-outs.

    Attributes:
        region: 1-based inclusive target coordinate label of the source core.
        start / end: 0-based start and exclusive end of the core on the target.
        loop_name: the grafted loop's registry name.
        orientation: strand order used.
        hairpin_rna: the transcribed hairpin, 5'->3' (incl. terminator as U's).
        insert_top / insert_bottom: the orderable DNA oligos, 5'->3'; annealed,
            they form the double-stranded insert (bottom = RC of top).
        fold: the most stable hairpin the stem+loop folds into (loop excluded
            terminator), or None if it folds into nothing.
        folds_as_designed: True when that most stable fold IS the intended
            stem (full core length) closing the intended loop - i.e. the
            construct does not prefer a competing structure.
        starts_with_g: the transcript's 5' base is G (the U6 Pol III +1
            preference; the graft-junction initiation requirement).
        premature_terminator_free: the hairpin body (stem + loop, before the 3'
            terminator) carries no >= 4 U run - including at the loop/stem
            junctions - that would terminate Pol III early.
        max_homopolymer_run: longest single-nucleotide run in the stem+loop
            (the terminator's own poly-T run is excluded by construction).
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
    starts_with_g: bool
    premature_terminator_free: bool
    max_homopolymer_run: int

    def __repr__(self) -> str:
        ok = "folds-as-designed" if self.folds_as_designed else "MISFOLD-RISK"
        return (
            f"ShrnaConstruct(region={self.region}, loop={self.loop_name}, "
            f"{ok}, 5'={'G' if self.starts_with_g else 'non-G'})"
        )


def convert(
    sense: BioSequence,
    antisense: BioSequence,
    *,
    region: str,
    start: int,
    end: int,
    options: ShrnaOptions | None = None,
) -> ShrnaConstruct:
    """Build a stem-loop shRNA construct from one duplex core.

    'sense' and 'antisense' are the paired-core strands (RNA, 5'->3'), exactly
    as siRNA design defines them: sense is the target site and antisense its
    reverse-complement. Overhangs are not part of the stem and must not be
    passed in. The core is used unchanged as the stem, so the guide is preserved.
    """
    opts = options or ShrnaOptions()
    sense_rna = sense.as_rna()
    antisense_rna = antisense.as_rna()

    loop_seq = opts.loop.seq
    stem_loop = sense_rna.seq + loop_seq + antisense_rna.seq
    terminator = "U" * opts.terminator_length if opts.include_terminator else ""
    hairpin_str = stem_loop + terminator

    hairpin_rna = BioSequence(
        hairpin_str, id=f"shRNA:{region}", moltype=MolType.RNA
    )
    insert_top = hairpin_rna.reverse_transcribe()
    insert_bottom = insert_top.reverse_complement()

    # Fold the stem+loop only: the poly-T terminator is a single-stranded tail,
    # not part of the hairpin, so it must not enter the fold search.
    fold = most_stable_hairpin(stem_loop)
    folds_as_designed = (
        fold is not None
        and fold.stem_length == len(sense_rna)
        and fold.loop == loop_seq
    )

    run, _ = longest_homopolymer_run(stem_loop)

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
        starts_with_g=hairpin_str[:1] == "G",
        premature_terminator_free=_POL3_TERMINATOR_RUN not in stem_loop,
        max_homopolymer_run=run,
    )


def convert_candidate(
    candidate: SiRNACandidate, options: ShrnaOptions | None = None
) -> ShrnaConstruct:
    """Convenience wrapper: convert a SiRNACandidate straight to a construct."""
    return convert(
        candidate.sense,
        candidate.antisense,
        region=candidate.region,
        start=candidate.position,
        end=candidate.end,
        options=options,
    )
