"""siRNA to shRNA conversion (pipeline stage 4): stem-loop scaffold.

Carries a validated siRNA candidate into a cloneable stem-loop shRNA. Both
strands are encoded in full - paired core plus its 3' overhang - so the guide
the design stage validated (seed-Tm, asymmetry, off-target) is exactly the
guide the hairpin expresses, overhang included. The construct is emitted in
the sense-loop-antisense (SA) orientation as the bare hairpin insert (top and
bottom DNA oligos), with a Pol III poly-T terminator but without any
vector-specific restriction tails: the wet-lab adds the sticky ends for
whichever vector they clone into.

Unlike the naive templating of earlier tools, each construct is scored: the
hairpin is folded (core.hairpin) to confirm it forms the intended stem-loop
rather than a competing structure, and a few legible QC flags (Pol III 5'-A/G
initiation, loop terminator-safety, homopolymer runs) are attached.

Rationale and references in docs/shrna-conversion.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from core.hairpin import HairpinFold, most_stable_hairpin
from core.selection import longest_homopolymer_run
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate
from shrna.loops import DEFAULT_LOOP, LOOPS, LoopScaffold

# A run of >= 4 U/T terminates Pol III (Bogenhagen & Brown 1981; Gao et al.
# 2018). The hairpin body must not contain one, including at the loop/stem
# junctions, which the upstream homopolymer gate cannot see.
_POL3_TERMINATOR_RUN = "UUUU"

# Minimum poly-T terminator length: T4 is the minimal Pol III terminator, with
# full efficiency at >= 6 T (Gao et al. 2018).
MIN_TERMINATOR_LENGTH = 4
DEFAULT_TERMINATOR_LENGTH = 6


def _as_encoded_overhang(overhang: str) -> str:
    """Render a candidate's 3' overhang as the hairpin encodes it.

    An expressed hairpin cannot carry the 2'-deoxythymidines of an ordered
    dTdT duplex, so those positions are encoded as ribo-U. The DNA insert
    oligo carries T either way.
    """
    return overhang.upper().replace("T", "U")


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
        drop_premature_terminators: withhold constructs whose body carries
            a premature Pol III terminator. On by default: grafting can
            build a U run across a junction that no upstream gate sees,
            and such a construct truncates before the hairpin completes.
    """

    loop: LoopScaffold = LOOPS[DEFAULT_LOOP]
    orientation: Orientation = Orientation.SA
    include_terminator: bool = True
    terminator_length: int = DEFAULT_TERMINATOR_LENGTH
    drop_premature_terminators: bool = True

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
        fold: the most stable hairpin the body (both arms plus loop, 3' tail
            excluded) folds into, or None if it folds into nothing.
        folds_as_designed: True when that most stable fold IS the intended
            stem (the paired core) closing the intended loop - i.e. the
            construct does not prefer a competing structure. The sense
            overhang is unpaired where it sits, so it counts as loop, not
            stem; a fold that pairs it up instead shifts the Dicer register
            and is flagged.
        starts_with_purine: the transcript's 5' base is A or G, the bases a
            human U6/7SK promoter initiates on efficiently at +1 (Gao et al.
            2017, Transcription 8:275). A C or U start is flagged.
        premature_terminator_free: the body carries no >= 4 U run - including
            at the stem/loop junctions - that would terminate Pol III early.
            The 3' tail is excluded: a U run where the guide overhang meets
            the terminator is the intended terminator, not a premature one.
        max_homopolymer_run: longest single-nucleotide run in the body (the
            3' tail's own poly-U run is excluded by construction).
        fold_structure: 'fold' in dot-bracket notation over the body, or
            empty when it folds into nothing. Paste-able into any RNA tool.
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
    fold_structure: str = ""

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
    """Build a stem-loop shRNA construct from one duplex core.

    'sense' and 'antisense' are the paired-core strands (RNA, 5'->3'), exactly
    as siRNA design defines them: sense is the target site and antisense its
    reverse-complement. The 3' overhangs come in separately so the converter
    knows where each core ends: both strands are encoded in full, but the
    sense overhang lands between stem and loop and the guide overhang in the
    3' tail, so the stem is the core alone.
    """
    opts = options or ShrnaOptions()
    sense_rna = sense.as_rna()
    antisense_rna = antisense.as_rna()
    sense_oh = _as_encoded_overhang(sense_overhang)
    guide_oh = _as_encoded_overhang(guide_overhang)

    # The sense overhang is unpaired where it sits, so it reads as loop; the
    # guide overhang is a single-stranded 3' end, so it joins the terminator.
    intended_loop = sense_oh + opts.loop.seq
    body = sense_rna.seq + intended_loop + antisense_rna.seq
    terminator = "U" * opts.terminator_length if opts.include_terminator else ""
    hairpin_str = body + guide_oh + terminator

    hairpin_rna = BioSequence(
        hairpin_str, id=f"shRNA:{region}", moltype=MolType.RNA
    )
    insert_top = hairpin_rna.reverse_transcribe()
    insert_bottom = insert_top.reverse_complement()

    # Fold the body only: the 3' tail (guide overhang plus terminator) is
    # single-stranded, so it must not enter the fold search.
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
        fold_structure="" if fold is None else fold.dot_bracket(len(body)),
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
        sense_overhang=candidate.sense_overhang,
        guide_overhang=candidate.guide_overhang,
        options=options,
    )
