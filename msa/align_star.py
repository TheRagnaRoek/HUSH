"""Star alignment: every panel member aligned pairwise against the target.

The target is the centre of the star. Each panel member is aligned to it on its
own, and the resulting pairwise alignment is projected onto ungapped target
coordinates. Because every member is aligned to the same reference, the rows
share one coordinate system without a multiple alignment ever being built - and
that coordinate system is the one siRNA candidates already live in.

This is Gusfield's center-star construction, which carries a factor-2 bound on
the sum-of-pairs cost when the underlying distance obeys the triangle
inequality. The classical method picks the centre to minimise the summed
pairwise cost; here the centre is fixed by the pipeline (candidates are
enumerated from the target), so the bound's premise does not apply and no
optimality is claimed. What the construction buys us is a consistent frame, not
a quality guarantee.

Scope and deferrals:

* Adequate while panel members are near-identical - HLA alleles of one locus run
  ~98-99% identity. Accuracy degrades as members diverge, because a single
  centre cannot represent the relationships within the panel the way a guide
  tree does. Divergent families want a progressive backend instead.
* Gaps are placed without regard to reading frame. For coding sequences a
  codon-aware aligner would not open frame-shifting gaps; the default gap costs
  below make indels expensive enough that this rarely bites on closely related
  members, but it is an approximation, not a constraint.
* Where several alignments score equally the first is taken. Among equals the
  choice is arbitrary, so a position reported as differing is not necessarily
  the only reading of the same indel.
* The scoring parameters are project defaults, not published constants. They are
  exposed so a caller can state them explicitly rather than inherit them.

Rationale and the validation record belong in docs/design-rules.md once this
unit is reviewed.
"""

from __future__ import annotations

from Bio import Align

from core.sequence import BioSequence, MolType, SequencePanel
from msa.alignment import GAP, AlignedRow, AlignmentBackend, PanelAlignment
from msa.msa_errors import AlignmentError

# --------------------------------------------------------------------------- #
# Scoring defaults
# --------------------------------------------------------------------------- #

# Project defaults, not published constants. Chosen so that on closely related
# coding sequences a substitution is always cheaper than opening a gap: an indel
# has to buy back at least four mismatches before it pays for itself. That keeps
# the aligner from scattering gaps to chase a few mismatches, which is the
# failure mode that would silently move which positions look polymorphic.
DEFAULT_MATCH = 2.0
DEFAULT_MISMATCH = -1.0
DEFAULT_OPEN_GAP = -10.0
DEFAULT_EXTEND_GAP = -0.5

# Terminal gaps are free by default: a member that is a partial CDS should hang
# off the ends rather than be stretched across the target to avoid a penalty.
# The uncovered positions surface as GAP in the projected row, which is exactly
# the "not covered" state downstream units need to see.
DEFAULT_END_GAP = 0.0


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #


def align_star(
    target: BioSequence,
    panel: SequencePanel,
    *,
    match: float = DEFAULT_MATCH,
    mismatch: float = DEFAULT_MISMATCH,
    open_gap: float = DEFAULT_OPEN_GAP,
    extend_gap: float = DEFAULT_EXTEND_GAP,
    end_gap: float = DEFAULT_END_GAP,
) -> PanelAlignment:
    """Align every panel member against the target and project onto its frame.

    Members are converted to the target's molecule type first, so a DNA target
    and an RNA member are compared base for base rather than mismatching on
    every T/U.

    Args:
        target: The target sequence, ungapped. Its coordinates are the frame
            every row is expressed in.
        panel: The related sequences to align against it.
        match: Score for a matching pair of residues.
        mismatch: Score for a differing pair.
        open_gap: Cost of opening a gap. Strongly negative by default so an
            indel is not used to save a handful of mismatches.
        extend_gap: Cost of each further gap position.
        end_gap: Cost of a gap at either end of either sequence. Zero by
            default, so a partial member is not forced to stretch.

    Returns:
        PanelAlignment with one row per panel member, in panel order, each
        projected onto ungapped target coordinates.

    Raises:
        AlignmentError: When the target or panel cannot be aligned - an empty
            panel, a panel of mixed molecule types, or a member the aligner
            rejects.
    """
    _validate_inputs(target, panel)
    aligner = _make_aligner(match, mismatch, open_gap, extend_gap, end_gap)

    rows = [
        _align_one(aligner, target, member)
        for member in panel
    ]
    return PanelAlignment.from_rows(target, rows, backend=AlignmentBackend.STAR)


# --------------------------------------------------------------------------- #
# Input checks
# --------------------------------------------------------------------------- #


def _validate_inputs(target: BioSequence, panel: SequencePanel) -> None:
    """Reject inputs that cannot produce a usable alignment."""
    if not len(target):
        raise AlignmentError("Cannot align against an empty target sequence.")
    if not len(panel):
        raise AlignmentError(
            f"Panel {panel.name!r} is empty; there is nothing to align."
        )
    if not panel.is_uniform:
        kinds = sorted(m.value for m in panel.moltypes)
        raise AlignmentError(
            f"Panel {panel.name!r} mixes molecule types {kinds}; align a panel "
            "of one kind."
        )
    if target.is_nucleotide != panel.sequences[0].is_nucleotide:
        raise AlignmentError(
            f"Target is {target.moltype.value} but the panel is "
            f"{panel.sequences[0].moltype.value}; these cannot be aligned."
        )


# --------------------------------------------------------------------------- #
# Aligner setup
# --------------------------------------------------------------------------- #


def _make_aligner(
    match: float,
    mismatch: float,
    open_gap: float,
    extend_gap: float,
    end_gap: float,
) -> Align.PairwiseAligner:
    """A global aligner configured with the given costs.

    Global mode is required: the target must be covered end to end, since every
    target position has to receive a residue or a gap.
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    # Order matters: open_gap_score/extend_gap_score set every gap score,
    # terminal ones included. The end-gap scores must be assigned afterwards to
    # take effect.
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    # Biopython renamed these to end_insertion_score / end_deletion_score and
    # warns on the old spellings. The old ones are kept because they are the
    # form guaranteed by the biopython>=1.83 floor in pyproject.toml; switching
    # means raising that floor, which is a project-wide call.
    aligner.target_end_gap_score = end_gap
    aligner.query_end_gap_score = end_gap
    return aligner


def _as_target_moltype(member: BioSequence, target: BioSequence) -> BioSequence:
    """Return 'member' in the target's molecule type, so T/U cannot mismatch."""
    if not target.is_nucleotide or member.moltype is target.moltype:
        return member
    return member.as_rna() if target.moltype is MolType.RNA else member.reverse_transcribe()


# --------------------------------------------------------------------------- #
# One member
# --------------------------------------------------------------------------- #


def _align_one(
    aligner: Align.PairwiseAligner, target: BioSequence, member: BioSequence
) -> AlignedRow:
    """Align one member against the target and project it onto target space."""
    prepared = _as_target_moltype(member, target)
    try:
        alignments = aligner.align(target.seq, prepared.seq)
        # Among equally scoring alignments the first is taken; see the module
        # docstring. Materialising all of them is not attempted - the count can
        # be astronomically large for repetitive sequences.
        best = alignments[0]
    except Exception as exc:  # pragma: no cover - Biopython raises several types
        raise AlignmentError(
            f"Could not align {member.label!r} against {target.label!r}: {exc}"
        ) from exc

    bases, insertions = project_to_target(best[0], best[1])
    return AlignedRow(id=member.id, bases=bases, insertions=insertions)


# --------------------------------------------------------------------------- #
# Projection - the part that has to be right
# --------------------------------------------------------------------------- #


def project_to_target(
    aligned_target: str, aligned_member: str
) -> tuple[str, dict[int, str]]:
    """Collapse a pairwise alignment onto ungapped target coordinates.

    Walks the alignment column by column. A column where the target carries a
    residue contributes one character to the projected row - the member's
    residue, or GAP where the member has none. A column where the *target* is
    gapped has no target coordinate at all; the member's residue there is set
    aside and attached to the next target position as an insertion.

    Deliberately free of Biopython types: it takes two equal-length strings, so
    it can be tested against hand-written cases and reused by any other backend
    that can produce a pairwise alignment.

    Args:
        aligned_target: The target row of a pairwise alignment, gaps included.
        aligned_member: The member row, same length.

    Returns:
        (bases, insertions) - the projected row and the insertions keyed by the
        target position they precede.

    Raises:
        AlignmentError: When the two rows differ in length.
    """
    if len(aligned_target) != len(aligned_member):
        raise AlignmentError(
            f"Alignment rows differ in length: {len(aligned_target)} vs "
            f"{len(aligned_member)}."
        )

    bases: list[str] = []
    insertions: dict[int, str] = {}
    pending: list[str] = []
    position = 0  # next ungapped target position to be filled

    for target_char, member_char in zip(aligned_target, aligned_member):
        if target_char == GAP:
            # No target coordinate here - hold the residue for the next one.
            if member_char != GAP:
                pending.append(member_char)
            continue
        if pending:
            insertions[position] = "".join(pending)
            pending = []
        bases.append(member_char)
        position += 1

    if pending:
        # Trailing residues, past the end of the target.
        insertions[position] = "".join(pending)

    return "".join(bases), insertions
