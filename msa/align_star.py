"""Star alignment: each panel member aligned pairwise to the target.

The target is the fixed centre, so every row is projected onto the target
coordinates the siRNA candidates use. Method, scores and limits:
docs/design-rules.md (Allele targeting, Alignment).
"""

from __future__ import annotations

from Bio import Align

from core.sequence import BioSequence, MolType, SequencePanel
from msa.alignment import GAP, AlignedRow, AlignmentBackend, PanelAlignment
from msa.msa_errors import AlignmentError

DEFAULT_MATCH = 2.0
DEFAULT_MISMATCH = -1.0
DEFAULT_OPEN_GAP = -10.0
DEFAULT_EXTEND_GAP = -0.5

# Free terminal gaps: a partial member hangs off the ends and shows as GAP.
DEFAULT_END_GAP = 0.0


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
    """Align every member to the target; rows follow panel order.

    Members are converted to the target's molecule type first, so T and U do
    not mismatch. Raises AlignmentError for an empty or mixed-type panel or a
    member the aligner rejects.
    """
    _validate_inputs(target, panel)
    aligner = _make_aligner(match, mismatch, open_gap, extend_gap, end_gap)

    rows = [
        _align_one(aligner, target, member)
        for member in panel
    ]
    return PanelAlignment.from_rows(target, rows, backend=AlignmentBackend.STAR)


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


def _make_aligner(
    match: float,
    mismatch: float,
    open_gap: float,
    extend_gap: float,
    end_gap: float,
) -> Align.PairwiseAligner:
    """Global aligner: the target must be covered end to end."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    # open/extend_gap_score also set the terminal gap scores, so the end-gap
    # scores must be assigned after them.
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    # Biopython renamed these and warns on the old spelling; the old one is
    # what the biopython>=1.83 floor in pyproject.toml guarantees.
    aligner.target_end_gap_score = end_gap
    aligner.query_end_gap_score = end_gap
    return aligner


def _as_target_moltype(member: BioSequence, target: BioSequence) -> BioSequence:
    """The member in the target's molecule type."""
    if not target.is_nucleotide or member.moltype is target.moltype:
        return member
    return member.as_rna() if target.moltype is MolType.RNA else member.reverse_transcribe()


def _align_one(
    aligner: Align.PairwiseAligner, target: BioSequence, member: BioSequence
) -> AlignedRow:
    """Align one member and project it onto target coordinates."""
    prepared = _as_target_moltype(member, target)
    try:
        alignments = aligner.align(target.seq, prepared.seq)
        # Take the first of the equal-scoring alignments; their number can be
        # huge for repetitive sequences.
        best = alignments[0]
    except Exception as exc:  # pragma: no cover - Biopython raises several types
        raise AlignmentError(
            f"Could not align {member.label!r} against {target.label!r}: {exc}"
        ) from exc

    bases, insertions = project_to_target(best[0], best[1])
    return AlignedRow(id=member.id, bases=bases, insertions=insertions)


def project_to_target(
    aligned_target: str, aligned_member: str
) -> tuple[str, dict[int, str]]:
    """Collapse a pairwise alignment onto ungapped target coordinates.

    Returns (bases, insertions). A column with a target residue yields the
    member's residue, or GAP. Member residues opposite a target gap are held
    back and attached as an insertion to the next target position. Takes plain
    strings, so it does not depend on Biopython. Raises AlignmentError when the
    rows differ in length.
    """
    if len(aligned_target) != len(aligned_member):
        raise AlignmentError(
            f"Alignment rows differ in length: {len(aligned_target)} vs "
            f"{len(aligned_member)}."
        )

    bases: list[str] = []
    insertions: dict[int, str] = {}
    pending: list[str] = []
    position = 0  # next ungapped target position to fill

    for target_char, member_char in zip(aligned_target, aligned_member):
        if target_char == GAP:
            # No target coordinate: hold the residue for the next one.
            if member_char != GAP:
                pending.append(member_char)
            continue
        if pending:
            insertions[position] = "".join(pending)
            pending = []
        bases.append(member_char)
        position += 1

    if pending:
        # Trailing residues past the end of the target.
        insertions[position] = "".join(pending)

    return "".join(bases), insertions
