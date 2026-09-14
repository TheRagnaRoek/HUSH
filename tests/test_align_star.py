"""Tests for the star alignment backend (msa.align_star)."""

import random

import pytest

from core.sequence import BioSequence, MolType, SequencePanel
from msa.align_star import align_star, project_to_target
from msa.alignment import GAP, AlignmentBackend
from msa.msa_errors import AlignmentError

# A target with no internal repeats, so the placement of an indel introduced
# below is unambiguous and the expected alignment can be written down by hand.
TARGET_SEQ = "ATGCGTACGTTAGCCATTGCACGTAACGGATC"  # 32 nt
CUT_AT, CUT_LEN = 12, 3  # TARGET_SEQ[12:15] == "GCC"


@pytest.fixture
def target() -> BioSequence:
    return BioSequence(TARGET_SEQ, id="TARGET", moltype=MolType.DNA)


def _panel(*members: BioSequence) -> SequencePanel:
    return SequencePanel("test-panel", list(members))


def _dna(seq: str, seq_id: str) -> BioSequence:
    return BioSequence(seq, id=seq_id, moltype=MolType.DNA)


# --------------------------------------------------------------------------- #
# The projection, tested without an aligner
# --------------------------------------------------------------------------- #


class TestProjection:
    """project_to_target takes two strings, so every case can be hand-checked."""

    @pytest.mark.parametrize(
        "aligned_target,aligned_member,expected_bases,expected_insertions",
        [
            # No gaps anywhere: the member row passes through unchanged.
            ("ATGC", "ATGC", "ATGC", {}),
            ("ATGC", "ATGA", "ATGA", {}),
            # Gap in the member = a deletion. It has a target coordinate, so it
            # stays in the row as GAP.
            ("ATGCGTA", "ATG--TA", "ATG--TA", {}),
            ("ATGGCCGTC", "ATGGCC---", "ATGGCC---", {}),
            # Gap in the target = an insertion in the member. No target
            # coordinate, so it moves into `insertions`, keyed by the position
            # it precedes.
            ("ATG---C", "ATGAAAC", "ATGC", {3: "AAA"}),
            ("--AT", "AAAT", "AT", {0: "AA"}),
            ("AT--", "ATAA", "AT", {2: "AA"}),
            # Both at once.
            ("AT--GC-A", "ATAAG--A", "ATG-A", {2: "AA"}),
            ("", "", "", {}),
        ],
    )
    def test_hand_checked_cases(
        self, aligned_target, aligned_member, expected_bases, expected_insertions
    ):
        bases, insertions = project_to_target(aligned_target, aligned_member)
        assert bases == expected_bases
        assert insertions == expected_insertions

    def test_rows_of_different_length_are_rejected(self):
        with pytest.raises(AlignmentError):
            project_to_target("ATG", "AT")


class TestProjectionInvariants:
    """Properties that must hold for any pairwise alignment, not just the
    hand-written ones. Random alignments are cheap and catch off-by-one errors
    that fixed examples miss."""

    @staticmethod
    def _random_alignment(rng: random.Random, length: int = 200) -> tuple[str, str]:
        target_row, member_row = [], []
        for _ in range(length):
            roll = rng.random()
            if roll < 0.12:  # insertion in the member
                target_row.append(GAP)
                member_row.append(rng.choice("ACGT"))
            elif roll < 0.24:  # deletion in the member
                target_row.append(rng.choice("ACGT"))
                member_row.append(GAP)
            else:
                base = rng.choice("ACGT")
                target_row.append(base)
                member_row.append(base if rng.random() < 0.9 else rng.choice("ACGT"))
        return "".join(target_row), "".join(member_row)

    @pytest.fixture
    def cases(self) -> list[tuple[str, str]]:
        rng = random.Random(7)
        return [self._random_alignment(rng) for _ in range(300)]

    def test_row_length_always_matches_the_ungapped_target(self, cases):
        for aligned_target, aligned_member in cases:
            bases, _ = project_to_target(aligned_target, aligned_member)
            assert len(bases) == len(aligned_target.replace(GAP, ""))

    def test_insertions_are_anchored_inside_the_target(self, cases):
        for aligned_target, aligned_member in cases:
            _, insertions = project_to_target(aligned_target, aligned_member)
            width = len(aligned_target.replace(GAP, ""))
            assert all(0 <= position <= width for position in insertions)

    def test_no_residue_is_lost(self, cases):
        # Everything the member carried must survive, either in the row or as an
        # insertion payload.
        for aligned_target, aligned_member in cases:
            bases, insertions = project_to_target(aligned_target, aligned_member)
            recovered = bases.replace(GAP, "") + "".join(
                insertions[position] for position in sorted(insertions)
            )
            assert sorted(recovered) == sorted(aligned_member.replace(GAP, ""))


# --------------------------------------------------------------------------- #
# End to end, through Biopython
# --------------------------------------------------------------------------- #


class TestAlignStar:
    def test_identical_member_reproduces_the_target(self, target):
        panel = _panel(_dna(TARGET_SEQ, "identical"))
        alignment = align_star(target, panel)
        assert alignment["identical"].bases == TARGET_SEQ
        assert alignment["identical"].insertions == {}

    def test_substitution_keeps_the_row_aligned(self, target):
        # One base changed; no indel, so the row is the member itself.
        changed = TARGET_SEQ[:20] + ("T" if TARGET_SEQ[20] != "T" else "A") + TARGET_SEQ[21:]
        alignment = align_star(target, _panel(_dna(changed, "substitution")))
        row = alignment["substitution"]
        assert row.bases == changed
        assert row.insertions == {}
        assert row.coverage == 1.0

    def test_deletion_becomes_gaps_at_the_deleted_positions(self, target):
        member_seq = TARGET_SEQ[:CUT_AT] + TARGET_SEQ[CUT_AT + CUT_LEN :]
        alignment = align_star(target, _panel(_dna(member_seq, "deletion")))
        row = alignment["deletion"]

        expected = TARGET_SEQ[:CUT_AT] + GAP * CUT_LEN + TARGET_SEQ[CUT_AT + CUT_LEN :]
        assert row.bases == expected
        assert row.insertions == {}
        assert row.covered == len(TARGET_SEQ) - CUT_LEN

    def test_insertion_is_recorded_off_the_row(self, target):
        member_seq = TARGET_SEQ[:CUT_AT] + "TTT" + TARGET_SEQ[CUT_AT:]
        alignment = align_star(target, _panel(_dna(member_seq, "insertion")))
        row = alignment["insertion"]

        assert row.bases == TARGET_SEQ  # nothing missing from the target's frame
        assert row.insertions == {CUT_AT: "TTT"}

    def test_every_row_matches_the_target_width(self, target):
        member_seq = TARGET_SEQ[:CUT_AT] + TARGET_SEQ[CUT_AT + CUT_LEN :]
        panel = _panel(
            _dna(TARGET_SEQ, "identical"),
            _dna(member_seq, "deletion"),
            _dna(TARGET_SEQ[:CUT_AT] + "TTT" + TARGET_SEQ[CUT_AT:], "insertion"),
        )
        alignment = align_star(target, panel)
        assert all(len(row) == len(TARGET_SEQ) for row in alignment)

    def test_panel_order_and_ids_are_preserved(self, target):
        panel = _panel(
            _dna(TARGET_SEQ, "first"),
            _dna(TARGET_SEQ, "second"),
            _dna(TARGET_SEQ, "third"),
        )
        assert align_star(target, panel).ids == ["first", "second", "third"]

    def test_backend_is_recorded(self, target):
        alignment = align_star(target, _panel(_dna(TARGET_SEQ, "x")))
        assert alignment.backend is AlignmentBackend.STAR

    def test_same_input_gives_the_same_result(self, target):
        # Among equally scoring alignments the first is taken; that choice has
        # to be reproducible or two runs would disagree.
        panel = _panel(_dna(TARGET_SEQ[:CUT_AT] + TARGET_SEQ[CUT_AT + CUT_LEN :], "d"))
        first = align_star(target, panel)["d"]
        second = align_star(target, panel)["d"]
        assert first.bases == second.bases
        assert first.insertions == second.insertions


class TestMoleculeTypes:
    def test_rna_member_is_compared_against_a_dna_target(self, target):
        # Without conversion every T/U pair would mismatch and the row would be
        # unusable. The projected row comes back in the target's alphabet.
        rna = BioSequence(TARGET_SEQ.replace("T", "U"), id="rna", moltype=MolType.RNA)
        row = align_star(target, _panel(rna))["rna"]
        assert row.bases == TARGET_SEQ
        assert row.coverage == 1.0


class TestValidation:
    def test_empty_target_is_rejected(self):
        empty = BioSequence("", id="empty")
        with pytest.raises(AlignmentError):
            align_star(empty, _panel(_dna(TARGET_SEQ, "x")))

    def test_empty_panel_is_rejected(self, target):
        with pytest.raises(AlignmentError):
            align_star(target, _panel())

    def test_mixed_molecule_types_are_rejected(self, target):
        panel = _panel(
            _dna(TARGET_SEQ, "dna"),
            BioSequence("MKLVWQAGYF", id="protein", moltype=MolType.PROTEIN),
        )
        with pytest.raises(AlignmentError):
            align_star(target, panel)

    def test_protein_panel_against_a_nucleotide_target_is_rejected(self, target):
        panel = _panel(BioSequence("MKLVWQAGYF", id="protein", moltype=MolType.PROTEIN))
        with pytest.raises(AlignmentError):
            align_star(target, panel)

    def test_duplicate_member_ids_are_rejected(self, target):
        panel = _panel(_dna(TARGET_SEQ, "same"), _dna(TARGET_SEQ, "same"))
        with pytest.raises(AlignmentError):
            align_star(target, panel)
