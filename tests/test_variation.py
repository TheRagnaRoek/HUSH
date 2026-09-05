"""Tests for the per-position comparison (msa.variation)."""

import pytest

from core.sequence import BioSequence, MolType
from msa.alignment import AlignedRow, PanelAlignment
from msa.msa_errors import AlignmentError
from msa.variation import PanelVariation, SiteState

TARGET_SEQ = "ATGCGTACGT"  # 10 nt


@pytest.fixture
def target() -> BioSequence:
    return BioSequence(TARGET_SEQ, id="TARGET", moltype=MolType.DNA)


@pytest.fixture
def alignment(target) -> PanelAlignment:
    """Four members, one per situation worth distinguishing."""
    return PanelAlignment.from_rows(
        target,
        [
            AlignedRow("same", "ATGCGTACGT"),  # identical
            AlignedRow("substituted", "ATGCGTACGA"),  # differs at position 9
            AlignedRow("deleted", "ATGC---CGT"),  # uncovered at 4, 5, 6
            AlignedRow("partial", "---CGTACGT"),  # uncovered at 0, 1, 2
        ],
    )


@pytest.fixture
def variation(alignment) -> PanelVariation:
    return PanelVariation.from_alignment(alignment)


class TestStateCharacters:
    def test_the_three_states_are_distinct(self):
        values = {state.value for state in SiteState}
        assert len(values) == 3

    def test_uncovered_uses_the_alignment_gap_character(self):
        from msa.alignment import GAP

        assert SiteState.UNCOVERED.value == GAP


class TestComparison:
    @pytest.mark.parametrize(
        "member_id,expected",
        [
            ("same", "=========="),
            ("substituted", "=========X"),
            ("deleted", "====---==="),
            ("partial", "---======="),
        ],
    )
    def test_state_rows(self, variation, member_id, expected):
        assert variation.states_for(member_id) == expected

    def test_every_row_is_target_width(self, variation):
        assert all(len(row) == variation.width for row in variation.states.values())

    def test_ids_and_order_follow_the_panel(self, variation, alignment):
        assert variation.ids == alignment.ids

    def test_target_is_carried_through(self, variation, target):
        assert variation.target is target


class TestUncoveredIsNotMismatch:
    """The distinction this module exists for: missing data is not evidence."""

    def test_deletion_reads_as_uncovered_not_mismatch(self, variation):
        assert variation.state_at("deleted", 5) is SiteState.UNCOVERED
        assert variation.state_at("deleted", 5) is not SiteState.MISMATCH

    def test_uncovered_positions_are_not_polymorphic(self, variation):
        # Positions 4-6 are uncovered for one member and matched by the rest.
        # None of them may count as a difference.
        assert set(variation.polymorphic_positions()).isdisjoint({4, 5, 6})

    def test_uncovered_positions_are_not_conserved_either(self, variation):
        assert set(variation.conserved_positions()).isdisjoint({4, 5, 6})


class TestLookups:
    def test_state_at(self, variation):
        assert variation.state_at("substituted", 9) is SiteState.MISMATCH
        assert variation.state_at("substituted", 0) is SiteState.MATCH

    def test_column_reports_every_member(self, variation):
        column = variation.column(9)
        assert column == {
            "same": SiteState.MATCH,
            "substituted": SiteState.MISMATCH,
            "deleted": SiteState.MATCH,
            "partial": SiteState.MATCH,
        }

    def test_membership_and_length(self, variation):
        assert "deleted" in variation
        assert "absent" not in variation
        assert len(variation) == 4

    def test_a_window_is_a_slice(self, variation):
        assert variation.states_for("deleted")[3:7] == "=---"

    def test_unknown_member_raises(self, variation):
        with pytest.raises(AlignmentError):
            variation.states_for("absent")

    @pytest.mark.parametrize("position", [-1, 10, 99])
    def test_position_outside_the_target_raises(self, variation, position):
        with pytest.raises(IndexError):
            variation.column(position)


class TestQueries:
    def test_polymorphic_positions(self, variation):
        # Only position 9 differs, and only for one member.
        assert list(variation.polymorphic_positions()) == [9]

    def test_conserved_positions(self, variation):
        # Every member must both cover and match: 0-2 are uncovered for
        # 'partial', 4-6 for 'deleted', and 9 differs for 'substituted'.
        assert list(variation.conserved_positions()) == [3, 7, 8]

    def test_a_fully_identical_panel_is_conserved_everywhere(self, target):
        alignment = PanelAlignment.from_rows(
            target, [AlignedRow("a", TARGET_SEQ), AlignedRow("b", TARGET_SEQ)]
        )
        variation = PanelVariation.from_alignment(alignment)
        assert list(variation.conserved_positions()) == list(range(len(TARGET_SEQ)))
        assert list(variation.polymorphic_positions()) == []


class TestValidation:
    def test_rows_must_cover_exactly_the_alignment_members(self, alignment):
        with pytest.raises(AlignmentError):
            PanelVariation(alignment, {"same": "=" * 10})

    def test_rows_must_be_target_width(self, alignment):
        with pytest.raises(AlignmentError):
            PanelVariation(alignment, {member: "==" for member in alignment.ids})

    def test_unknown_state_characters_are_rejected(self, alignment):
        with pytest.raises(AlignmentError):
            PanelVariation(alignment, {member: "?" * 10 for member in alignment.ids})

    def test_states_are_copied_on_construction(self, alignment):
        source = {member: "=" * 10 for member in alignment.ids}
        variation = PanelVariation(alignment, source)
        source["same"] = "XXXXXXXXXX"
        assert variation.states_for("same") == "=========="
