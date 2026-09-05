"""Tests for the window evaluation (msa.windows)."""

import pytest

from core.sequence import BioSequence, MolType, SequencePanel
from msa.alignment import AlignedRow, PanelAlignment
from msa.msa_errors import MsaError
from msa.targeting import AlleleTargeting
from msa.variation import PanelVariation
from msa.windows import (
    WindowCriteria,
    WindowEvaluation,
    evaluate_all,
    evaluate_window,
    guide_position,
    regions,
    select,
)

#            0    5    10   15   20
TARGET_SEQ = "ATGCGTACGTTAGCCATTGCA"  # 21 nt
SUB_A2, SUB_A3 = 8, 10  # where the spared alleles differ
DEL_FROM, DEL_LEN = 15, 3  # where A4 is uncovered
LENGTH = 6


def _sub(position: int, base: str) -> str:
    replacement = base if TARGET_SEQ[position] != base else "C"
    return TARGET_SEQ[:position] + replacement + TARGET_SEQ[position + 1 :]


@pytest.fixture
def target() -> BioSequence:
    return BioSequence(TARGET_SEQ, id="TARGET", moltype=MolType.DNA)


@pytest.fixture
def panel(target) -> SequencePanel:
    """Ids only - targeting addresses members by id, not by sequence."""
    return SequencePanel(
        "panel",
        [
            BioSequence(TARGET_SEQ, id=name, moltype=MolType.DNA)
            for name in ("A1", "A2", "A3")
        ],
    )


@pytest.fixture
def variation(target) -> PanelVariation:
    """A1 identical, A2 differs at 8, A3 at 10, A4 uncovered at 15-17."""
    alignment = PanelAlignment.from_rows(
        target,
        [
            AlignedRow("A1", TARGET_SEQ),
            AlignedRow("A2", _sub(SUB_A2, "A")),
            AlignedRow("A3", _sub(SUB_A3, "G")),
            AlignedRow(
                "A4",
                TARGET_SEQ[:DEL_FROM] + "-" * DEL_LEN + TARGET_SEQ[DEL_FROM + DEL_LEN :],
            ),
        ],
    )
    return PanelVariation.from_alignment(alignment)


def _criteria(hit, spare=(), length=LENGTH, min_spare=1) -> WindowCriteria:
    return WindowCriteria(
        AlleleTargeting(hit=hit, spare=spare),
        length=length,
        min_spare_mismatches=min_spare,
    )


# --------------------------------------------------------------------------- #


class TestCriteria:
    def test_default_length_matches_the_sirna_core(self):
        from core.sirna import DEFAULT_LENGTH

        assert WindowCriteria(AlleleTargeting(hit={"A1"})).length == DEFAULT_LENGTH

    @pytest.mark.parametrize("length", [0, -1])
    def test_non_positive_length_is_rejected(self, length):
        with pytest.raises(MsaError):
            WindowCriteria(AlleleTargeting(hit={"A1"}), length=length)

    def test_requiring_zero_differences_is_rejected(self):
        # Zero would credit every window with discriminating against everything.
        with pytest.raises(MsaError):
            WindowCriteria(AlleleTargeting(hit={"A1"}), min_spare_mismatches=0)


class TestGuidePosition:
    """The one place the coordinate flip can be got backwards."""

    def test_the_two_ends_swap(self):
        assert guide_position(0, 19) == 18
        assert guide_position(18, 19) == 0

    def test_conversion_is_its_own_inverse(self):
        for position in range(19):
            assert guide_position(guide_position(position, 19), 19) == position

    def test_the_guide_seed_sits_at_the_three_prime_end_of_the_window(self):
        # Seed = guide positions 1..7 (0-based), i.e. the published 2-8.
        seed_in_window = {guide_position(g, 19) for g in range(1, 8)}
        assert seed_in_window == set(range(11, 18))
        # Which is the far end of the window, not the near one.
        assert min(seed_in_window) > 19 // 2

    @pytest.mark.parametrize("position", [-1, 19, 99])
    def test_positions_outside_the_window_are_rejected(self, position):
        with pytest.raises(IndexError):
            guide_position(position, 19)


class TestHitGate:
    def test_an_identical_allele_passes_everywhere(self, variation):
        criteria = _criteria(hit={"A1"})
        assert all(e.hit_pass for e in evaluate_all(variation, criteria))

    def test_a_differing_allele_fails_where_it_differs(self, variation):
        criteria = _criteria(hit={"A2"})
        failing = [e.position for e in evaluate_all(variation, criteria) if not e.hit_pass]
        # Every window covering position 8 must fail.
        assert set(failing) == set(range(SUB_A2 - LENGTH + 1, SUB_A2 + 1))

    def test_max_hit_mismatches_counts_the_worst_allele(self, variation):
        criteria = _criteria(hit={"A1", "A2"})
        window = evaluate_window(variation, criteria, SUB_A2)
        assert window.max_hit_mismatches == 1

    def test_an_uncovered_target_allele_fails(self, variation):
        criteria = _criteria(hit={"A4"})
        window = evaluate_window(variation, criteria, DEL_FROM)
        assert not window.hit_pass
        assert window.uncovered["A4"]
        # Uncovered, not differing - the distinction has to survive to here.
        assert window.mismatches["A4"] == ()
        assert any("not covered" in reason for reason in window.reasons)


class TestSpareGate:
    def test_a_window_needs_a_difference_from_every_spared_allele(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2", "A3"})
        chosen = [e.region for e in select(variation, criteria)]
        # Only windows covering both 8 and 10 can discriminate against both.
        assert chosen == ["6-11", "7-12", "8-13", "9-14"]

    def test_min_spare_mismatches_reports_the_weakest_allele(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2", "A3"})
        window = evaluate_window(variation, criteria, SUB_A2)
        assert window.min_spare_mismatches == 1

    def test_raising_the_threshold_rejects_single_differences(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2", "A3"}, min_spare=2)
        assert select(variation, criteria) == []

    def test_an_uncovered_spared_allele_cannot_be_discriminated_against(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A4"})
        window = evaluate_window(variation, criteria, DEL_FROM)
        assert window.spare_pass is False
        assert any("not covered" in reason for reason in window.reasons)

    def test_pan_allele_has_no_spare_gate(self, variation):
        criteria = _criteria(hit={"A1", "A2", "A3"})
        for window in evaluate_all(variation, criteria):
            assert window.spare_pass is None
            assert window.min_spare_mismatches is None

    def test_pan_allele_selects_windows_free_of_any_difference(self, variation):
        criteria = _criteria(hit={"A1", "A2", "A3"})
        chosen = {e.position for e in select(variation, criteria)}
        covering_a_difference = set(range(SUB_A2 - LENGTH + 1, SUB_A3 + 1))
        assert chosen.isdisjoint(covering_a_difference)


class TestBreakdown:
    def test_positions_are_window_relative(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2"})
        window = evaluate_window(variation, criteria, SUB_A2 - 2)
        assert window.mismatches["A2"] == (2,)

    def test_only_constrained_alleles_are_examined(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2"})
        window = evaluate_window(variation, criteria, 0)
        assert set(window.mismatches) == {"A1", "A2"}

    def test_region_label_is_one_based_and_inclusive(self, variation):
        window = evaluate_window(variation, _criteria(hit={"A1"}), 0)
        assert window.region == f"1-{LENGTH}"
        assert window.end == LENGTH

    def test_a_satisfying_window_carries_no_reasons(self, variation):
        window = evaluate_window(variation, _criteria(hit={"A1"}), 0)
        assert window.satisfies
        assert window.reasons == ()


class TestSliding:
    def test_every_window_of_the_target_is_evaluated(self, variation):
        criteria = _criteria(hit={"A1"})
        assert len(list(evaluate_all(variation, criteria))) == len(TARGET_SEQ) - LENGTH + 1

    def test_windows_come_in_target_order(self, variation):
        positions = [e.position for e in evaluate_all(variation, _criteria(hit={"A1"}))]
        assert positions == sorted(positions)

    def test_a_window_longer_than_the_target_is_rejected(self, variation):
        criteria = _criteria(hit={"A1"}, length=len(TARGET_SEQ) + 1)
        with pytest.raises(MsaError):
            list(evaluate_all(variation, criteria))

    def test_a_window_outside_the_target_is_rejected(self, variation):
        with pytest.raises(IndexError):
            evaluate_window(variation, _criteria(hit={"A1"}), len(TARGET_SEQ))

    def test_an_allele_absent_from_the_alignment_is_rejected(self, variation):
        with pytest.raises(MsaError):
            evaluate_window(variation, _criteria(hit={"absent"}), 0)


class TestRegions:
    def test_overlapping_windows_merge(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2", "A3"})
        assert list(regions(evaluate_all(variation, criteria))) == [(5, 14)]

    def test_no_satisfying_window_gives_no_region(self, variation):
        criteria = _criteria(hit={"A1"}, spare={"A2", "A3"}, min_spare=2)
        assert list(regions(evaluate_all(variation, criteria))) == []

    def test_separated_windows_stay_separate(self):
        def window(position: int, satisfies: bool) -> WindowEvaluation:
            return WindowEvaluation(
                position=position,
                length=2,
                mismatches={},
                uncovered={},
                insertions={},
                max_hit_mismatches=0,
                min_spare_mismatches=None,
                hit_pass=satisfies,
                spare_pass=None,
                satisfies=satisfies,
            )

        stream = [window(0, True), window(1, False), window(5, True)]
        assert list(regions(stream)) == [(0, 2), (5, 7)]
