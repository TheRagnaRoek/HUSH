"""Tests for the mismatch discrimination model (msa.discrimination)."""

import pytest

from core.seed import SEED_END, SEED_START
from msa.discrimination import (
    CONSISTENT_POSITIONS,
    HUANG_RECOMMENDED_POSITIONS,
    MULTI_MISMATCH,
    STRONGEST_POSITIONS,
    GuideRegion,
    MismatchClass,
    assess_mismatch,
    best_mismatch,
    classify_mismatch,
    discrimination_key,
    guide_position_1based,
    guide_region,
    huang_recommends,
)
from msa.msa_errors import MsaError
from msa.windows import guide_position

LENGTH = 19


def _at_guide(guide_position_1: int, target_base: str, spared_base: str):
    """An assessment for a mismatch at a given 1-based guide position."""
    return assess_mismatch(
        LENGTH - guide_position_1, LENGTH, target_base=target_base, spared_base=spared_base
    )


# --------------------------------------------------------------------------- #
# Level 1: the published table
# --------------------------------------------------------------------------- #


class TestHuangTable:
    """Table 1 of Huang et al. 2009, reproduced row by row."""

    PUBLISHED = {
        ("A", "C"): {8, 9, 10, 11, 16},
        ("A", "G"): {10, 11, 12, 16},
        ("A", "T"): {9, 10, 14, 15, 16},
        ("C", "A"): {11, 12, 13, 14, 16},
        ("C", "G"): {9, 10, 12, 13, 16},
        ("C", "T"): {9, 10, 11, 12, 14},
        ("G", "A"): {10, 13, 15, 16},
        ("G", "C"): {10, 11, 12, 13},
        ("G", "T"): {9, 10, 11, 12, 14},
        ("T", "A"): {9, 11, 12, 15},
        ("T", "C"): {10, 11, 12, 13},
        ("T", "G"): {10, 11, 12, 13, 16},
    }

    def test_every_base_change_has_a_row(self):
        # Twelve rows: four bases changing into three others each.
        assert len(HUANG_RECOMMENDED_POSITIONS) == 12
        assert set(HUANG_RECOMMENDED_POSITIONS) == set(self.PUBLISHED)

    @pytest.mark.parametrize("key", sorted(PUBLISHED))
    def test_row_matches_the_publication(self, key):
        assert HUANG_RECOMMENDED_POSITIONS[key] == self.PUBLISHED[key]

    @pytest.mark.parametrize("key", sorted(PUBLISHED))
    def test_lookup_agrees_with_the_table(self, key):
        spared, target = key
        listed = {p for p in range(1, LENGTH + 1) if huang_recommends(spared, target, p)}
        assert listed == self.PUBLISHED[key]

    def test_p10_and_p16_dominate_the_table(self):
        # The convergence with Schwarz that makes those two positions credible.
        rows_with_10 = sum(1 for ps in HUANG_RECOMMENDED_POSITIONS.values() if 10 in ps)
        rows_with_16 = sum(1 for ps in HUANG_RECOMMENDED_POSITIONS.values() if 16 in ps)
        assert rows_with_10 == 10
        assert rows_with_16 == 7

    def test_rna_and_dna_spellings_agree(self):
        assert huang_recommends("T", "C", 10) == huang_recommends("U", "C", 10)

    def test_identical_bases_are_not_a_row(self):
        # Not a mismatch, so not a recommendation either way.
        assert huang_recommends("A", "A", 10) is None

    @pytest.mark.parametrize("base", ["N", "R", "Y"])
    def test_ambiguity_codes_yield_none(self, base):
        assert huang_recommends(base, "C", 10) is None
        assert huang_recommends("A", base, 10) is None


# --------------------------------------------------------------------------- #
# Coordinates
# --------------------------------------------------------------------------- #


class TestGuidePosition:
    def test_the_two_ends_swap(self):
        assert guide_position_1based(0, LENGTH) == LENGTH
        assert guide_position_1based(LENGTH - 1, LENGTH) == 1

    def test_stays_in_step_with_the_zero_based_conversion(self):
        # msa.windows counts from 0 for internal use, this module from 1 because
        # both papers do. If either changes, this fails rather than drifting.
        for position in range(LENGTH):
            assert (
                guide_position_1based(position, LENGTH)
                == guide_position(position, LENGTH) + 1
            )

    @pytest.mark.parametrize("position", [-1, LENGTH, 99])
    def test_positions_outside_the_window_are_rejected(self, position):
        with pytest.raises(MsaError):
            guide_position_1based(position, LENGTH)


class TestGuideRegion:
    def test_position_one_is_the_anchor(self):
        # Held in Argonaute's MID pocket; it does not pair with the target.
        assert guide_region(1) is GuideRegion.ANCHOR

    @pytest.mark.parametrize("position", range(SEED_START, SEED_END + 1))
    def test_the_seed_follows_core_seed(self, position):
        assert guide_region(position) is GuideRegion.SEED

    @pytest.mark.parametrize("position", [SEED_END + 1, 10, 16, 19])
    def test_everything_past_the_seed(self, position):
        assert guide_region(position) is GuideRegion.POST_SEED

    def test_zero_is_rejected(self):
        with pytest.raises(MsaError):
            guide_region(0)


# --------------------------------------------------------------------------- #
# Level 2: mismatch class
# --------------------------------------------------------------------------- #


class TestMismatchClass:
    @pytest.mark.parametrize(
        "guide,other,expected",
        [
            ("A", "A", MismatchClass.PURINE_PURINE),
            ("A", "G", MismatchClass.PURINE_PURINE),
            ("G", "A", MismatchClass.PURINE_PURINE),
            ("G", "G", MismatchClass.PURINE_PURINE),
            ("A", "C", MismatchClass.PURINE_PYRIMIDINE),
            ("G", "U", MismatchClass.PURINE_PYRIMIDINE),
            ("C", "A", MismatchClass.PYRIMIDINE_PURINE),
            ("U", "G", MismatchClass.PYRIMIDINE_PURINE),
            ("C", "C", MismatchClass.PYRIMIDINE_PYRIMIDINE),
            ("U", "U", MismatchClass.PYRIMIDINE_PYRIMIDINE),
        ],
    )
    def test_classification(self, guide, other, expected):
        assert classify_mismatch(guide, other) is expected

    def test_t_and_u_are_the_same_pyrimidine(self):
        assert classify_mismatch("T", "G") is classify_mismatch("U", "G")

    @pytest.mark.parametrize("base", ["N", "R", "Y"])
    def test_ambiguity_yields_none(self, base):
        assert classify_mismatch(base, "A") is None
        assert classify_mismatch("A", base) is None


class TestAssessment:
    def test_the_guide_carries_the_complement_of_the_target(self):
        # The siRNA is designed against the target, so it pairs with it and
        # mismatches the spared allele.
        assessment = _at_guide(10, target_base="C", spared_base="A")
        assert assessment.guide_base == "G"
        assert assessment.spared_base == "A"
        assert assessment.mismatch_class is MismatchClass.PURINE_PURINE

    def test_mutation_type_is_spared_to_target(self):
        assert _at_guide(10, target_base="C", spared_base="A").mutation_type == "A>C"

    def test_mutation_type_uses_dna_spelling(self):
        assert _at_guide(10, target_base="C", spared_base="U").mutation_type == "T>C"

    def test_window_and_guide_positions_are_both_kept(self):
        assessment = _at_guide(16, target_base="C", spared_base="A")
        assert assessment.guide_position == 16
        assert assessment.window_position == LENGTH - 16


class TestSourcesAgree:
    def test_disagreement_is_reported_not_resolved(self):
        # Huang lists guide position 8 for A>C; position 8 is inside the seed,
        # which Schwarz found to be the weak region.
        assessment = _at_guide(8, target_base="C", spared_base="A")
        assert assessment.huang_recommended is True
        assert assessment.region is GuideRegion.SEED
        assert assessment.sources_agree is False

    def test_agreement_where_both_point_the_same_way(self):
        assessment = _at_guide(16, target_base="C", spared_base="A")
        assert assessment.huang_recommended is True
        assert assessment.region is GuideRegion.POST_SEED
        assert assessment.sources_agree is True

    def test_unknown_never_counts_as_agreement(self):
        assert _at_guide(10, target_base="N", spared_base="A").sources_agree is False


# --------------------------------------------------------------------------- #
# Ranking
# --------------------------------------------------------------------------- #


class TestSchwarzPositionSets:
    def test_strongest_are_the_two_singled_out(self):
        assert STRONGEST_POSITIONS == {10, 16}

    def test_consistent_holds_the_remaining_named_positions(self):
        assert CONSISTENT_POSITIONS == {5, 9, 12, 13, 14}

    def test_the_two_sets_are_disjoint(self):
        assert not STRONGEST_POSITIONS & CONSISTENT_POSITIONS

    def test_p5_is_kept_although_it_lies_in_the_seed(self):
        # Reproduced as published; region is checked first, so it still ranks
        # with the seed.
        assert 5 in CONSISTENT_POSITIONS
        assert guide_region(5) is GuideRegion.SEED


class TestDiscriminationKey:
    def test_no_mismatch_sorts_lowest(self):
        assert discrimination_key([]) < discrimination_key(
            [_at_guide(4, target_base="C", spared_base="A")]
        )

    def test_two_mismatches_outrank_one_however_well_placed(self):
        two_in_seed = [
            _at_guide(4, target_base="C", spared_base="A"),
            _at_guide(6, target_base="C", spared_base="A"),
        ]
        one_ideal = [_at_guide(16, target_base="C", spared_base="A")]
        assert discrimination_key(two_in_seed) > discrimination_key(one_ideal)

    def test_the_threshold_is_two(self):
        assert MULTI_MISMATCH == 2
        one = [_at_guide(16, target_base="C", spared_base="A")]
        assert discrimination_key(one)[0] is False
        two = one + [_at_guide(12, target_base="C", spared_base="A")]
        assert discrimination_key(two)[0] is True

    def test_huang_outranks_every_schwarz_criterion(self):
        # Two cases that agree on every Schwarz criterion - same guide position,
        # same region, both purine:purine - and differ only in whether Huang's
        # table lists that position for the base change.
        listed = [_at_guide(16, target_base="C", spared_base="A")]  # A>C lists 16
        unlisted = [_at_guide(16, target_base="C", spared_base="G")]  # G>C does not
        assert huang_recommends("A", "C", 16) is True
        assert huang_recommends("G", "C", 16) is False
        assert listed[0].region is unlisted[0].region
        assert listed[0].mismatch_class is unlisted[0].mismatch_class
        assert discrimination_key(listed) > discrimination_key(unlisted)

    def test_post_seed_outranks_the_seed(self):
        post = [_at_guide(16, target_base="C", spared_base="A")]
        seed = [_at_guide(8, target_base="C", spared_base="A")]
        assert discrimination_key(post) > discrimination_key(seed)

    def test_strongest_positions_outrank_the_merely_consistent(self):
        strongest = [_at_guide(16, target_base="G", spared_base="C")]
        consistent = [_at_guide(13, target_base="G", spared_base="C")]
        assert discrimination_key(strongest) > discrimination_key(consistent)

    def test_more_mismatches_break_an_otherwise_equal_tie(self):
        two = [
            _at_guide(16, target_base="C", spared_base="A"),
            _at_guide(4, target_base="C", spared_base="A"),
        ]
        three = two + [_at_guide(5, target_base="C", spared_base="A")]
        assert discrimination_key(three) > discrimination_key(two)

    def test_the_key_has_the_documented_shape(self):
        key = discrimination_key([_at_guide(16, target_base="C", spared_base="A")])
        assert len(key) == 7
        assert all(isinstance(entry, bool) for entry in key[:6])
        assert isinstance(key[6], int)


class TestBestMismatch:
    def test_picks_the_best_placed_of_several(self):
        weak = _at_guide(4, target_base="C", spared_base="A")
        strong = _at_guide(16, target_base="C", spared_base="A")
        assert best_mismatch([weak, strong]) is strong
        assert best_mismatch([strong, weak]) is strong

    def test_none_for_an_empty_sequence(self):
        assert best_mismatch([]) is None
