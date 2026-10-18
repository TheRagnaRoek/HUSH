"""Tests for the mismatch discrimination model (msa.discrimination)."""

import pytest

from core.seed import SEED_END, SEED_START
from msa.discrimination import (
    CONSISTENT_POSITIONS,
    DU_LOW_TOLERANCE,
    HUANG_RECOMMENDED_POSITIONS,
    LOW_TOLERANCE_WINDOW,
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
    ordered_by_placement,
)
from msa.msa_errors import MsaError
from msa.windows import guide_position

LENGTH = 19


def _at_guide(guide_position_1: int, target_base: str, spared_base: str):
    """An assessment for a mismatch at a given 1-based guide position."""
    return assess_mismatch(
        LENGTH - guide_position_1, LENGTH, target_base=target_base, spared_base=spared_base
    )


# Level 1: the published table


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


# Coordinates


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


# Level 2: mismatch class


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


# Ranking


class TestSchwarzPositionSets:
    def test_strongest_are_the_positions_replicated_in_every_assay(self):
        # P13 and P16 appear in all five of Schwarz's reported sets, P10 in four.
        assert STRONGEST_POSITIONS == {10, 13, 16}

    def test_consistent_holds_the_remaining_named_positions(self):
        assert CONSISTENT_POSITIONS == {9, 12, 14}

    def test_the_two_sets_are_disjoint(self):
        assert not STRONGEST_POSITIONS & CONSISTENT_POSITIONS

    def test_the_seed_positions_are_not_used(self):
        # Schwarz saw discrimination at P3, P4 and P5, but from one siRNA whose
        # seed pairing they call the most stable of their set, and which loaded
        # predominantly its sense strand. Neither set reproduces them.
        for position in (3, 4, 5):
            assert position not in STRONGEST_POSITIONS
            assert position not in CONSISTENT_POSITIONS
        assert guide_region(5) is GuideRegion.SEED


class TestLowToleranceWindow:
    """The window the three sources converge on."""

    def test_du_region_is_p9_to_p15(self):
        # Du's target-site positions 5-11, converted by guide = 20 - target.
        assert DU_LOW_TOLERANCE == {9, 10, 11, 12, 13, 14, 15}

    def test_the_conversion_matches_the_published_sequences(self):
        # Du print both strands; the guide is the exact reverse complement of
        # the target site, which is what fixes guide = 20 - target.
        target = "CTTATTGGAGAGAGCACGA"
        guide = "TCGTGCTCTCTCCAATAAG"
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        assert "".join(complement[b] for b in reversed(target)) == guide
        assert len(target) == 19
        for target_position in (5, 11):  # the bounds Du reports
            assert 20 - target_position in DU_LOW_TOLERANCE

    def test_window_is_du_plus_p16(self):
        # P16 sits in Du's medium-tolerance region but is Schwarz's
        # best-replicated position and is listed in 7 of Huang's 12 rows.
        assert LOW_TOLERANCE_WINDOW == DU_LOW_TOLERANCE | {16}
        assert 16 not in DU_LOW_TOLERANCE

    def test_huang_table_sits_almost_entirely_inside_the_window(self):
        entries = [p for ps in HUANG_RECOMMENDED_POSITIONS.values() for p in ps]
        assert len(entries) == 55
        assert sum(1 for p in entries if p in LOW_TOLERANCE_WINDOW) == 54

    def test_every_strongest_position_lies_inside_the_window(self):
        assert STRONGEST_POSITIONS <= LOW_TOLERANCE_WINDOW
        assert CONSISTENT_POSITIONS <= LOW_TOLERANCE_WINDOW


class TestDiscriminationKey:
    """Each test isolates one criterion by holding every other one equal."""

    def test_no_mismatch_sorts_lowest(self):
        assert discrimination_key([]) < discrimination_key(
            [_at_guide(4, target_base="C", spared_base="A")]
        )

    def test_huang_outranks_every_other_criterion(self):
        # Same position, same region, same class; only the table listing differs.
        listed = [_at_guide(16, target_base="C", spared_base="A")]  # A>C lists 16
        unlisted = [_at_guide(16, target_base="C", spared_base="G")]  # G>C does not
        assert huang_recommends("A", "C", 16) is True
        assert huang_recommends("G", "C", 16) is False
        assert listed[0].region is unlisted[0].region
        assert listed[0].mismatch_class is unlisted[0].mismatch_class
        assert discrimination_key(listed) > discrimination_key(unlisted)

    def test_the_window_outranks_a_position_merely_3prime_of_the_seed(self):
        # P14 and P17 are both post-seed, neither is a strongest position,
        # G>C lists neither, and both are purine:purine. Only the window differs.
        inside = [_at_guide(14, target_base="C", spared_base="G")]
        outside = [_at_guide(17, target_base="C", spared_base="G")]
        assert inside[0].region is outside[0].region
        assert inside[0].mismatch_class is outside[0].mismatch_class
        assert 14 in LOW_TOLERANCE_WINDOW and 17 not in LOW_TOLERANCE_WINDOW
        assert discrimination_key(inside) > discrimination_key(outside)

    def test_strongest_positions_outrank_the_rest_of_the_window(self):
        # T>A lists neither P13 nor P14; both post-seed, both in the window,
        # both pyrimidine:pyrimidine.
        strongest = [_at_guide(13, target_base="A", spared_base="T")]
        merely_inside = [_at_guide(14, target_base="A", spared_base="T")]
        assert huang_recommends("T", "A", 13) is False
        assert huang_recommends("T", "A", 14) is False
        assert discrimination_key(strongest) > discrimination_key(merely_inside)

    def test_post_seed_outranks_the_seed_outside_the_window(self):
        # Both outside the window and unlisted, both purine:purine - so the
        # region is the only thing left to separate them.
        post = [_at_guide(17, target_base="C", spared_base="G")]
        seed = [_at_guide(7, target_base="C", spared_base="G")]
        assert post[0].mismatch_class is seed[0].mismatch_class
        assert discrimination_key(post) > discrimination_key(seed)

    def test_purine_purine_decides_last(self):
        # Same seed position, same region, neither listed; only the class differs.
        purine = [_at_guide(7, target_base="C", spared_base="G")]
        pyrimidine = [_at_guide(7, target_base="A", spared_base="C")]
        assert purine[0].mismatch_class is MismatchClass.PURINE_PURINE
        assert pyrimidine[0].mismatch_class is MismatchClass.PYRIMIDINE_PYRIMIDINE
        assert purine[0].region is pyrimidine[0].region
        assert discrimination_key(purine) > discrimination_key(pyrimidine)

    def test_p5_is_no_longer_credited_above_its_neighbours(self):
        # P5 used to carry a tier of its own. It now ranks exactly like any other
        # seed position with the same class.
        p5 = [_at_guide(5, target_base="C", spared_base="G")]
        p4 = [_at_guide(4, target_base="C", spared_base="G")]
        assert discrimination_key(p5) == discrimination_key(p4)


class TestCountIsNotRanked:
    """The number of differences must never enter the key - no source measures it
    against placement. This is the property the whole design note defends."""

    def test_one_well_placed_difference_outranks_three_poorly_placed_ones(self):
        # Not an endorsement of the one-difference window: the ranking simply
        # does not look at the count, and the renderers report it separately.
        one = [_at_guide(16, target_base="C", spared_base="A")]
        three = [
            _at_guide(11, target_base="C", spared_base="G"),
            _at_guide(15, target_base="C", spared_base="G"),
            _at_guide(17, target_base="C", spared_base="G"),
        ]
        assert discrimination_key(one) > discrimination_key(three)

    def test_adding_differences_never_changes_the_key(self):
        one = [_at_guide(16, target_base="C", spared_base="A")]
        many = one + [
            _at_guide(4, target_base="C", spared_base="A"),
            _at_guide(5, target_base="C", spared_base="A"),
        ]
        assert discrimination_key(many) == discrimination_key(one)

    def test_the_key_has_the_documented_shape(self):
        key = discrimination_key([_at_guide(16, target_base="C", spared_base="A")])
        assert len(key) == 5
        assert all(isinstance(entry, bool) for entry in key)

    def test_multi_mismatch_is_reporting_only(self):
        # Still exported as the reporting flag (the two-mismatch convention of
        # Naito et al. 2009), but it must not appear anywhere in the key.
        assert MULTI_MISMATCH == 2
        one = [_at_guide(16, target_base="C", spared_base="A")]
        two = one + [_at_guide(12, target_base="C", spared_base="A")]
        assert discrimination_key(one) == discrimination_key(two)


class TestOrderedByPlacement:
    def test_best_placed_comes_first(self):
        weak = _at_guide(4, target_base="C", spared_base="A")
        strong = _at_guide(16, target_base="C", spared_base="A")
        assert ordered_by_placement([weak, strong])[0] is strong

    def test_every_assessment_is_kept(self):
        given = [
            _at_guide(4, target_base="C", spared_base="A"),
            _at_guide(16, target_base="C", spared_base="A"),
            _at_guide(11, target_base="C", spared_base="G"),
        ]
        assert sorted(a.guide_position for a in ordered_by_placement(given)) == [
            4,
            11,
            16,
        ]

    def test_it_leads_with_the_mismatch_the_key_rests_on(self):
        given = [
            _at_guide(4, target_base="C", spared_base="A"),
            _at_guide(16, target_base="C", spared_base="A"),
        ]
        assert ordered_by_placement(given)[0] is best_mismatch(given)

    def test_empty_stays_empty(self):
        assert ordered_by_placement([]) == ()


class TestBestMismatch:
    def test_picks_the_best_placed_of_several(self):
        weak = _at_guide(4, target_base="C", spared_base="A")
        strong = _at_guide(16, target_base="C", spared_base="A")
        assert best_mismatch([weak, strong]) is strong
        assert best_mismatch([strong, weak]) is strong

    def test_none_for_an_empty_sequence(self):
        assert best_mismatch([]) is None
