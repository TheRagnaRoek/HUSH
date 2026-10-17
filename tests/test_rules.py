"""Tests for the classical Ui-Tei and Reynolds rule filters (core.rules)."""

import pytest

from core.rules import (
    AMARZGUIOUI_DEFAULT_THRESHOLD,
    AMARZGUIOUI_DIFFERENTIAL_MAX,
    AMARZGUIOUI_DIFFERENTIAL_MIN,
    REYNOLDS_DEFAULT_THRESHOLD,
    UI_TEI_MAX_GC_RUN,
    _longest_run,
    amarzguioui,
    reynolds,
    ui_tei,
)
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate, enumerate_candidates


def make_candidate(sense: str, antisense: str) -> SiRNACandidate:
    """Candidate with strands set directly (not necessarily complementary), to isolate one rule."""
    return SiRNACandidate(
        position=0,
        length=len(sense),
        sense=BioSequence(sense, id="s", moltype=MolType.RNA),
        antisense=BioSequence(antisense, id="a", moltype=MolType.RNA),
    )


# An all-passing baseline: guide 5' = A (i), sense 5' = G (ii), guide[:7] all
# A/U (iii), no GC run > 9 (iv).
ALL_PASS_SENSE = "GAUAUAUAUAUAUAUAUAU"
ALL_PASS_GUIDE = "AUAUAUAUAUAUAUAUAUA"


class TestLongestRun:
    def test_empty(self):
        assert _longest_run("", frozenset("GC")) == 0

    def test_no_match(self):
        assert _longest_run("AUAUAU", frozenset("GC")) == 0

    def test_run_is_maximal_not_total(self):
        # Two separate GC runs of 3 and 2; longest is 3, not the total 5.
        assert _longest_run("GCGAUGC", frozenset("GC")) == 3


class TestAllRulesPass:
    def test_baseline_satisfied(self):
        r = ui_tei(make_candidate(ALL_PASS_SENSE, ALL_PASS_GUIDE))
        assert r.satisfied is True
        assert r.score == 4


class TestRuleIsolation:
    def test_rule_i_guide_5p_not_au_fails(self):
        r = ui_tei(make_candidate(ALL_PASS_SENSE, "GUAUAUAUAUAUAUAUAUA"))
        assert r.guide_5p_au is False
        assert r.satisfied is False
        assert r.score == 3

    def test_rule_ii_sense_5p_not_gc_fails(self):
        r = ui_tei(make_candidate("AAUAUAUAUAUAUAUAUAU", ALL_PASS_GUIDE))
        assert r.sense_5p_gc is False
        assert r.satisfied is False

    def test_rule_iii_guide_window_not_au_rich_fails(self):
        # guide[0] = A (passes i) but only 1 A/U in the first 7 (fails iii).
        r = ui_tei(make_candidate(ALL_PASS_SENSE, "AGCGCGCAUAUAUAUAUAU"))
        assert r.guide_5p_au is True
        assert r.guide_5p_au_rich is False
        assert r.satisfied is False

    def test_rule_iv_long_gc_run_fails(self):
        # sense[0] = G (passes ii) but a GC run of 10 (fails iv).
        r = ui_tei(make_candidate("GCGCGCGCGCAUAUAUAUA", ALL_PASS_GUIDE))
        assert r.sense_5p_gc is True
        assert r.no_long_gc_run is False
        assert r.satisfied is False


class TestRuleIVBoundary:
    def test_gc_run_of_exactly_max_passes(self):
        sense = "GC" * (UI_TEI_MAX_GC_RUN // 2) + "G" + "AUAUAUAUA"  # run of 9
        sense = sense[:19]
        r = ui_tei(make_candidate(sense, ALL_PASS_GUIDE))
        assert _longest_run(sense, frozenset("GC")) == UI_TEI_MAX_GC_RUN
        assert r.no_long_gc_run is True

    def test_gc_run_one_over_max_fails(self):
        sense = "GCGCGCGCGC" + "AUAUAUAUA"  # run of 10
        r = ui_tei(make_candidate(sense, ALL_PASS_GUIDE))
        assert _longest_run(sense, frozenset("GC")) == UI_TEI_MAX_GC_RUN + 1
        assert r.no_long_gc_run is False


class TestAmbiguousBases:
    def test_ambiguous_terminus_fails_base_identity_rules(self):
        r = ui_tei(make_candidate("NAUAUAUAUAUAUAUAUAU", "NUAUAUAUAUAUAUAUAUA"))
        assert r.sense_5p_gc is False   # N is not G/C
        assert r.guide_5p_au is False   # N is not A/U


class TestIntegrationWithEnumerator:
    def test_ui_tei_accepts_real_candidates(self):
        target = BioSequence("ATGGCTACGTTAGCCATTGCACGTAAC", id="t",
                             moltype=MolType.DNA)
        results = [ui_tei(c) for c in enumerate_candidates(target, length=19)]
        assert len(results) == 9
        assert all(0 <= r.score <= 4 for r in results)


# Reynolds

# A 19-mer that maximises the assessable Reynolds score (internal repeat is
# deferred). By 1-based position: pos3='A', pos10='U', pos13='A' (no penalty),
# pos19='A' (+1 and A/U, no G/C penalty), positions 15-19 all A/U (+5), and
# GC = 8/19 = 42% (in range). Positive points 5+1+1+1+1 = 9, no penalty.
REY_ALL_PASS = "CAAGCGCGCUUCAAAUAUA"


def rey_candidate(sense: str):
    return make_candidate(sense, "A" * len(sense))


def mutate(seq: str, index: int, base: str) -> str:
    chars = list(seq)
    chars[index] = base
    return "".join(chars)


class TestReynoldsBaseline:
    # The arithmetic tests fix the internal-repeat term (False) so they isolate
    # the motif/penalty scoring from the hairpin model, which is exercised
    # separately below and in test_hairpin.py.
    def test_all_positive_no_penalty(self):
        r = reynolds(rey_candidate(REY_ALL_PASS), no_internal_repeat=False)
        assert r.au_count_3p == 5
        assert (r.gc_in_range, r.a_at_pos19, r.a_at_pos3, r.u_at_pos10) == (
            True,
            True,
            True,
            True,
        )
        assert (r.g_at_pos13, r.gc_at_pos19) == (False, False)
        assert r.penalty == 0
        assert r.score == 9  # internal-repeat term fixed False
        assert r.max_possible == 10
        assert r.meets_threshold(REYNOLDS_DEFAULT_THRESHOLD) is True


class TestReynoldsInternalRepeat:
    def test_computed_by_default(self):
        # The baseline can only fold a marginal stem (Tm well below 20 C), so
        # the internal-repeat point is earned and the score reaches its ceiling.
        r = reynolds(rey_candidate(REY_ALL_PASS))
        assert r.no_internal_repeat is True
        assert r.score == 10
        assert r.max_possible == 10

    def test_override_true_adds_a_point_over_false(self):
        no = reynolds(rey_candidate(REY_ALL_PASS), no_internal_repeat=False)
        yes = reynolds(rey_candidate(REY_ALL_PASS), no_internal_repeat=True)
        assert yes.score == no.score + 1

    def test_strong_internal_hairpin_loses_the_point(self):
        # A 7-bp G/C stem closing a 3-nt loop melts far above 20 C, so this
        # candidate has internal structure and forfeits the point by default.
        r = reynolds(rey_candidate("GGGGGGGAAACCCCCCCAU"))
        assert r.no_internal_repeat is False


class TestReynoldsAdditive3pWindow:
    """Positions 15-19 contribute +1 PER A/U base, not a single threshold point."""

    def test_no_au_in_tail_is_zero(self):
        assert reynolds(rey_candidate("C" * 14 + "GCGCG")).au_count_3p == 0

    def test_all_five_au(self):
        assert reynolds(rey_candidate("C" * 14 + "AUAUA")).au_count_3p == 5

    def test_partial_count(self):
        # last five 'GAUCG' -> A and U only.
        assert reynolds(rey_candidate("C" * 14 + "GAUCG")).au_count_3p == 2


class TestReynoldsPosition19Gradient:
    """Position 19 scores A +2, U +1, G or C -1 (A counts for A-at-19 and the A/U window)."""

    def test_a_at_19_scores_two_over_g(self):
        a19 = reynolds(rey_candidate(REY_ALL_PASS), no_internal_repeat=False)
        g19 = reynolds(
            rey_candidate(mutate(REY_ALL_PASS, 18, "G")), no_internal_repeat=False
        )
        assert a19.a_at_pos19 is True and a19.gc_at_pos19 is False
        assert g19.gc_at_pos19 is True and g19.a_at_pos19 is False
        assert a19.au_count_3p == 5 and g19.au_count_3p == 4
        # A -> 9, G -> 6: the promised three-point spread at this position.
        assert a19.score == 9
        assert g19.score == 6

    def test_u_at_19_is_one_below_a(self):
        u19 = reynolds(
            rey_candidate(mutate(REY_ALL_PASS, 18, "U")), no_internal_repeat=False
        )
        assert (u19.a_at_pos19, u19.gc_at_pos19) == (False, False)
        assert u19.au_count_3p == 5  # U still counts toward the window
        assert u19.score == 8


class TestReynoldsPenalties:
    def test_g_at_pos13_deducts_a_point(self):
        r = reynolds(
            rey_candidate(mutate(REY_ALL_PASS, 12, "G")), no_internal_repeat=False
        )
        assert r.g_at_pos13 is True
        assert r.gc_in_range is True  # baseline has GC headroom for the swap
        assert r.penalty == 1
        assert r.score == 8

    def test_gc_at_pos19_deducts_a_point(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 18, "C")))
        assert r.gc_at_pos19 is True
        assert r.a_at_pos19 is False


class TestReynoldsCriterionIsolation:
    def test_gc_below_range_fails(self):
        r = reynolds(rey_candidate("UAAUAAUAAUAAUAAUAAU"))  # 0% GC
        assert r.gc_in_range is False

    def test_gc_above_range_fails(self):
        r = reynolds(rey_candidate("GCGCGCGCGCGCGCGCGAU"))  # ~89% GC
        assert r.gc_in_range is False

    def test_non_a_at_pos3_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 2, "C")))
        assert r.a_at_pos3 is False

    def test_non_u_at_pos10_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 9, "A")))
        assert r.u_at_pos10 is False


class TestReynoldsLongerCoreAnchoring:
    def test_3p_rules_anchor_to_actual_terminus_not_position_19(self):
        # 24-mer: index 18 is 'A' but the true 3' end (index 23) is 'G'.
        r = reynolds(rey_candidate("A" * 23 + "G"))
        assert r.a_at_pos19 is False
        assert r.gc_at_pos19 is True

    def test_24mer_ending_in_a_scores_the_terminus(self):
        r = reynolds(rey_candidate("A" * 24))
        assert r.a_at_pos19 is True
        assert r.gc_at_pos19 is False


class TestReynoldsLengthGuard:
    def test_short_core_raises(self):
        with pytest.raises(ValueError):
            reynolds(rey_candidate("A" * 18))


class TestReynoldsIntegration:
    def test_reynolds_over_real_candidates(self):
        target = BioSequence("ATGGCTACGTTAGCCATTGCACGTAACGGCATATGCC", id="t",
                             moltype=MolType.DNA)
        results = [reynolds(c) for c in enumerate_candidates(target, length=19)]
        assert results  # non-empty
        for r in results:
            assert -2 <= r.score <= r.max_possible == 10
            assert isinstance(r.no_internal_repeat, bool)


# Amarzguioui

# A 19-mer hitting all three positive motifs and no negatives:
#   pos1 = G  -> S1 (+1), not U1
#   pos6 = A  -> A6 (+1)
#   pos19 = A -> W19 (+1), not G19
AMARZ_ALL_POSITIVE = "GCCCCACCCCCCCCCCCCA"


def amarz_candidate(sense: str):
    return make_candidate(sense, "A" * len(sense))


class TestAmarzguiouiMotifs:
    def test_all_positive_motifs_score_three(self):
        r = amarzguioui(amarz_candidate(AMARZ_ALL_POSITIVE))
        assert (r.s1_gc, r.a6_a, r.w19_au) == (True, True, True)
        assert (r.u1, r.g19) == (False, False)
        assert r.motif_score == 3
        # Motif points alone already reach the threshold of 3.
        assert r.meets_threshold(AMARZGUIOUI_DEFAULT_THRESHOLD) is True

    def test_negative_motifs_subtract(self):
        # pos1 = U -> U1 (-1) and not S1; pos19 = G -> G19 (-1) and not W19.
        r = amarzguioui(amarz_candidate("UCCCCCCCCCCCCCCCCCG"))
        assert (r.u1, r.g19) == (True, True)
        assert (r.s1_gc, r.w19_au) == (False, False)
        assert r.motif_score == -2

    def test_s1_and_u1_are_mutually_exclusive_at_pos1(self):
        r_gc = amarzguioui(amarz_candidate(AMARZ_ALL_POSITIVE))
        assert r_gc.s1_gc is True and r_gc.u1 is False
        # pos1 = A: neither S1 nor U1.
        r_a = amarzguioui(amarz_candidate("ACCCCACCCCCCCCCCCCA"))
        assert r_a.s1_gc is False and r_a.u1 is False


class TestAmarzguiouiEndDifferential:
    """DT3 is measured 3' end MINUS 5' end - the sign is easy to invert."""

    def test_au_rich_3p_end_is_positive(self):
        # 5' triplet GGG (0 A/U), 3' triplet AAA (3) -> differential +3. This is
        # the functional orientation: the sense 3' end is the guide's 5' end.
        r = amarzguioui(amarz_candidate("GGG" + "C" * 13 + "AAA"))
        assert r.end_au_differential == 3

    def test_au_rich_5p_end_is_negative(self):
        # 5' triplet AAA (3 A/U), 3' triplet GGG (0) -> differential -3.
        r = amarzguioui(amarz_candidate("AAA" + "C" * 13 + "GGG"))
        assert r.end_au_differential == -3

    def test_symmetric_ends_give_zero(self):
        r = amarzguioui(amarz_candidate("AAC" + "C" * 13 + "CAA"))
        assert r.end_au_differential == 0

    def test_differential_enters_score_directly(self):
        r = amarzguioui(amarz_candidate("GGG" + "C" * 13 + "AAA"))
        assert r.score == r.motif_score + r.end_au_differential

    def test_stays_within_published_bounds(self):
        target = BioSequence(
            "ATGGCTACGTTAGCCATTGCACGTAACGGCATATGCC", id="t", moltype=MolType.DNA
        )
        for r in (amarzguioui(c) for c in enumerate_candidates(target, length=19)):
            assert (
                AMARZGUIOUI_DIFFERENTIAL_MIN
                <= r.end_au_differential
                <= AMARZGUIOUI_DIFFERENTIAL_MAX
            )


class TestAmarzguiouiPublishedTable3:
    """Reproduce all 46 scores of Amarzguioui and Prydz 2004 Table 3; this fixes the DT3 sign."""

    def test_every_row_reproduces(self, amarzguioui_table3):
        assert len(amarzguioui_table3) == 46
        for row in amarzguioui_table3:
            r = amarzguioui(amarz_candidate(row["sense"]))
            assert r.end_au_differential == row["dt3"], row["name"]
            assert r.score == row["score"], row["name"]

    def test_class_means_match_the_paper(self, amarzguioui_table3):
        # The paper reports mean scores of 2.6 (functional) and 0 (not).
        def mean(functional: bool) -> float:
            scores = [
                amarzguioui(amarz_candidate(r["sense"])).score
                for r in amarzguioui_table3
                if (r["functional"] == "yes") is functional
            ]
            return sum(scores) / len(scores)

        assert mean(True) == pytest.approx(2.6, abs=0.05)
        assert mean(False) == pytest.approx(0.0, abs=0.05)

    def test_threshold_enriches_for_functional_sirnas(self, amarzguioui_table3):
        selected = [
            r
            for r in amarzguioui_table3
            if amarzguioui(amarz_candidate(r["sense"])).meets_threshold()
        ]
        functional = [r for r in selected if r["functional"] == "yes"]
        # The paper: 11 of 14 selected from the test set are functional.
        assert len(selected) == 14
        assert len(functional) == 11


class TestAmarzguiouiLengthGuard:
    def test_short_core_raises(self):
        with pytest.raises(ValueError):
            amarzguioui(amarz_candidate("A" * 18))


class TestAmarzguiouiIntegration:
    def test_over_real_candidates(self):
        target = BioSequence("ATGGCTACGTTAGCCATTGCACGTAACGGCATATGCC", id="t",
                             moltype=MolType.DNA)
        results = [amarzguioui(c) for c in enumerate_candidates(target, length=19)]
        assert results
        for r in results:
            assert -2 <= r.motif_score <= 3
            assert -5 <= r.score <= 6
