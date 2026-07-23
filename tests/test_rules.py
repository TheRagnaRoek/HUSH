"""Tests for the classical Ui-Tei and Reynolds rule filters (core.rules)."""

import pytest

from core.rules import (
    REYNOLDS_DEFAULT_THRESHOLD,
    UI_TEI_MAX_GC_RUN,
    _longest_run,
    reynolds,
    ui_tei,
)
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate, enumerate_candidates


def make_candidate(sense: str, antisense: str) -> SiRNACandidate:
    """Build a candidate with strands set directly.

    The Ui-Tei rules read the two strands independently, so tests set sense and
    antisense by hand (they need not be true complements) to isolate one rule at
    a time without reasoning through reverse-complementation.
    """
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


# --------------------------------------------------------------------------- #
# Reynolds
# --------------------------------------------------------------------------- #

# A 19-mer satisfying all seven assessable Reynolds criteria (III is deferred):
#   idx2='A' (V), idx9='U' (VI), idx12='C' (VIII not G), idx18='A' (IV, VII),
#   3'-terminal 5 = UAUGA has 4 A/U (II), GC content 6/19 = 32% (I).
REY_ALL_PASS = "CAAUGACAUUAGCAUAUGA"


def rey_candidate(sense: str):
    return make_candidate(sense, "A" * len(sense))


def mutate(seq: str, index: int, base: str) -> str:
    chars = list(seq)
    chars[index] = base
    return "".join(chars)


class TestReynoldsBaseline:
    def test_all_assessable_pass_and_iii_pending(self):
        r = reynolds(rey_candidate(REY_ALL_PASS))
        assert r.score == 7
        assert r.max_possible == 7
        assert r.low_internal_structure is None
        assert r.pending == ["low_internal_structure"]
        assert r.meets_threshold(REYNOLDS_DEFAULT_THRESHOLD) is True


class TestReynoldsCriterionIIIDeferral:
    def test_supplying_true_adds_the_point(self):
        r = reynolds(rey_candidate(REY_ALL_PASS), low_internal_structure=True)
        assert r.max_possible == 8
        assert r.score == 8
        assert r.pending == []

    def test_supplying_false_counts_as_assessed_but_no_point(self):
        r = reynolds(rey_candidate(REY_ALL_PASS), low_internal_structure=False)
        assert r.max_possible == 8
        assert r.score == 7
        assert r.pending == []


class TestReynoldsCriterionIsolation:
    def test_i_gc_below_range_fails(self):
        r = reynolds(rey_candidate("UAAUAAUAAUAAUAAUAAU"))  # 0% GC
        assert r.gc_in_range is False

    def test_i_gc_above_range_fails(self):
        r = reynolds(rey_candidate("GCGCGCGCGCGCGCGCGAU"))  # ~89% GC
        assert r.gc_in_range is False

    def test_ii_au_poor_tail_fails(self):
        r = reynolds(rey_candidate(REY_ALL_PASS[:14] + "GCGCG"))  # tail 0 A/U
        assert r.au_rich_3p is False

    def test_iv_non_a_3p_end_fails_but_vii_independent(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 18, "U")))
        assert r.a_at_3p_end is False
        assert r.not_gc_at_3p_end is True  # U is not G/C

    def test_v_non_a_at_pos3_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 2, "C")))
        assert r.a_at_pos3 is False

    def test_vi_non_u_at_pos10_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 9, "A")))
        assert r.u_at_pos10 is False

    def test_vii_gc_3p_end_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 18, "C")))
        assert r.not_gc_at_3p_end is False
        assert r.a_at_3p_end is False

    def test_viii_g_at_pos13_fails(self):
        r = reynolds(rey_candidate(mutate(REY_ALL_PASS, 12, "G")))
        assert r.not_g_at_pos13 is False


class TestReynoldsLongerCoreAnchoring:
    def test_3p_rules_anchor_to_actual_terminus_not_position_19(self):
        # 24-mer: index 18 is 'A' but the true 3' end (index 23) is 'G'.
        ends_in_g = "A" * 23 + "G"
        r = reynolds(rey_candidate(ends_in_g))
        assert r.a_at_3p_end is False
        assert r.not_gc_at_3p_end is False

    def test_24mer_ending_in_a_passes_3p_rules(self):
        r = reynolds(rey_candidate("A" * 24))
        assert r.a_at_3p_end is True
        assert r.not_gc_at_3p_end is True


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
            assert 0 <= r.score <= r.max_possible == 7
            assert r.pending == ["low_internal_structure"]
