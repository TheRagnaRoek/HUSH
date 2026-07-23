"""Tests for the classical Ui-Tei rule filter (core.rules)."""

from core.rules import (
    UI_TEI_MAX_GC_RUN,
    _longest_run,
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
