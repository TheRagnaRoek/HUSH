"""Tests for the candidate selection filter (core.selection)."""

import pytest

from core.inputs import load_sequences
from core.rules import amarzguioui, reynolds, ui_tei
from core.seed import seed_tm
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates
from core.selection import (
    CandidateEvaluation,
    RuleMode,
    SelectionCriteria,
    evaluate,
    evaluate_all,
    select,
)


def one_candidate(sense: str):
    core = BioSequence(sense, id="c", moltype=MolType.RNA)
    return next(iter(enumerate_candidates(core, length=len(sense))))


@pytest.fixture
def hla_candidates(hla_a01_fasta):
    target = load_sequences(hla_a01_fasta)[0]
    return list(enumerate_candidates(target))


class TestCriteriaValidation:
    def test_rejects_non_rulemode(self):
        with pytest.raises(TypeError, match="rule_mode"):
            SelectionCriteria(rule_mode="all")  # type: ignore[arg-type]

    def test_rejects_inverted_gc_window(self):
        with pytest.raises(ValueError, match="gc_min <= gc_max"):
            SelectionCriteria(use_gc_filter=True, gc_min=0.6, gc_max=0.3)

    def test_enabled_rule_names_tracks_toggles(self):
        c = SelectionCriteria(
            use_ui_tei=True, use_reynolds=False, use_amarzguioui=True
        )
        assert c.enabled_rule_names == ("ui_tei", "amarzguioui")


class TestOnlyEnabledCriteriaAreComputed:
    def test_disabled_criteria_are_none(self):
        criteria = SelectionCriteria(
            use_ui_tei=True,
            use_reynolds=False,
            use_amarzguioui=False,
            use_seed_filter=False,
        )
        ev = evaluate(one_candidate("GAUAUAUAUAUAUAUAUAU"), criteria)
        assert ev.ui_tei is not None
        assert ev.reynolds is None
        assert ev.amarzguioui is None
        assert ev.seed is None
        assert ev.asymmetry is None

    def test_no_rules_and_no_filters_selects_everything(self):
        criteria = SelectionCriteria(
            use_ui_tei=False,
            use_reynolds=False,
            use_amarzguioui=False,
            use_seed_filter=False,
        )
        ev = evaluate(one_candidate("GCGCGCGCGCGCGCGCGCG"), criteria)
        assert ev.rule_pass is True  # no rules -> no constraint
        assert ev.selected is True


class TestRuleModes:
    # A candidate contrived so the three rule sets disagree, checked against the
    # rule functions directly so the mode logic is what is under test.
    SENSE = "GAUACAUUGGCAUAUAUGA"

    def _votes(self):
        c = one_candidate(self.SENSE)
        return (
            ui_tei(c).satisfied,
            reynolds(c).meets_threshold(),
            amarzguioui(c).meets_threshold(),
        )

    def test_all_mode_needs_every_enabled_rule(self):
        criteria = SelectionCriteria(rule_mode=RuleMode.ALL, use_seed_filter=False)
        ev = evaluate(one_candidate(self.SENSE), criteria)
        assert ev.rule_pass is all(self._votes())

    def test_any_mode_needs_one_enabled_rule(self):
        criteria = SelectionCriteria(rule_mode=RuleMode.ANY, use_seed_filter=False)
        ev = evaluate(one_candidate(self.SENSE), criteria)
        assert ev.rule_pass is any(self._votes())

    def test_any_is_never_stricter_than_all(self, hla_candidates):
        all_c = SelectionCriteria(rule_mode=RuleMode.ALL, use_seed_filter=False)
        any_c = SelectionCriteria(rule_mode=RuleMode.ANY, use_seed_filter=False)
        selected_all = {e.candidate.position for e in select(hla_candidates, all_c)}
        selected_any = {e.candidate.position for e in select(hla_candidates, any_c)}
        assert selected_all <= selected_any

    def test_single_algorithm_matches_that_algorithm(self, hla_candidates):
        # With only Reynolds enabled, selection is exactly Reynolds' own verdict.
        criteria = SelectionCriteria(
            use_ui_tei=False,
            use_reynolds=True,
            use_amarzguioui=False,
            use_seed_filter=False,
        )
        selected = {e.candidate.position for e in select(hla_candidates, criteria)}
        expected = {
            c.position for c in hla_candidates if reynolds(c).meets_threshold()
        }
        assert selected == expected


class TestIndependentFilters:
    def test_gc_filter_is_a_separate_gate(self):
        # A 100% GC core fails a tight GC window regardless of the rules.
        criteria = SelectionCriteria(
            use_ui_tei=False,
            use_reynolds=False,
            use_amarzguioui=False,
            use_seed_filter=False,
            use_gc_filter=True,
            gc_min=0.30,
            gc_max=0.52,
        )
        ev = evaluate(one_candidate("GCGCGCGCGCGCGCGCGCG"), criteria)
        assert ev.gc_pass is False
        assert ev.selected is False
        assert any("GC" in r for r in ev.reasons)

    def test_gc_filter_default_off_does_not_constrain(self):
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False,
        )
        ev = evaluate(one_candidate("GCGCGCGCGCGCGCGCGCG"), criteria)
        assert ev.gc_pass is None

    def test_seed_filter_matches_seed_module(self, hla_candidates):
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=True,
        )
        for c in hla_candidates[:40]:
            ev = evaluate(c, criteria)
            assert ev.seed_pass is seed_tm(c).passes

    def test_seed_threshold_is_customisable(self):
        # A stricter threshold can only reject more (or the same).
        strong = one_candidate("GCGCGCGCGCGCGCGCGCG")
        lax = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=True, seed_tm_threshold=99.0,
        )
        strict = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=True, seed_tm_threshold=0.0,
        )
        assert evaluate(strong, lax).seed_pass is True
        assert evaluate(strong, strict).seed_pass is False

    def test_filters_and_on_top_of_rules(self):
        # Passes the (disabled) rule step but fails the seed filter -> rejected.
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=True, seed_tm_threshold=0.0,
        )
        ev = evaluate(one_candidate("GCGCGCGCGCGCGCGCGCG"), criteria)
        assert ev.rule_pass is True
        assert ev.seed_pass is False
        assert ev.selected is False


class TestSelectPreservesOrderAndDropsNothingSilently:
    def test_select_returns_subset_in_position_order(self, hla_candidates):
        criteria = SelectionCriteria(rule_mode=RuleMode.ANY)
        selected = select(hla_candidates, criteria)
        positions = [e.candidate.position for e in selected]
        assert positions == sorted(positions)
        assert all(isinstance(e, CandidateEvaluation) for e in selected)
        assert all(e.selected for e in selected)

    def test_evaluate_all_covers_every_candidate(self, hla_candidates):
        criteria = SelectionCriteria()
        evaluations = list(evaluate_all(hla_candidates, criteria))
        assert len(evaluations) == len(hla_candidates)
        selected = [e for e in evaluations if e.selected]
        assert len(select(hla_candidates, criteria)) == len(selected)


class TestReasons:
    def test_selected_has_no_reasons(self):
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False,
        )
        ev = evaluate(one_candidate("GAUAUAUAUAUAUAUAUAU"), criteria)
        assert ev.selected is True
        assert ev.reasons == ()

    def test_rejected_lists_each_failed_gate(self):
        criteria = SelectionCriteria(
            use_ui_tei=True, use_reynolds=False, use_amarzguioui=False,
            rule_mode=RuleMode.ALL,
            use_gc_filter=True, gc_min=0.0, gc_max=0.1,
            use_seed_filter=True, seed_tm_threshold=0.0,
        )
        # A GC-rich core: fails Ui-Tei, the GC window and the seed filter at once.
        ev = evaluate(one_candidate("GCGCGCGCGCGCGCGCGCG"), criteria)
        assert ev.selected is False
        assert len(ev.reasons) == 3


class TestRealHlaTarget:
    def test_default_criteria_are_strict_but_nonempty_reasons_visible(
        self, hla_candidates
    ):
        # Default = all three rules (ALL) + seed filter. On this GC-rich CDS that
        # is very strict; the point is that every rejection carries a reason.
        criteria = SelectionCriteria()
        evaluations = list(evaluate_all(hla_candidates, criteria))
        rejected = [e for e in evaluations if not e.selected]
        assert rejected  # this target does reject most candidates
        assert all(e.reasons for e in rejected)

    def test_relaxing_to_any_rule_without_seed_selects_more(self, hla_candidates):
        strict = SelectionCriteria()  # ALL rules + seed
        relaxed = SelectionCriteria(rule_mode=RuleMode.ANY, use_seed_filter=False)
        assert len(select(hla_candidates, relaxed)) >= len(
            select(hla_candidates, strict)
        )
