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
    longest_homopolymer_run,
    select,
)
from core.thermo import ThermoParams


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

    def test_rejects_non_thermoparams(self):
        with pytest.raises(TypeError, match="thermo_params"):
            SelectionCriteria(thermo_params="xia")  # type: ignore[arg-type]

    def test_rejects_non_positive_max_run(self):
        with pytest.raises(ValueError, match="max_homopolymer_run"):
            SelectionCriteria(use_homopolymer_filter=True, max_homopolymer_run=0)


class TestLongestHomopolymerRun:
    def test_counts_longest_run_and_reports_base(self):
        assert longest_homopolymer_run("ACGGGGUA") == (4, "G")

    def test_single_run_spanning_whole_string(self):
        assert longest_homopolymer_run("AAAA") == (4, "A")

    def test_no_repeats_gives_run_of_one(self):
        # All runs are length 1; the first base seen wins the tie.
        assert longest_homopolymer_run("ACGU") == (1, "A")

    def test_empty_string_is_zero(self):
        assert longest_homopolymer_run("") == (0, "")


class TestHomopolymerFilter:
    def _criteria(self, **kw):
        # Isolate the homopolymer gate from every other criterion.
        base = dict(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False, use_homopolymer_filter=True,
        )
        base.update(kw)
        return SelectionCriteria(**base)

    def test_disabled_leaves_fields_none(self):
        criteria = SelectionCriteria(
            use_seed_filter=False, use_homopolymer_filter=False
        )
        ev = evaluate(one_candidate("ACGUACGUACGUACGUACG"), criteria)
        assert ev.homopolymer_pass is None
        assert ev.homopolymer_run is None

    def test_run_over_limit_is_rejected(self):
        # A GGGG run (4) exceeds the default max of 3.
        ev = evaluate(one_candidate("ACGGGGUACGUACGUACGU"), self._criteria())
        assert ev.homopolymer_run == 4
        assert ev.homopolymer_pass is False
        assert ev.selected is False
        assert any("homopolymer" in r for r in ev.reasons)

    def test_run_at_limit_passes(self):
        # A GGG run (3) is exactly at the default max and is allowed.
        ev = evaluate(one_candidate("ACGGGUACGUACGUACGUA"), self._criteria())
        assert ev.homopolymer_run == 3
        assert ev.homopolymer_pass is True
        assert ev.selected is True

    def test_max_run_is_tunable(self):
        seq = "ACGGGUACGUACGUACGUA"  # longest run 3
        assert evaluate(seq_c := one_candidate(seq), self._criteria(max_homopolymer_run=2)).homopolymer_pass is False
        assert evaluate(seq_c, self._criteria(max_homopolymer_run=3)).homopolymer_pass is True


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
            use_homopolymer_filter=False,
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


class TestHairpinFilter:
    # A sense strand whose two GC arms fold into a stable stem; and one with no
    # complementary stem, so nothing folds.
    SELF_FOLDER = "GCGCGCAUUUUGCGCGCUU"
    NON_FOLDER = "A" * 19

    def _base(self, **kw):
        return SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False, use_homopolymer_filter=False, **kw,
        )

    def test_default_off_does_not_constrain(self):
        ev = evaluate(one_candidate(self.SELF_FOLDER), self._base())
        assert ev.hairpin_pass is None
        assert ev.hairpin_tm is None

    def test_rejects_stable_internal_structure(self):
        ev = evaluate(
            one_candidate(self.SELF_FOLDER),
            self._base(use_hairpin_filter=True, hairpin_tm_max=20.0),
        )
        assert ev.hairpin_tm is not None and ev.hairpin_tm > 20.0
        assert ev.hairpin_pass is False
        assert ev.selected is False
        assert any("internal structure" in r for r in ev.reasons)

    def test_passes_when_nothing_folds(self):
        ev = evaluate(
            one_candidate(self.NON_FOLDER),
            self._base(use_hairpin_filter=True),
        )
        assert ev.hairpin_tm is None
        assert ev.hairpin_pass is True
        assert ev.selected is True

    def test_threshold_is_customisable(self):
        c = one_candidate(self.SELF_FOLDER)
        lax = evaluate(c, self._base(use_hairpin_filter=True, hairpin_tm_max=200.0))
        strict = evaluate(c, self._base(use_hairpin_filter=True, hairpin_tm_max=0.0))
        assert lax.hairpin_pass is True
        assert strict.hairpin_pass is False


class TestAdvancedThermoKnobs:
    SENSE = "GCAUGCAUGCAUGCAUGCA"

    def _base(self, **kw):
        return SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False, **kw
        )

    def test_params_change_asymmetry_but_not_seed(self):
        # The parameter set is scoped to asymmetry/hairpin; the seed filter stays
        # pinned to its own (Freier) internals regardless of this choice.
        c = one_candidate(self.SENSE)
        xia = evaluate(c, self._base(
            use_seed_filter=True, use_asymmetry_filter=True,
            thermo_params=ThermoParams.XIA_1998,
        ))
        freier = evaluate(c, self._base(
            use_seed_filter=True, use_asymmetry_filter=True,
            thermo_params=ThermoParams.FREIER_1986,
        ))
        assert xia.asymmetry.differential != pytest.approx(
            freier.asymmetry.differential
        )
        assert xia.seed.max_tm == pytest.approx(freier.seed.max_tm)

    def test_asymmetry_window_flows_through(self):
        c = one_candidate(self.SENSE)
        w4 = evaluate(c, self._base(
            use_seed_filter=False, use_asymmetry_filter=True, asymmetry_window=4
        ))
        w6 = evaluate(c, self._base(
            use_seed_filter=False, use_asymmetry_filter=True, asymmetry_window=6
        ))
        assert w4.asymmetry.window == 4
        assert w6.asymmetry.window == 6

    def test_hairpin_uses_selected_parameter_set(self):
        c = one_candidate("GCGCGCAUUUUGCGCGCUU")
        xia = evaluate(c, self._base(
            use_seed_filter=False, use_hairpin_filter=True,
            hairpin_tm_max=500.0, thermo_params=ThermoParams.XIA_1998,
        ))
        freier = evaluate(c, self._base(
            use_seed_filter=False, use_hairpin_filter=True,
            hairpin_tm_max=500.0, thermo_params=ThermoParams.FREIER_1986,
        ))
        assert xia.hairpin_tm != pytest.approx(freier.hairpin_tm)


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
        # Default = rule union (ANY) + seed filter + homopolymer gate. The seed
        # filter is very strict on this GC-rich CDS; the point is that every
        # rejection still carries a reason.
        criteria = SelectionCriteria()
        evaluations = list(evaluate_all(hla_candidates, criteria))
        rejected = [e for e in evaluations if not e.selected]
        assert rejected  # this target does reject most candidates
        assert all(e.reasons for e in rejected)

    def test_dropping_the_seed_filter_selects_more(self, hla_candidates):
        default = SelectionCriteria()  # ANY rules + seed + homopolymer
        relaxed = SelectionCriteria(use_seed_filter=False)
        assert len(select(hla_candidates, relaxed)) >= len(
            select(hla_candidates, default)
        )
