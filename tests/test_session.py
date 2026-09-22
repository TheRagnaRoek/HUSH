"""Tests for the service layer (service.session.HushSession).

These exercise the orchestration seam directly - the same pipeline the CLI and
the future HTTP API drive - independently of any renderer. The CLI smoke tests
(test_cli) cover the terminal-facing behaviour on top of this.
"""

from types import SimpleNamespace

import pytest

import service.session as session_mod
from core.errors import HushError
from core.sequence import MolType
from msa.targeting import TargetingMode
from service import AlleleTargetingSpec, DesignParams, HushSession
from service.results import ConvertResult, DesignResult, InspectResult
from service.session import sort_selection

_ALLELE_1 = "HLA-DEMO*01:01"  # a member id in hla_panel_example.fasta
_ALLELE_2 = "HLA-DEMO*02:01"


# --------------------------------------------------------------------------- #
# loading and inspection
# --------------------------------------------------------------------------- #


def test_load_target_makes_it_current(target_fasta):
    session = HushSession()
    target = session.load_target(target_fasta)
    assert session.target is target
    assert target.is_nucleotide
    assert not session.has_panel


def test_load_panel_sets_has_panel(target_fasta, panel_fasta):
    session = HushSession()
    session.load_target(target_fasta)
    panel = session.load_panel(panel_fasta)
    assert session.has_panel
    assert len(panel) == 3


def test_inspect_returns_design_input_preview(target_fasta, panel_fasta):
    session = HushSession()
    session.load_target(target_fasta)
    session.load_panel(panel_fasta)
    result = session.inspect()
    assert isinstance(result, InspectResult)
    assert result.has_panel
    assert result.target.id == session.target.id


def test_design_without_target_raises(target_fasta):
    with pytest.raises(HushError):
        HushSession().design()


# --------------------------------------------------------------------------- #
# single-target design
# --------------------------------------------------------------------------- #


def test_design_bare_default_selects(demo_lowgc_fasta):
    # The shipped low-GC demo selects under the full default paradigm (seed on).
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    result = session.design()
    assert isinstance(result, DesignResult)
    assert result.evaluations  # candidates were enumerated
    assert result.selected  # and some cleared every gate
    # Selected is a subset of all evaluations, all with a passing verdict.
    assert all(e.selected for e in result.selected)


def test_design_offtarget_skipped_without_index(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    result = session.design()  # off-target on by default, no db configured
    assert result.offtarget_skipped is True
    assert result.criteria.use_offtarget_filter is False


def test_design_selected_is_position_ordered_by_default(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    positions = [e.candidate.position for e in session.design().selected]
    assert positions == sorted(positions)


def test_design_protein_target_raises_hush_error(mini_pdb):
    session = HushSession()
    session.load_target(mini_pdb)
    with pytest.raises(HushError):
        session.design()


def test_design_inverted_gc_window_raises_valueerror(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    params = DesignParams(use_gc_filter=True, gc_min=0.60, gc_max=0.30)
    with pytest.raises(ValueError):
        session.design(params)


# --------------------------------------------------------------------------- #
# allele-aware design (panel targeting)
# --------------------------------------------------------------------------- #


def _panel_session(target_fasta, panel_fasta):
    session = HushSession()
    session.load_target(target_fasta)
    session.load_panel(panel_fasta)
    return session


def test_design_bare_panel_is_pan_allele(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    result = session.design(DesignParams(use_seed_filter=False))
    assert result.targeting is not None
    assert result.targeting.mode is TargetingMode.PAN_ALLELE
    assert result.targeting.discriminating is False
    assert result.window_by_pos is not None


def test_design_allele_specific_discriminates(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(
        use_seed_filter=False,
        targeting=AlleleTargetingSpec(allele_specific=_ALLELE_1),
    )
    result = session.design(params)
    assert result.targeting.mode is TargetingMode.ALLELE_SPECIFIC
    assert result.targeting.discriminating is True
    assert result.targeting.hit == frozenset({_ALLELE_1})


def test_design_explicit_hit_and_spare(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(
        use_seed_filter=False,
        targeting=AlleleTargetingSpec(hit=(_ALLELE_1,), spare=(_ALLELE_2,)),
    )
    targeting = session.design(params).targeting
    assert targeting.hit == frozenset({_ALLELE_1})
    assert targeting.spare == frozenset({_ALLELE_2})


def test_design_target_specific_spares_the_whole_panel(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(
        use_seed_filter=False,
        targeting=AlleleTargetingSpec(target_specific=True),
    )
    targeting = session.design(params).targeting
    assert targeting.mode is TargetingMode.TARGET_SPECIFIC
    assert targeting.hit == frozenset()
    assert targeting.discriminating is True
    # Every panel member is spared.
    assert len(targeting.spare) == 3


def test_design_unknown_allele_raises_hush_error(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(targeting=AlleleTargetingSpec(allele_specific="NOT_IN_PANEL"))
    with pytest.raises(HushError):
        session.design(params)


def test_design_mutually_exclusive_targeting_raises_valueerror(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(
        targeting=AlleleTargetingSpec(pan_allele=True, allele_specific=_ALLELE_1)
    )
    with pytest.raises(ValueError):
        session.design(params)


def test_design_spare_without_hit_raises_valueerror(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    params = DesignParams(targeting=AlleleTargetingSpec(spare=(_ALLELE_1,)))
    with pytest.raises(ValueError):
        session.design(params)


def test_design_unknown_sort_raises_valueerror(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    with pytest.raises(ValueError):
        session.design(DesignParams(sort="nonsense"))


def test_design_sort_discrimination_needs_something_spared(target_fasta, panel_fasta):
    session = _panel_session(target_fasta, panel_fasta)
    # Pan-allele targeting spares nothing, so there is no discrimination to rank.
    params = DesignParams(use_seed_filter=False, sort="discrimination")
    with pytest.raises(ValueError):
        session.design(params)


def test_design_default_sort_is_auto(target_fasta, panel_fasta):
    # 'auto' must not raise for a pan-allele run (it falls back to position),
    # where explicit 'discrimination' would.
    session = _panel_session(target_fasta, panel_fasta)
    result = session.design(DesignParams(use_seed_filter=False))
    assert result.targeting.discriminating is False  # pan-allele default
    assert result.selected == sorted(
        result.selected, key=lambda e: e.candidate.position
    )


class TestSortSelectionRouting:
    """sort_selection's 'auto' routing, isolated from the pipeline."""

    @staticmethod
    def _cand(position):
        return SimpleNamespace(candidate=SimpleNamespace(position=position))

    def test_auto_ranks_by_discrimination_when_something_is_spared(self):
        a, b = self._cand(0), self._cand(1)
        window_by_pos = {
            0: SimpleNamespace(discrimination_key=(1,)),
            1: SimpleNamespace(discrimination_key=(2,)),
        }
        targeting = SimpleNamespace(discriminating=True)
        out = sort_selection([a, b], "auto", window_by_pos, targeting)
        # Higher key first: position 1 (key 2) outranks position 0 (key 1).
        assert [e.candidate.position for e in out] == [1, 0]

    def test_auto_falls_back_to_position_when_nothing_is_spared(self):
        a, b = self._cand(0), self._cand(1)
        targeting = SimpleNamespace(discriminating=False)
        out = sort_selection([b, a], "auto", None, targeting)
        assert out == [b, a]  # position sort leaves the input order untouched

    def test_explicit_position_keeps_target_order_even_when_discriminating(self):
        a, b = self._cand(0), self._cand(1)
        window_by_pos = {
            0: SimpleNamespace(discrimination_key=(1,)),
            1: SimpleNamespace(discrimination_key=(2,)),
        }
        targeting = SimpleNamespace(discriminating=True)
        out = sort_selection([a, b], "position", window_by_pos, targeting)
        assert out == [a, b]


# --------------------------------------------------------------------------- #
# alignment caching
# --------------------------------------------------------------------------- #


def test_alignment_computed_once_across_designs(target_fasta, panel_fasta, monkeypatch):
    session = _panel_session(target_fasta, panel_fasta)
    calls = {"n": 0}
    real_align = session_mod.align_star

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real_align(*args, **kwargs)

    monkeypatch.setattr(session_mod, "align_star", counting)
    session.design(DesignParams(use_seed_filter=False))
    session.design(DesignParams(use_seed_filter=False, length=21))
    assert calls["n"] == 1  # the second design reused the cached alignment


def test_reloading_panel_invalidates_cache(target_fasta, panel_fasta, monkeypatch):
    session = _panel_session(target_fasta, panel_fasta)
    calls = {"n": 0}
    real_align = session_mod.align_star

    def counting(*args, **kwargs):
        calls["n"] += 1
        return real_align(*args, **kwargs)

    monkeypatch.setattr(session_mod, "align_star", counting)
    session.design(DesignParams(use_seed_filter=False))
    session.load_panel(panel_fasta)  # a fresh panel must drop the cached alignment
    session.design(DesignParams(use_seed_filter=False))
    assert calls["n"] == 2


# --------------------------------------------------------------------------- #
# convert (in-memory, no TSV round-trip)
# --------------------------------------------------------------------------- #


def test_convert_uses_last_design_selection(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    result = session.design()
    converted = session.convert()
    assert isinstance(converted, ConvertResult)
    assert len(converted.constructs) == len(result.selected)
    assert all(c.hairpin_rna.moltype is MolType.RNA for c in converted.constructs)


def test_convert_without_design_raises(target_fasta):
    session = HushSession()
    session.load_target(target_fasta)
    with pytest.raises(HushError):
        session.convert()
