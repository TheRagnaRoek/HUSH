"""Tests for service.session.HushSession, independent of any renderer."""

from types import SimpleNamespace

import pytest

import service.session as session_mod
from core.errors import HushError
from core.offtarget import BruteForceSearcher
from core.sequence import MolType
from msa.targeting import TargetingMode
from service import AlleleTargetingSpec, DesignParams, HushSession
from service.results import ConvertResult, DesignResult, InspectResult
from shrna.convert import ShrnaOptions
from service.session import sort_selection

_ALLELE_1 = "HLA-DEMO*01:01"  # a member id in hla_panel_example.fasta
_ALLELE_2 = "HLA-DEMO*02:01"


# loading and inspection


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


# single-target design


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


# allele-aware design (panel targeting)


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


# alignment caching


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


# convert (in-memory, no TSV round-trip)


def test_convert_uses_last_design_selection(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    result = session.design()
    converted = session.convert()
    assert isinstance(converted, ConvertResult)
    # Every selected candidate is accounted for: converted, or withheld by a
    # conversion filter. Nothing disappears between the two stages.
    assert converted.total_converted == len(result.selected)
    assert len(converted.constructs) + len(converted.dropped) == len(result.selected)
    assert all(c.hairpin_rna.moltype is MolType.RNA for c in converted.constructs)
    assert all(c.premature_terminator_free for c in converted.constructs)
    assert not any(c.premature_terminator_free for c in converted.dropped)


def test_convert_without_design_raises(target_fasta):
    session = HushSession()
    session.load_target(target_fasta)
    with pytest.raises(HushError):
        session.convert()


# conversion filters


def test_convert_withholds_premature_terminators_by_default(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    session.design()
    filtered = session.convert()
    kept = session.convert(ShrnaOptions(drop_premature_terminators=False))
    # Opting out returns everything; the default is a strict subset of it.
    assert not kept.dropped
    assert len(kept.constructs) == filtered.total_converted
    assert len(filtered.constructs) < len(kept.constructs)
    assert all(c.premature_terminator_free for c in filtered.constructs)
    assert not any(c.premature_terminator_free for c in filtered.dropped)


def test_apply_convert_filters_is_a_partition():
    options = ShrnaOptions()
    constructs = [
        SimpleNamespace(premature_terminator_free=True),
        SimpleNamespace(premature_terminator_free=False),
        SimpleNamespace(premature_terminator_free=True),
    ]
    result = session_mod.apply_convert_filters(constructs, options)
    assert len(result.constructs) == 2
    assert len(result.dropped) == 1
    assert result.total_converted == 3


# sort_selection guards


def test_unknown_sort_mode_is_rejected(demo_lowgc_fasta):
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    with pytest.raises(ValueError):
        session.design(DesignParams(sort="nonesuch"))


def test_discrimination_sort_without_a_panel_is_rejected(demo_lowgc_fasta):
    # Nothing is spared without a panel, so there is no discrimination to rank.
    session = HushSession()
    session.load_target(demo_lowgc_fasta)
    with pytest.raises(ValueError):
        session.design(DesignParams(sort="discrimination"))


# off-target self resolution

# An IMGT-named target, as the bundled HLA files are: no RefSeq accession sits
# in the header, so nothing but homology can tell the gate what the design owns.
_IMGT_TARGET = "ATGCTTAGCCGATTACAGAGTCTAGGCATTACCGTAAGGCTTCAGAGGTCATTACGGCAT"
_IMGT_ALLELE = _IMGT_TARGET[:30] + "A" + _IMGT_TARGET[31:]  # one allelic difference
_SELF_TRANSCRIPT = "NM_900001.1"


class _StubBackend(BruteForceSearcher):
    """A brute-force searcher standing in for bowtie, batched interface included."""

    def precompute(self, reads):
        return self


@pytest.fixture
def stub_index(monkeypatch, tmp_path):
    """Point the session's bowtie backend at an in-memory two-transcript reference."""
    refs = [
        (_SELF_TRANSCRIPT, _IMGT_TARGET),  # the design's own transcript
        ("NM_900002.1", "TTTAAACCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTT"),
    ]
    monkeypatch.setattr(
        session_mod,
        "BowtieSearcher",
        lambda prefix, max_mismatch=3, bowtie_path="bowtie": _StubBackend(
            refs, max_mismatch=max_mismatch
        ),
    )
    return tmp_path / "index"


def _imgt_session(tmp_path, with_panel=True):
    path = tmp_path / "allele.fasta"
    path.write_text(f">HLA-A02:01:01:01\n{_IMGT_TARGET}\n", encoding="ascii")
    session = HushSession()
    session.load_target(path)
    if with_panel:
        panel = tmp_path / "panel.fasta"
        panel.write_text(f">HLA-A24:02:01:01\n{_IMGT_ALLELE}\n", encoding="ascii")
        session.load_panel(panel)
    return session


def _bare_params(stub_index, **overrides):
    # Only the off-target gate is under test here.
    return DesignParams(
        use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
        use_seed_filter=False, use_homopolymer_filter=False,
        offtarget_db=stub_index, **overrides,
    )


def test_imgt_target_resolves_its_own_transcript(stub_index, tmp_path):
    # The regression: an IMGT header yields no accession, so before resolution
    # the gate counted the target's own transcript against every candidate.
    session = _imgt_session(tmp_path)
    result = session.design(_bare_params(stub_index))
    own = result.criteria.offtarget_self_transcripts
    assert _SELF_TRANSCRIPT in own
    assert own.resolved  # found by homology, not by the header
    assert not own.declared
    assert result.selected


def test_panel_members_seed_resolution(stub_index, tmp_path):
    # A panel member is an intended target, so its transcripts are owned too.
    with_panel = _imgt_session(tmp_path).design(_bare_params(stub_index))
    assert _SELF_TRANSCRIPT in with_panel.criteria.offtarget_self_transcripts


def test_resolution_can_be_turned_off(stub_index, tmp_path):
    # Without resolution the gate is back to measuring on-target homology: the
    # target's own transcript is a 0-mismatch hit, so nothing survives.
    session = _imgt_session(tmp_path)
    result = session.design(
        _bare_params(stub_index, resolve_offtarget_self=False)
    )
    assert not result.criteria.offtarget_self_transcripts
    assert not result.selected


def test_declared_accessions_are_a_list(stub_index, tmp_path):
    session = _imgt_session(tmp_path)
    result = session.design(
        _bare_params(
            stub_index,
            offtarget_self_accessions=("NM_900001", "NM_900002.3"),
            resolve_offtarget_self=False,
        )
    )
    own = result.criteria.offtarget_self_transcripts
    assert len(own) == 2
    assert "NM_900002" in own  # version-insensitive
    assert result.selected


def test_self_fasta_declares_a_family_by_sequence(stub_index, tmp_path):
    # A self FASTA covers sequences that are neither the target nor the panel.
    family = tmp_path / "family.fasta"
    family.write_text(f">HLA-A01:01:01:01\n{_IMGT_TARGET}\n", encoding="ascii")
    session = _imgt_session(tmp_path, with_panel=False)
    result = session.design(
        _bare_params(stub_index, offtarget_self_fasta=family)
    )
    assert _SELF_TRANSCRIPT in result.criteria.offtarget_self_transcripts


def test_coverage_threshold_is_carried_on_the_criteria(stub_index, tmp_path):
    session = _imgt_session(tmp_path)
    result = session.design(_bare_params(stub_index, offtarget_self_coverage=0.4))
    assert result.criteria.offtarget_self_transcripts.min_coverage == 0.4


def test_refseq_header_still_contributes(stub_index, tmp_path):
    # A RefSeq-named target keeps the old header route, now alongside resolution.
    path = tmp_path / "refseq.fasta"
    path.write_text(
        f">NM_900001.8 Homo sapiens demo transcript, mRNA\n{_IMGT_TARGET}\n",
        encoding="ascii",
    )
    session = HushSession()
    session.load_target(path)
    result = session.design(
        _bare_params(stub_index, resolve_offtarget_self=False)
    )
    assert result.criteria.offtarget_self_transcripts.declared == {"nm_900001"}
