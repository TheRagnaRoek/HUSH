"""Tests for the MSA / allele module (msa.targeting, msa.msa_errors)."""

from dataclasses import fields

import pytest

from core.errors import HushError
from core.sequence import BioSequence, MolType, SequencePanel
from msa.msa_errors import MsaError, PanelError, TargetingError
from msa.targeting import AlleleTargeting, TargetingMode

A2 = "HLA-A*02:01"
A1 = "HLA-A*01:01"
A3 = "HLA-A*03:01"


@pytest.fixture
def panel() -> SequencePanel:
    """A three-member allele panel differing only in the last base."""
    return SequencePanel(
        "hla-a",
        [
            BioSequence("ATGGCCGTCATGGCGCCCCGA", id=A2, moltype=MolType.DNA),
            BioSequence("ATGGCCGTCATGGCGCCCCGT", id=A1, moltype=MolType.DNA),
            BioSequence("ATGGCCGTCATGGCGCCCCGC", id=A3, moltype=MolType.DNA),
        ],
    )


class TestNamedConstructors:
    def test_allele_specific_spares_every_other_member(self, panel):
        t = AlleleTargeting.allele_specific(A2, panel)
        assert t.hit == {A2}
        assert t.spare == {A1, A3}

    def test_pan_allele_hits_everything_and_spares_nothing(self, panel):
        t = AlleleTargeting.pan_allele(panel)
        assert t.hit == {A1, A2, A3}
        assert t.spare == frozenset()

    def test_subset_spares_the_remainder(self, panel):
        t = AlleleTargeting.subset([A2, A1], panel)
        assert t.hit == {A2, A1}
        assert t.spare == {A3}

    def test_target_specific_spares_the_whole_panel(self, panel):
        # The loaded target is silenced by construction, so specificity means
        # hitting nothing in the panel and sparing all of it.
        t = AlleleTargeting.target_specific(panel)
        assert t.hit == frozenset()
        assert t.spare == {A1, A2, A3}
        assert t.mode is TargetingMode.TARGET_SPECIFIC
        assert t.discriminating

    def test_target_specific_is_the_opposite_of_pan_allele(self, panel):
        target = AlleleTargeting.target_specific(panel)
        pan = AlleleTargeting.pan_allele(panel)
        assert target.hit == pan.spare
        assert target.spare == pan.hit

    def test_constructors_assign_every_panel_member(self, panel):
        # Nothing is left unconstrained when a request is built from a panel.
        for t in (
            AlleleTargeting.allele_specific(A2, panel),
            AlleleTargeting.pan_allele(panel),
            AlleleTargeting.subset([A2], panel),
            AlleleTargeting.target_specific(panel),
        ):
            assert t.unconstrained(panel) == frozenset()

    def test_subset_of_the_whole_panel_is_pan_allele(self, panel):
        assert AlleleTargeting.subset(panel.ids, panel) == AlleleTargeting.pan_allele(panel)


class TestModeIsDerived:
    @pytest.mark.parametrize(
        "hit,spare,expected",
        [
            ({A2}, {A1}, TargetingMode.ALLELE_SPECIFIC),
            ({A2}, set(), TargetingMode.PAN_ALLELE),
            ({A2, A1}, set(), TargetingMode.PAN_ALLELE),
            ({A2, A1}, {A3}, TargetingMode.SUBSET),
            (set(), {A1}, TargetingMode.TARGET_SPECIFIC),
            (set(), {A1, A2, A3}, TargetingMode.TARGET_SPECIFIC),
        ],
    )
    def test_mode(self, hit, spare, expected):
        assert AlleleTargeting(hit=hit, spare=spare).mode is expected

    def test_mode_is_not_a_stored_field(self):
        # mode is a property computed from hit/spare, so it can never disagree
        # with them. If it ever becomes a field, this test should fail loudly.
        assert "mode" not in {f.name for f in fields(AlleleTargeting)}

    def test_same_hit_different_spare_gives_different_modes(self):
        assert AlleleTargeting(hit={A2, A1}).mode is TargetingMode.PAN_ALLELE
        assert AlleleTargeting(hit={A2, A1}, spare={A3}).mode is TargetingMode.SUBSET

    def test_discriminating_matches_the_presence_of_spared_alleles(self):
        assert not AlleleTargeting(hit={A2}).discriminating
        assert AlleleTargeting(hit={A2}, spare={A1}).discriminating


class TestValidation:
    def test_empty_hit_and_spare_raises(self):
        # A request that constrains nothing has no meaning; an empty hit is only
        # valid when something is spared (target-specific).
        with pytest.raises(TargetingError):
            AlleleTargeting(hit=set())

    def test_empty_hit_with_spare_is_allowed(self):
        t = AlleleTargeting(hit=set(), spare={A1})
        assert t.hit == frozenset()
        assert t.spare == {A1}
        assert t.mode is TargetingMode.TARGET_SPECIFIC

    def test_allele_in_both_sets_raises(self):
        with pytest.raises(TargetingError):
            AlleleTargeting(hit={A2}, spare={A2, A1})

    def test_requested_allele_absent_from_panel_raises(self, panel):
        with pytest.raises(TargetingError):
            AlleleTargeting.allele_specific("HLA-B*07:02", panel)

    def test_subset_reports_every_missing_allele(self, panel):
        with pytest.raises(TargetingError) as excinfo:
            AlleleTargeting.subset([A2, "HLA-B*07:02"], panel)
        assert "HLA-B*07:02" in str(excinfo.value)

    def test_duplicate_panel_ids_raise(self):
        dup = SequencePanel(
            "dup",
            [
                BioSequence("ATGC", id=A2, moltype=MolType.DNA),
                BioSequence("ATGG", id=A2, moltype=MolType.DNA),
            ],
        )
        with pytest.raises(PanelError):
            AlleleTargeting.pan_allele(dup)

    def test_blank_panel_ids_raise(self):
        blank = SequencePanel("blank", [BioSequence("ATGC", moltype=MolType.DNA)])
        with pytest.raises(PanelError):
            AlleleTargeting.pan_allele(blank)


class TestNormalisation:
    @pytest.mark.parametrize("ids", [[A2, A1], (A2, A1), {A2, A1}, iter([A2, A1])])
    def test_any_iterable_of_ids_is_accepted(self, ids):
        assert AlleleTargeting(hit=ids).hit == {A2, A1}

    def test_sets_are_stored_as_frozensets(self):
        t = AlleleTargeting(hit=[A2], spare=[A1])
        assert isinstance(t.hit, frozenset)
        assert isinstance(t.spare, frozenset)

    def test_duplicates_in_the_input_collapse(self):
        assert AlleleTargeting(hit=[A2, A2, A1]).hit == {A2, A1}


class TestImmutability:
    def test_fields_cannot_be_reassigned(self):
        t = AlleleTargeting(hit={A2})
        with pytest.raises(Exception):  # dataclasses.FrozenInstanceError
            t.hit = frozenset({A1})

    def test_equal_requests_compare_and_hash_equal(self):
        a = AlleleTargeting(hit=[A2, A1], spare=[A3])
        b = AlleleTargeting(hit=[A1, A2], spare=[A3])
        assert a == b
        assert hash(a) == hash(b)
        assert len({a, b}) == 1


class TestErrorHierarchy:
    @pytest.mark.parametrize("exc", [PanelError, TargetingError])
    def test_module_errors_are_catchable_as_msa_and_hush_errors(self, exc):
        assert issubclass(exc, MsaError)
        assert issubclass(exc, HushError)

    def test_targeting_failure_is_caught_by_except_hush_error(self, panel):
        with pytest.raises(HushError):
            AlleleTargeting.allele_specific("HLA-B*07:02", panel)


class TestScaling:
    """An IPD-IMGT/HLA locus can hold thousands of alleles; panel checks must
    stay linear in panel size."""

    def test_large_panel_is_handled(self):
        big = SequencePanel(
            "big",
            [
                BioSequence("ATGGCCGTCATGGCGCCCCGA", id=f"HLA-A*{i:05d}", moltype=MolType.DNA)
                for i in range(2000)
            ],
        )
        t = AlleleTargeting.allele_specific("HLA-A*00042", big)
        assert t.hit == {"HLA-A*00042"}
        assert len(t.spare) == 1999

    def test_error_message_is_truncated_for_large_id_lists(self):
        big = SequencePanel(
            "big",
            [
                BioSequence("ATGC", id=f"HLA-A*{i:05d}", moltype=MolType.DNA)
                for i in range(2000)
            ],
        )
        with pytest.raises(TargetingError) as excinfo:
            AlleleTargeting.subset([f"MISSING-{i}" for i in range(25)], big)
        assert "more" in str(excinfo.value)
