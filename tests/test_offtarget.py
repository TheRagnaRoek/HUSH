"""Tests for the transcriptome homology off-target filter (core.offtarget)."""

import shutil
import subprocess

import pytest

from core.offtarget import (
    DEFAULT_SELF_MIN_COVERAGE,
    BowtieSearcher,
    BruteForceSearcher,
    OffTargetHit,
    OffTargetProfile,
    PrecomputedSearcher,
    SelfTranscripts,
    _parse_bowtie_line,
    evaluate_offtarget,
    resolve_self_transcripts,
)
from core.selection import SelectionCriteria, evaluate
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates

# The planted fixture (tests/data/offtarget_mini.fasta):
#   NM_000001.1  the self transcript, exact sense site (0 mismatch)
#   NM_000002.1  a guide off-target: the sense site with 2 mismatches
#   NM_000003.1  a passenger off-target: the antisense site with 1 mismatch
#   NM_000004.1  a clean, unrelated transcript
_REFS = {
    "NM_000001.1": "TTAAGGCCTAATGCTTAGCCGATTACAGAGATCGATCGA",
    "NM_000002.1": "TTAAGGCCTAATGCGTAGCCGAGTACAGAGATCGATCGA",
    "NM_000003.1": "TTAAGGCCTATCTGTAATCTGCTAAGCATGATCGATCGA",
    "NM_000004.1": "TTTAAACCCGGGAAATTTCCCGGGAAATTTCCCGGGAAA",
}
_SENSE_CORE = "ATGCTTAGCCGATTACAGA"  # 19 nt; the guide targets this site
_SELF = "NM_000001.1"


def make_candidate(sense: str):
    """A candidate whose sense strand is 'sense' and antisense its revcomp."""
    core = BioSequence(sense, id="core", moltype=MolType.DNA)
    return next(iter(enumerate_candidates(core, length=len(sense))))


@pytest.fixture
def brute() -> BruteForceSearcher:
    return BruteForceSearcher(list(_REFS.items()), max_mismatch=3)


@pytest.fixture
def candidate():
    return make_candidate(_SENSE_CORE)


class TestBruteForceSearch:
    def test_finds_exact_and_mismatch_sites(self, brute):
        hits = brute.search(_SENSE_CORE)
        by_acc = {h.accession: h.mismatches for h in hits}
        assert by_acc["NM_000001.1"] == 0  # self, exact
        assert by_acc["NM_000002.1"] == 2  # planted 2-mismatch off-target
        assert "NM_000004.1" not in by_acc  # unrelated, no hit within 3

    def test_hits_carry_position(self, brute):
        hit = next(h for h in brute.search(_SENSE_CORE) if h.accession == _SELF)
        assert hit.position == 10  # after the 10-nt left flank

    def test_rna_read_is_matched_as_dna(self, brute):
        # U-spelled reads must match the DNA reference.
        assert brute.search("AUGCUUAGCCGAUUACAGA") == brute.search(_SENSE_CORE)

    def test_depth_limits_reported_hits(self):
        shallow = BruteForceSearcher(list(_REFS.items()), max_mismatch=1)
        # The 2-mismatch off-target is invisible at depth 1.
        accs = {h.accession for h in shallow.search(_SENSE_CORE)}
        assert accs == {"NM_000001.1"}

    def test_from_fasta_matches_in_memory(self, offtarget_mini_fasta, brute):
        from_file = BruteForceSearcher.from_fasta(offtarget_mini_fasta, max_mismatch=3)
        assert from_file.search(_SENSE_CORE) == brute.search(_SENSE_CORE)


class TestEvaluateOffTarget:
    def test_both_strands_are_searched(self, brute, candidate):
        prof = evaluate_offtarget(candidate, brute, threshold=3, self_transcripts=_SELF)
        assert prof.guide_min_mismatch == 2       # from the sense-site off-target
        assert prof.passenger_min_mismatch == 1   # from the antisense-site off-target
        assert prof.min_mismatch == 1

    def test_on_target_is_excluded(self, brute, candidate):
        excluded = evaluate_offtarget(candidate, brute, self_transcripts=_SELF)
        assert excluded.guide_min_mismatch == 2   # the 0-mismatch self hit is gone
        kept = evaluate_offtarget(candidate, brute, self_transcripts=None)
        assert kept.guide_min_mismatch == 0       # without exclusion, self dominates

    def test_exclusion_ignores_accession_version(self, brute, candidate):
        # A version-less self accession still matches the versioned transcript.
        prof = evaluate_offtarget(candidate, brute, self_transcripts="NM_000001")
        assert prof.guide_min_mismatch == 2

    def test_counts_by_mismatch_level(self, brute, candidate):
        prof = evaluate_offtarget(candidate, brute, self_transcripts=_SELF)
        assert prof.guide_counts == (0, 0, 1, 0)      # one hit at 2 mismatches
        assert prof.passenger_counts == (0, 1, 0, 0)  # one hit at 1 mismatch

    def test_clean_candidate_has_no_offtarget(self, brute):
        # A sense core absent from the reference yields no hits -> passes.
        clean = make_candidate("GAGAGAGAGAGAGAGAGAG")
        prof = evaluate_offtarget(clean, brute, threshold=3)
        assert prof.min_mismatch is None
        assert prof.passes is True


class TestSelfTranscripts:
    """The owned-transcript set: a gene family is never one accession."""

    def test_of_normalises_versions_and_case(self):
        own = SelfTranscripts.of(declared=["NM_002116.8", "nr_001434"])
        assert "NM_002116" in own
        assert "NM_002116.3" in own  # any version of the same transcript
        assert "NR_001434.4" in own
        assert "NM_002117" not in own

    def test_declared_and_resolved_both_count(self):
        own = SelfTranscripts.of(declared=["NM_1"], resolved=["XM_2", "XM_3"])
        assert len(own) == 3
        assert own.accessions == {"nm_1", "xm_2", "xm_3"}

    def test_empty_set_is_falsey(self):
        assert not SelfTranscripts.of()
        assert SelfTranscripts.of(declared=["NM_1"])

    def test_blank_accessions_are_dropped(self):
        assert len(SelfTranscripts.of(declared=["NM_1", "", "  "])) == 1


class TestSetExclusion:
    """Exclusion takes a set, which is what HLA needs (see docs/design-rules.md)."""

    def test_several_transcripts_are_excluded_at_once(self, brute, candidate):
        # Both the self transcript and the planted guide off-target are owned.
        own = SelfTranscripts.of(declared=["NM_000001", "NM_000002"])
        prof = evaluate_offtarget(candidate, brute, self_transcripts=own)
        assert prof.guide_min_mismatch is None  # nothing left on the guide side
        assert prof.passenger_min_mismatch == 1  # the other strand is untouched

    def test_one_accession_is_not_enough(self, brute, candidate):
        # The single-accession behaviour HUSH used to have: a sibling transcript
        # of the same family still counts against the candidate.
        prof = evaluate_offtarget(candidate, brute, self_transcripts="NM_000001")
        assert prof.guide_min_mismatch == 2

    def test_iterable_of_accessions_is_accepted(self, brute, candidate):
        prof = evaluate_offtarget(
            candidate, brute, self_transcripts=["NM_000001", "NM_000002"]
        )
        assert prof.guide_min_mismatch is None


class TestResolveSelfTranscripts:
    """Resolving the owned set by exact k-mer coverage of the design's sequences."""

    # A 60-nt target, and transcripts carrying progressively less of it.
    TARGET = (
        "ATGCTTAGCCGATTACAGAGTCTAGGCATTACCGTAAGGCTTCAGAGGTCATTACGGCAT"
    )
    PANEL = "GGTTCAGTACCGATTTCAGGGTCAT"  # 25 nt, shares no 19-mer with TARGET
    REFS = [
        ("NM_100001.1", TARGET),                       # whole target: coverage 1.0
        ("NM_100002.1", TARGET[:30] + "TTTTTTTTTT"),   # 12 of 42 tiles: 0.29
        ("NM_100003.1", "CCCC" + TARGET[:19] + "CCCC"),  # 1 of 42 tiles: 0.02
        ("NM_100004.1", "TTTAAACCCGGGAAATTTCCCGGGAAATTTCCCGGG"),  # nothing
        ("NM_100005.1", PANEL),                        # the panel member in full
    ]

    @pytest.fixture
    def refs(self) -> BruteForceSearcher:
        return BruteForceSearcher(self.REFS, max_mismatch=0)

    def test_high_coverage_transcript_is_resolved(self, refs):
        own = resolve_self_transcripts(refs, [self.TARGET], length=19, min_coverage=0.5)
        assert own == {"nm_100001"}

    def test_threshold_widens_the_family(self, refs):
        own = resolve_self_transcripts(refs, [self.TARGET], length=19, min_coverage=0.25)
        assert own == {"nm_100001", "nm_100002"}

    def test_incidental_single_kmer_is_not_self(self, refs):
        # One shared 19-mer is a genuine off-target, not family membership.
        own = resolve_self_transcripts(
            refs, [self.TARGET], length=19, min_coverage=DEFAULT_SELF_MIN_COVERAGE
        )
        assert "nm_100003" not in own
        assert "nm_100004" not in own

    def test_coverage_is_judged_per_input_not_pooled(self, refs):
        # The panel member is 25 nt against a 60-nt target; pooling the tiles
        # would dilute its own transcript below any useful threshold.
        own = resolve_self_transcripts(
            refs, [self.TARGET, self.PANEL], length=19, min_coverage=0.5
        )
        assert own == {"nm_100001", "nm_100005"}

    def test_batched_backend_is_used_when_available(self, refs):
        calls = []

        class Batched(BruteForceSearcher):
            def precompute(self, reads):
                calls.append(sorted(reads))
                return self

        batched = Batched(self.REFS, max_mismatch=0)
        own = resolve_self_transcripts(batched, [self.TARGET], length=19)
        assert len(calls) == 1  # one pass over the pooled tiles, not one per tile
        assert "nm_100001" in own

    def test_no_sequences_resolves_to_nothing(self, refs):
        assert resolve_self_transcripts(refs, [], length=19) == frozenset()
        # A sequence shorter than one tile has no k-mers to judge.
        assert resolve_self_transcripts(refs, ["ACGT"], length=19) == frozenset()

    @pytest.mark.parametrize("bad", [0, -0.1, 1.5])
    def test_coverage_outside_the_unit_range_rejected(self, refs, bad):
        with pytest.raises(ValueError, match="min_coverage"):
            resolve_self_transcripts(refs, [self.TARGET], length=19, min_coverage=bad)

    def test_length_must_be_positive(self, refs):
        with pytest.raises(ValueError, match="length"):
            resolve_self_transcripts(refs, [self.TARGET], length=0)


class TestThresholdBoundary:
    def test_passes_when_threshold_at_or_below_nearest(self, brute, candidate):
        # Nearest off-target is at 1 mismatch.
        assert evaluate_offtarget(candidate, brute, 1, self_transcripts=_SELF).passes
        assert not evaluate_offtarget(candidate, brute, 2, self_transcripts=_SELF).passes

    def test_per_strand_pass_flags(self, brute, candidate):
        prof = evaluate_offtarget(candidate, brute, threshold=2, self_transcripts=_SELF)
        assert prof.guide_passes is True       # guide nearest is 2 (>= 2)
        assert prof.passenger_passes is False  # passenger nearest is 1 (< 2)
        assert prof.passes is False            # the stricter strand wins


class TestProfileModel:
    def test_min_mismatch_handles_missing_strand(self):
        prof = OffTargetProfile(
            guide_min_mismatch=None, guide_counts=(0, 0, 0, 0),
            passenger_min_mismatch=2, passenger_counts=(0, 0, 1, 0),
            max_mismatch=3, threshold=3,
        )
        assert prof.min_mismatch == 2
        assert prof.passes is False

    def test_absent_offtarget_passes(self):
        prof = OffTargetProfile(
            guide_min_mismatch=None, guide_counts=(0,),
            passenger_min_mismatch=None, passenger_counts=(0,),
            max_mismatch=0, threshold=1,
        )
        assert prof.min_mismatch is None
        assert prof.passes is True


class TestSelectionIntegration:
    def test_rejected_with_reason(self, brute, candidate):
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False,
            use_offtarget_filter=True, offtarget_min_mismatch=3,
            offtarget_searcher=brute,
            offtarget_self_transcripts=SelfTranscripts.of(declared=_SELF),
        )
        ev = evaluate(candidate, criteria)
        assert ev.selected is False
        assert ev.offtarget_pass is False
        assert any("off-target" in r for r in ev.reasons)

    def test_selected_when_threshold_low(self, brute, candidate):
        criteria = SelectionCriteria(
            use_ui_tei=False, use_reynolds=False, use_amarzguioui=False,
            use_seed_filter=False,
            use_offtarget_filter=True, offtarget_min_mismatch=1,
            offtarget_searcher=brute,
            offtarget_self_transcripts=SelfTranscripts.of(declared=_SELF),
        )
        ev = evaluate(candidate, criteria)
        assert ev.selected is True
        assert ev.offtarget_pass is True


class TestCriteriaValidation:
    def test_filter_without_searcher_rejected(self):
        with pytest.raises(ValueError, match="requires an offtarget_searcher"):
            SelectionCriteria(use_offtarget_filter=True, offtarget_searcher=None)

    def test_threshold_must_be_positive(self, brute):
        with pytest.raises(ValueError, match="must be >= 1"):
            SelectionCriteria(
                use_offtarget_filter=True, offtarget_searcher=brute,
                offtarget_min_mismatch=0,
            )

    def test_search_depth_must_reach_threshold(self):
        shallow = BruteForceSearcher(list(_REFS.items()), max_mismatch=1)
        with pytest.raises(ValueError, match="too shallow"):
            SelectionCriteria(
                use_offtarget_filter=True, offtarget_searcher=shallow,
                offtarget_min_mismatch=3,
            )


class TestBowtieSearcher:
    def test_depth_above_bowtie_cap_rejected(self):
        # bowtie -v allows 0-3; a deeper request is refused before any binary check.
        with pytest.raises(ValueError, match="0-3 mismatches"):
            BowtieSearcher("nonexistent_index", max_mismatch=4)

    @pytest.mark.parametrize(
        "descriptors, expected_mm",
        [("", 0), ("5:A>G", 1), ("8:A>G,12:C>T", 2), ("1:A>C,4:G>T,9:C>A", 3)],
    )
    def test_parse_bowtie_line_counts_mismatches(self, descriptors, expected_mm):
        # bowtie default output: name, strand, ref, 0-based offset, read, quals,
        # reserved, mismatch descriptors.
        line = "\t".join(
            ["ATGCTTAGCCGATTACAGA", "+", "NM_000002.1", "10",
             "ATGCTTAGCCGATTACAGA", "IIIIIIIIIIIIIIIIIII", "0", descriptors]
        )
        name, hit = _parse_bowtie_line(line)
        assert name == "ATGCTTAGCCGATTACAGA"
        assert hit.accession == "NM_000002.1"
        assert hit.position == 10
        assert hit.mismatches == expected_mm

    def test_parse_bowtie_line_without_trailing_field(self):
        # A perfect hit whose trailing (empty) descriptor field is absent.
        line = "\t".join(["r", "+", "NM_1", "0", "AC", "II", "0"])
        _, hit = _parse_bowtie_line(line)
        assert hit.mismatches == 0

    def test_precomputed_lookup(self):
        hits = {"ATGC": [OffTargetHit("NM_1", 0, 0)]}
        pre = PrecomputedSearcher(hits, max_mismatch=3)
        assert pre.search("ATGC") == hits["ATGC"]
        assert pre.search("AUGC") == hits["ATGC"]  # RNA spelling normalised
        assert pre.search("GGGG") == []


# bowtie parity - only when bowtie is installed (CI), skipped otherwise.

_HAS_BOWTIE = shutil.which("bowtie") is not None and shutil.which("bowtie-build") is not None


@pytest.mark.skipif(not _HAS_BOWTIE, reason="bowtie / bowtie-build not on PATH")
class TestBowtieParity:
    @pytest.fixture
    def bowtie_index(self, offtarget_mini_fasta, tmp_path):
        prefix = tmp_path / "mini"
        subprocess.run(
            ["bowtie-build", str(offtarget_mini_fasta), str(prefix)],
            check=True, capture_output=True, text=True,
        )
        return prefix

    def test_matches_brute_force(self, bowtie_index, brute, candidate):
        bowtie = BowtieSearcher(bowtie_index, max_mismatch=3)
        reads = [candidate.sense.seq, candidate.antisense.seq]
        precomputed = bowtie.precompute(reads)

        bowtie_prof = evaluate_offtarget(candidate, precomputed, self_transcripts=_SELF)
        brute_prof = evaluate_offtarget(candidate, brute, self_transcripts=_SELF)
        assert bowtie_prof.guide_min_mismatch == brute_prof.guide_min_mismatch
        assert bowtie_prof.passenger_min_mismatch == brute_prof.passenger_min_mismatch
        assert bowtie_prof.guide_counts == brute_prof.guide_counts
        assert bowtie_prof.passenger_counts == brute_prof.passenger_counts
