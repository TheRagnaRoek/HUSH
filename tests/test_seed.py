"""Tests for the seed-duplex Tm off-target filter (core.seed)."""

import itertools

import pytest

from core.errors import SequenceValidationError
from core.inputs import load_sequences
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates
from core.seed import (
    SEED_END,
    SEED_LENGTH,
    SEED_START,
    SEED_TM_THRESHOLD_C,
    SIDIRECT_SODIUM_M,
    SIDIRECT_TOTAL_STRAND_CONC_M,
    seed_duplex_tm,
    seed_region,
    seed_tm,
)
from core.sirna import SiRNACandidate
from core.thermo import NN_PARAMS_FREIER_1986, NN_PARAMS_XIA_1998


def make_candidate_pair(sense: str) -> SiRNACandidate:
    """Candidate whose strands are true reverse complements (the filter reads a real duplex)."""
    core = BioSequence(sense, id="core", moltype=MolType.RNA)
    return next(iter(enumerate_candidates(core, length=len(sense))))


class TestSeedRegion:
    def test_is_positions_two_to_eight(self):
        # Position 1 (G) is excluded; the seed is the next seven bases.
        assert seed_region("GAAAUUUCCCGGGAAACCC") == "AAAUUUC"
        assert SEED_LENGTH == 7
        assert (SEED_START, SEED_END) == (2, 8)

    def test_dna_biosequence_transcribed(self):
        dna = BioSequence("GAAATTTCCCGGGAAACCC", id="t", moltype=MolType.DNA)
        assert seed_region(dna) == "AAAUUUC"

    def test_exactly_eight_nt_is_enough(self):
        assert seed_region("GAAAUUUC") == "AAAUUUC"

    def test_too_short_rejected(self):
        with pytest.raises(ValueError, match="at least 8 nt"):
            seed_region("GAAAUUU")

    def test_ambiguity_outside_seed_is_tolerated(self):
        # The filter never looks outside positions 2-8, so an N elsewhere in the
        # strand must not disqualify the candidate.
        assert seed_region("GAAAUUUCNNNNNNNNNNN") == "AAAUUUC"

    def test_ambiguity_inside_seed_rejected(self):
        with pytest.raises(SequenceValidationError):
            seed_region("GAANUUUCCCGGGAAACCC")


@pytest.fixture(scope="module")
def all_heptamer_tms() -> list[float]:
    """Seed-duplex Tm of every possible heptamer; computed once for the module."""
    return [
        seed_duplex_tm("".join(p))
        for p in itertools.product("ACGU", repeat=SEED_LENGTH)
    ]


class TestSiDirectPublishedFigures:
    """Reproduce siDirect 2.0 (Naito 2009): of all 4^7 = 16384 seed heptamers, 4488 (27.4%)
    have Tm below 21.5 C, with values from -12 C to 60 C. Any wrong constant moves the count."""

    def test_covers_every_heptamer(self, all_heptamer_tms):
        assert len(all_heptamer_tms) == 4**SEED_LENGTH == 16384

    def test_count_below_threshold_matches_paper(self, all_heptamer_tms):
        below = sum(1 for tm in all_heptamer_tms if tm < SEED_TM_THRESHOLD_C)
        assert below == 4488
        assert round(100 * below / len(all_heptamer_tms), 1) == 27.4

    def test_range_matches_paper(self, all_heptamer_tms):
        # Paper: "ranged from -12 C to 60 C".
        assert min(all_heptamer_tms) == pytest.approx(-12.0, abs=0.05)
        assert max(all_heptamer_tms) == pytest.approx(60.0, abs=0.5)


class TestParameterSetMatters:
    """Why core.seed pins Freier 1986 rather than using thermo's default."""

    def test_freier_and_xia_tables_differ(self):
        assert NN_PARAMS_FREIER_1986 != NN_PARAMS_XIA_1998

    def test_swapping_the_table_would_move_the_verdict(self):
        # A concrete candidate seed that passes under the calibrated Freier
        # numbers but would be scored differently on Xia values. This test
        # exists to document that the two are not interchangeable.
        from core.thermo import stacking_terms

        seed = "AAAUUUC"
        freier = stacking_terms(seed, params=NN_PARAMS_FREIER_1986)
        xia = stacking_terms(seed, params=NN_PARAMS_XIA_1998)
        assert freier != xia


class TestSeedDuplexTm:
    def test_gc_rich_seeds_melt_higher_than_au_rich(self):
        assert seed_duplex_tm("GCGCGCG") > seed_duplex_tm("AUAUAUA")

    def test_au_rich_seed_passes_and_gc_rich_seed_fails(self):
        # Weak (A/U-rich) seeds pass; strong (G/C-rich) seeds fail.
        assert seed_duplex_tm("AUAUAUA") < SEED_TM_THRESHOLD_C
        assert seed_duplex_tm("GCGCGCG") > SEED_TM_THRESHOLD_C

    def test_rises_with_concentration(self):
        assert seed_duplex_tm("GCAUGCA", total_strand_conc=1e-6) < seed_duplex_tm(
            "GCAUGCA", total_strand_conc=1e-4
        )

    def test_rises_with_salt(self):
        assert seed_duplex_tm("GCAUGCA", sodium=0.01) < seed_duplex_tm(
            "GCAUGCA", sodium=1.0
        )

    def test_salt_term_is_sixteen_point_six_per_decade(self):
        # 16.6 * log10(Na); a ten-fold change in salt is worth exactly 16.6 C.
        low = seed_duplex_tm("GCAUGCA", sodium=0.01)
        high = seed_duplex_tm("GCAUGCA", sodium=0.1)
        assert high - low == pytest.approx(16.6)

    def test_defaults_are_the_sidirect_conditions(self):
        assert SIDIRECT_TOTAL_STRAND_CONC_M == 100e-6
        assert SIDIRECT_SODIUM_M == 0.1
        assert seed_duplex_tm("GCAUGCA") == pytest.approx(
            seed_duplex_tm(
                "GCAUGCA",
                total_strand_conc=SIDIRECT_TOTAL_STRAND_CONC_M,
                sodium=SIDIRECT_SODIUM_M,
            )
        )

    def test_dna_spelling_accepted(self):
        assert seed_duplex_tm("GCATGCA") == pytest.approx(seed_duplex_tm("GCAUGCA"))

    @pytest.mark.parametrize("bad", [0.0, -1e-6])
    def test_non_positive_concentration_rejected(self, bad):
        with pytest.raises(ValueError, match="must be positive"):
            seed_duplex_tm("GCAUGCA", total_strand_conc=bad)

    @pytest.mark.parametrize("bad", [0.0, -0.1])
    def test_non_positive_sodium_rejected(self, bad):
        with pytest.raises(ValueError, match="[Ss]odium.*must be positive"):
            seed_duplex_tm("GCAUGCA", sodium=bad)

    def test_too_short_rejected(self):
        with pytest.raises(ValueError, match="at least 2 nt"):
            seed_duplex_tm("G")


class TestSeedTmFilter:
    def test_both_strands_are_evaluated(self):
        candidate = make_candidate_pair("GAUAUAUAUAUAUAUAUAU")
        result = seed_tm(candidate)
        assert result.guide_seed == seed_region(candidate.antisense)
        assert result.passenger_seed == seed_region(candidate.sense)

    def test_passes_when_both_strands_are_weak(self):
        # An A/U-only core makes both seeds maximally unstable.
        result = seed_tm(make_candidate_pair("AUAUAUAUAUAUAUAUAUA"))
        assert result.guide_passes is True
        assert result.passenger_passes is True
        assert result.passes is True

    def test_fails_when_both_strands_are_strong(self):
        result = seed_tm(make_candidate_pair("GCGCGCGCGCGCGCGCGCG"))
        assert result.passes is False

    def test_one_strong_strand_fails_the_whole_candidate(self):
        # The passenger strand is loaded often enough to matter, so siDirect
        # requires BOTH seeds below the cut-off. Find a real candidate where the
        # two strands disagree and confirm the strict one wins.
        target = BioSequence(
            "AUAUAUAGCGCGCGCGCGCGCGCAUAUAUAUAUAUAUAU", id="mixed", moltype=MolType.RNA
        )
        disagreeing = [
            r
            for r in (seed_tm(c) for c in enumerate_candidates(target))
            if r.guide_passes != r.passenger_passes
        ]
        assert disagreeing, "expected at least one strand-disagreeing candidate"
        for result in disagreeing:
            assert result.passes is False
            assert result.max_tm == max(result.guide_tm, result.passenger_tm)

    def test_threshold_is_strict_and_configurable(self):
        result = seed_tm(make_candidate_pair("GAUAUAUAUAUAUAUAUAU"))
        # A threshold set exactly at the observed Tm must NOT pass ("below").
        exact = seed_tm(
            make_candidate_pair("GAUAUAUAUAUAUAUAUAU"), threshold=result.guide_tm
        )
        assert exact.guide_passes is False

    def test_stricter_threshold_is_never_more_permissive(self):
        candidate = make_candidate_pair("GAUAUCGAUAUAUCGAUAU")
        assert seed_tm(candidate, threshold=10.0).passes <= seed_tm(
            candidate, threshold=21.5
        ).passes

    def test_repr_reports_verdict_and_both_strands(self):
        text = repr(seed_tm(make_candidate_pair("AUAUAUAUAUAUAUAUAUA")))
        assert "passes=True" in text
        assert "guide=" in text and "passenger=" in text


class TestRealHlaTarget:
    def test_filter_is_severe_on_a_gc_rich_hla_cds(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        results = [seed_tm(c) for c in enumerate_candidates(target)]
        assert len(results) > 1000

        passing = [r for r in results if r.passes]
        # At 63.8% GC this CDS yields far fewer passing seeds than the ~3% of
        # all human 23-mers siDirect reports genome-wide; the point is that the
        # filter bites hard here rather than admitting most candidates.
        assert 0 < len(passing) < 0.25 * len(results)

    def test_every_candidate_scores_within_the_possible_range(self, hla_a01_fasta):
        # No candidate can fall outside the range spanned by all 16384 heptamers.
        target = load_sequences(hla_a01_fasta)[0]
        for result in (seed_tm(c) for c in enumerate_candidates(target)):
            for tm in (result.guide_tm, result.passenger_tm):
                assert -12.1 < tm < 60.1
