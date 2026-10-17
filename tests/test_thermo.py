"""Tests for the nearest-neighbor duplex thermodynamics (core.thermo)."""

import math

import pytest

from core.errors import MolTypeError, SequenceValidationError
from core.inputs import load_sequences
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates
from core.thermo import (
    INIT_DH,
    INIT_DS,
    NN_PARAMS,
    REFERENCE_TEMPERATURE_C,
    SYMMETRY_DS,
    TERMINAL_AU_DH,
    TERMINAL_AU_DS,
    ZERO_CELSIUS_K,
    duplex_thermo,
    gibbs_free_energy,
    is_self_complementary,
    melting_temperature,
    stacking_terms,
)

RC = {"A": "U", "U": "A", "G": "C", "C": "G"}


def revcomp(seq: str) -> str:
    return "".join(RC[b] for b in reversed(seq))


class TestParameterTable:
    def test_all_sixteen_steps_present(self):
        steps = {a + b for a in "ACGU" for b in "ACGU"}
        assert set(NN_PARAMS) == steps

    def test_reverse_complement_steps_agree(self):
        # The defining symmetry: a step read from the other strand is its
        # reverse-complement and must carry identical parameters.
        for step, values in NN_PARAMS.items():
            assert NN_PARAMS[revcomp(step)] == values

    def test_spot_values_match_xia_1998(self):
        # Two published entries, checked against the paper's table as typed.
        assert NN_PARAMS["GC"] == (-14.88, -36.9)  # GC/GC, the most stable step
        assert NN_PARAMS["AA"] == (-6.82, -19.0)  # AA/UU, among the weakest

    def test_stacking_is_always_stabilising(self):
        for dh, ds in NN_PARAMS.values():
            assert dh < 0
            assert ds < 0


class TestStackingTerms:
    def test_single_base_has_no_steps(self):
        assert stacking_terms("A") == (0.0, 0.0)

    def test_two_bases_are_one_step(self):
        assert stacking_terms("GC") == NN_PARAMS["GC"]

    def test_steps_are_summed_over_overlapping_windows(self):
        # GCAU is three overlapping steps: GC, CA, AU.
        dh, ds = stacking_terms("GCAU")
        expected_dh = sum(NN_PARAMS[s][0] for s in ("GC", "CA", "AU"))
        expected_ds = sum(NN_PARAMS[s][1] for s in ("GC", "CA", "AU"))
        assert dh == pytest.approx(expected_dh)
        assert ds == pytest.approx(expected_ds)

    def test_excludes_initiation_and_end_terms(self):
        # duplex_thermo adds only initiation and the end corrections to the stack.
        stack_dh, stack_ds = stacking_terms("GCAU")
        duplex = duplex_thermo("GCAU")
        assert duplex.dh - stack_dh == pytest.approx(INIT_DH + TERMINAL_AU_DH)
        assert duplex.ds - stack_ds == pytest.approx(INIT_DS + TERMINAL_AU_DS)


class TestHandComputedDuplexes:
    """Full duplexes worked out by hand from the published parameters."""

    def test_gc_self_complementary_duplex(self):
        # Steps: GC. Both ends G/C, so no terminal penalty. GC is its own
        # reverse-complement, so the symmetry correction applies.
        t = duplex_thermo("GC")
        assert t.dh == pytest.approx(-14.88 + INIT_DH)
        assert t.ds == pytest.approx(-36.9 + INIT_DS + SYMMETRY_DS)
        assert t.terminal_au == 0
        assert t.self_complementary is True

    def test_gcau_duplex_with_one_au_end(self):
        # Steps GC, CA, AU; 5' end is G (no penalty), 3' end is U (one penalty);
        # revcomp(GCAU) = AUGC, so not self-complementary.
        t = duplex_thermo("GCAU")
        assert t.dh == pytest.approx(-34.70 + INIT_DH + TERMINAL_AU_DH)
        assert t.ds == pytest.approx(-90.5 + INIT_DS + TERMINAL_AU_DS)
        assert t.terminal_au == 1
        assert t.self_complementary is False
        assert t.dg37 == pytest.approx(-2.093, abs=1e-3)


class TestTerminalAuPenalty:
    @pytest.mark.parametrize(
        "seq, expected",
        [
            ("GCGCGC", 0),  # both ends G/C
            ("ACGCGC", 1),  # 5' A only
            ("GCGCGU", 1),  # 3' U only
            ("ACGCGU", 2),  # both ends A/U
        ],
    )
    def test_ends_counted(self, seq, expected):
        assert duplex_thermo(seq).terminal_au == expected

    def test_penalty_destabilises(self):
        # Same length and same GC-rich interior, but A/U ends: the A/U-ended
        # duplex must come out less stable (less negative dG).
        gc_ended = duplex_thermo("GCGCGC")
        au_ended = duplex_thermo("ACGCGU")
        assert au_ended.dg37 > gc_ended.dg37


class TestSelfComplementarity:
    @pytest.mark.parametrize("seq", ["GC", "AU", "GCAUAUGC", "AAGCUU"])
    def test_detected(self, seq):
        assert is_self_complementary(seq) is True
        assert duplex_thermo(seq).self_complementary is True

    @pytest.mark.parametrize("seq", ["GCAU", "AAAA", "ACGUACGUAAG"])
    def test_not_detected(self, seq):
        assert is_self_complementary(seq) is False
        assert duplex_thermo(seq).self_complementary is False

    def test_symmetry_correction_applied_once(self):
        palindrome = duplex_thermo("GCAUAUGC")
        stack_ds = stacking_terms("GCAUAUGC")[1]
        assert palindrome.ds - stack_ds == pytest.approx(INIT_DS + SYMMETRY_DS)


class TestFreeEnergy:
    def test_gibbs_matches_definition(self):
        # dS is in cal and dH in kcal, so dS is scaled by 1000.
        assert gibbs_free_energy(-10.0, -25.0, 37.0) == pytest.approx(
            -10.0 - (37.0 + ZERO_CELSIUS_K) * -25.0 / 1000.0
        )

    def test_dg37_is_dg_at_reference_temperature(self):
        t = duplex_thermo("GCAUGCAUGC")
        assert t.dg37 == t.dg(REFERENCE_TEMPERATURE_C)

    def test_stability_falls_as_temperature_rises(self):
        t = duplex_thermo("GCAUGCAUGC")
        assert t.dg(25.0) < t.dg(37.0) < t.dg(50.0)

    def test_gc_rich_more_stable_than_au_rich(self):
        gc_rich = duplex_thermo("GCGCGCGCGC")
        au_rich = duplex_thermo("AUAUAUAUAU")
        assert gc_rich.dg37 < au_rich.dg37


class TestMeltingTemperature:
    def test_rises_with_concentration(self):
        t = duplex_thermo("GCAUGCAUGCAUGCAUGCA")
        assert t.tm(1e-7) < t.tm(1e-6) < t.tm(1e-5)

    def test_gc_rich_melts_higher(self):
        assert duplex_thermo("GCGCGCGCGC").tm() > duplex_thermo("AUAUAUAUAU").tm()

    def test_self_complementary_uses_undivided_concentration(self):
        # x = 1 rather than 4 for a self-complementary strand, so at the same
        # total concentration its effective term is four-fold higher.
        dh, ds = -100.0, -280.0
        assert melting_temperature(
            dh, ds, 1e-6, self_complementary=True
        ) == pytest.approx(melting_temperature(dh, ds, 4e-6, self_complementary=False))

    def test_non_positive_concentration_rejected(self):
        with pytest.raises(ValueError, match="must be positive"):
            melting_temperature(-100.0, -280.0, 0.0)


class TestBiopythonCrossCheck:
    """Check against the Xia 1998 table of Biopython (RNA_NN2), an independent implementation.
    Its gas constant is 1.987 against 1.9872 here, hence the tolerance."""

    TOLERANCE_C = 0.01

    @pytest.mark.parametrize(
        "seq",
        [
            "GCGCGCGCGC",
            "AUAUAUAUAU",
            "CUGAUCGAUCGUAGCUA",
            "AAGCUU",
            "ACGUACGUAAG",
            "GGGCCCAAAUUU",
            "UUUUUUUUUU",
        ],
    )
    def test_tm_matches_biopython(self, seq):
        mt = pytest.importorskip("Bio.SeqUtils.MeltingTemp")

        total_conc = 1e-6
        t = duplex_thermo(seq)
        # Biopython's tables are spelled in DNA letters, takes concentrations in
        # nM, and needs self-complementarity declared rather than detected. For
        # two distinct strands it uses dnac1 - dnac2/2, i.e. half of each equals
        # a quarter of the total; for a self-complementary strand it takes
        # dnac1 as the total.
        if t.self_complementary:
            expected = mt.Tm_NN(
                seq.replace("U", "T"),
                nn_table=mt.RNA_NN2,
                saltcorr=0,
                selfcomp=True,
                dnac1=total_conc * 1e9,
                dnac2=total_conc * 1e9,
            )
        else:
            expected = mt.Tm_NN(
                seq.replace("U", "T"),
                nn_table=mt.RNA_NN2,
                saltcorr=0,
                dnac1=total_conc / 2 * 1e9,
                dnac2=total_conc / 2 * 1e9,
            )
        assert t.tm(total_conc) == pytest.approx(expected, abs=self.TOLERANCE_C)


class TestInputHandling:
    def test_dna_string_read_as_rna(self):
        assert duplex_thermo("GCAT").dh == duplex_thermo("GCAU").dh

    def test_biosequence_dna_target_is_transcribed(self):
        dna = BioSequence("GCATGCAT", id="t", moltype=MolType.DNA)
        assert duplex_thermo(dna).sequence == "GCAUGCAU"

    def test_whitespace_and_case_normalised(self):
        assert duplex_thermo(" gc au \n").sequence == "GCAU"

    def test_ambiguous_bases_rejected(self):
        with pytest.raises(SequenceValidationError, match="unambiguous RNA"):
            duplex_thermo("GCNAU")

    def test_protein_biosequence_rejected(self):
        protein = BioSequence("MKWVTFISLL", id="p", moltype=MolType.PROTEIN)
        with pytest.raises(MolTypeError):
            duplex_thermo(protein)

    @pytest.mark.parametrize("seq", ["", "G"])
    def test_too_short_for_a_duplex(self, seq):
        with pytest.raises(ValueError, match="at least 2 nt"):
            duplex_thermo(seq)


class TestRealHlaTarget:
    """End-to-end over a real HLA-A*01:01 coding sequence."""

    def test_every_candidate_scores_finite_and_stable(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        candidates = list(enumerate_candidates(target))
        assert len(candidates) > 1000  # ~1.1 kb of CDS, 19-mer windows

        energies = [duplex_thermo(c.sense).dg37 for c in candidates]
        assert all(math.isfinite(dg) for dg in energies)
        # A 19-mer perfect duplex is strongly stable under any base composition.
        assert max(energies) < -10.0

    def test_gc_content_tracks_stability(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        candidates = list(enumerate_candidates(target))

        most_stable = min(candidates, key=lambda c: duplex_thermo(c.sense).dg37)
        least_stable = max(candidates, key=lambda c: duplex_thermo(c.sense).dg37)
        assert most_stable.gc_content > least_stable.gc_content

    def test_both_strands_of_a_candidate_agree(self, hla_a01_fasta):
        # A perfect duplex is the same molecule read from either side, so the
        # sense and antisense strands must give identical thermodynamics.
        target = load_sequences(hla_a01_fasta)[0]
        for candidate in list(enumerate_candidates(target))[:50]:
            sense = duplex_thermo(candidate.sense)
            antisense = duplex_thermo(candidate.antisense)
            assert sense.dh == pytest.approx(antisense.dh)
            assert sense.ds == pytest.approx(antisense.ds)
