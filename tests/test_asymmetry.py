"""Tests for duplex-end thermodynamic asymmetry (core.asymmetry)."""

import pytest

from core.asymmetry import (
    DEFAULT_WINDOW_BP,
    MIN_WINDOW_BP,
    EndAsymmetry,
    end_asymmetry,
    terminal_stacking_dg,
)
from core.errors import SequenceValidationError
from core.inputs import load_sequences
from core.rules import amarzguioui
from core.sequence import BioSequence, MolType
from core.sirna import SiRNACandidate, enumerate_candidates
from core.thermo import (
    NN_PARAMS_FREIER_1986,
    NN_PARAMS_XIA_1998,
    gibbs_free_energy,
    stacking_terms,
)


def make_candidate(sense: str) -> SiRNACandidate:
    """Build a candidate whose strands are true reverse-complements."""
    core = BioSequence(sense, id="core", moltype=MolType.RNA)
    return next(iter(enumerate_candidates(core, length=len(sense))))


# A duplex that is G/C-rich at the sense 5' end and A/U-rich at the sense 3'
# end. The sense 3' end is the same physical end as the guide's 5' end, so the
# guide's 5' end is the loose one - the functional orientation.
GUIDE_FAVOURED = "GCGCGCGCGCGCGCAUAUA"
# The same duplex with the strand roles swapped: taking the reverse-complement
# as the sense strand describes the same molecule read from the other side, so
# the loose end now sits at the sense 5' end and the passenger would be loaded.
PASSENGER_FAVOURED = (
    BioSequence(GUIDE_FAVOURED, moltype=MolType.RNA).reverse_complement().seq
)
# Equal 4-bp windows at both ends, so neither end is looser than the other.
SYMMETRIC = "GCGC" + "A" * 11 + "GCGC"


class TestTerminalStackingDg:
    def test_uses_stacking_only_no_initiation(self):
        # Must equal the bare stacking sum of the window, with no duplex
        # initiation or terminal A-U penalty folded in.
        strand = "GCGCGCGCGCGCGCGCGCG"
        expected = gibbs_free_energy(*stacking_terms(strand[:4]))
        assert terminal_stacking_dg(strand, window=4) == pytest.approx(expected)

    def test_gc_window_more_stable_than_au_window(self):
        assert terminal_stacking_dg("GCGCGCGCGCGCGCGCGCG") < terminal_stacking_dg(
            "AUAUAUAUAUAUAUAUAUA"
        )

    def test_reads_from_the_five_prime_end(self):
        # Only the first 'window' bases matter.
        assert terminal_stacking_dg("GCGC" + "A" * 15) == terminal_stacking_dg(
            "GCGC" + "U" * 15
        )

    def test_larger_window_is_more_stable(self):
        strand = "GCGCGCGCGCGCGCGCGCG"
        assert terminal_stacking_dg(strand, window=6) < terminal_stacking_dg(
            strand, window=4
        )

    def test_parameter_set_defaults_to_xia_and_freier_differs(self):
        strand = "GCAUGCAUGCAUGCAUGCA"
        default = terminal_stacking_dg(strand, window=4)
        xia = terminal_stacking_dg(strand, window=4, params=NN_PARAMS_XIA_1998)
        freier = terminal_stacking_dg(strand, window=4, params=NN_PARAMS_FREIER_1986)
        assert default == pytest.approx(xia)   # default is Xia
        assert freier != pytest.approx(xia)    # the choice actually bites


class TestSignConvention:
    """The differential's sign carries the whole mechanistic claim."""

    def test_loose_guide_end_gives_positive_differential(self):
        result = end_asymmetry(make_candidate(GUIDE_FAVOURED))
        assert result.differential > 0
        assert result.favours_guide is True
        assert result.preferred_strand == "guide"

    def test_loose_sense_end_gives_negative_differential(self):
        result = end_asymmetry(make_candidate(PASSENGER_FAVOURED))
        assert result.differential < 0
        assert result.favours_guide is False
        assert result.preferred_strand == "passenger"

    def test_differential_is_guide_minus_sense(self):
        result = end_asymmetry(make_candidate(GUIDE_FAVOURED))
        assert result.differential == pytest.approx(
            result.guide_5p_dg - result.sense_5p_dg
        )

    def test_swapping_strand_roles_flips_the_sign(self):
        # Same molecule, read from the other side: the differential must negate
        # exactly. A duplex cannot prefer both of its strands.
        forward = end_asymmetry(make_candidate(GUIDE_FAVOURED))
        swapped = end_asymmetry(make_candidate(PASSENGER_FAVOURED))
        assert forward.differential == pytest.approx(-swapped.differential)

    def test_energetically_symmetric_duplex_prefers_neither(self):
        result = end_asymmetry(make_candidate(SYMMETRIC))
        assert result.differential == pytest.approx(0.0)
        assert result.favours_guide is False
        assert result.preferred_strand == "none"

    def test_reduces_to_the_sense_strand_alone(self):
        # The guide's 5' window is the reverse-complement of the sense strand's
        # 3' window, and nearest-neighbor energies are invariant under
        # reverse-complementation - so the differential is computable from the
        # sense strand's two ends. This is a useful cross-check on the strand
        # bookkeeping in end_asymmetry().
        sense = GUIDE_FAVOURED
        expected = gibbs_free_energy(
            *stacking_terms(sense[-4:])
        ) - gibbs_free_energy(*stacking_terms(sense[:4]))
        assert end_asymmetry(make_candidate(sense)).differential == pytest.approx(
            expected
        )


class TestAgreementWithAmarzguioui:
    """The thermodynamic differential and the A/U proxy (Amarzguioui DT3) should agree in sign;
    a disagreement would mean one is inverted."""

    def test_signs_agree_on_constructed_duplexes(self):
        for sense, expected_positive in (
            (GUIDE_FAVOURED, True),
            (PASSENGER_FAVOURED, False),
        ):
            candidate = make_candidate(sense)
            thermo = end_asymmetry(candidate).differential
            proxy = amarzguioui(candidate).end_au_differential
            assert (thermo > 0) is expected_positive
            assert (proxy > 0) is expected_positive

    def test_correlate_across_a_real_transcript(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        candidates = list(enumerate_candidates(target))

        pairs = [
            (end_asymmetry(c).differential, amarzguioui(c).end_au_differential)
            for c in candidates
        ]
        # Among candidates where the A/U proxy is unambiguous, the thermodynamic
        # differential should agree in sign the large majority of the time.
        decisive = [(t, p) for t, p in pairs if abs(p) >= 2]
        assert len(decisive) > 50
        agree = sum(1 for t, p in decisive if (t > 0) == (p > 0))
        assert agree / len(decisive) > 0.9


class TestWindowValidation:
    def test_default_window_is_four(self):
        assert DEFAULT_WINDOW_BP == 4
        assert end_asymmetry(make_candidate(GUIDE_FAVOURED)).window == 4

    def test_window_below_minimum_rejected(self):
        with pytest.raises(ValueError, match="at least 2 base pairs"):
            end_asymmetry(make_candidate(GUIDE_FAVOURED), window=MIN_WINDOW_BP - 1)

    def test_overlapping_windows_rejected(self):
        # A 19-bp core cannot support two 10-bp end windows.
        with pytest.raises(ValueError, match="overlap"):
            end_asymmetry(make_candidate(GUIDE_FAVOURED), window=10)

    def test_largest_non_overlapping_window_allowed(self):
        result = end_asymmetry(make_candidate(GUIDE_FAVOURED), window=9)
        assert result.window == 9

    def test_window_changes_the_measurement(self):
        candidate = make_candidate("GCGCGCGCGCGCGCAUAUA")
        assert end_asymmetry(candidate, window=2).differential != pytest.approx(
            end_asymmetry(candidate, window=6).differential
        )


class TestTemperature:
    def test_default_is_thirty_seven(self):
        assert end_asymmetry(make_candidate(GUIDE_FAVOURED)).temperature_c == 37.0

    def test_temperature_changes_the_free_energies(self):
        cold = end_asymmetry(make_candidate(GUIDE_FAVOURED), temperature_c=25.0)
        warm = end_asymmetry(make_candidate(GUIDE_FAVOURED), temperature_c=50.0)
        assert cold.guide_5p_dg < warm.guide_5p_dg


class TestParameterSet:
    def test_params_flow_through_and_default_is_xia(self):
        c = make_candidate("GCAUGCAUGCAUGCAUGCA")
        default = end_asymmetry(c).differential
        xia = end_asymmetry(c, params=NN_PARAMS_XIA_1998).differential
        freier = end_asymmetry(c, params=NN_PARAMS_FREIER_1986).differential
        assert default == pytest.approx(xia)
        assert freier != pytest.approx(xia)


class TestInputValidation:
    def test_ambiguous_residue_rejected(self):
        core = BioSequence("GCGCGCGCGCGCGCAUANA", id="c", moltype=MolType.RNA)
        candidate = SiRNACandidate(0, 19, core, core.reverse_complement())
        with pytest.raises(SequenceValidationError):
            end_asymmetry(candidate)


class TestRepr:
    def test_reports_differential_and_preference(self):
        text = repr(end_asymmetry(make_candidate(GUIDE_FAVOURED)))
        assert "preferred=guide" in text
        assert "kcal/mol" in text


class TestRealHlaTarget:
    def test_measures_every_candidate(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        results = [end_asymmetry(c) for c in enumerate_candidates(target)]
        assert len(results) > 1000
        assert all(isinstance(r, EndAsymmetry) for r in results)
        # Both orientations occur in a real transcript.
        assert any(r.favours_guide for r in results)
        assert any(not r.favours_guide for r in results)
