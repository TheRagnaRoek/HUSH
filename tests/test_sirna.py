"""Tests for the siRNA candidate model and enumerator (core.sirna)."""

import pytest

from core.errors import MolTypeError
from core.sequence import BioSequence, MolType
from core.sirna import (
    DEFAULT_LENGTH,
    SiRNACandidate,
    enumerate_candidates,
)


class TestEnumeration:
    def test_candidate_count_is_length_minus_window_plus_one(self):
        target = BioSequence("AUGCAUGCAUGCAUGCAUGCA", id="t", moltype=MolType.RNA)
        # 21-nt target, 19-nt window -> 21 - 19 + 1 = 3 candidates.
        candidates = list(enumerate_candidates(target, length=19))
        assert len(candidates) == 3

    def test_positions_are_contiguous_from_zero(self):
        target = BioSequence("A" * 25, id="t", moltype=MolType.RNA)
        positions = [c.position for c in enumerate_candidates(target, length=19)]
        assert positions == [0, 1, 2, 3, 4, 5, 6]

    def test_default_length_is_used_when_unspecified(self):
        target = BioSequence("ACGU" * 6, id="t", moltype=MolType.RNA)  # 24 nt
        first = next(enumerate_candidates(target))
        assert first.length == DEFAULT_LENGTH

    def test_dna_target_is_transcribed_to_rna(self):
        target = BioSequence("ATGC" * 5, id="t", moltype=MolType.DNA)  # 20 nt
        first = next(enumerate_candidates(target, length=19))
        assert first.sense.moltype is MolType.RNA
        assert "T" not in first.sense.seq

    def test_input_validated_eagerly_before_iteration(self):
        # Bad length must raise on the call, not lazily during iteration.
        target = BioSequence("A" * 25, id="t", moltype=MolType.RNA)
        with pytest.raises(ValueError):
            enumerate_candidates(target, length=30)


class TestStrandCorrectness:
    def test_sense_matches_target_site(self):
        target = BioSequence("AUGGCUACGUUAGCCAUUGCA", id="t", moltype=MolType.RNA)
        first = next(enumerate_candidates(target, length=19))
        assert first.sense.seq == target.seq[:19]

    def test_antisense_is_reverse_complement_of_sense(self):
        target = BioSequence("AUGGCUACGUUAGCCAUUGCA", id="t", moltype=MolType.RNA)
        first = next(enumerate_candidates(target, length=19))
        assert first.antisense.seq == first.sense.reverse_complement().seq

    def test_strands_are_rna(self):
        target = BioSequence("AUGGCUACGUUAGCCAUUGCA", id="t", moltype=MolType.RNA)
        first = next(enumerate_candidates(target, length=19))
        assert first.sense.moltype is MolType.RNA
        assert first.antisense.moltype is MolType.RNA


class TestCandidateProperties:
    def _first(self):
        target = BioSequence("GCGCGCGCGCAUAUAUAUAU", id="t", moltype=MolType.RNA)
        return next(enumerate_candidates(target, length=19))

    def test_end_and_region_labels(self):
        c = self._first()
        assert c.end == 19
        assert c.region == "1-19"

    def test_guide_and_passenger_aliases(self):
        c = self._first()
        assert c.guide is c.antisense
        assert c.passenger is c.sense

    def test_gc_content_matches_sense(self):
        c = self._first()
        assert c.gc_content == c.sense.gc_content


class TestValidation:
    def test_length_below_range_raises(self):
        target = BioSequence("A" * 25, id="t", moltype=MolType.RNA)
        with pytest.raises(ValueError):
            enumerate_candidates(target, length=18)

    def test_length_above_range_raises(self):
        target = BioSequence("A" * 30, id="t", moltype=MolType.RNA)
        with pytest.raises(ValueError):
            enumerate_candidates(target, length=25)

    def test_protein_target_raises_moltype_error(self):
        target = BioSequence("MKAstop".replace("stop", ""), id="p",
                             moltype=MolType.PROTEIN)
        with pytest.raises(MolTypeError):
            enumerate_candidates(target, length=19)

    def test_target_shorter_than_window_raises(self):
        target = BioSequence("ACGUACGU", id="t", moltype=MolType.RNA)  # 8 nt
        with pytest.raises(ValueError):
            enumerate_candidates(target, length=19)
