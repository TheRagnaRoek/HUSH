"""Tests for the sequence domain model (core.sequence)."""

import pytest

from core.errors import MolTypeError, SequenceValidationError
from core.sequence import (
    BioSequence,
    DesignInput,
    MolType,
    SequencePanel,
    detect_moltype,
)


class TestNormalisationAndValidation:
    def test_seq_is_uppercased_and_whitespace_stripped(self):
        s = BioSequence("  aug c\nauc ", id="x")
        assert s.seq == "AUGCAUC"

    def test_empty_sequence_is_unknown(self):
        assert BioSequence("").moltype is MolType.UNKNOWN

    def test_invalid_residue_for_explicit_moltype_raises(self):
        with pytest.raises(SequenceValidationError):
            BioSequence("ACGTXQ", moltype=MolType.DNA)

    def test_gaps_are_allowed_for_aligned_sequences(self):
        s = BioSequence("ACG-T.A", moltype=MolType.DNA)
        assert s.seq == "ACG-T.A"


class TestMolTypeDetection:
    @pytest.mark.parametrize(
        "seq,expected",
        [
            ("ATGCGTACGT", MolType.DNA),
            ("AUGCGUACGU", MolType.RNA),
            ("MKLVWQAGYF", MolType.PROTEIN),
            ("", MolType.UNKNOWN),
        ],
    )
    def test_detect_moltype(self, seq, expected):
        assert detect_moltype(seq) == expected

    def test_mixed_t_and_u_is_not_nucleotide(self):
        # A sequence containing both T and U is ambiguous -> protein bucket.
        assert detect_moltype("ATUGATUG") is MolType.PROTEIN

    def test_explicit_moltype_overrides_detection(self):
        s = BioSequence("ACGACG", moltype=MolType.PROTEIN)
        assert s.moltype is MolType.PROTEIN


class TestNucleotideOperations:
    def test_transcribe_dna_to_rna(self):
        assert BioSequence("ATGC", moltype=MolType.DNA).transcribe().seq == "AUGC"

    def test_transcribe_rna_is_identity(self):
        rna = BioSequence("AUGC", moltype=MolType.RNA)
        assert rna.transcribe() is rna

    def test_reverse_transcribe_rna_to_dna(self):
        assert BioSequence("AUGC").reverse_transcribe().seq == "ATGC"

    def test_dna_reverse_complement(self):
        assert BioSequence("ATGCGT").reverse_complement().seq == "ACGCAT"

    def test_rna_reverse_complement(self):
        rc = BioSequence("AUGC", moltype=MolType.RNA).reverse_complement()
        assert rc.seq == "GCAU"
        assert rc.moltype is MolType.RNA

    def test_complement_preserves_ambiguity_codes(self):
        assert BioSequence("RYSWKM", moltype=MolType.DNA).complement().seq == "YRSWMK"

    def test_gc_content(self):
        assert BioSequence("GGCCAT", moltype=MolType.DNA).gc_content == pytest.approx(4 / 6)

    def test_nucleotide_op_on_protein_raises(self):
        prot = BioSequence("MKLVWQ", moltype=MolType.PROTEIN)
        with pytest.raises(MolTypeError):
            prot.transcribe()
        with pytest.raises(MolTypeError):
            _ = prot.gc_content


class TestSubseq:
    def test_subseq_slices_and_annotates_id(self):
        s = BioSequence("AUGCAUGC", id="t", moltype=MolType.RNA)
        sub = s.subseq(2, 5)
        assert sub.seq == "GCA"
        assert sub.id == "t:3-5"

    @pytest.mark.parametrize("start,end", [(-1, 3), (0, 99), (5, 2)])
    def test_invalid_slice_raises(self, start, end):
        with pytest.raises(ValueError):
            BioSequence("AUGCAUGC").subseq(start, end)


class TestSequencePanel:
    def _panel(self):
        return SequencePanel(
            "fam",
            [
                BioSequence("AUGC", id="a", moltype=MolType.RNA),
                BioSequence("AUGG", id="b", moltype=MolType.RNA),
            ],
        )

    def test_len_iter_ids(self):
        p = self._panel()
        assert len(p) == 2
        assert p.ids == ["a", "b"]
        assert [s.id for s in p] == ["a", "b"]

    def test_get_by_id(self):
        p = self._panel()
        assert p.get("b").seq == "AUGG"
        assert p.get("missing") is None

    def test_is_uniform(self):
        assert self._panel().is_uniform
        mixed = SequencePanel(
            "m",
            [
                BioSequence("AUGC", moltype=MolType.RNA),
                BioSequence("MKLV", moltype=MolType.PROTEIN),
            ],
        )
        assert not mixed.is_uniform


class TestDesignInput:
    def test_has_panel(self):
        target = BioSequence("AUGC", moltype=MolType.RNA)
        assert not DesignInput(target).has_panel
        panel = SequencePanel("p", [BioSequence("AUGC", moltype=MolType.RNA)])
        assert DesignInput(target, panel).has_panel

    def test_require_nucleotide_target_ok(self):
        target = BioSequence("AUGC", moltype=MolType.RNA)
        assert DesignInput(target).require_nucleotide_target() is target

    def test_require_nucleotide_target_rejects_protein(self):
        target = BioSequence("MKLVWQ", moltype=MolType.PROTEIN)
        with pytest.raises(MolTypeError):
            DesignInput(target).require_nucleotide_target()
