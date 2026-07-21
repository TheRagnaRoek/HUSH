"""Tests for the Input Layer (core.inputs)."""

import pytest

from core.errors import InputParseError
from core.inputs import (
    build_design_input,
    detect_format,
    load_panel,
    load_sequences,
    load_target,
    read_fasta,
    read_pdb,
)
from core.sequence import MolType


class TestFormatDetection:
    @pytest.mark.parametrize(
        "name,fmt",
        [
            ("x.fasta", "fasta"),
            ("x.fa", "fasta"),
            ("x.fna", "fasta"),
            ("x.pdb", "pdb"),
            ("x.ent", "pdb"),
            ("x.cif", "cif"),
        ],
    )
    def test_detect_format(self, tmp_path, name, fmt):
        p = tmp_path / name
        p.write_text(">a\nACGT\n")
        assert detect_format(p) == fmt

    def test_unknown_extension_raises(self, tmp_path):
        p = tmp_path / "x.txt"
        p.write_text("ACGT")
        with pytest.raises(InputParseError):
            detect_format(p)


class TestReadFasta:
    def test_reads_target_example(self, target_fasta):
        seqs = read_fasta(target_fasta)
        assert len(seqs) == 1
        assert seqs[0].id == "HUSH_demo_target"
        assert seqs[0].moltype is MolType.RNA
        assert "synthetic example" in seqs[0].description

    def test_reads_multi_record_panel(self, panel_fasta):
        seqs = read_fasta(panel_fasta)
        assert len(seqs) == 3
        assert seqs[0].id == "HLA-DEMO*01:01"
        assert all(s.moltype is MolType.RNA for s in seqs)

    def test_metadata_records_source(self, target_fasta):
        seq = read_fasta(target_fasta)[0]
        assert seq.metadata["source_format"] == "fasta"
        assert seq.metadata["source"].endswith("target_example.fasta")

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(InputParseError):
            read_fasta(tmp_path / "nope.fasta")

    def test_empty_file_raises(self, tmp_path):
        p = tmp_path / "empty.fasta"
        p.write_text("")
        with pytest.raises(InputParseError):
            read_fasta(p)


class TestReadPdb:
    def test_reads_protein_seqres(self, mini_pdb):
        seqs = read_pdb(mini_pdb)
        assert len(seqs) == 1
        assert seqs[0].moltype is MolType.PROTEIN
        assert seqs[0].seq == "MKLVWQAG"
        assert seqs[0].metadata.get("chain") == "A"

    def test_missing_structure_raises(self, tmp_path):
        with pytest.raises(InputParseError):
            read_pdb(tmp_path / "nope.pdb")


class TestLoaders:
    def test_load_target_single(self, target_fasta):
        target = load_target(target_fasta)
        assert target.id == "HUSH_demo_target"
        assert "selected_from" not in target.metadata

    def test_load_target_from_multi_records_first_and_flags(self, panel_fasta):
        target = load_target(panel_fasta)
        assert target.id == "HLA-DEMO*01:01"
        assert target.metadata["selected_from"] == 3

    def test_load_panel(self, panel_fasta):
        panel = load_panel(panel_fasta, name="HLA-DEMO")
        assert panel.name == "HLA-DEMO"
        assert len(panel) == 3
        assert panel.is_uniform

    def test_load_panel_default_name_is_stem(self, panel_fasta):
        panel = load_panel(panel_fasta)
        assert panel.name == "hla_panel_example"

    def test_load_sequences_dispatch_pdb(self, mini_pdb):
        seqs = load_sequences(mini_pdb)
        assert seqs[0].moltype is MolType.PROTEIN


class TestBuildDesignInput:
    def test_target_only(self, target_fasta):
        design = build_design_input(target_fasta)
        assert design.target.id == "HUSH_demo_target"
        assert not design.has_panel
        # target is RNA -> downstream requirement satisfied
        assert design.require_nucleotide_target() is design.target

    def test_target_with_panel(self, target_fasta, panel_fasta):
        design = build_design_input(target_fasta, panel_fasta, panel_name="HLA-DEMO")
        assert design.has_panel
        assert design.panel.name == "HLA-DEMO"
        assert len(design.panel) == 3
