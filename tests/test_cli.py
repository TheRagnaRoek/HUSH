"""Smoke tests for the CLI (app.cli)."""

from typer.testing import CliRunner

from app.cli import app

runner = CliRunner()


def test_input_target_only(target_fasta):
    result = runner.invoke(app, ["input", str(target_fasta)])
    assert result.exit_code == 0
    assert "HUSH_demo_target" in result.stdout
    assert "ready for siRNA design" in result.stdout


def test_input_with_panel(target_fasta, panel_fasta):
    result = runner.invoke(
        app, ["input", str(target_fasta), "--panel", str(panel_fasta)]
    )
    assert result.exit_code == 0
    assert "Panel" in result.stdout
    assert "3 sequences" in result.stdout


def test_input_protein_pdb_warns(mini_pdb):
    result = runner.invoke(app, ["input", str(mini_pdb)])
    assert result.exit_code == 0
    # Protein target should trigger the nucleotide-required warning, not crash.
    assert "needs the coding transcript" in result.stdout


def test_inspect_lists_sequences(panel_fasta):
    result = runner.invoke(app, ["inspect", str(panel_fasta)])
    assert result.exit_code == 0
    assert "3 sequence(s)" in result.stdout


def test_input_bad_file_exits_nonzero(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("not a sequence file")
    result = runner.invoke(app, ["input", str(bad)])
    assert result.exit_code != 0
