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


# --------------------------------------------------------------------------- #
# design command
# --------------------------------------------------------------------------- #


def test_design_default_runs(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta)])
    assert result.exit_code == 0
    assert "Criteria" in result.stdout
    assert "candidate 19-mer" in result.stdout
    assert "selected:" in result.stdout


def test_design_relaxed_selects_and_lists(target_fasta):
    result = runner.invoke(
        app, ["design", str(target_fasta), "--mode", "any", "--no-seed", "-n", "3"]
    )
    assert result.exit_code == 0
    assert "mode: ANY" in result.stdout
    assert "guide 5->3" in result.stdout  # the candidate table rendered


def test_design_single_rule(target_fasta):
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--rules", "reynolds", "--no-seed", "-n", "1"],
    )
    assert result.exit_code == 0
    assert "reynolds  (mode:" in result.stdout


def test_design_unknown_rule_errors(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta), "--rules", "bogus"])
    assert result.exit_code != 0


def test_design_bad_mode_errors(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta), "--mode", "xor"])
    assert result.exit_code != 0


def test_design_inverted_gc_window_errors(target_fasta):
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--gc", "--gc-min", "60", "--gc-max", "30"],
    )
    assert result.exit_code == 1
    assert "Design error" in result.output


def test_design_protein_target_errors(mini_pdb):
    result = runner.invoke(app, ["design", str(mini_pdb)])
    assert result.exit_code == 1
    assert "Design error" in result.output


def test_design_writes_tsv(target_fasta, tmp_path):
    out = tmp_path / "selected.tsv"
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--mode", "any", "--no-seed",
            "-o", str(out),
        ],
    )
    assert result.exit_code == 0
    assert out.exists()
    lines = out.read_text(encoding="ascii").splitlines()
    assert lines[0].startswith("region\tstart\tend\tguide_5to3")
    assert "guide_overhang" in lines[0]
    assert "homopolymer_max_run" in lines[0]
    assert len(lines) > 1  # header plus at least one selected candidate


def test_design_reports_overhang_policy(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta)])
    assert result.exit_code == 0
    assert "natural 3' overhang" in result.stdout


def test_design_dtdt_overhang_renders_lowercase_tt(target_fasta):
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--mode", "any", "--no-seed",
         "--overhang", "dtdt", "-n", "3"],
    )
    assert result.exit_code == 0
    assert "dtdt 3' overhang" in result.stdout
    assert "tt " in result.stdout  # strands end in the lowercase dTdT overhang


def test_design_unknown_overhang_errors(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta), "--overhang", "bogus"])
    assert result.exit_code != 0


def test_design_homopolymer_filter_shows_funnel_row(target_fasta):
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--mode", "any", "--no-seed", "--homopolymer"],
    )
    assert result.exit_code == 0
    assert "homopolymer:" in result.stdout
    assert "reject runs > 3 nt" in result.stdout
