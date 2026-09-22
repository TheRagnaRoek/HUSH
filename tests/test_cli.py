"""Smoke tests for the CLI (app.cli)."""

from typer.testing import CliRunner

from app.cli import app

runner = CliRunner()


def test_inspect_with_panel_target_readiness(target_fasta, panel_fasta):
    # --panel switches inspect into design-input preview mode.
    result = runner.invoke(
        app, ["inspect", str(target_fasta), "--panel", str(panel_fasta)]
    )
    assert result.exit_code == 0
    assert "HUSH_demo_target" in result.stdout
    assert "ready for siRNA design" in result.stdout
    assert "Panel" in result.stdout
    assert "3 sequences" in result.stdout


def test_inspect_panel_protein_pdb_warns(mini_pdb, panel_fasta):
    result = runner.invoke(
        app, ["inspect", str(mini_pdb), "--panel", str(panel_fasta)]
    )
    assert result.exit_code == 0
    # Protein target should trigger the nucleotide-required warning, not crash.
    assert "needs the coding transcript" in result.stdout


def test_inspect_lists_sequences(panel_fasta):
    result = runner.invoke(app, ["inspect", str(panel_fasta)])
    assert result.exit_code == 0
    assert "3 sequence(s)" in result.stdout


def test_inspect_listing_shows_readiness(target_fasta):
    result = runner.invoke(app, ["inspect", str(target_fasta)])
    assert result.exit_code == 0
    assert "ready for siRNA design" in result.stdout


def test_inspect_listing_protein_flags_transcript_needed(mini_pdb):
    result = runner.invoke(app, ["inspect", str(mini_pdb)])
    assert result.exit_code == 0
    assert "coding transcript" in result.stdout


def test_inspect_bad_file_exits_nonzero(tmp_path):
    bad = tmp_path / "bad.txt"
    bad.write_text("not a sequence file")
    result = runner.invoke(app, ["inspect", str(bad)])
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


def test_design_lowgc_demo_selects_under_bare_default(demo_lowgc_fasta):
    # The shipped low-GC demo must yield candidates under the full default
    # paradigm (seed on), so a newcomer's first bare run is non-empty.
    result = runner.invoke(app, ["design", str(demo_lowgc_fasta)])
    assert result.exit_code == 0
    assert "seed Tm:" in result.stdout  # seed filter is on by default
    # The selected count sits on its own funnel line; it must be non-zero.
    selected_line = next(
        line for line in result.stdout.splitlines() if "selected:" in line
    )
    assert int(selected_line.split("selected:")[1].strip()) > 0


def test_design_default_mode_is_union(target_fasta):
    # The best-practice default combines the rule sets by ANY (union).
    result = runner.invoke(app, ["design", str(target_fasta)])
    assert result.exit_code == 0
    assert "mode: ANY" in result.stdout


def test_design_offtarget_on_by_default_warns_without_index(target_fasta):
    # Off-target scanning defaults on but degrades to a warning with no index.
    result = runner.invoke(app, ["design", str(target_fasta)])
    assert result.exit_code == 0
    assert "off-target scanning is on by default" in result.stdout


def test_design_no_offtarget_silences_the_reminder(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta), "--no-offtarget"])
    assert result.exit_code == 0
    assert "off-target scanning is on by default" not in result.stdout


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


def test_design_homopolymer_filter_on_by_default(target_fasta):
    # The homopolymer gate is part of the default paradigm.
    result = runner.invoke(app, ["design", str(target_fasta), "--no-seed"])
    assert result.exit_code == 0
    assert "homopolymer:" in result.stdout
    assert "reject runs > 3 nt" in result.stdout


def test_design_allow_homopolymers_disables_the_gate(target_fasta):
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--no-seed", "--allow-homopolymers"],
    )
    assert result.exit_code == 0
    # The criteria block reports the gate as off, and no funnel row for it.
    assert "homopolymer:off" in result.stdout.replace(" ", "")


# --------------------------------------------------------------------------- #
# design command - allele-aware panel targeting (MSA stage)
# --------------------------------------------------------------------------- #

_ALLELE_1 = "HLA-DEMO*01:01"  # a member id in hla_panel_example.fasta


def test_design_panel_pan_allele_gates_and_reports(target_fasta, panel_fasta):
    # A bare --panel defaults to pan-allele targeting.
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--mode", "any", "--no-seed",
        ],
    )
    assert result.exit_code == 0
    assert "pan-allele" in result.stdout
    # The allele gate appears in both the criteria block and the funnel.
    assert "allele:" in result.stdout


def test_design_allele_specific_shows_spare_margin(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", _ALLELE_1, "--mode", "any", "--no-seed",
        ],
    )
    assert result.exit_code == 0
    assert "allele-specific" in result.stdout
    # Discriminating targeting adds the two per-candidate allele columns.
    assert "minDiff" in result.stdout
    assert "disc" in result.stdout


def test_design_explicit_hit_and_spare(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--hit", _ALLELE_1, "--spare", "HLA-DEMO*02:01",
            "--mode", "any", "--no-seed",
        ],
    )
    assert result.exit_code == 0
    assert "minDiff" in result.stdout


def _design_sorted(target_fasta, panel_fasta, sort: str):
    return runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", _ALLELE_1, "--mode", "any", "--no-seed",
            "--sort", sort,
        ],
    )


def test_design_sort_discrimination(target_fasta, panel_fasta):
    result = _design_sorted(target_fasta, panel_fasta, "discrimination")
    assert result.exit_code == 0
    # The legend has to say which allele the two columns describe, since both
    # report the weakest one rather than an average.
    assert "weakest spared allele" in result.stdout


def test_design_default_sort_is_auto_discrimination_when_discriminating(
    target_fasta, panel_fasta
):
    # With no --sort, a discriminating request ranks by discrimination (auto),
    # so the bare run matches an explicit --sort discrimination.
    implicit = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", _ALLELE_1, "--mode", "any", "--no-seed",
        ],
    )
    discrimination = _design_sorted(target_fasta, panel_fasta, "discrimination")
    assert implicit.exit_code == 0
    assert implicit.stdout == discrimination.stdout


def test_design_sort_position_is_still_available(target_fasta, panel_fasta):
    result = _design_sorted(target_fasta, panel_fasta, "position")
    assert result.exit_code == 0


def test_design_target_specific_reports_mode(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--target-specific", "--mode", "any", "--no-seed",
        ],
    )
    assert result.exit_code == 0
    assert "target-specific" in result.stdout
    assert "spare all" in result.stdout


def test_design_explain_shows_per_allele_breakdown(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", _ALLELE_1, "--mode", "any", "--no-seed",
            "--explain",
        ],
    )
    assert result.exit_code == 0
    assert "Per-allele discrimination" in result.stdout


def test_design_sort_discrimination_needs_something_spared(target_fasta, panel_fasta):
    # Pan-allele targeting asks for no discrimination, so there is none to rank.
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--pan-allele", "--sort", "discrimination",
        ],
    )
    assert result.exit_code != 0


def test_design_unknown_sort_errors(target_fasta, panel_fasta):
    assert _design_sorted(target_fasta, panel_fasta, "nonsense").exit_code != 0


def test_design_targeting_without_panel_errors(target_fasta):
    result = runner.invoke(app, ["design", str(target_fasta), "--pan-allele"])
    assert result.exit_code != 0


def test_design_unknown_allele_errors(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", "NOT_IN_PANEL",
        ],
    )
    assert result.exit_code == 1
    assert "Design error" in result.output


def test_design_spare_without_hit_errors(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--spare", _ALLELE_1,
        ],
    )
    assert result.exit_code != 0


def test_design_mutually_exclusive_targeting_errors(target_fasta, panel_fasta):
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--pan-allele", "--allele-specific", _ALLELE_1,
        ],
    )
    assert result.exit_code != 0


def test_design_panel_tsv_has_allele_columns(target_fasta, panel_fasta, tmp_path):
    out = tmp_path / "allele.tsv"
    result = runner.invoke(
        app,
        [
            "design", str(target_fasta), "--panel", str(panel_fasta),
            "--allele-specific", _ALLELE_1, "--mode", "any", "--no-seed",
            "-o", str(out),
        ],
    )
    assert result.exit_code == 0
    header = out.read_text(encoding="ascii").splitlines()[0]
    assert "allele_hit_max_mismatch" in header
    assert "allele_spare_min_mismatch" in header
    # The full deciding-mismatch criteria, not just position, reach the TSV.
    assert "allele_spare_count" in header
    assert "disc_strongest_position" in header
    assert "disc_multi_mismatch" in header


# --------------------------------------------------------------------------- #
# convert command (stage 4: siRNA -> shRNA)
# --------------------------------------------------------------------------- #


def _design_tsv(target_fasta, tmp_path):
    """Produce a selected-candidate TSV to feed the converter."""
    out = tmp_path / "selected.tsv"
    result = runner.invoke(
        app,
        ["design", str(target_fasta), "--no-seed", "--no-offtarget", "-o", str(out)],
    )
    assert result.exit_code == 0, result.output
    return out


def test_convert_runs_and_reports_scaffold(target_fasta, tmp_path):
    sel = _design_tsv(target_fasta, tmp_path)
    result = runner.invoke(app, ["convert", str(sel), "--loop", "jensen"])
    assert result.exit_code == 0, result.output
    assert "stem-loop (SA)" in result.stdout
    assert "UGUGCUU" in result.stdout  # the jensen loop is echoed
    assert "fold as designed" in result.stdout


def test_convert_writes_tsv_with_expected_columns_and_rows(target_fasta, tmp_path):
    sel = _design_tsv(target_fasta, tmp_path)
    out = tmp_path / "shrna.tsv"
    result = runner.invoke(app, ["convert", str(sel), "-o", str(out)])
    assert result.exit_code == 0, result.output
    lines = out.read_text(encoding="ascii").splitlines()
    assert lines[0].startswith("region\tstart\tend\tloop_name\tloop_seq")
    assert "insert_top_dna_5to3" in lines[0]
    assert "folds_as_designed" in lines[0]
    # One construct per selected candidate in the input TSV.
    assert len(lines) - 1 == len(sel.read_text(encoding="ascii").splitlines()) - 1


def test_convert_bottom_oligo_is_reverse_complement_of_top(target_fasta, tmp_path):
    sel = _design_tsv(target_fasta, tmp_path)
    out = tmp_path / "shrna.tsv"
    runner.invoke(app, ["convert", str(sel), "-o", str(out)])
    header, first = out.read_text(encoding="ascii").splitlines()[:2]
    cols = header.split("\t")
    row = dict(zip(cols, first.split("\t")))
    top, bottom = row["insert_top_dna_5to3"], row["insert_bottom_dna_5to3"]
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    assert bottom == "".join(complement[b] for b in reversed(top))


def test_convert_unknown_loop_errors(target_fasta, tmp_path):
    sel = _design_tsv(target_fasta, tmp_path)
    result = runner.invoke(app, ["convert", str(sel), "--loop", "bogus"])
    assert result.exit_code != 0


def test_convert_rejects_non_design_tsv(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("foo\tbar\n1\t2\n", encoding="ascii")
    result = runner.invoke(app, ["convert", str(bad)])
    assert result.exit_code == 1
    assert "Convert error" in result.output
