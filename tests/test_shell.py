"""Tests for the interactive shell (app.shell.HushShell).

The command loop is driven directly through ``dispatch`` - one line at a time -
so state and output can be checked without a live stdin. The rendering itself is
the CLI's (covered by test_cli); these tests cover the shell's state, its flag
handling, and that a bad command never tears down the session.
"""

import pytest

from app.shell import HushShell, apply_convert_flags, apply_design_flags
from service.params import AlleleTargetingSpec, DesignParams

_ALLELE_1 = "HLA-DEMO*01:01"


# --------------------------------------------------------------------------- #
# apply_design_flags (unit)
# --------------------------------------------------------------------------- #


def test_apply_flags_only_changes_named_options():
    base = DesignParams()
    updated = apply_design_flags(base, ["--no-seed", "--length", "21"])
    assert updated.use_seed_filter is False
    assert updated.length == 21
    # Everything else is inherited unchanged.
    assert updated.use_ui_tei == base.use_ui_tei
    assert updated.rule_mode == base.rule_mode


def test_apply_flags_gc_percent_becomes_fraction():
    updated = apply_design_flags(DesignParams(), ["--gc", "--gc-min", "40", "--gc-max", "60"])
    assert updated.use_gc_filter is True
    assert updated.gc_min == pytest.approx(0.40)
    assert updated.gc_max == pytest.approx(0.60)


def test_apply_flags_allow_homopolymers_disables_gate():
    updated = apply_design_flags(DesignParams(), ["--allow-homopolymers"])
    assert updated.use_homopolymer_filter is False


def test_apply_flags_targeting_builds_spec():
    updated = apply_design_flags(DesignParams(), ["--allele-specific", _ALLELE_1])
    assert isinstance(updated.targeting, AlleleTargetingSpec)
    assert updated.targeting.allele_specific == _ALLELE_1


def test_apply_flags_unknown_rule_raises():
    with pytest.raises(ValueError):
        apply_design_flags(DesignParams(), ["--rules", "bogus"])


def test_apply_flags_unknown_option_raises():
    with pytest.raises(ValueError):
        apply_design_flags(DesignParams(), ["--frobnicate"])


def test_apply_convert_flags_builds_options():
    options = apply_convert_flags(["--loop", "psuper", "--no-terminator"])
    assert options.loop.name == "psuper"
    assert options.include_terminator is False


# --------------------------------------------------------------------------- #
# loading and state
# --------------------------------------------------------------------------- #


def test_load_sets_target(target_fasta, capsys):
    shell = HushShell()
    assert shell.dispatch(f"load {target_fasta}") is True
    assert shell.session.target is not None
    assert "loaded target" in capsys.readouterr().out


def test_panel_sets_panel(target_fasta, panel_fasta):
    shell = HushShell()
    shell.dispatch(f"load {target_fasta}")
    shell.dispatch(f"panel {panel_fasta}")
    assert shell.session.has_panel


def test_set_persists_criteria(target_fasta):
    shell = HushShell()
    shell.dispatch("set --no-seed --gc")
    assert shell.params.use_seed_filter is False
    assert shell.params.use_gc_filter is True


def test_set_no_args_prints_criteria(capsys):
    HushShell().dispatch("set")
    assert "Criteria" in capsys.readouterr().out


def test_reset_restores_defaults():
    shell = HushShell()
    shell.dispatch("set --no-seed")
    shell.dispatch("reset")
    assert shell.params.use_seed_filter is DesignParams().use_seed_filter


def test_limit_updates_display_cap():
    shell = HushShell()
    shell.dispatch("limit 5")
    assert shell.limit == 5


# --------------------------------------------------------------------------- #
# design / convert
# --------------------------------------------------------------------------- #


def test_design_without_target_reports_error(capsys):
    shell = HushShell()
    assert shell.dispatch("design") is True  # the shell survives the error
    assert "error:" in capsys.readouterr().out
    assert shell.session.last_design is None


def test_design_runs_and_records_last_design(demo_lowgc_fasta, capsys):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    capsys.readouterr()
    shell.dispatch("design")
    assert shell.session.last_design is not None
    assert "selected:" in capsys.readouterr().out


def test_design_inline_flags_do_not_persist(demo_lowgc_fasta):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    shell.dispatch("design --no-seed")  # one-off override
    assert shell.params.use_seed_filter is True  # persistent criteria untouched


def test_convert_after_design(demo_lowgc_fasta, capsys):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    shell.dispatch("design")
    capsys.readouterr()
    shell.dispatch("convert --loop jensen")
    assert "stem-loop (SA)" in capsys.readouterr().out


def test_convert_without_design_reports_error(target_fasta, capsys):
    shell = HushShell()
    shell.dispatch(f"load {target_fasta}")
    capsys.readouterr()
    assert shell.dispatch("convert") is True
    assert "error:" in capsys.readouterr().out


# --------------------------------------------------------------------------- #
# save, control flow, robustness
# --------------------------------------------------------------------------- #


def test_save_writes_both_tsvs_after_convert(demo_lowgc_fasta, tmp_path):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    shell.dispatch("design")
    shell.dispatch("convert")
    out = tmp_path / "run1"
    shell.dispatch(f"save {out}")
    assert (out / "design.tsv").exists()
    assert (out / "shrna.tsv").exists()


def test_save_writes_design_only_without_convert(demo_lowgc_fasta, tmp_path):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    shell.dispatch("design")
    out = tmp_path / "run1"
    shell.dispatch(f"save {out}")
    assert (out / "design.tsv").exists()
    assert not (out / "shrna.tsv").exists()


def test_new_design_clears_stale_convert(demo_lowgc_fasta):
    shell = HushShell()
    shell.dispatch(f"load {demo_lowgc_fasta}")
    shell.dispatch("design")
    shell.dispatch("convert")
    assert shell.last_convert is not None
    shell.dispatch("design")  # a fresh design supersedes the old convert
    assert shell.last_convert is None


def test_save_without_design_reports_error(tmp_path, capsys):
    shell = HushShell()
    assert shell.dispatch(f"save {tmp_path / 'x'}") is True
    assert "error:" in capsys.readouterr().out


def test_quit_ends_loop_but_help_does_not():
    shell = HushShell()
    assert shell.dispatch("help") is True
    assert shell.dispatch("quit") is False


def test_unknown_command_is_reported(capsys):
    shell = HushShell()
    assert shell.dispatch("frobnicate x") is True
    assert "unknown command" in capsys.readouterr().out


def test_bad_flag_does_not_tear_down_shell(target_fasta, capsys):
    shell = HushShell()
    shell.dispatch(f"load {target_fasta}")
    capsys.readouterr()
    assert shell.dispatch("set --bogus") is True
    assert "error:" in capsys.readouterr().out
    # The session is intact and still usable afterwards.
    assert shell.session.target is not None
