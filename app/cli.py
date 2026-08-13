"""HUSH command-line interface.

CLI-first entry point for the pipeline. Currently exposes the Input Layer via
'hush input'; siRNA design, MSA/allele, and shRNA conversion commands are
added as those stages come online.

Run 'hush --help' (after 'pip install -e .') or 'python -m app.cli --help'.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from core.asymmetry import DEFAULT_WINDOW_BP
from core.errors import HushError
from core.hairpin import DEFAULT_INTERNAL_TM_THRESHOLD_C
from core.inputs import build_design_input, load_sequences
from core.rules import AMARZGUIOUI_DEFAULT_THRESHOLD, REYNOLDS_DEFAULT_THRESHOLD
from core.seed import SEED_TM_THRESHOLD_C
from core.selection import RuleMode, SelectionCriteria, evaluate_all
from core.sequence import BioSequence, DesignInput
from core.sirna import enumerate_candidates
from core.thermo import REFERENCE_TEMPERATURE_C, ThermoParams

# Help-panel label that groups the expert thermodynamic knobs in `design --help`.
_ADVANCED = "Advanced (thermodynamic internals)"

app = typer.Typer(
    add_completion=False,
    help="HUSH - Human genome Universal Silencing Hairpin-designer.",
    no_args_is_help=True,
)


def _echo_sequence(seq: BioSequence, index: int | None = None) -> None:
    prefix = f"  [{index}] " if index is not None else "  "
    header = f"{prefix}{seq.id or '<unnamed>'}"
    typer.echo(typer.style(header, bold=True))
    if seq.description:
        typer.echo(f"      {seq.description}")
    line = f"      type={seq.moltype.value}  length={seq.length}"
    if seq.is_nucleotide and seq.length:
        line += f"  GC={seq.gc_content:.1%}"
    typer.echo(line)
    preview = seq.seq if seq.length <= 60 else f"{seq.seq[:57]}..."
    typer.echo(f"      {preview}")


def _summarize_design_input(design: DesignInput) -> None:
    typer.echo(typer.style("Target", fg=typer.colors.CYAN, bold=True))
    _echo_sequence(design.target)

    if design.has_panel:
        assert design.panel is not None
        typer.echo("")
        typer.echo(
            typer.style(
                f"Panel: {design.panel.name}  ({len(design.panel)} sequences)",
                fg=typer.colors.CYAN,
                bold=True,
            )
        )
        uniform = "uniform" if design.panel.is_uniform else "MIXED types"
        typer.echo(
            f"  molecule types: "
            f"{', '.join(sorted(t.value for t in design.panel.moltypes))} "
            f"({uniform})"
        )
        for i, seq in enumerate(design.panel):
            _echo_sequence(seq, index=i)

    # Foreshadow the downstream requirement without failing the inspection.
    typer.echo("")
    if design.target.is_nucleotide:
        typer.echo(
            typer.style(
                "[OK] Target is a nucleotide sequence - ready for siRNA design.",
                fg=typer.colors.GREEN,
            )
        )
    else:
        typer.echo(
            typer.style(
                f"[!] Target is {design.target.moltype.value}; siRNA/shRNA design "
                "needs the coding transcript (mRNA/CDS). Supply the nucleotide "
                "sequence before designing.",
                fg=typer.colors.YELLOW,
            )
        )


@app.command("input")
def input_cmd(
    target: Path = typer.Argument(
        ...,
        exists=True,
        dir_okay=False,
        readable=True,
        help="Target sequence file (FASTA or PDB).",
    ),
    panel: Optional[Path] = typer.Option(
        None,
        "--panel",
        "-p",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Optional related-sequence panel (e.g. an HLA allele family FASTA).",
    ),
    target_format: Optional[str] = typer.Option(
        None, "--target-format", help="Override target format: fasta | pdb | cif."
    ),
    panel_format: Optional[str] = typer.Option(
        None, "--panel-format", help="Override panel format: fasta | pdb | cif."
    ),
) -> None:
    """Load and inspect target (and optional panel) sequences.

    This exercises the Input Layer: it parses the file(s), validates residues,
    auto-detects molecule type, and reports what the downstream pipeline will
    receive.
    """
    try:
        design = build_design_input(
            target,
            panel,
            target_fmt=target_format,
            panel_fmt=panel_format,
        )
    except HushError as exc:
        typer.secho(f"Input error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    _summarize_design_input(design)


@app.command("inspect")
def inspect_cmd(
    path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Any FASTA/PDB file to list sequences from.",
    ),
    fmt: Optional[str] = typer.Option(
        None, "--format", "-f", help="Override format: fasta | pdb | cif."
    ),
) -> None:
    """List every sequence in a file (no target/panel roles assigned)."""
    try:
        seqs = load_sequences(path, fmt=fmt)
    except HushError as exc:
        typer.secho(f"Input error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    typer.echo(f"{len(seqs)} sequence(s) in {path.name}:")
    for i, seq in enumerate(seqs):
        _echo_sequence(seq, index=i)


# --------------------------------------------------------------------------- #
# design command
# --------------------------------------------------------------------------- #

_RULE_ALIASES = {
    "ui-tei": "ui_tei",
    "uitei": "ui_tei",
    "ui_tei": "ui_tei",
    "reynolds": "reynolds",
    "amarzguioui": "amarzguioui",
    "amarz": "amarzguioui",
}


def _parse_rules(spec: str) -> tuple[bool, bool, bool]:
    """Turn a --rules spec into (use_ui_tei, use_reynolds, use_amarzguioui)."""
    text = spec.strip().lower()
    if text in ("", "all"):
        return True, True, True
    if text == "none":
        return False, False, False
    flags = {"ui_tei": False, "reynolds": False, "amarzguioui": False}
    for token in (t.strip() for t in text.split(",")):
        if not token:
            continue
        key = _RULE_ALIASES.get(token)
        if key is None:
            raise typer.BadParameter(
                f"unknown rule {token!r}; choose from ui-tei, reynolds, "
                "amarzguioui (comma-separated), or 'all'/'none'.",
                param_hint="--rules",
            )
        flags[key] = True
    return flags["ui_tei"], flags["reynolds"], flags["amarzguioui"]


def _parse_mode(spec: str) -> RuleMode:
    try:
        return RuleMode(spec.strip().lower())
    except ValueError:
        raise typer.BadParameter(
            f"unknown mode {spec!r}; choose 'all' or 'any'.", param_hint="--mode"
        )


_PARAMS_ALIASES = {
    "xia": ThermoParams.XIA_1998,
    "xia1998": ThermoParams.XIA_1998,
    "freier": ThermoParams.FREIER_1986,
    "freier1986": ThermoParams.FREIER_1986,
}


def _parse_params(spec: str) -> ThermoParams:
    key = _PARAMS_ALIASES.get(spec.strip().lower())
    if key is None:
        raise typer.BadParameter(
            f"unknown parameter set {spec!r}; choose 'xia' or 'freier'.",
            param_hint="--params",
        )
    return key


def _rules_code(evaluation, criteria: SelectionCriteria) -> str:
    """Compact per-rule pass code, e.g. 'UR-' (Ui-Tei + Reynolds pass, Amarz fails)."""
    parts = []
    if evaluation.ui_tei is not None:
        parts.append("U" if evaluation.ui_tei.satisfied else "-")
    if evaluation.reynolds is not None:
        passed = evaluation.reynolds.meets_threshold(criteria.reynolds_threshold)
        parts.append("R" if passed else "-")
    if evaluation.amarzguioui is not None:
        passed = evaluation.amarzguioui.meets_threshold(criteria.amarzguioui_threshold)
        parts.append("A" if passed else "-")
    return "".join(parts) or "-"


def _print_criteria(criteria: SelectionCriteria) -> None:
    typer.echo(typer.style("Criteria", bold=True))
    if criteria.enabled_rule_names:
        names = ", ".join(n.replace("_", "-") for n in criteria.enabled_rule_names)
        typer.echo(f"  rules:      {names}  (mode: {criteria.rule_mode.value.upper()})")
    else:
        typer.echo("  rules:      none")
    seed = (
        f"< {criteria.seed_tm_threshold:.1f} C (both strands)"
        if criteria.use_seed_filter
        else "off"
    )
    typer.echo(f"  seed Tm:    {seed}")
    gc = (
        f"{criteria.gc_min * 100:.0f}-{criteria.gc_max * 100:.0f}%"
        if criteria.use_gc_filter
        else "off"
    )
    typer.echo(f"  GC filter:  {gc}")
    hairpin = (
        f"Tm <= {criteria.hairpin_tm_max:.1f} C"
        if criteria.use_hairpin_filter
        else "off"
    )
    typer.echo(f"  hairpin:    {hairpin}")
    asym = (
        f"favour guide ({criteria.asymmetry_window} bp @ {criteria.asymmetry_temp_c:.0f} C)"
        if criteria.use_asymmetry_filter
        else "off"
    )
    typer.echo(f"  asymmetry:  {asym}")
    if criteria.use_asymmetry_filter or criteria.use_hairpin_filter:
        typer.echo(f"  params:     {criteria.thermo_params.value}")


def _print_funnel(criteria: SelectionCriteria, evaluations: list, selected: list) -> None:
    def row(label: str, count: int) -> str:
        return f"  {label:<16}{count}"

    typer.echo(typer.style("Passing each enabled gate", bold=True))
    typer.echo(
        row(f"rules ({criteria.rule_mode.value}):", sum(1 for e in evaluations if e.rule_pass))
    )
    if criteria.use_gc_filter:
        typer.echo(row("GC window:", sum(1 for e in evaluations if e.gc_pass)))
    if criteria.use_seed_filter:
        typer.echo(row("seed Tm:", sum(1 for e in evaluations if e.seed_pass)))
    if criteria.use_hairpin_filter:
        typer.echo(row("hairpin:", sum(1 for e in evaluations if e.hairpin_pass)))
    if criteria.use_asymmetry_filter:
        typer.echo(row("asymmetry:", sum(1 for e in evaluations if e.asymmetry_pass)))
    typer.echo(typer.style(row("selected:", len(selected)), fg=typer.colors.GREEN, bold=True))


def _print_table(
    selected: list, criteria: SelectionCriteria, length: int, limit: int
) -> None:
    seq_w = max(length + 1, 14)
    header = (
        f"  {'region':<11} {'guide 5->3':<{seq_w}} {'passenger 5->3':<{seq_w}} {'GC%':>4}"
    )
    if criteria.enabled_rule_names:
        header += f" {'rules':<5}"
    if criteria.use_seed_filter:
        header += f" {'seedTm':>7}"
    if criteria.use_hairpin_filter:
        header += f" {'hpTm':>6}"
    typer.echo(typer.style(header, bold=True))

    shown = selected if limit == 0 else selected[:limit]
    for e in shown:
        c = e.candidate
        row = (
            f"  {c.region:<11} {c.guide.seq:<{seq_w}} {c.passenger.seq:<{seq_w}} "
            f"{c.gc_content * 100:>4.0f}"
        )
        if criteria.enabled_rule_names:
            row += f" {_rules_code(e, criteria):<5}"
        if criteria.use_seed_filter:
            row += f" {e.seed.max_tm:>7.1f}"
        if criteria.use_hairpin_filter:
            row += f" {('  none' if e.hairpin_tm is None else f'{e.hairpin_tm:.1f}'):>6}"
        typer.echo(row)
    if limit and len(selected) > limit:
        typer.echo(
            f"  ... {len(selected) - limit} more "
            "(use --limit 0 to show all, or -o file.tsv)"
        )


def _write_design_tsv(path: Path, selected: list, criteria: SelectionCriteria) -> None:
    lines = [
        "region\tstart\tend\tguide_5to3\tpassenger_5to3\tgc_percent\t"
        "rules_passed\tseed_max_tm_c\thairpin_tm_c"
    ]
    for e in selected:
        c = e.candidate
        seed_tm_val = f"{e.seed.max_tm:.2f}" if e.seed is not None else ""
        hairpin_val = "" if e.hairpin_tm is None else f"{e.hairpin_tm:.2f}"
        lines.append(
            "\t".join(
                [
                    c.region,
                    str(c.position + 1),
                    str(c.end),
                    c.guide.seq,
                    c.passenger.seq,
                    f"{c.gc_content * 100:.1f}",
                    _rules_code(e, criteria),
                    seed_tm_val,
                    hairpin_val,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


@app.command("design")
def design_cmd(
    target: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Target mRNA/CDS file (FASTA or PDB).",
    ),
    length: int = typer.Option(
        19, "--length", "-l", help="siRNA core length (19-24)."
    ),
    rules: str = typer.Option(
        "all", "--rules",
        help="Rule sets: comma list of ui-tei,reynolds,amarzguioui, or all/none.",
    ),
    mode: str = typer.Option(
        "all", "--mode",
        help="Combine enabled rules by 'all' (intersection) or 'any' (union).",
    ),
    reynolds_threshold: int = typer.Option(
        REYNOLDS_DEFAULT_THRESHOLD, "--reynolds-threshold",
        help="Reynolds score cut-off.",
    ),
    amarzguioui_threshold: int = typer.Option(
        AMARZGUIOUI_DEFAULT_THRESHOLD, "--amarzguioui-threshold",
        help="Amarzguioui score cut-off.",
    ),
    seed: bool = typer.Option(
        True, "--seed/--no-seed",
        help="Apply the seed-duplex off-target Tm filter (both strands).",
    ),
    seed_tm: float = typer.Option(
        SEED_TM_THRESHOLD_C, "--seed-tm", help="Seed-duplex Tm cut-off in C."
    ),
    gc: bool = typer.Option(
        False, "--gc/--no-gc", help="Apply a standalone GC-content window filter."
    ),
    gc_min: float = typer.Option(
        30.0, "--gc-min", help="GC window lower bound (percent)."
    ),
    gc_max: float = typer.Option(
        52.0, "--gc-max", help="GC window upper bound (percent)."
    ),
    asymmetry: bool = typer.Option(
        False, "--asymmetry/--no-asymmetry",
        help="Require duplex-end asymmetry to favour the guide strand.",
    ),
    hairpin: bool = typer.Option(
        False, "--hairpin/--no-hairpin",
        help="Reject candidates with stable internal (self-hairpin) structure.",
        rich_help_panel=_ADVANCED,
    ),
    hairpin_tm: float = typer.Option(
        DEFAULT_INTERNAL_TM_THRESHOLD_C, "--hairpin-tm",
        help="Internal-structure Tm cut-off in C (reject if the most stable "
        "self-hairpin melts above this).",
        rich_help_panel=_ADVANCED,
    ),
    asym_window: int = typer.Option(
        DEFAULT_WINDOW_BP, "--asym-window",
        help="Terminal base pairs compared at each duplex end for asymmetry.",
        rich_help_panel=_ADVANCED,
    ),
    asym_temp: float = typer.Option(
        REFERENCE_TEMPERATURE_C, "--asym-temp",
        help="Temperature in C for the asymmetry free energy.",
        rich_help_panel=_ADVANCED,
    ),
    params: str = typer.Option(
        "xia", "--params",
        help="Nearest-neighbor parameter set for the asymmetry and hairpin "
        "readouts only (never the seed filter): xia | freier.",
        rich_help_panel=_ADVANCED,
    ),
    limit: int = typer.Option(
        20, "--limit", "-n", help="Max candidates to print (0 = all)."
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Write all selected candidates to a TSV file."
    ),
) -> None:
    """Enumerate siRNA candidates for a target and filter them by chosen criteria.

    The three classical rule sets toggle independently and combine by --mode;
    the seed-Tm, GC, hairpin and asymmetry filters are applied on top. Nothing is
    ranked: selected candidates are reported in target-position order. Expert
    thermodynamic knobs live in the Advanced option group (see --help).
    """
    use_ut, use_rey, use_amz = _parse_rules(rules)
    rule_mode = _parse_mode(mode)
    thermo_params = _parse_params(params)
    try:
        criteria = SelectionCriteria(
            use_ui_tei=use_ut,
            use_reynolds=use_rey,
            use_amarzguioui=use_amz,
            rule_mode=rule_mode,
            reynolds_threshold=reynolds_threshold,
            amarzguioui_threshold=amarzguioui_threshold,
            use_seed_filter=seed,
            seed_tm_threshold=seed_tm,
            use_gc_filter=gc,
            gc_min=gc_min / 100.0,
            gc_max=gc_max / 100.0,
            use_hairpin_filter=hairpin,
            hairpin_tm_max=hairpin_tm,
            use_asymmetry_filter=asymmetry,
            asymmetry_window=asym_window,
            asymmetry_temp_c=asym_temp,
            thermo_params=thermo_params,
        )
        design_input = build_design_input(target)
        mrna = design_input.require_nucleotide_target()
        candidates = list(enumerate_candidates(mrna, length=length))
        evaluations = list(evaluate_all(candidates, criteria))
    except (HushError, ValueError) as exc:
        typer.secho(f"Design error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    selected = [e for e in evaluations if e.selected]

    typer.echo(typer.style(f"Target: {mrna.id or '<unnamed>'}", fg=typer.colors.CYAN, bold=True))
    typer.echo(f"  {len(mrna)} nt  |  {len(evaluations)} candidate {length}-mer(s)")
    typer.echo("")
    _print_criteria(criteria)
    typer.echo("")
    _print_funnel(criteria, evaluations, selected)
    typer.echo("")

    if not selected:
        typer.secho(
            "No candidates selected. Try relaxing: --mode any, --no-seed, "
            "or a wider --gc range.",
            fg=typer.colors.YELLOW,
        )
    else:
        _print_table(selected, criteria, length, limit)

    if output is not None:
        _write_design_tsv(output, selected, criteria)
        typer.echo(f"\nWrote {len(selected)} selected candidate(s) to {output}")


if __name__ == "__main__":
    app()
