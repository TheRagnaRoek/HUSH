"""HUSH command-line interface.

CLI-first entry point for the pipeline. Exposes the Input Layer via
'hush inspect' and siRNA candidate design via 'hush design'. The MSA/allele
stage is engaged from within 'hush design' by supplying '--panel': the panel is
aligned to the target and each candidate is additionally gated on its allele
targeting. The shRNA conversion command is added as that stage comes online.

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
from core.offtarget import (
    DEFAULT_OFFTARGET_MAX_MISMATCH,
    DEFAULT_OFFTARGET_MIN_MISMATCH,
)
from msa.discrimination import (
    LOW_TOLERANCE_WINDOW,
    MULTI_MISMATCH,
    STRONGEST_POSITIONS,
    GuideRegion,
    MismatchAssessment,
    MismatchClass,
    ordered_by_placement,
)
from msa.targeting import AlleleTargeting, TargetingMode
from msa.windows import DEFAULT_MIN_SPARE_MISMATCHES, WindowEvaluation
from core.rules import AMARZGUIOUI_DEFAULT_THRESHOLD, REYNOLDS_DEFAULT_THRESHOLD
from core.seed import SEED_TM_THRESHOLD_C
from core.selection import DEFAULT_MAX_HOMOPOLYMER_RUN, RuleMode, SelectionCriteria
from core.sequence import BioSequence, DesignInput, MolType
from core.sirna import OverhangPolicy
from core.thermo import REFERENCE_TEMPERATURE_C, ThermoParams
from app.parsing import (
    parse_mode,
    parse_overhang,
    parse_params,
    parse_rules,
    split_ids,
)
from service.params import AlleleTargetingSpec, DesignParams
from service.results import ConvertResult, DesignResult
from service.session import HushSession
from shrna.convert import DEFAULT_TERMINATOR_LENGTH, ShrnaOptions, convert
from shrna.loops import DEFAULT_LOOP, LOOPS, LoopScaffold, get_loop

# Help-panel label that groups the expert thermodynamic knobs in `design --help`.
_ADVANCED = "Advanced (thermodynamic internals)"

# Help-panel label grouping the allele-aware (MSA panel) options in `design --help`.
_ALLELE = "Allele-aware (MSA panel targeting)"

app = typer.Typer(
    add_completion=False,
    help="HUSH - Human genome Universal Silencing Hairpin-designer.",
    no_args_is_help=True,
)


def _echo_sequence(
    seq: BioSequence, index: int | None = None, show_readiness: bool = False
) -> None:
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
    if show_readiness:
        if seq.is_nucleotide:
            typer.echo(
                typer.style(
                    "      [OK] nucleotide - ready for siRNA design.",
                    fg=typer.colors.GREEN,
                )
            )
        else:
            typer.echo(
                typer.style(
                    f"      [!] {seq.moltype.value}; siRNA/shRNA design needs the "
                    "coding transcript (mRNA/CDS) first.",
                    fg=typer.colors.YELLOW,
                )
            )


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


@app.command("inspect")
def inspect_cmd(
    path: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Sequence file to inspect (FASTA or PDB/mmCIF).",
    ),
    panel: Optional[Path] = typer.Option(
        None,
        "--panel",
        "-p",
        exists=True,
        dir_okay=False,
        readable=True,
        help="Optional related-sequence panel (e.g. an HLA allele family FASTA). "
        "When given, PATH is treated as the design target.",
    ),
    fmt: Optional[str] = typer.Option(
        None, "--format", "-f", help="Override PATH format: fasta | pdb | cif."
    ),
    panel_format: Optional[str] = typer.Option(
        None, "--panel-format", help="Override panel format: fasta | pdb | cif."
    ),
) -> None:
    """Inspect sequences in a file, or preview the design input.

    Without --panel, lists every sequence in the file (no roles assigned) and
    flags each as ready for siRNA design or needing a coding transcript first.
    With --panel, treats PATH as the design target: it parses both files,
    validates residues, auto-detects molecule type, assigns target/panel roles,
    and reports what the downstream siRNA pipeline will receive.
    """
    try:
        if panel is None:
            seqs = load_sequences(path, fmt=fmt)
            typer.echo(f"{len(seqs)} sequence(s) in {path.name}:")
            for i, seq in enumerate(seqs):
                _echo_sequence(seq, index=i, show_readiness=True)
        else:
            design = build_design_input(
                path,
                panel,
                target_fmt=fmt,
                panel_fmt=panel_format,
            )
            _summarize_design_input(design)
    except HushError as exc:
        typer.secho(f"Input error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)


# --------------------------------------------------------------------------- #
# design command
# --------------------------------------------------------------------------- #

# The flag-vocabulary conversions live in app.parsing (shared with the shell)
# and raise ValueError; here they are re-wrapped as Typer usage errors so a bad
# --flag still reports with its hint and the usage exit code.


def _parse_rules(spec: str) -> tuple[bool, bool, bool]:
    try:
        return parse_rules(spec)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--rules")


def _parse_mode(spec: str) -> RuleMode:
    try:
        return parse_mode(spec)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--mode")


def _parse_params(spec: str) -> ThermoParams:
    try:
        return parse_params(spec)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--params")


def _parse_overhang(spec: str) -> OverhangPolicy:
    try:
        return parse_overhang(spec)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--overhang")


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


def _print_criteria(
    criteria: SelectionCriteria,
    targeting: AlleleTargeting | None = None,
    min_spare: int = DEFAULT_MIN_SPARE_MISMATCHES,
) -> None:
    typer.echo(typer.style("Criteria", bold=True))
    if targeting is not None:
        mode = targeting.mode.value.replace("_", "-")
        if targeting.mode is TargetingMode.TARGET_SPECIFIC:
            allele = (
                f"{mode}: loaded target only, spare all {len(targeting.spare)} "
                f"panel allele(s) (>= {min_spare} diff each)"
            )
        elif targeting.discriminating:
            allele = (
                f"{mode}: hit {sorted(targeting.hit)}, spare "
                f"{len(targeting.spare)} allele(s) (>= {min_spare} diff each)"
            )
        else:
            allele = f"{mode}: hit {len(targeting.hit)} allele(s), spare none"
        typer.echo(f"  allele:     {allele}")
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
    if criteria.use_offtarget_filter:
        off = f"reject nearest < {criteria.offtarget_min_mismatch} mismatch(es)"
        if criteria.offtarget_self_accession:
            off += f" (self: {criteria.offtarget_self_accession})"
    else:
        off = "off"
    typer.echo(f"  off-target: {off}")
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
    homopolymer = (
        f"reject runs > {criteria.max_homopolymer_run} nt"
        if criteria.use_homopolymer_filter
        else "off"
    )
    typer.echo(f"  homopolymer:{homopolymer}")
    asym = (
        f"favour guide ({criteria.asymmetry_window} bp @ {criteria.asymmetry_temp_c:.0f} C)"
        if criteria.use_asymmetry_filter
        else "off"
    )
    typer.echo(f"  asymmetry:  {asym}")
    if criteria.use_asymmetry_filter or criteria.use_hairpin_filter:
        typer.echo(f"  params:     {criteria.thermo_params.value}")


def _print_funnel(
    criteria: SelectionCriteria,
    evaluations: list,
    selected: list,
    window_by_pos: dict[int, WindowEvaluation] | None = None,
) -> None:
    def row(label: str, count: int) -> str:
        return f"  {label:<16}{count}"

    typer.echo(typer.style("Passing each enabled gate", bold=True))
    typer.echo(
        row(f"rules ({criteria.rule_mode.value}):", sum(1 for e in evaluations if e.rule_pass))
    )
    if window_by_pos is not None:
        typer.echo(
            row(
                "allele:",
                sum(
                    1
                    for e in evaluations
                    if window_by_pos[e.candidate.position].satisfies
                ),
            )
        )
    if criteria.use_gc_filter:
        typer.echo(row("GC window:", sum(1 for e in evaluations if e.gc_pass)))
    if criteria.use_seed_filter:
        typer.echo(row("seed Tm:", sum(1 for e in evaluations if e.seed_pass)))
    if criteria.use_offtarget_filter:
        typer.echo(row("off-target:", sum(1 for e in evaluations if e.offtarget_pass)))
    if criteria.use_hairpin_filter:
        typer.echo(row("hairpin:", sum(1 for e in evaluations if e.hairpin_pass)))
    if criteria.use_homopolymer_filter:
        typer.echo(row("homopolymer:", sum(1 for e in evaluations if e.homopolymer_pass)))
    if criteria.use_asymmetry_filter:
        typer.echo(row("asymmetry:", sum(1 for e in evaluations if e.asymmetry_pass)))
    typer.echo(typer.style(row("selected:", len(selected)), fg=typer.colors.GREEN, bold=True))


# Short codes for the criteria that set a mismatch's discrimination rank, so a
# whole assessment reads in one cell instead of just its position.
_REGION_CODE = {
    GuideRegion.POST_SEED: "post",
    GuideRegion.SEED: "seed",
    GuideRegion.ANCHOR: "anch",
}
_CLASS_CODE = {
    MismatchClass.PURINE_PURINE: "pp",
    MismatchClass.PURINE_PYRIMIDINE: "pu",
    MismatchClass.PYRIMIDINE_PURINE: "up",
    MismatchClass.PYRIMIDINE_PYRIMIDINE: "uu",
}


def _position_tier(guide_position: int) -> str:
    """S for a strongest position (P10/P13/P16), W inside the window, else '-'."""
    if guide_position in STRONGEST_POSITIONS:
        return "S"
    if guide_position in LOW_TOLERANCE_WINDOW:
        return "W"
    return "-"


def _mismatch_code(d: MismatchAssessment) -> str:
    """Every ranking criterion of one mismatch in one cell.

    'P16 post H+ S pp' reads: guide position 16; 3' of the seed; listed in
    Huang's table for this base change (H+, vs H- not listed, H? ambiguous); a
    strongest position (S; W = inside the P9-P16 window, - = outside);
    purine:purine.
    """
    huang = "H?" if d.huang_recommended is None else ("H+" if d.huang_recommended else "H-")
    cls = "?" if d.mismatch_class is None else _CLASS_CODE[d.mismatch_class]
    return f"P{d.guide_position} {_REGION_CODE[d.region]} {huang} {_position_tier(d.guide_position)} {cls}"


def _other_positions(window: WindowEvaluation) -> str:
    """The weakest allele's remaining differences, best-placed first.

    The rank rests on one mismatch, but the others are what let a user judge
    count against placement - the trade-off msa.discrimination refuses to make.
    Collapsing the row to the deciding mismatch alone would hide them.
    """
    rest = window.deciding_mismatches[1:]
    if not rest:
        return ""
    return " +" + ",".join(f"P{a.guide_position}" for a in rest)


def _discrimination_cell(window: WindowEvaluation) -> str:
    """Full-criteria readout of the mismatch a window's discrimination rests on,
    followed by the weakest allele's remaining difference positions.

    See _mismatch_code for the encoding. 'P16 post H+ S pp +P11,P4' is a window
    whose rank rests on P16 but which differs from that allele at three places.
    """
    deciding = window.deciding_mismatch
    if deciding is None:
        return "none"
    return _mismatch_code(deciding) + _other_positions(window)


def _discrimination_fields(window: WindowEvaluation) -> list[str]:
    """The deciding mismatch as ten TSV fields. Empty ones where it has none.

    The last field carries every difference to the weakest spared allele, in
    ranking order, so a downstream reader is not left with the best one only.
    """
    deciding = window.deciding_mismatch
    if deciding is None:
        return [""] * 10
    huang = "" if deciding.huang_recommended is None else str(deciding.huang_recommended)
    count = window.min_spare_mismatches or 0
    return [
        str(deciding.guide_position),
        deciding.mutation_type,
        "" if deciding.mismatch_class is None else deciding.mismatch_class.value,
        huang,
        deciding.region.value,
        str(deciding.sources_agree),
        str(deciding.guide_position in STRONGEST_POSITIONS),
        str(deciding.guide_position in LOW_TOLERANCE_WINDOW),
        str(count >= MULTI_MISMATCH),
        ",".join(f"P{a.guide_position}" for a in window.deciding_mismatches),
    ]


def _print_window_discrimination(window: WindowEvaluation) -> None:
    """Per-spared-allele breakdown for one window, weakest allele first.

    Where the table row collapses to the deciding (weakest) allele, this shows
    every spared allele: its difference count, all its mismatch positions, and
    its own best mismatch's full criteria - so an siRNA can be judged against the
    whole panel, not just the allele it distinguishes least well.
    """
    rows = window.spared_alleles_by_discrimination
    if not rows:
        typer.echo("      (nothing spared - no discrimination to assess)")
        return
    for index, (allele, assessments) in enumerate(rows):
        # Ranking order, not numeric order: the leftmost position is the one the
        # allele's rank rests on, and the rest are what the ranking ignores.
        ordered = ordered_by_placement(assessments)
        best = ordered[0] if ordered else None
        positions = ",".join(f"P{a.guide_position}" for a in ordered)
        code = _mismatch_code(best) if best is not None else "none"
        agree = "yes" if (best is not None and best.sources_agree) else "no"
        mutation = best.mutation_type if best is not None else ""
        marker = "  [weakest]" if index == 0 else ""
        typer.echo(
            f"      {allele:<26} diff={len(assessments):<2} [{positions}]  "
            f"{code:<17} {mutation:<4} agree={agree}{marker}"
        )


def render_discrimination_detail(
    selected: list, window_by_pos: dict[int, WindowEvaluation], limit: int
) -> None:
    """Print the full per-allele discrimination breakdown for the shown candidates."""
    shown = selected if limit == 0 else selected[:limit]
    typer.echo("")
    typer.echo(typer.style("Per-allele discrimination (weakest allele first)", bold=True))
    for e in shown:
        c = e.candidate
        typer.echo(f"  {c.region}  guide 5->3 {c.guide_full}")
        _print_window_discrimination(window_by_pos[c.position])
    typer.secho(
        "  code: P<pos> region(post/seed/anch) Huang(H+/H-/H?) "
        "position(S=P10/P13/P16, W=inside P9-P16, -=outside) class(pp/pu/up/uu); "
        "agree = Huang and region concur. Positions are listed in ranking order, "
        "the first one being what the allele's rank rests on.",
        fg=typer.colors.BRIGHT_BLACK,
    )
    if limit and len(selected) > limit:
        typer.echo(
            f"  ... {len(selected) - limit} more (use --limit 0 to show all)"
        )


def _print_table(
    selected: list,
    criteria: SelectionCriteria,
    length: int,
    limit: int,
    window_by_pos: dict[int, WindowEvaluation] | None = None,
    targeting: AlleleTargeting | None = None,
) -> None:
    # The minDiff/disc columns are shown only when the request discriminates.
    # minDiff is the weakest per-spared-allele difference count; disc places that
    # allele's best-positioned difference against the published evidence.
    show_spare = (
        window_by_pos is not None
        and targeting is not None
        and targeting.discriminating
    )
    # Strands are shown as ordered: paired core plus the 2-nt 3' overhang.
    seq_w = max(length + 3, 14)
    header = (
        f"  {'region':<11} {'guide 5->3':<{seq_w}} {'passenger 5->3':<{seq_w}} {'GC%':>4}"
    )
    if criteria.enabled_rule_names:
        header += f" {'rules':<5}"
    if show_spare:
        header += f" {'minDiff':>7} {'placement':<26}"
    if criteria.use_seed_filter:
        header += f" {'seedTm':>7}"
    if criteria.use_offtarget_filter:
        header += f" {'offMM':>6}"
    if criteria.use_hairpin_filter:
        header += f" {'hpTm':>6}"
    if criteria.use_homopolymer_filter:
        header += f" {'run':>4}"
    typer.echo(typer.style(header, bold=True))

    shown = selected if limit == 0 else selected[:limit]
    for e in shown:
        c = e.candidate
        row = (
            f"  {c.region:<11} {c.guide_full:<{seq_w}} {c.passenger_full:<{seq_w}} "
            f"{c.gc_content * 100:>4.0f}"
        )
        if criteria.enabled_rule_names:
            row += f" {_rules_code(e, criteria):<5}"
        if show_spare:
            window = window_by_pos[c.position]
            margin = window.min_spare_mismatches
            row += f" {('none' if margin is None else str(margin)):>7}"
            row += f" {_discrimination_cell(window):<26}"
        if criteria.use_seed_filter:
            row += f" {e.seed.max_tm:>7.1f}"
        if criteria.use_offtarget_filter:
            mm = None if e.offtarget is None else e.offtarget.min_mismatch
            row += f" {('none' if mm is None else str(mm)):>6}"
        if criteria.use_hairpin_filter:
            row += f" {('  none' if e.hairpin_tm is None else f'{e.hairpin_tm:.1f}'):>6}"
        if criteria.use_homopolymer_filter:
            row += f" {e.homopolymer_run:>4}"
        typer.echo(row)
    if show_spare:
        # Both allele columns report the weakest spared allele, not an average:
        # an siRNA has to spare every one of them, so the one it distinguishes
        # least well is the one that decides.
        typer.secho(
            "  minDiff and placement describe the weakest spared allele - the one "
            "this window distinguishes least well.",
            fg=typer.colors.BRIGHT_BLACK,
        )
        typer.secho(
            "  placement: P<pos> region(post/seed/anch) Huang(H+/H-/H?) "
            "position(S=P10/P13/P16, W=inside P9-P16, -=outside) class(pp/pu/up/uu), "
            "then '+' the allele's further differences.",
            fg=typer.colors.BRIGHT_BLACK,
        )
        typer.secho(
            "  Rows are ranked on the leading mismatch only. minDiff and the '+' "
            "positions are reported, never weighed against placement - no source "
            "measures that trade-off (docs/msa_ranking_note.md).",
            fg=typer.colors.BRIGHT_BLACK,
        )
        typer.secho(
            "  For every spared allele (not just the weakest), pass --explain, or "
            "'assess <region>' in the shell.",
            fg=typer.colors.BRIGHT_BLACK,
        )
    if limit and len(selected) > limit:
        typer.echo(
            f"  ... {len(selected) - limit} more "
            "(use --limit 0 to show all, or -o file.tsv)"
        )


def _combined_offtarget_counts(profile) -> list[int]:
    """Guide + passenger off-target hit counts at mismatch levels 0..3."""
    out = [0, 0, 0, 0]
    for counts in (profile.guide_counts, profile.passenger_counts):
        for i, n in enumerate(counts):
            if i < len(out):
                out[i] += n
    return out


def _write_design_tsv(
    path: Path,
    selected: list,
    criteria: SelectionCriteria,
    window_by_pos: dict[int, WindowEvaluation] | None = None,
) -> None:
    # guide_5to3 / passenger_5to3 are the full ordered strands (core + 3'
    # overhang); the *_overhang columns name the trailing overhang explicitly.
    # start/end are the paired-core coordinates on the target.
    header = (
        "region\tstart\tend\tguide_5to3\tguide_overhang\tpassenger_5to3\t"
        "passenger_overhang\tgc_percent\trules_passed\tseed_max_tm_c\t"
        "hairpin_tm_c\thomopolymer_max_run\tofftarget_min_mismatch\t"
        "offtarget_hits_0\tofftarget_hits_1\tofftarget_hits_2\tofftarget_hits_3"
    )
    if window_by_pos is not None:
        # Allele-window read-outs at this candidate's position: the worst
        # difference among hit alleles (0 for a selected candidate) and the
        # weakest discrimination among spared alleles (empty for pan-allele).
        header += (
            "\tallele_hit_max_mismatch\tallele_spare_min_mismatch"
            "\tallele_spare_count"
            # The single mismatch the window's discrimination rank rests on:
            # the best-placed difference from the weakest spared allele, with
            # every criterion that sets its rank. disc_all_positions carries
            # that allele's differences in full, in ranking order, because the
            # rank itself ignores everything after the first.
            "\tdisc_guide_position\tdisc_mutation\tdisc_mismatch_class"
            "\tdisc_huang_recommended\tdisc_region\tdisc_sources_agree"
            "\tdisc_strongest_position\tdisc_in_window\tdisc_multi_mismatch"
            "\tdisc_all_positions"
        )
    lines = [header]
    for e in selected:
        c = e.candidate
        seed_tm_val = f"{e.seed.max_tm:.2f}" if e.seed is not None else ""
        hairpin_val = "" if e.hairpin_tm is None else f"{e.hairpin_tm:.2f}"
        run_val = "" if e.homopolymer_run is None else str(e.homopolymer_run)
        if e.offtarget is None:
            off_min, off_counts = "", ["", "", "", ""]
        else:
            m = e.offtarget.min_mismatch
            off_min = "" if m is None else str(m)
            off_counts = [str(n) for n in _combined_offtarget_counts(e.offtarget)]
        fields = [
            c.region,
            str(c.position + 1),
            str(c.end),
            c.guide_full,
            c.guide_overhang,
            c.passenger_full,
            c.sense_overhang,
            f"{c.gc_content * 100:.1f}",
            _rules_code(e, criteria),
            seed_tm_val,
            hairpin_val,
            run_val,
            off_min,
            *off_counts,
        ]
        if window_by_pos is not None:
            window = window_by_pos[c.position]
            spare_min = window.min_spare_mismatches
            fields += [
                str(window.max_hit_mismatches),
                "" if spare_min is None else str(spare_min),
                str(len(window.mismatch_assessments)),
                *_discrimination_fields(window),
            ]
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def render_design_result(
    result: DesignResult, limit: int, explain: bool = False
) -> None:
    """Print a design run: header, criteria, funnel, and the selected table.

    The shared renderer behind both 'hush design' and the shell's 'design', so
    the two stay identical. 'limit' caps the printed rows (0 = all). With
    'explain' and a discriminating request, the full per-allele breakdown is
    appended for the shown candidates.
    """
    mrna = result.target
    typer.echo(
        typer.style(f"Target: {mrna.id or '<unnamed>'}", fg=typer.colors.CYAN, bold=True)
    )
    typer.echo(
        f"  {len(mrna)} nt  |  {len(result.evaluations)} candidate "
        f"{result.length}-mer(s)  |  {result.overhang.value} 3' overhang"
    )
    if result.offtarget_skipped:
        typer.secho(
            "  [!] off-target scanning is on by default but no --offtarget-db was "
            "given; skipping it. Build an index with tools/build_offtarget_db.py, "
            "or pass --no-offtarget to silence this.",
            fg=typer.colors.YELLOW,
        )
    typer.echo("")
    _print_criteria(result.criteria, result.targeting, result.min_spare_mismatches)
    typer.echo("")
    _print_funnel(
        result.criteria, result.evaluations, result.selected, result.window_by_pos
    )
    typer.echo("")

    if not result.selected:
        hint = "--no-seed (strict on GC-rich targets), or --allow-homopolymers"
        if result.window_by_pos is not None:
            hint += " (or a targeting request the panel can satisfy)"
        typer.secho(
            f"No candidates selected. Try relaxing: {hint}.", fg=typer.colors.YELLOW
        )
    else:
        _print_table(
            result.selected,
            result.criteria,
            result.length,
            limit,
            result.window_by_pos,
            result.targeting,
        )
        if (
            explain
            and result.window_by_pos is not None
            and result.targeting is not None
            and result.targeting.discriminating
        ):
            render_discrimination_detail(
                result.selected, result.window_by_pos, limit
            )


@app.command("design")
def design_cmd(
    target: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Target mRNA/CDS file (FASTA or PDB).",
    ),
    length: int = typer.Option(
        19, "--length", "-l", help="siRNA core length (19-24)."
    ),
    overhang: str = typer.Option(
        "natural", "--overhang",
        help="3' overhang on each strand: natural (mRNA-derived) | dtdt | uu.",
    ),
    rules: str = typer.Option(
        "all", "--rules",
        help="Rule sets: comma list of ui-tei,reynolds,amarzguioui, or all/none.",
    ),
    mode: str = typer.Option(
        "any", "--mode",
        help="Combine enabled rules by 'any' (union, the default - each rule "
        "weighted equally) or 'all' (intersection).",
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
    offtarget: bool = typer.Option(
        True, "--offtarget/--no-offtarget",
        help="Apply the transcriptome homology off-target filter (on by default; "
        "needs a bowtie index via --offtarget-db, else it is skipped with a "
        "warning). Pass --no-offtarget to silence the reminder.",
    ),
    offtarget_db: Optional[Path] = typer.Option(
        None, "--offtarget-db",
        help="bowtie index prefix for the reference transcriptome "
        "(build with tools/build_offtarget_db.py).",
    ),
    offtarget_min_mismatch: int = typer.Option(
        DEFAULT_OFFTARGET_MIN_MISMATCH, "--offtarget-min-mismatch",
        help="Reject a candidate whose nearest off-target has fewer than this many "
        "mismatches.",
    ),
    offtarget_max_mismatch: int = typer.Option(
        DEFAULT_OFFTARGET_MAX_MISMATCH, "--offtarget-max-mismatch",
        help="Mismatch depth the search enumerates (bowtie -v allows 0-3).",
    ),
    offtarget_self_accession: Optional[str] = typer.Option(
        None, "--offtarget-self-accession",
        help="Accession of the target's own transcript to exclude from off-target "
        "hits (auto-detected from the FASTA header when omitted).",
    ),
    bowtie_path: str = typer.Option(
        "bowtie", "--bowtie-path", help="Path to the bowtie executable.",
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
    allow_homopolymers: bool = typer.Option(
        False, "--allow-homopolymers/--forbid-homopolymers",
        help="Whether to keep candidates with a single-nucleotide run over "
        "--max-run. Forbidden by default (poly-runs hurt siRNA stability and "
        "cause polymerase read-through); pass --allow-homopolymers to keep them.",
    ),
    max_run: int = typer.Option(
        DEFAULT_MAX_HOMOPOLYMER_RUN, "--max-run",
        help="Longest tolerated single-nucleotide run in the core (e.g. 3).",
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
    panel: Optional[Path] = typer.Option(
        None, "--panel", "-p",
        exists=True, dir_okay=False, readable=True,
        help="Allele-family FASTA (e.g. an HLA locus). Turns on the MSA/allele "
        "stage: candidates are gated to windows that hit the requested alleles "
        "and discriminate against the spared ones.",
        rich_help_panel=_ALLELE,
    ),
    panel_format: Optional[str] = typer.Option(
        None, "--panel-format",
        help="Override panel format: fasta | pdb | cif.",
        rich_help_panel=_ALLELE,
    ),
    pan_allele: bool = typer.Option(
        False, "--pan-allele",
        help="Hit every panel member and spare none (conserved-region targeting). "
        "The default when --panel is given without another targeting flag.",
        rich_help_panel=_ALLELE,
    ),
    allele_specific: Optional[str] = typer.Option(
        None, "--allele-specific",
        help="Hit exactly this panel member id and spare the rest "
        "(allele-specific targeting).",
        rich_help_panel=_ALLELE,
    ),
    target_specific: bool = typer.Option(
        False, "--target-specific",
        help="Spare every panel member and hit none of them: design specific to "
        "the loaded target (silenced by construction) against the whole reference "
        "panel. Use when the target is loaded apart from the panel.",
        rich_help_panel=_ALLELE,
    ),
    hit: Optional[str] = typer.Option(
        None, "--hit",
        help="Comma list of panel member ids to hit. With --spare the two sets "
        "are explicit; without it, the remaining members are spared.",
        rich_help_panel=_ALLELE,
    ),
    spare: Optional[str] = typer.Option(
        None, "--spare",
        help="Comma list of panel member ids to discriminate against "
        "(refines --hit).",
        rich_help_panel=_ALLELE,
    ),
    min_spare_mismatches: int = typer.Option(
        DEFAULT_MIN_SPARE_MISMATCHES, "--min-spare-mismatches",
        help="Differences a spared allele must carry inside a window before the "
        "window is credited with discriminating against it.",
        rich_help_panel=_ALLELE,
    ),
    sort: str = typer.Option(
        "auto", "--sort",
        help="Report order: auto (default - discrimination when the request "
        "spares something, else position), position (nothing ranked), or "
        "discrimination (best first, judged by each window's weakest spared "
        "allele; needs a request that spares something).",
        rich_help_panel=_ALLELE,
    ),
    explain: bool = typer.Option(
        False, "--explain",
        help="After the table, print the full per-allele discrimination breakdown "
        "(every spared allele, not just the weakest) for the shown candidates.",
        rich_help_panel=_ALLELE,
    ),
    limit: int = typer.Option(
        20, "--limit", "-n", help="Max candidates to print (0 = all)."
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Write all selected candidates to a TSV file."
    ),
) -> None:
    """Enumerate siRNA candidates for a target and filter them by chosen criteria.

    The best-practice default applies all three rule sets combined by ANY (union),
    the seed-duplex off-target Tm filter, the homopolymer gate, and transcriptome
    off-target scanning when a --offtarget-db is configured. The GC, hairpin and
    asymmetry gates are opt-in. Ordering is 'auto' by default: a request that
    spares alleles is ranked best-discrimination first, and everything else stays
    in target-position order ('--sort position' forces target order). Expert
    thermodynamic knobs live in the Advanced option group (see --help).

    Supply --panel to make the run allele-aware: the panel is aligned to the
    target (MSA/allele stage) and each candidate is additionally gated on whether
    its window hits the requested alleles and discriminates against the spared
    ones. The targeting intent is pan-allele by default, or --allele-specific,
    --target-specific (spare the whole panel), or --hit (+ --spare). This allele
    gate is ANDed with the other filters; add --explain for the full per-allele
    breakdown.
    """
    use_ut, use_rey, use_amz = _parse_rules(rules)
    rule_mode = _parse_mode(mode)
    thermo_params = _parse_params(params)
    overhang_policy = _parse_overhang(overhang)
    if panel is None and (
        pan_allele or allele_specific or target_specific or hit or spare
    ):
        raise typer.BadParameter(
            "allele targeting needs a panel; pass --panel with an allele-family "
            "FASTA.",
            param_hint="--panel",
        )

    # Turn the parsed flags into a framework-agnostic request, then let the
    # session wire the pipeline together and hand back a lossless result; the
    # printers below render it. The allele-targeting intent is only assembled
    # when a panel is given (the guard above rejects it otherwise).
    targeting_spec = None
    if panel is not None:
        targeting_spec = AlleleTargetingSpec(
            pan_allele=pan_allele,
            allele_specific=allele_specific,
            target_specific=target_specific,
            hit=tuple(split_ids(hit)) if hit is not None else None,
            spare=tuple(split_ids(spare)) if spare is not None else None,
        )
    design_params = DesignParams(
        length=length,
        overhang=overhang_policy,
        use_ui_tei=use_ut,
        use_reynolds=use_rey,
        use_amarzguioui=use_amz,
        rule_mode=rule_mode,
        reynolds_threshold=reynolds_threshold,
        amarzguioui_threshold=amarzguioui_threshold,
        use_seed_filter=seed,
        seed_tm_threshold=seed_tm,
        use_offtarget_filter=offtarget,
        offtarget_db=offtarget_db,
        offtarget_min_mismatch=offtarget_min_mismatch,
        offtarget_max_mismatch=offtarget_max_mismatch,
        offtarget_self_accession=offtarget_self_accession,
        bowtie_path=bowtie_path,
        use_gc_filter=gc,
        gc_min=gc_min / 100.0,
        gc_max=gc_max / 100.0,
        use_homopolymer_filter=not allow_homopolymers,
        max_homopolymer_run=max_run,
        use_asymmetry_filter=asymmetry,
        asymmetry_window=asym_window,
        asymmetry_temp_c=asym_temp,
        use_hairpin_filter=hairpin,
        hairpin_tm_max=hairpin_tm,
        thermo_params=thermo_params,
        targeting=targeting_spec,
        min_spare_mismatches=min_spare_mismatches,
        sort=sort,
    )

    session = HushSession()
    try:
        session.load_target(target)
        if panel is not None:
            session.load_panel(panel, fmt=panel_format)
        result = session.design(design_params)
    except (HushError, ValueError) as exc:
        typer.secho(f"Design error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    render_design_result(result, limit, explain=explain)

    if output is not None:
        _write_design_tsv(output, result.selected, result.criteria, result.window_by_pos)
        typer.echo(f"\nWrote {len(result.selected)} selected candidate(s) to {output}")


# --------------------------------------------------------------------------- #
# convert command (stage 4: siRNA -> shRNA)
# --------------------------------------------------------------------------- #


def _parse_loop(spec: str) -> LoopScaffold:
    try:
        return get_loop(spec)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--loop")


# Columns the converter needs from a 'hush design' TSV to rebuild the core.
_DESIGN_TSV_REQUIRED = {"region", "start", "end", "passenger_5to3", "passenger_overhang"}


def _read_design_tsv(path: Path) -> list[dict[str, str]]:
    """Read a 'hush design' TSV into a list of column-keyed row dicts."""
    lines = path.read_text(encoding="ascii").splitlines()
    if not lines:
        raise ValueError(f"{path} is empty.")
    header = lines[0].split("\t")
    missing = _DESIGN_TSV_REQUIRED - set(header)
    if missing:
        raise ValueError(
            f"{path} is missing column(s) {sorted(missing)}; expected a TSV from "
            "'hush design -o'."
        )
    idx = {name: i for i, name in enumerate(header)}
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        rows.append({name: fields[i] for name, i in idx.items()})
    return rows


def _core_from_row(row: dict[str, str]) -> tuple[BioSequence, BioSequence]:
    """Rebuild the (sense, antisense) duplex core from one design-TSV row.

    The passenger strand as ordered is the core plus its 3' overhang; stripping
    the recorded overhang leaves the paired core (the mRNA-matching sense
    strand), and the guide is its reverse-complement.
    """
    passenger_full = row["passenger_5to3"]
    overhang = row["passenger_overhang"]
    core_str = passenger_full[: len(passenger_full) - len(overhang)] if overhang else passenger_full
    sense = BioSequence(core_str, moltype=MolType.RNA)
    return sense, sense.reverse_complement()


def _print_convert_table(constructs: list, options: ShrnaOptions, limit: int) -> None:
    hp_w = max(len(c.hairpin_rna.seq) for c in constructs)
    header = (
        f"  {'region':<11} {'hairpin 5->3':<{hp_w}} {'fold':>5} {'5AG':>3} "
        f"{'termOK':>6} {'run':>4}"
    )
    typer.echo(typer.style(header, bold=True))
    shown = constructs if limit == 0 else constructs[:limit]
    for c in shown:
        fold = "ok" if c.folds_as_designed else "RISK"
        row = (
            f"  {c.region:<11} {c.hairpin_rna.seq:<{hp_w}} {fold:>5} "
            f"{(c.hairpin_rna.seq[:1] if c.starts_with_purine else '-'):>3} "
            f"{('yes' if c.premature_terminator_free else 'NO'):>6} "
            f"{c.max_homopolymer_run:>4}"
        )
        typer.echo(row)
    typer.secho(
        "  fold: 'ok' = folds into the intended stem-loop; 'RISK' = a competing "
        "structure is more stable.",
        fg=typer.colors.BRIGHT_BLACK,
    )
    typer.secho(
        "  Order the top/bottom oligos from the TSV (-o); anneal them to form the "
        "insert, then add your vector's sticky-end tails.",
        fg=typer.colors.BRIGHT_BLACK,
    )
    if limit and len(constructs) > limit:
        typer.echo(
            f"  ... {len(constructs) - limit} more "
            "(use --limit 0 to show all, or -o file.tsv)"
        )


def _write_convert_tsv(path: Path, constructs: list) -> None:
    header = (
        "region\tstart\tend\tloop_name\tloop_seq\torientation\thairpin_rna_5to3\t"
        "insert_top_dna_5to3\tinsert_bottom_dna_5to3\tfolds_as_designed\t"
        "fold_dg37\tfold_tm_c\tstarts_with_purine\tpremature_terminator_free\t"
        "max_homopolymer_run"
    )
    lines = [header]
    for c in constructs:
        dg37 = "" if c.fold is None else f"{c.fold.dg37:.2f}"
        tm = "" if c.fold is None else f"{c.fold.tm:.2f}"
        fields = [
            c.region,
            str(c.start + 1),
            str(c.end),
            c.loop_name,
            LOOPS[c.loop_name].seq,
            c.orientation.value,
            c.hairpin_rna.seq,
            c.insert_top.seq,
            c.insert_bottom.seq,
            str(c.folds_as_designed),
            dg37,
            tm,
            str(c.starts_with_purine),
            str(c.premature_terminator_free),
            str(c.max_homopolymer_run),
        ]
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def render_convert_result(result: ConvertResult, limit: int) -> None:
    """Print a convert run: scaffold header, summary, and the construct table.

    The shared renderer behind both 'hush convert' and the shell's 'convert'.
    """
    options = result.options
    loop = options.loop
    term = f"{options.terminator_length}x T" if options.include_terminator else "none"
    typer.echo(
        typer.style(
            f"Scaffold: stem-loop (SA)  |  loop {loop.name} {loop.seq}  |  "
            f"terminator {term}",
            fg=typer.colors.CYAN,
            bold=True,
        )
    )
    typer.echo(f"  {loop.citation}")
    typer.echo("")

    constructs = result.constructs
    if not constructs:
        typer.secho(
            "No constructs to show - nothing was converted.", fg=typer.colors.YELLOW
        )
        return

    fold_ok = sum(1 for c in constructs if c.folds_as_designed)
    purine_start = sum(1 for c in constructs if c.starts_with_purine)
    typer.echo(typer.style("Summary", bold=True))
    typer.echo(f"  constructs:        {len(constructs)}")
    typer.echo(f"  fold as designed:  {fold_ok}")
    typer.echo(f"  5'-A/G (U6 Pol III): {purine_start}")
    typer.echo("")
    _print_convert_table(constructs, options, limit)


@app.command("convert")
def convert_cmd(
    selected: Path = typer.Argument(
        ..., exists=True, dir_okay=False, readable=True,
        help="Selected-candidate TSV from 'hush design -o'.",
    ),
    loop: str = typer.Option(
        DEFAULT_LOOP, "--loop",
        help="Loop scaffold: jensen (UGUGCUU, default) | psuper (UUCAAGAGA) | "
        "trc (CUCGAG).",
    ),
    terminator: bool = typer.Option(
        True, "--terminator/--no-terminator",
        help="Append the Pol III poly-T terminator to the hairpin insert.",
    ),
    terminator_length: int = typer.Option(
        DEFAULT_TERMINATOR_LENGTH, "--terminator-length",
        help="Length of the poly-T terminator run (>= 4 when included).",
    ),
    limit: int = typer.Option(
        20, "--limit", "-n", help="Max constructs to print (0 = all)."
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Write all constructs to a TSV file."
    ),
) -> None:
    """Convert selected siRNA candidates into stem-loop shRNA constructs.

    Reads a 'hush design -o' TSV and, for each selected candidate, grafts its
    19-bp duplex core (unchanged) onto the chosen loop as a sense-loop-antisense
    hairpin. The bare insert is emitted as top/bottom DNA oligos plus a Pol III
    poly-T terminator - vector-agnostic, so add your vector's restriction
    sticky-end tails before cloning. Each construct is folded to confirm it forms
    the intended stem-loop, and QC flags (Pol III 5'-A/G start, loop
    terminator-safety, homopolymer runs) are reported.
    """
    loop_scaffold = _parse_loop(loop)
    try:
        options = ShrnaOptions(
            loop=loop_scaffold,
            include_terminator=terminator,
            terminator_length=terminator_length,
        )
        rows = _read_design_tsv(selected)
        constructs = []
        for row in rows:
            sense, antisense = _core_from_row(row)
            constructs.append(
                convert(
                    sense,
                    antisense,
                    region=row["region"],
                    start=int(row["start"]) - 1,
                    end=int(row["end"]),
                    options=options,
                )
            )
    except (HushError, ValueError) as exc:
        typer.secho(f"Convert error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    render_convert_result(ConvertResult(constructs=constructs, options=options), limit)

    if output is not None:
        _write_convert_tsv(output, constructs)
        typer.echo(f"\nWrote {len(constructs)} construct(s) to {output}")


@app.command("shell")
def shell_cmd() -> None:
    """Start an interactive session: load a target and panel once, then design.

    The one-shot commands re-parse and re-load on every call; the shell keeps a
    live instance so you load the target (and panel) once, tune the criteria with
    'set', and re-run 'design'/'convert' as often as you like - the panel
    alignment is computed once and reused. Type 'help' inside for the commands.
    """
    # Imported here so the command definitions above load without pulling in the
    # shell (which imports this module for its renderers).
    from app.shell import HushShell

    HushShell().run()


if __name__ == "__main__":
    app()
