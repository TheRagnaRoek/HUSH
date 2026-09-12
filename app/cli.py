"""HUSH command-line interface.

CLI-first entry point for the pipeline. Exposes the Input Layer via
'hush inspect' and siRNA candidate design via 'hush design'. The MSA/allele
stage is engaged from within 'hush design' by supplying '--panel': the panel is
aligned to the target and each candidate is additionally gated on its allele
targeting. The shRNA conversion command is added as that stage comes online.

Run 'hush --help' (after 'pip install -e .') or 'python -m app.cli --help'.
"""

from __future__ import annotations

import re
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
    BowtieSearcher,
)
from msa.align_star import align_star
from msa.discrimination import MismatchClass
from msa.targeting import AlleleTargeting
from msa.variation import PanelVariation
from msa.windows import (
    DEFAULT_MIN_SPARE_MISMATCHES,
    WindowCriteria,
    WindowEvaluation,
)
from msa.windows import evaluate_windows
from core.rules import AMARZGUIOUI_DEFAULT_THRESHOLD, REYNOLDS_DEFAULT_THRESHOLD
from core.seed import SEED_TM_THRESHOLD_C
from core.selection import (
    DEFAULT_MAX_HOMOPOLYMER_RUN,
    RuleMode,
    SelectionCriteria,
    evaluate_all,
)
from core.sequence import BioSequence, DesignInput
from core.sirna import OverhangPolicy, enumerate_candidates
from core.thermo import REFERENCE_TEMPERATURE_C, ThermoParams

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


_OVERHANG_ALIASES = {
    "natural": OverhangPolicy.NATURAL,
    "mrna": OverhangPolicy.NATURAL,
    "dtdt": OverhangPolicy.DTDT,
    "tt": OverhangPolicy.DTDT,
    "uu": OverhangPolicy.UU,
}


def _parse_overhang(spec: str) -> OverhangPolicy:
    key = _OVERHANG_ALIASES.get(spec.strip().lower())
    if key is None:
        raise typer.BadParameter(
            f"unknown overhang {spec!r}; choose 'natural', 'dtdt', or 'uu'.",
            param_hint="--overhang",
        )
    return key


def _split_ids(spec: str) -> list[str]:
    """Split a comma-separated id list, dropping blanks and surrounding space."""
    return [token.strip() for token in spec.split(",") if token.strip()]


def _build_targeting(
    panel,
    *,
    pan_allele: bool,
    allele_specific: str | None,
    hit: str | None,
    spare: str | None,
) -> AlleleTargeting:
    """Turn the allele-targeting flags into one AlleleTargeting request.

    The four flags name three mutually exclusive intents (pan-allele, allele-
    specific, explicit hit/spare); pan-allele is the default when a panel is
    given without another flag. --spare only refines an explicit --hit.
    """
    intents = [
        name
        for name, chosen in (
            ("--pan-allele", pan_allele),
            ("--allele-specific", allele_specific is not None),
            ("--hit", hit is not None),
        )
        if chosen
    ]
    if len(intents) > 1:
        raise typer.BadParameter(
            f"{' and '.join(intents)} are mutually exclusive; choose one "
            "targeting intent.",
            param_hint="--pan-allele/--allele-specific/--hit",
        )
    if spare is not None and hit is None:
        raise typer.BadParameter(
            "--spare only refines an explicit --hit; name the alleles to hit too.",
            param_hint="--spare",
        )

    if allele_specific is not None:
        return AlleleTargeting.allele_specific(allele_specific, panel)
    if hit is not None:
        hits = _split_ids(hit)
        if spare is not None:
            return AlleleTargeting(hit=hits, spare=_split_ids(spare))
        return AlleleTargeting.subset(hits, panel)
    # No explicit intent: broad, conserved-region targeting across the family.
    return AlleleTargeting.pan_allele(panel)


# RefSeq accession, e.g. NM_005018.3 / NR_ / XM_ / XR_, with optional version.
_REFSEQ_ACCESSION = re.compile(r"\b([NX][MR]_\d+(?:\.\d+)?)\b")


def _detect_accession(seq: BioSequence) -> str | None:
    """Best-effort RefSeq accession from a target's id/description, for self-exclusion.

    The off-target search must drop hits to the candidate's own transcript. The
    accession usually sits in the FASTA header (e.g. 'NM_005018.3' in the PDCD1
    fixture); pass --offtarget-self-accession to override when it does not.
    """
    match = _REFSEQ_ACCESSION.search(f"{seq.id} {seq.description}")
    return match.group(1) if match else None


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
        if targeting.discriminating:
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


def _discrimination_cell(window: WindowEvaluation) -> str:
    """Compact readout of the mismatch a window's discrimination rests on.

    'P16*pp' reads: guide position 16, listed in Huang's table for this base
    change (*), purine:purine (pp). A missing mark means the opposite; '?' means
    the bases were ambiguous and neither source applies.
    """
    deciding = window.deciding_mismatch
    if deciding is None:
        return "none"
    if deciding.huang_recommended is None:
        listed = "?"
    else:
        listed = "*" if deciding.huang_recommended else ""
    purine = "pp" if deciding.mismatch_class is MismatchClass.PURINE_PURINE else ""
    return f"P{deciding.guide_position}{listed}{purine}"


def _apply_sort(
    selected: list,
    sort: str,
    window_by_pos: dict[int, WindowEvaluation] | None,
    targeting: AlleleTargeting | None,
) -> list:
    """Order the selected candidates as the user asked.

    Position order is the default and the project's convention: nothing is ranked
    unless asked for. Discrimination order is opt-in because it imposes a
    priority between published findings that the measurements do not fully
    settle - see msa.discrimination.DISCRIMINATION_ORDER_NOTE.
    """
    if sort == "position":
        return selected
    if sort != "discrimination":
        raise typer.BadParameter(
            f"unknown sort {sort!r}; choose 'position' or 'discrimination'.",
            param_hint="--sort",
        )
    if window_by_pos is None or targeting is None or not targeting.discriminating:
        raise typer.BadParameter(
            "--sort discrimination needs a request that spares something; "
            "pan-allele targeting asks for no discrimination to rank.",
            param_hint="--sort",
        )
    # Position breaks ties so two runs never disagree.
    return sorted(
        selected,
        key=lambda e: (
            window_by_pos[e.candidate.position].discrimination_key,
            -e.candidate.position,
        ),
        reverse=True,
    )


def _discrimination_fields(window: WindowEvaluation) -> list[str]:
    """The deciding mismatch, as six TSV fields. Empty ones where it has none."""
    deciding = window.deciding_mismatch
    if deciding is None:
        return [""] * 6
    huang = "" if deciding.huang_recommended is None else str(deciding.huang_recommended)
    return [
        str(deciding.guide_position),
        deciding.mutation_type,
        "" if deciding.mismatch_class is None else deciding.mismatch_class.value,
        huang,
        deciding.region.value,
        str(deciding.sources_agree),
    ]


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
        header += f" {'minDiff':>7} {'disc':>7}"
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
            row += f" {_discrimination_cell(window):>7}"
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
            "  minDiff and disc describe the weakest spared allele - the one this "
            "window distinguishes least well.",
            fg=typer.colors.BRIGHT_BLACK,
        )
        typer.secho(
            "  disc: guide position, * = listed by Huang et al. for this base "
            "change, pp = purine:purine.",
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
            # The single mismatch the window's discrimination rank rests on:
            # the best-placed difference from the weakest spared allele.
            "\tdisc_guide_position\tdisc_mutation\tdisc_mismatch_class"
            "\tdisc_huang_recommended\tdisc_region\tdisc_sources_agree"
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
                *_discrimination_fields(window),
            ]
        lines.append("\t".join(fields))
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
        "position", "--sort",
        help="Report order: position (default, nothing ranked) or discrimination "
        "(best first, judged by each window's weakest spared allele; needs a "
        "request that spares something).",
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
    asymmetry gates are opt-in. Nothing is ranked by default: selected candidates
    are reported in target-position order, and '--sort discrimination' is the one
    way to ask for an ordering. Expert thermodynamic knobs live in the Advanced
    option group (see --help).

    Supply --panel to make the run allele-aware: the panel is aligned to the
    target (MSA/allele stage) and each candidate is additionally gated on whether
    its window hits the requested alleles and discriminates against the spared
    ones. The targeting intent is pan-allele by default, or --allele-specific /
    --hit (+ --spare). This allele gate is ANDed with the other filters.
    """
    use_ut, use_rey, use_amz = _parse_rules(rules)
    rule_mode = _parse_mode(mode)
    thermo_params = _parse_params(params)
    overhang_policy = _parse_overhang(overhang)
    if panel is None and (pan_allele or allele_specific or hit or spare):
        raise typer.BadParameter(
            "allele targeting needs a panel; pass --panel with an allele-family "
            "FASTA.",
            param_hint="--panel",
        )
    try:
        design_input = build_design_input(target, panel, panel_fmt=panel_format)
        mrna = design_input.require_nucleotide_target()
        candidates = list(
            enumerate_candidates(mrna, length=length, overhang=overhang_policy)
        )

        # Allele gate (MSA / stage 3): align the panel to the target and evaluate
        # each candidate window's targeting. Windows and candidates share the
        # target's coordinates one-for-one, so a candidate is looked up by its
        # start position.
        targeting = None
        window_by_pos = None
        if design_input.has_panel:
            targeting = _build_targeting(
                design_input.panel,
                pan_allele=pan_allele,
                allele_specific=allele_specific,
                hit=hit,
                spare=spare,
            )
            variation = PanelVariation.from_alignment(
                align_star(mrna, design_input.panel)
            )
            window_criteria = WindowCriteria(
                targeting=targeting,
                length=length,
                min_spare_mismatches=min_spare_mismatches,
            )
            window_by_pos = {
                ev.position: ev
                for ev in evaluate_windows(variation, window_criteria)
            }

        # Off-target scanning is on by default, but it needs a reference index.
        # With none configured we skip it (a warning is printed below) rather than
        # error, so a bare run still works; --no-offtarget silences the reminder.
        offtarget_searcher = None
        self_accession = None
        offtarget_active = offtarget and offtarget_db is not None
        offtarget_skipped = offtarget and offtarget_db is None
        if offtarget_active:
            self_accession = offtarget_self_accession or _detect_accession(mrna)
            # One batched bowtie run over every strand of every candidate, then a
            # cheap dict lookup per candidate inside the selection loop.
            backend = BowtieSearcher(
                offtarget_db,
                max_mismatch=offtarget_max_mismatch,
                bowtie_path=bowtie_path,
            )
            reads = [c.sense.seq for c in candidates]
            reads += [c.antisense.seq for c in candidates]
            offtarget_searcher = backend.precompute(reads)

        criteria = SelectionCriteria(
            use_ui_tei=use_ut,
            use_reynolds=use_rey,
            use_amarzguioui=use_amz,
            rule_mode=rule_mode,
            reynolds_threshold=reynolds_threshold,
            amarzguioui_threshold=amarzguioui_threshold,
            use_seed_filter=seed,
            seed_tm_threshold=seed_tm,
            use_offtarget_filter=offtarget_active,
            offtarget_min_mismatch=offtarget_min_mismatch,
            offtarget_searcher=offtarget_searcher,
            offtarget_self_accession=self_accession,
            use_gc_filter=gc,
            gc_min=gc_min / 100.0,
            gc_max=gc_max / 100.0,
            use_hairpin_filter=hairpin,
            hairpin_tm_max=hairpin_tm,
            use_homopolymer_filter=not allow_homopolymers,
            max_homopolymer_run=max_run,
            use_asymmetry_filter=asymmetry,
            asymmetry_window=asym_window,
            asymmetry_temp_c=asym_temp,
            thermo_params=thermo_params,
        )
        evaluations = list(evaluate_all(candidates, criteria))
    except (HushError, ValueError) as exc:
        typer.secho(f"Design error: {exc}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1)

    def _allele_ok(evaluation) -> bool:
        return (
            window_by_pos is None
            or window_by_pos[evaluation.candidate.position].satisfies
        )

    selected = [e for e in evaluations if e.selected and _allele_ok(e)]
    selected = _apply_sort(selected, sort, window_by_pos, targeting)

    typer.echo(typer.style(f"Target: {mrna.id or '<unnamed>'}", fg=typer.colors.CYAN, bold=True))
    typer.echo(
        f"  {len(mrna)} nt  |  {len(evaluations)} candidate {length}-mer(s)  "
        f"|  {overhang_policy.value} 3' overhang"
    )
    if offtarget_skipped:
        typer.secho(
            "  [!] off-target scanning is on by default but no --offtarget-db was "
            "given; skipping it. Build an index with tools/build_offtarget_db.py, "
            "or pass --no-offtarget to silence this.",
            fg=typer.colors.YELLOW,
        )
    typer.echo("")
    _print_criteria(criteria, targeting, min_spare_mismatches)
    typer.echo("")
    _print_funnel(criteria, evaluations, selected, window_by_pos)
    typer.echo("")

    if not selected:
        hint = "--no-seed (strict on GC-rich targets), or --allow-homopolymers"
        if window_by_pos is not None:
            hint += " (or a targeting request the panel can satisfy)"
        typer.secho(f"No candidates selected. Try relaxing: {hint}.", fg=typer.colors.YELLOW)
    else:
        _print_table(selected, criteria, length, limit, window_by_pos, targeting)

    if output is not None:
        _write_design_tsv(output, selected, criteria, window_by_pos)
        typer.echo(f"\nWrote {len(selected)} selected candidate(s) to {output}")


if __name__ == "__main__":
    app()
