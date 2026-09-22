"""Build API response models from service-layer results.

The counterpart to the CLI's printers: where app.cli renders a DesignResult to a
terminal table, this renders it to api.schemas models (the same columns, as
JSON). Kept separate from the schemas so the wire shape and the mapping evolve
independently.
"""

from __future__ import annotations

from core.selection import CandidateEvaluation, SelectionCriteria
from core.sequence import BioSequence, DesignInput, SequencePanel
from msa.targeting import AlleleTargeting
from msa.windows import WindowEvaluation
from service.results import ConvertResult, DesignResult
from shrna.convert import ShrnaConstruct
from shrna.loops import LOOPS

from api import schemas

# The nucleotide-target requirement, surfaced by inspect for a protein target.
_PROTEIN_NOTE = (
    "siRNA/shRNA design needs the coding transcript (mRNA/CDS); resolve the "
    "protein to its transcript first."
)


# --------------------------------------------------------------------------- #
# sequences
# --------------------------------------------------------------------------- #


def sequence_summary(seq: BioSequence) -> schemas.SequenceSummary:
    return schemas.SequenceSummary(
        id=seq.id,
        description=seq.description,
        length=seq.length,
        moltype=seq.moltype.value,
        is_nucleotide=seq.is_nucleotide,
        gc_percent=seq.gc_content * 100 if seq.is_nucleotide and seq.length else None,
    )


def panel_summary(panel: SequencePanel) -> schemas.PanelSummary:
    return schemas.PanelSummary(
        name=panel.name,
        count=len(panel),
        uniform=panel.is_uniform,
        members=[sequence_summary(s) for s in panel],
    )


def inspect_response(design_input: DesignInput) -> schemas.InspectResponse:
    ready = design_input.target.is_nucleotide
    return schemas.InspectResponse(
        target=sequence_summary(design_input.target),
        target_ready=ready,
        note=None if ready else _PROTEIN_NOTE,
        panel=panel_summary(design_input.panel) if design_input.has_panel else None,
    )


# --------------------------------------------------------------------------- #
# design
# --------------------------------------------------------------------------- #


def _rules_code(e: CandidateEvaluation, criteria: SelectionCriteria) -> str:
    """Compact per-rule pass code, e.g. 'UR-' (Ui-Tei + Reynolds pass, Amarz fails)."""
    parts = []
    if e.ui_tei is not None:
        parts.append("U" if e.ui_tei.satisfied else "-")
    if e.reynolds is not None:
        parts.append("R" if e.reynolds.meets_threshold(criteria.reynolds_threshold) else "-")
    if e.amarzguioui is not None:
        parts.append(
            "A" if e.amarzguioui.meets_threshold(criteria.amarzguioui_threshold) else "-"
        )
    return "".join(parts) or "-"


def _offtarget_hits(profile) -> list[int]:
    """Guide + passenger off-target hit counts at mismatch levels 0..3."""
    out = [0, 0, 0, 0]
    for counts in (profile.guide_counts, profile.passenger_counts):
        for i, n in enumerate(counts):
            if i < len(out):
                out[i] += n
    return out


def _discrimination(window: WindowEvaluation) -> schemas.DiscriminationResult | None:
    deciding = window.deciding_mismatch
    if deciding is None:
        return None
    return schemas.DiscriminationResult(
        guide_position=deciding.guide_position,
        mutation=deciding.mutation_type,
        mismatch_class=None if deciding.mismatch_class is None else deciding.mismatch_class.value,
        huang_recommended=deciding.huang_recommended,
        region=deciding.region.value,
        sources_agree=deciding.sources_agree,
    )


def criteria_summary(criteria: SelectionCriteria) -> schemas.CriteriaSummary:
    return schemas.CriteriaSummary(
        rules=[n.replace("_", "-") for n in criteria.enabled_rule_names],
        rule_mode=criteria.rule_mode.value,
        seed_filter=criteria.use_seed_filter,
        seed_tm_threshold_c=criteria.seed_tm_threshold if criteria.use_seed_filter else None,
        offtarget_filter=criteria.use_offtarget_filter,
        offtarget_min_mismatch=(
            criteria.offtarget_min_mismatch if criteria.use_offtarget_filter else None
        ),
        offtarget_self_accession=criteria.offtarget_self_accession,
        gc_filter=criteria.use_gc_filter,
        gc_min=criteria.gc_min if criteria.use_gc_filter else None,
        gc_max=criteria.gc_max if criteria.use_gc_filter else None,
        hairpin_filter=criteria.use_hairpin_filter,
        hairpin_tm_max_c=criteria.hairpin_tm_max if criteria.use_hairpin_filter else None,
        homopolymer_filter=criteria.use_homopolymer_filter,
        max_homopolymer_run=(
            criteria.max_homopolymer_run if criteria.use_homopolymer_filter else None
        ),
        asymmetry_filter=criteria.use_asymmetry_filter,
    )


def targeting_summary(targeting: AlleleTargeting) -> schemas.TargetingSummary:
    return schemas.TargetingSummary(
        mode=targeting.mode.value,
        hit=sorted(targeting.hit),
        spare=sorted(targeting.spare),
        discriminating=targeting.discriminating,
    )


def _funnel(result: DesignResult) -> dict[str, int]:
    c = result.criteria
    evals = result.evaluations
    funnel = {"rules": sum(1 for e in evals if e.rule_pass)}
    if result.window_by_pos is not None:
        funnel["allele"] = sum(
            1 for e in evals if result.window_by_pos[e.candidate.position].satisfies
        )
    if c.use_gc_filter:
        funnel["gc"] = sum(1 for e in evals if e.gc_pass)
    if c.use_seed_filter:
        funnel["seed"] = sum(1 for e in evals if e.seed_pass)
    if c.use_offtarget_filter:
        funnel["offtarget"] = sum(1 for e in evals if e.offtarget_pass)
    if c.use_hairpin_filter:
        funnel["hairpin"] = sum(1 for e in evals if e.hairpin_pass)
    if c.use_homopolymer_filter:
        funnel["homopolymer"] = sum(1 for e in evals if e.homopolymer_pass)
    if c.use_asymmetry_filter:
        funnel["asymmetry"] = sum(1 for e in evals if e.asymmetry_pass)
    funnel["selected"] = len(result.selected)
    return funnel


def candidate_result(
    e: CandidateEvaluation,
    criteria: SelectionCriteria,
    window: WindowEvaluation | None,
) -> schemas.CandidateResult:
    c = e.candidate
    allele = None
    if window is not None:
        allele = schemas.AlleleResult(
            hit_max_mismatch=window.max_hit_mismatches,
            spare_min_mismatch=window.min_spare_mismatches,
            discrimination=_discrimination(window),
        )
    return schemas.CandidateResult(
        region=c.region,
        start=c.position + 1,
        end=c.end,
        guide_5to3=c.guide_full,
        guide_overhang=c.guide_overhang,
        passenger_5to3=c.passenger_full,
        passenger_overhang=c.sense_overhang,
        gc_percent=c.gc_content * 100,
        rules_passed=_rules_code(e, criteria),
        seed_max_tm_c=e.seed.max_tm if e.seed is not None else None,
        offtarget_min_mismatch=e.offtarget.min_mismatch if e.offtarget is not None else None,
        offtarget_hits=_offtarget_hits(e.offtarget) if e.offtarget is not None else None,
        hairpin_tm_c=e.hairpin_tm,
        homopolymer_max_run=e.homopolymer_run,
        allele=allele,
    )


def design_response(result: DesignResult) -> schemas.DesignResponse:
    window_by_pos = result.window_by_pos
    selected = [
        candidate_result(
            e,
            result.criteria,
            None if window_by_pos is None else window_by_pos.get(e.candidate.position),
        )
        for e in result.selected
    ]
    return schemas.DesignResponse(
        target=sequence_summary(result.target),
        length=result.length,
        overhang=result.overhang.value,
        candidate_count=len(result.evaluations),
        offtarget_skipped=result.offtarget_skipped,
        criteria=criteria_summary(result.criteria),
        funnel=_funnel(result),
        targeting=targeting_summary(result.targeting) if result.targeting else None,
        selected=selected,
    )


# --------------------------------------------------------------------------- #
# convert
# --------------------------------------------------------------------------- #


def construct_result(c: ShrnaConstruct) -> schemas.ConstructResult:
    return schemas.ConstructResult(
        region=c.region,
        start=c.start + 1,
        end=c.end,
        loop_name=c.loop_name,
        loop_seq=LOOPS[c.loop_name].seq,
        orientation=c.orientation.value,
        hairpin_rna_5to3=c.hairpin_rna.seq,
        insert_top_dna_5to3=c.insert_top.seq,
        insert_bottom_dna_5to3=c.insert_bottom.seq,
        folds_as_designed=c.folds_as_designed,
        fold_dg37=None if c.fold is None else c.fold.dg37,
        fold_tm_c=None if c.fold is None else c.fold.tm,
        starts_with_g=c.starts_with_g,
        premature_terminator_free=c.premature_terminator_free,
        max_homopolymer_run=c.max_homopolymer_run,
    )


def convert_response(result: ConvertResult) -> schemas.ConvertResponse:
    options = result.options
    constructs = result.constructs
    return schemas.ConvertResponse(
        loop_name=options.loop.name,
        loop_seq=options.loop.seq,
        terminator_length=(
            options.terminator_length if options.include_terminator else None
        ),
        construct_count=len(constructs),
        fold_ok_count=sum(1 for c in constructs if c.folds_as_designed),
        g_start_count=sum(1 for c in constructs if c.starts_with_g),
        constructs=[construct_result(c) for c in constructs],
    )
