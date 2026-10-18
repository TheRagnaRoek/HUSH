"""Pydantic request and response models: the JSON contract and OpenAPI schema.

Request models convert to service.params with to_domain(); response models are
built by api.serialize. Enums travel as plain strings. See docs/api.md.
"""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field

from core.asymmetry import DEFAULT_WINDOW_BP
from core.hairpin import DEFAULT_INTERNAL_TM_THRESHOLD_C
from core.offtarget import (
    DEFAULT_OFFTARGET_MAX_MISMATCH,
    DEFAULT_OFFTARGET_MIN_MISMATCH,
    DEFAULT_SELF_MIN_COVERAGE,
)
from core.rules import AMARZGUIOUI_DEFAULT_THRESHOLD, REYNOLDS_DEFAULT_THRESHOLD
from core.seed import SEED_TM_THRESHOLD_C
from core.selection import DEFAULT_MAX_HOMOPOLYMER_RUN, RuleMode
from core.sirna import OverhangPolicy
from core.thermo import REFERENCE_TEMPERATURE_C, ThermoParams
from msa.windows import DEFAULT_MIN_SPARE_MISMATCHES
from service.params import AlleleTargetingSpec, DesignParams
from shrna.convert import (
    DEFAULT_TERMINATOR_LENGTH,
    ShrnaOptions,
    StemOverhang,
    get_stem_overhang,
)
from shrna.loops import DEFAULT_LOOP, get_loop


def _to_enum(cls, value: str, what: str):
    try:
        return cls(value)
    except ValueError:
        allowed = ", ".join(m.value for m in cls)
        raise ValueError(f"unknown {what} {value!r}; choose from {allowed}.")


class TargetingSchema(BaseModel):
    """A targeting intent by panel member id."""

    pan_allele: bool = False
    allele_specific: Optional[str] = None
    target_specific: bool = False
    hit: Optional[list[str]] = None
    spare: Optional[list[str]] = None

    def to_domain(self) -> AlleleTargetingSpec:
        return AlleleTargetingSpec(
            pan_allele=self.pan_allele,
            allele_specific=self.allele_specific,
            target_specific=self.target_specific,
            hit=tuple(self.hit) if self.hit is not None else None,
            spare=tuple(self.spare) if self.spare is not None else None,
        )


class DesignParamsSchema(BaseModel):
    """Design parameters; a bare {} gives the defaults. GC bounds are fractions."""

    length: int = 19
    overhang: str = OverhangPolicy.NATURAL.value

    use_ui_tei: bool = True
    use_reynolds: bool = True
    use_amarzguioui: bool = True
    rule_mode: str = RuleMode.ANY.value
    reynolds_threshold: int = REYNOLDS_DEFAULT_THRESHOLD
    amarzguioui_threshold: int = AMARZGUIOUI_DEFAULT_THRESHOLD

    use_seed_filter: bool = True
    seed_tm_threshold: float = SEED_TM_THRESHOLD_C

    use_offtarget_filter: bool = True
    offtarget_db: Optional[str] = None
    offtarget_min_mismatch: int = DEFAULT_OFFTARGET_MIN_MISMATCH
    offtarget_max_mismatch: int = DEFAULT_OFFTARGET_MAX_MISMATCH
    offtarget_self_accessions: list[str] = Field(default_factory=list)
    resolve_offtarget_self: bool = True
    offtarget_self_coverage: float = DEFAULT_SELF_MIN_COVERAGE
    # bowtie_path and offtarget_self_fasta are absent on purpose; see docs/api.md.

    use_gc_filter: bool = False
    gc_min: float = 0.30
    gc_max: float = 0.52

    use_homopolymer_filter: bool = True
    max_homopolymer_run: int = DEFAULT_MAX_HOMOPOLYMER_RUN

    use_asymmetry_filter: bool = False
    asymmetry_window: int = DEFAULT_WINDOW_BP
    asymmetry_temp_c: float = REFERENCE_TEMPERATURE_C

    use_hairpin_filter: bool = False
    hairpin_tm_max: float = DEFAULT_INTERNAL_TM_THRESHOLD_C

    thermo_params: str = ThermoParams.XIA_1998.value

    targeting: Optional[TargetingSchema] = None
    min_spare_mismatches: int = DEFAULT_MIN_SPARE_MISMATCHES
    sort: str = "auto"

    def to_domain(self) -> DesignParams:
        from pathlib import Path

        return DesignParams(
            length=self.length,
            overhang=_to_enum(OverhangPolicy, self.overhang, "overhang"),
            use_ui_tei=self.use_ui_tei,
            use_reynolds=self.use_reynolds,
            use_amarzguioui=self.use_amarzguioui,
            rule_mode=_to_enum(RuleMode, self.rule_mode, "mode"),
            reynolds_threshold=self.reynolds_threshold,
            amarzguioui_threshold=self.amarzguioui_threshold,
            use_seed_filter=self.use_seed_filter,
            seed_tm_threshold=self.seed_tm_threshold,
            use_offtarget_filter=self.use_offtarget_filter,
            offtarget_db=Path(self.offtarget_db) if self.offtarget_db else None,
            offtarget_min_mismatch=self.offtarget_min_mismatch,
            offtarget_max_mismatch=self.offtarget_max_mismatch,
            offtarget_self_accessions=tuple(self.offtarget_self_accessions),
            resolve_offtarget_self=self.resolve_offtarget_self,
            offtarget_self_coverage=self.offtarget_self_coverage,
            use_gc_filter=self.use_gc_filter,
            gc_min=self.gc_min,
            gc_max=self.gc_max,
            use_homopolymer_filter=self.use_homopolymer_filter,
            max_homopolymer_run=self.max_homopolymer_run,
            use_asymmetry_filter=self.use_asymmetry_filter,
            asymmetry_window=self.asymmetry_window,
            asymmetry_temp_c=self.asymmetry_temp_c,
            use_hairpin_filter=self.use_hairpin_filter,
            hairpin_tm_max=self.hairpin_tm_max,
            thermo_params=_to_enum(ThermoParams, self.thermo_params, "parameter set"),
            targeting=self.targeting.to_domain() if self.targeting else None,
            min_spare_mismatches=self.min_spare_mismatches,
            sort=self.sort,
        )


class ConvertOptionsSchema(BaseModel):
    """shRNA conversion options."""

    loop: str = DEFAULT_LOOP
    stem_overhang: str = StemOverhang.KEEP.value
    include_terminator: bool = True
    terminator_length: int = DEFAULT_TERMINATOR_LENGTH
    drop_premature_terminators: bool = True

    def to_domain(self) -> ShrnaOptions:
        return ShrnaOptions(
            loop=get_loop(self.loop),
            stem_overhang=get_stem_overhang(self.stem_overhang),
            include_terminator=self.include_terminator,
            terminator_length=self.terminator_length,
            drop_premature_terminators=self.drop_premature_terminators,
        )


class InspectRequest(BaseModel):
    target_fasta: str = Field(..., description="Target as FASTA or a bare sequence.")
    panel_fasta: Optional[str] = Field(None, description="Optional panel FASTA.")


class DesignRequest(BaseModel):
    target_fasta: str
    panel_fasta: Optional[str] = None
    params: DesignParamsSchema = Field(default_factory=DesignParamsSchema)


class CreateSessionRequest(BaseModel):
    target_fasta: str
    panel_fasta: Optional[str] = None


class SequenceSummary(BaseModel):
    id: str
    description: str
    length: int
    moltype: str
    is_nucleotide: bool
    gc_percent: Optional[float] = None


class PanelSummary(BaseModel):
    name: str
    count: int
    uniform: bool
    members: list[SequenceSummary]


class InspectResponse(BaseModel):
    target: SequenceSummary
    target_ready: bool
    note: Optional[str] = None
    panel: Optional[PanelSummary] = None


class CriteriaSummary(BaseModel):
    rules: list[str]
    rule_mode: str
    seed_filter: bool
    seed_tm_threshold_c: Optional[float] = None
    offtarget_filter: bool
    offtarget_min_mismatch: Optional[int] = None
    # Transcripts the gate treated as the design's own; empty on an active
    # filter means every hit counted (see docs/api.md).
    offtarget_self_accessions: list[str] = Field(default_factory=list)
    offtarget_self_declared: int = 0
    offtarget_self_resolved: int = 0
    offtarget_self_coverage: Optional[float] = None
    gc_filter: bool
    gc_min: Optional[float] = None
    gc_max: Optional[float] = None
    hairpin_filter: bool
    hairpin_tm_max_c: Optional[float] = None
    homopolymer_filter: bool
    max_homopolymer_run: Optional[int] = None
    asymmetry_filter: bool


class TargetingSummary(BaseModel):
    mode: str
    hit: list[str]
    spare: list[str]
    discriminating: bool


class DiscriminationResult(BaseModel):
    """The best-placed difference to one spared allele, which sets its rank.

    all_guide_positions lists every difference to that allele in ranking order.
    """

    guide_position: int
    mutation: str
    mismatch_class: Optional[str] = None
    region: str
    sources_agree: bool
    huang_recommended: Optional[bool] = None
    strongest_position: bool
    # Membership of the P9-P16 window; the name is kept for older clients.
    consistent_position: bool
    # Two or more differences; reported only, not ranked.
    multi_mismatch: bool
    # Defaulted so older clients still validate.
    all_guide_positions: list[int] = []


class PerAlleleDiscrimination(BaseModel):
    """One spared allele: its difference count and best-placed difference."""

    allele: str
    diff_count: int
    discrimination: Optional[DiscriminationResult] = None


class AlleleResult(BaseModel):
    hit_max_mismatch: int
    spare_min_mismatch: Optional[int] = None
    spare_count: int = 0
    # The weakest allele's discrimination, then all spared alleles weakest first.
    discrimination: Optional[DiscriminationResult] = None
    per_allele: list[PerAlleleDiscrimination] = []


class CandidateResult(BaseModel):
    region: str
    start: int
    end: int
    guide_5to3: str
    guide_overhang: str
    passenger_5to3: str
    passenger_overhang: str
    gc_percent: float
    rules_passed: str
    seed_max_tm_c: Optional[float] = None
    offtarget_min_mismatch: Optional[int] = None
    offtarget_hits: Optional[list[int]] = None
    hairpin_tm_c: Optional[float] = None
    homopolymer_max_run: Optional[int] = None
    allele: Optional[AlleleResult] = None


class DesignResponse(BaseModel):
    target: SequenceSummary
    length: int
    overhang: str
    candidate_count: int
    offtarget_skipped: bool
    criteria: CriteriaSummary
    funnel: dict[str, int]
    targeting: Optional[TargetingSummary] = None
    selected: list[CandidateResult]


class ConstructResult(BaseModel):
    region: str
    start: int
    end: int
    loop_name: str
    loop_seq: str
    orientation: str
    hairpin_rna_5to3: str
    insert_top_dna_5to3: str
    insert_bottom_dna_5to3: str
    folds_as_designed: bool
    fold_dg37: Optional[float] = None
    fold_tm_c: Optional[float] = None
    starts_with_purine: bool
    premature_terminator_free: bool
    max_homopolymer_run: int
    sense_overhang: str
    guide_overhang: str


class ConvertResponse(BaseModel):
    loop_name: str
    loop_seq: str
    stem_overhang: str = StemOverhang.KEEP.value
    terminator_length: Optional[int] = None
    construct_count: int
    fold_ok_count: int
    purine_start_count: int
    converted_count: int
    withheld_count: int
    constructs: list[ConstructResult]
    withheld: list[ConstructResult] = []


class SessionResponse(BaseModel):
    id: str
    target: SequenceSummary
    panel: Optional[PanelSummary] = None
    has_last_design: bool


class HealthResponse(BaseModel):
    """Liveness, plus what the off-target check needs on this host."""

    status: str
    service: str
    version: str
    # bowtie1 found on PATH (searched as the default bowtie_path).
    bowtie_available: bool
    # The default transcriptome index exists; 'hush reference build' makes it.
    offtarget_reference: bool


class ErrorResponse(BaseModel):
    detail: str
