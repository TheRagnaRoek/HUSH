"""Candidate selection: every enabled criterion is an independent gate, ANDed.

There is no combined score. Each candidate keeps its per-criterion breakdown.
Selected candidates keep input order. See docs/design-rules.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Iterator

from core.asymmetry import DEFAULT_WINDOW_BP, EndAsymmetry, end_asymmetry
from core.hairpin import DEFAULT_INTERNAL_TM_THRESHOLD_C, internal_structure_tm
from core.offtarget import (
    DEFAULT_OFFTARGET_MIN_MISMATCH,
    OffTargetProfile,
    OffTargetSearcher,
    SelfTranscripts,
    evaluate_offtarget,
)
from core.rules import (
    AMARZGUIOUI_DEFAULT_THRESHOLD,
    REYNOLDS_DEFAULT_THRESHOLD,
    AmarzguiouiResult,
    ReynoldsResult,
    UiTeiResult,
    amarzguioui,
    reynolds,
    ui_tei,
)
from core.seed import SEED_TM_THRESHOLD_C, SeedTmResult, seed_tm
from core.sirna import SiRNACandidate
from core.thermo import REFERENCE_TEMPERATURE_C, ThermoParams


class RuleMode(str, Enum):
    """How enabled rule sets combine: ALL (intersection) or ANY (union)."""

    ALL = "all"
    ANY = "any"


# Standalone GC window, as fractions.
DEFAULT_GC_MIN = 0.30
DEFAULT_GC_MAX = 0.52

# A run longer than this is rejected.
DEFAULT_MAX_HOMOPOLYMER_RUN = 3


def longest_homopolymer_run(seq: str) -> tuple[int, str]:
    """(length, base) of the longest single-nucleotide run."""
    best_len, best_base = 0, ""
    run_len, prev = 0, ""
    for base in seq:
        run_len = run_len + 1 if base == prev else 1
        prev = base
        if run_len > best_len:
            best_len, best_base = run_len, base
    return best_len, best_base


@dataclass(frozen=True)
class SelectionCriteria:
    """Which gates to apply, and their values.

    The defaults are the paradigm in docs/cli.md. The off-target gate is off
    here because it needs a searcher; the service layer enables it. The
    self-transcript set is resolved by the service layer.
    """

    use_ui_tei: bool = True
    use_reynolds: bool = True
    use_amarzguioui: bool = True
    rule_mode: RuleMode = RuleMode.ANY
    reynolds_threshold: int = REYNOLDS_DEFAULT_THRESHOLD
    amarzguioui_threshold: int = AMARZGUIOUI_DEFAULT_THRESHOLD

    use_gc_filter: bool = False
    gc_min: float = DEFAULT_GC_MIN
    gc_max: float = DEFAULT_GC_MAX

    use_seed_filter: bool = True
    seed_tm_threshold: float = SEED_TM_THRESHOLD_C

    use_offtarget_filter: bool = False
    offtarget_min_mismatch: int = DEFAULT_OFFTARGET_MIN_MISMATCH
    offtarget_searcher: OffTargetSearcher | None = None
    offtarget_self_transcripts: SelfTranscripts | None = None

    use_hairpin_filter: bool = False
    hairpin_tm_max: float = DEFAULT_INTERNAL_TM_THRESHOLD_C

    # Judged on the core only.
    use_homopolymer_filter: bool = True
    max_homopolymer_run: int = DEFAULT_MAX_HOMOPOLYMER_RUN

    use_asymmetry_filter: bool = False

    asymmetry_window: int = DEFAULT_WINDOW_BP
    asymmetry_temp_c: float = REFERENCE_TEMPERATURE_C
    # Applies to the asymmetry and hairpin read-outs only, never to the seed
    # filter or the rule sets.
    thermo_params: ThermoParams = ThermoParams.XIA_1998

    def __post_init__(self) -> None:
        if not isinstance(self.rule_mode, RuleMode):
            raise TypeError(
                f"rule_mode must be a RuleMode, got {type(self.rule_mode).__name__}"
            )
        if not isinstance(self.thermo_params, ThermoParams):
            raise TypeError(
                "thermo_params must be a ThermoParams, got "
                f"{type(self.thermo_params).__name__}"
            )
        if self.use_gc_filter and not (0.0 <= self.gc_min <= self.gc_max <= 1.0):
            raise ValueError(
                f"GC window must satisfy 0 <= gc_min <= gc_max <= 1; got "
                f"[{self.gc_min}, {self.gc_max}]."
            )
        if self.use_homopolymer_filter and self.max_homopolymer_run < 1:
            raise ValueError(
                f"max_homopolymer_run must be >= 1; got {self.max_homopolymer_run}."
            )
        if self.use_offtarget_filter:
            if self.offtarget_searcher is None:
                raise ValueError(
                    "use_offtarget_filter requires an offtarget_searcher "
                    "(a reference database handle)."
                )
            if self.offtarget_min_mismatch < 1:
                raise ValueError(
                    "offtarget_min_mismatch must be >= 1; got "
                    f"{self.offtarget_min_mismatch}."
                )
            # The search must reach threshold-1 mismatches to see every
            # disqualifying hit.
            if self.offtarget_searcher.max_mismatch < self.offtarget_min_mismatch - 1:
                raise ValueError(
                    f"offtarget searcher depth {self.offtarget_searcher.max_mismatch} "
                    f"is too shallow for offtarget_min_mismatch "
                    f"{self.offtarget_min_mismatch}; it must be >= "
                    f"{self.offtarget_min_mismatch - 1}."
                )

    @property
    def enabled_rule_names(self) -> tuple[str, ...]:
        """Names of the enabled rule sets."""
        names = []
        if self.use_ui_tei:
            names.append("ui_tei")
        if self.use_reynolds:
            names.append("reynolds")
        if self.use_amarzguioui:
            names.append("amarzguioui")
        return tuple(names)


@dataclass(frozen=True)
class CandidateEvaluation:
    """One candidate's verdict under a SelectionCriteria.

    Results and *_pass flags are None for gates that were not enabled.
    'reasons' names the gates a rejected candidate failed.
    """

    candidate: SiRNACandidate
    ui_tei: UiTeiResult | None
    reynolds: ReynoldsResult | None
    amarzguioui: AmarzguiouiResult | None
    seed: SeedTmResult | None
    offtarget: OffTargetProfile | None
    asymmetry: EndAsymmetry | None
    hairpin_tm: float | None
    homopolymer_run: int | None
    gc_content: float
    rule_pass: bool
    gc_pass: bool | None
    seed_pass: bool | None
    offtarget_pass: bool | None
    hairpin_pass: bool | None
    homopolymer_pass: bool | None
    asymmetry_pass: bool | None
    selected: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)

    def __repr__(self) -> str:
        verdict = "selected" if self.selected else "rejected"
        return (
            f"CandidateEvaluation(region={self.candidate.region}, {verdict}"
            + ("" if self.selected else f", reasons={list(self.reasons)}")
            + ")"
        )


def _rule_pass(votes: list[bool], mode: RuleMode) -> bool:
    """Combine the enabled rule verdicts; no enabled rule means no constraint."""
    if not votes:
        return True
    return all(votes) if mode is RuleMode.ALL else any(votes)


def evaluate(
    candidate: SiRNACandidate, criteria: SelectionCriteria
) -> CandidateEvaluation:
    """Evaluate one candidate; only the enabled gates are computed."""
    ut = ui_tei(candidate) if criteria.use_ui_tei else None
    rey = reynolds(candidate) if criteria.use_reynolds else None
    amz = amarzguioui(candidate) if criteria.use_amarzguioui else None

    votes: list[bool] = []
    if ut is not None:
        votes.append(ut.satisfied)
    if rey is not None:
        votes.append(rey.meets_threshold(criteria.reynolds_threshold))
    if amz is not None:
        votes.append(amz.meets_threshold(criteria.amarzguioui_threshold))
    rule_pass = _rule_pass(votes, criteria.rule_mode)

    gc = candidate.gc_content
    gc_pass = (
        criteria.gc_min <= gc <= criteria.gc_max if criteria.use_gc_filter else None
    )

    seed = (
        seed_tm(candidate, threshold=criteria.seed_tm_threshold)
        if criteria.use_seed_filter
        else None
    )
    seed_pass = seed.passes if seed is not None else None

    if criteria.use_offtarget_filter:
        assert criteria.offtarget_searcher is not None  # guarded in __post_init__
        offtarget = evaluate_offtarget(
            candidate,
            criteria.offtarget_searcher,
            threshold=criteria.offtarget_min_mismatch,
            self_transcripts=criteria.offtarget_self_transcripts,
        )
        offtarget_pass = offtarget.passes
    else:
        offtarget = None
        offtarget_pass = None

    if criteria.use_hairpin_filter:
        hairpin_tm = internal_structure_tm(
            candidate.sense, params=criteria.thermo_params.table
        )
        hairpin_pass = hairpin_tm is None or hairpin_tm <= criteria.hairpin_tm_max
    else:
        hairpin_tm = None
        hairpin_pass = None

    if criteria.use_homopolymer_filter:
        homopolymer_run, homopolymer_base = longest_homopolymer_run(candidate.sense.seq)
        homopolymer_pass = homopolymer_run <= criteria.max_homopolymer_run
    else:
        homopolymer_run, homopolymer_base = None, ""
        homopolymer_pass = None

    asym = (
        end_asymmetry(
            candidate,
            window=criteria.asymmetry_window,
            temperature_c=criteria.asymmetry_temp_c,
            params=criteria.thermo_params.table,
        )
        if criteria.use_asymmetry_filter
        else None
    )
    asymmetry_pass = asym.favours_guide if asym is not None else None

    gates = [rule_pass]
    gates += [
        p
        for p in (
            gc_pass, seed_pass, offtarget_pass,
            hairpin_pass, homopolymer_pass, asymmetry_pass,
        )
        if p is not None
    ]
    selected = all(gates)

    reasons = () if selected else _reasons(
        criteria, rule_pass, gc, gc_pass, seed, seed_pass,
        offtarget, offtarget_pass,
        hairpin_tm, hairpin_pass, homopolymer_run, homopolymer_base,
        homopolymer_pass, asym, asymmetry_pass,
    )

    return CandidateEvaluation(
        candidate=candidate,
        ui_tei=ut,
        reynolds=rey,
        amarzguioui=amz,
        seed=seed,
        offtarget=offtarget,
        asymmetry=asym,
        hairpin_tm=hairpin_tm,
        homopolymer_run=homopolymer_run,
        gc_content=gc,
        rule_pass=rule_pass,
        gc_pass=gc_pass,
        seed_pass=seed_pass,
        offtarget_pass=offtarget_pass,
        hairpin_pass=hairpin_pass,
        homopolymer_pass=homopolymer_pass,
        asymmetry_pass=asymmetry_pass,
        selected=selected,
        reasons=reasons,
    )


def _reasons(
    criteria: SelectionCriteria,
    rule_pass: bool,
    gc: float,
    gc_pass: bool | None,
    seed: SeedTmResult | None,
    seed_pass: bool | None,
    offtarget: OffTargetProfile | None,
    offtarget_pass: bool | None,
    hairpin_tm: float | None,
    hairpin_pass: bool | None,
    homopolymer_run: int | None,
    homopolymer_base: str,
    homopolymer_pass: bool | None,
    asym: EndAsymmetry | None,
    asymmetry_pass: bool | None,
) -> tuple[str, ...]:
    """Plain-language list of the gates a rejected candidate failed."""
    out: list[str] = []
    if not rule_pass:
        joined = "/".join(criteria.enabled_rule_names)
        out.append(f"rule selection not met ({criteria.rule_mode.value} of {joined})")
    if gc_pass is False:
        out.append(
            f"GC {gc:.0%} outside [{criteria.gc_min:.0%}, {criteria.gc_max:.0%}]"
        )
    if seed_pass is False and seed is not None:
        out.append(
            f"seed off-target Tm {seed.max_tm:.1f}C >= {seed.threshold:.1f}C"
        )
    if offtarget_pass is False and offtarget is not None:
        out.append(
            f"off-target at {offtarget.min_mismatch} mismatch(es) "
            f"< {offtarget.threshold}"
        )
    if hairpin_pass is False and hairpin_tm is not None:
        out.append(
            f"internal structure Tm {hairpin_tm:.1f}C > {criteria.hairpin_tm_max:.1f}C"
        )
    if homopolymer_pass is False and homopolymer_run is not None:
        out.append(
            f"homopolymer run of {homopolymer_run} {homopolymer_base} "
            f"(> {criteria.max_homopolymer_run})"
        )
    if asymmetry_pass is False and asym is not None:
        out.append(f"{asym.preferred_strand}-favoured loading")
    return tuple(out)


def select(
    candidates: Iterable[SiRNACandidate], criteria: SelectionCriteria
) -> list[CandidateEvaluation]:
    """The selected evaluations, in input order."""
    return [
        evaluation
        for evaluation in evaluate_all(candidates, criteria)
        if evaluation.selected
    ]


def evaluate_all(
    candidates: Iterable[SiRNACandidate], criteria: SelectionCriteria
) -> Iterator[CandidateEvaluation]:
    """Lazily evaluate all candidates, selected or not."""
    for candidate in candidates:
        yield evaluate(candidate, criteria)
