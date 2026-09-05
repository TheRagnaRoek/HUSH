"""Candidate selection - the design pipeline's option-rich filter step.

Applies a user's design choices to the enumerated candidates. Deliberately NOT a
ranking: the signals live on incommensurable scales, so any single combined score
would need weights HUSH cannot justify without efficacy data. Instead every
criterion is an explicit, optional gate (in the spirit of siDirect), and each
candidate keeps its full per-criterion breakdown so any verdict is legible.

Two groups of criteria:

* The three rule sets (Ui-Tei, Reynolds, Amarzguioui), toggled independently and
  combined by RuleMode.ALL (must satisfy every enabled set) or RuleMode.ANY (kept
  if any enabled set passes). No rule enabled -> no rule constraint.
* Independent filters, each optional, ANDed on top: a GC window, the seed-duplex
  off-target Tm, and the loading asymmetry.

Customisable values (GC range, seed Tm, rule thresholds) are options, not edits to
the rules. The standalone GC filter is separate from Reynolds' internal GC
criterion (fixed at the published 30-52%). Nothing is ranked: selected candidates
are returned in target-position order.

Rationale in docs/design_notes.md (selection section).
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

# --------------------------------------------------------------------------- #
# Rule combination mode
# --------------------------------------------------------------------------- #


class RuleMode(str, Enum):
    """How the enabled rule sets combine.

    ALL: intersection - a candidate must satisfy every enabled rule set.
    ANY: union - a candidate is kept if any enabled rule set passes it.
    """

    ALL = "all"
    ANY = "any"


# --------------------------------------------------------------------------- #
# GC-content window default
# --------------------------------------------------------------------------- #

# The standalone GC filter's default window, as a fraction. It mirrors Reynolds'
# published 30-52%, but here it is a tunable option applied independently, not
# the rule criterion. Off by default (GC is already judged inside Reynolds).
DEFAULT_GC_MIN = 0.30
DEFAULT_GC_MAX = 0.52

# Longest run of one nucleotide the homopolymer filter tolerates; a run longer
# than this (>= 4 by default) is rejected. Guards against synthesis-troublesome
# and G-quadruplex-prone stretches (e.g. GGGG). Off by default.
DEFAULT_MAX_HOMOPOLYMER_RUN = 3


def longest_homopolymer_run(seq: str) -> tuple[int, str]:
    """Return (length, base) of the longest single-nucleotide run in a strand."""
    best_len, best_base = 0, ""
    run_len, prev = 0, ""
    for base in seq:
        run_len = run_len + 1 if base == prev else 1
        prev = base
        if run_len > best_len:
            best_len, best_base = run_len, base
    return best_len, best_base


# --------------------------------------------------------------------------- #
# Criteria
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SelectionCriteria:
    """A user's selection options: which criteria to apply and with what values.

    Rule sets:
        use_ui_tei / use_reynolds / use_amarzguioui: enable each rule set.
        rule_mode: ALL or ANY combination of the enabled rule sets.
        reynolds_threshold / amarzguioui_threshold: the graded rule sets' score
            cut-offs (Ui-Tei is an all-or-nothing filter and has no threshold).

    Independent filters (each ANDed on top of the rule step when enabled):
        use_gc_filter, gc_min, gc_max: a standalone GC-content window. Distinct
            from Reynolds criterion I (which stays at the published 30-52%).
        use_seed_filter, seed_tm_threshold: the siDirect seed-duplex off-target
            Tm filter, applied to both strands.
        use_offtarget_filter, offtarget_min_mismatch, offtarget_searcher,
            offtarget_self_accession: the transcriptome homology off-target
            filter. Rejects a candidate whose nearest off-target (on either
            strand, excluding its own source transcript) has fewer than
            offtarget_min_mismatch mismatches. Needs a searcher (a reference
            database handle); distinct from the seed filter above, which is a
            sequence-intrinsic proxy rather than a database search.
        use_hairpin_filter, hairpin_tm_max: a standalone internal-structure gate -
            reject a candidate whose most stable self-hairpin melts above
            hairpin_tm_max. Distinct from Reynolds criterion (its internal-repeat
            sub-criterion stays Xia-pinned at 20 C); this is the tunable one.
        use_homopolymer_filter, max_homopolymer_run: reject a candidate whose
            core contains a single-nucleotide run longer than max_homopolymer_run
            (e.g. GGGG). Judged on the core only; overhangs are not counted.
        use_asymmetry_filter: require the duplex end asymmetry to favour loading
            the guide strand (sign only - no magnitude threshold).

    Advanced thermodynamic knobs (safe to leave at their defaults):
        asymmetry_window: terminal base pairs compared at each duplex end.
        asymmetry_temp_c: temperature the asymmetry dG is evaluated at.
        thermo_params: nearest-neighbor parameter set for the asymmetry and
            hairpin readouts ONLY. It never touches the seed filter (pinned to
            Freier, on which its 21.5 C threshold is calibrated) or the rule sets.
    """

    use_ui_tei: bool = True
    use_reynolds: bool = True
    use_amarzguioui: bool = True
    rule_mode: RuleMode = RuleMode.ALL
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
    offtarget_self_accession: str | None = None

    use_hairpin_filter: bool = False
    hairpin_tm_max: float = DEFAULT_INTERNAL_TM_THRESHOLD_C

    use_homopolymer_filter: bool = False
    max_homopolymer_run: int = DEFAULT_MAX_HOMOPOLYMER_RUN

    use_asymmetry_filter: bool = False

    asymmetry_window: int = DEFAULT_WINDOW_BP
    asymmetry_temp_c: float = REFERENCE_TEMPERATURE_C
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
            # The search must reach at least threshold-1 mismatches to see every
            # disqualifying near-match.
            if self.offtarget_searcher.max_mismatch < self.offtarget_min_mismatch - 1:
                raise ValueError(
                    f"offtarget searcher depth {self.offtarget_searcher.max_mismatch} "
                    f"is too shallow for offtarget_min_mismatch "
                    f"{self.offtarget_min_mismatch}; it must be >= "
                    f"{self.offtarget_min_mismatch - 1}."
                )

    @property
    def enabled_rule_names(self) -> tuple[str, ...]:
        """Names of the rule sets that are switched on, in a stable order."""
        names = []
        if self.use_ui_tei:
            names.append("ui_tei")
        if self.use_reynolds:
            names.append("reynolds")
        if self.use_amarzguioui:
            names.append("amarzguioui")
        return tuple(names)


# --------------------------------------------------------------------------- #
# Evaluation result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CandidateEvaluation:
    """The outcome of applying one SelectionCriteria to one candidate.

    Component results are present only for the criteria that were enabled; a
    disabled criterion is None (not computed, not a black box - simply not asked
    for). Each *_pass flag is True/False when its criterion was applied and None
    when it was not.

    'selected' is the overall verdict; 'reasons' lists, in plain language, the
    criteria a rejected candidate failed (empty when selected).
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


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def _rule_pass(votes: list[bool], mode: RuleMode) -> bool:
    """Combine the enabled rule verdicts. No rules enabled -> no constraint."""
    if not votes:
        return True
    return all(votes) if mode is RuleMode.ALL else any(votes)


def evaluate(
    candidate: SiRNACandidate, criteria: SelectionCriteria
) -> CandidateEvaluation:
    """Apply 'criteria' to one candidate and return the full breakdown.

    Only the enabled criteria are computed. The rule sets are combined by
    'criteria.rule_mode'; the enabled independent filters are then ANDed on top.
    """
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
            self_accession=criteria.offtarget_self_accession,
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
    """Human-readable list of the criteria a rejected candidate failed."""
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
    """Evaluate every candidate and return the selected ones, in input order.

    For the full set of evaluations (selected and rejected alike), map evaluate()
    over the candidates directly.
    """
    return [
        evaluation
        for evaluation in evaluate_all(candidates, criteria)
        if evaluation.selected
    ]


def evaluate_all(
    candidates: Iterable[SiRNACandidate], criteria: SelectionCriteria
) -> Iterator[CandidateEvaluation]:
    """Lazily evaluate every candidate against 'criteria' (selected or not)."""
    for candidate in candidates:
        yield evaluate(candidate, criteria)
