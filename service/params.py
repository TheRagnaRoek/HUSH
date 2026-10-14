"""Request objects for the HUSH service layer.

Framework-agnostic, already-parsed design/convert requests. A CLI turns flag
strings into these; a future HTTP layer turns JSON into these; the service
itself speaks only this vocabulary. Fields mirror the tunable pipeline knobs and
default to the best-practice paradigm (see 'hush design --help' and
docs/cli.md). GC bounds are fractions here, not percent.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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


@dataclass(frozen=True)
class AlleleTargetingSpec:
    """A targeting intent by panel-member id, before it is resolved to sets.

    The four named intents (pan-allele, allele-specific, target-specific,
    explicit hit/spare) are mutually exclusive. Resolution against a panel -
    which fills in the spared set and validates the ids - happens in the session,
    where the panel is known. The all-default spec is the pan-allele intent (hit
    everything, spare nothing).
    """

    pan_allele: bool = False
    allele_specific: str | None = None
    target_specific: bool = False
    hit: tuple[str, ...] | None = None
    spare: tuple[str, ...] | None = None


@dataclass(frozen=True)
class DesignParams:
    """Everything the design stage needs beyond the loaded target and panel."""

    length: int = 19
    overhang: OverhangPolicy = OverhangPolicy.NATURAL

    # Classical rules.
    use_ui_tei: bool = True
    use_reynolds: bool = True
    use_amarzguioui: bool = True
    rule_mode: RuleMode = RuleMode.ANY
    reynolds_threshold: int = REYNOLDS_DEFAULT_THRESHOLD
    amarzguioui_threshold: int = AMARZGUIOUI_DEFAULT_THRESHOLD

    # Seed-duplex off-target Tm.
    use_seed_filter: bool = True
    seed_tm_threshold: float = SEED_TM_THRESHOLD_C

    # Transcriptome off-target. use_offtarget_filter is the intent; with no
    # offtarget_db the session falls back to the reference built by
    # 'hush reference build', and only skips (loudly) if that is absent too.
    use_offtarget_filter: bool = True
    offtarget_db: Path | None = None
    offtarget_min_mismatch: int = DEFAULT_OFFTARGET_MIN_MISMATCH
    offtarget_max_mismatch: int = DEFAULT_OFFTARGET_MAX_MISMATCH
    bowtie_path: str = "bowtie"

    # What the off-target gate treats as the design's own rather than a hit.
    # Accessions name transcripts outright; resolve_offtarget_self finds the rest
    # by homology to the target, the panel and offtarget_self_fasta, which is the
    # only workable route for a gene family (one accession never covers one).
    offtarget_self_accessions: tuple[str, ...] = ()
    offtarget_self_fasta: Path | None = None
    resolve_offtarget_self: bool = True
    offtarget_self_coverage: float = DEFAULT_SELF_MIN_COVERAGE

    # Standalone GC window (fractions in [0, 1], not percent).
    use_gc_filter: bool = False
    gc_min: float = 0.30
    gc_max: float = 0.52

    # Homopolymer run.
    use_homopolymer_filter: bool = True
    max_homopolymer_run: int = DEFAULT_MAX_HOMOPOLYMER_RUN

    # Duplex-end loading asymmetry.
    use_asymmetry_filter: bool = False
    asymmetry_window: int = DEFAULT_WINDOW_BP
    asymmetry_temp_c: float = REFERENCE_TEMPERATURE_C

    # Internal (self-hairpin) structure.
    use_hairpin_filter: bool = False
    hairpin_tm_max: float = DEFAULT_INTERNAL_TM_THRESHOLD_C

    # Nearest-neighbour parameter set for the asymmetry and hairpin readouts
    # only (never the seed filter).
    thermo_params: ThermoParams = ThermoParams.XIA_1998

    # Allele-aware targeting; consulted only when a panel is loaded.
    targeting: AlleleTargetingSpec | None = None
    min_spare_mismatches: int = DEFAULT_MIN_SPARE_MISMATCHES

    # Report ordering: "auto" (discrimination when the request spares something,
    # else position), "position" (nothing ranked), or "discrimination".
    sort: str = "auto"
