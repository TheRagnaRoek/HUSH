"""HushSession - the stateful orchestration layer.

One session is one working instance of HUSH: a loaded target (the template
mRNA/CDS) and an optional related-sequence panel, plus the MSA alignment cached
so repeated design runs with different parameters do not re-align. Both the CLI
and the future HTTP API drive the pipeline through this object rather than
reaching into core/msa/shrna directly.

Every method returns a structured result (service.results) and raises the domain
exceptions (core.errors.HushError, or ValueError for an incoherent request). It
prints nothing and knows nothing about terminals or HTTP - turning a result into
a table or JSON, and a ValueError into an exit code or a 4xx, is the caller's job.
"""

from __future__ import annotations

import re
from pathlib import Path

from core.errors import HushError
from core.inputs import load_panel as _load_panel_file
from core.inputs import load_target as _load_target_file
from core.offtarget import BowtieSearcher
from core.selection import SelectionCriteria, evaluate_all
from core.sequence import BioSequence, DesignInput, SequencePanel
from core.sirna import SiRNACandidate, enumerate_candidates
from msa.align_star import align_star
from msa.targeting import AlleleTargeting
from msa.variation import PanelVariation
from msa.windows import WindowCriteria, evaluate_windows
from shrna.convert import ShrnaOptions, convert_candidate

from service.params import AlleleTargetingSpec, DesignParams
from service.results import ConvertResult, DesignResult, InspectResult

# RefSeq accession, e.g. NM_005018.3 / NR_ / XM_ / XR_, with optional version.
_REFSEQ_ACCESSION = re.compile(r"\b([NX][MR]_\d+(?:\.\d+)?)\b")


def _detect_accession(seq: BioSequence) -> str | None:
    """Best-effort RefSeq accession from a target's id/description.

    The off-target search must drop hits to the candidate's own transcript. The
    accession usually sits in the FASTA header; pass an explicit self-accession
    to override when it does not.
    """
    match = _REFSEQ_ACCESSION.search(f"{seq.id} {seq.description}")
    return match.group(1) if match else None


def resolve_targeting(
    spec: AlleleTargetingSpec | None, panel: SequencePanel
) -> AlleleTargeting:
    """Resolve a targeting intent into an AlleleTargeting against a panel.

    Raises ValueError for an incoherent intent (mutually exclusive intents, or
    spare without hit), and the msa TargetingError/PanelError (both HushError)
    when a named id is not in the panel. A None spec is the default: pan-allele.
    """
    spec = spec or AlleleTargetingSpec()
    intents = [
        name
        for name, chosen in (
            ("pan_allele", spec.pan_allele),
            ("allele_specific", spec.allele_specific is not None),
            ("hit", spec.hit is not None),
        )
        if chosen
    ]
    if len(intents) > 1:
        raise ValueError(
            f"{' and '.join(intents)} are mutually exclusive; choose one "
            "targeting intent."
        )
    if spec.spare is not None and spec.hit is None:
        raise ValueError(
            "spare only refines an explicit hit; name the alleles to hit too."
        )

    if spec.allele_specific is not None:
        return AlleleTargeting.allele_specific(spec.allele_specific, panel)
    if spec.hit is not None:
        if spec.spare is not None:
            return AlleleTargeting(hit=spec.hit, spare=spec.spare)
        return AlleleTargeting.subset(spec.hit, panel)
    return AlleleTargeting.pan_allele(panel)


def sort_selection(selected, sort, window_by_pos, targeting):
    """Order selected candidates. Position (the default) ranks nothing.

    Discrimination order is opt-in because it imposes a priority between
    published findings the measurements do not fully settle; it needs a request
    that spares something to have anything to rank. Position breaks ties so two
    runs never disagree.
    """
    if sort == "position":
        return selected
    if sort != "discrimination":
        raise ValueError(
            f"unknown sort {sort!r}; choose 'position' or 'discrimination'."
        )
    if window_by_pos is None or targeting is None or not targeting.discriminating:
        raise ValueError(
            "sort discrimination needs a request that spares something; "
            "pan-allele targeting asks for no discrimination to rank."
        )
    return sorted(
        selected,
        key=lambda e: (
            window_by_pos[e.candidate.position].discrimination_key,
            -e.candidate.position,
        ),
        reverse=True,
    )


class HushSession:
    """A loaded target + optional panel, and the operations over them.

    Construct empty and load inputs with load_target / load_panel (or pass them
    in), then call inspect / design / convert. The alignment underlying the
    allele stage is computed once per (target, panel) and reused across design
    calls; changing either input invalidates it.
    """

    def __init__(
        self,
        target: BioSequence | None = None,
        panel: SequencePanel | None = None,
    ) -> None:
        self._target = target
        self._panel = panel
        self._variation: PanelVariation | None = None
        self._last_design: DesignResult | None = None

    # -- loading ----------------------------------------------------------- #

    def load_target(self, path: str | Path, fmt: str | None = None) -> BioSequence:
        """Load the target (template) sequence from a file and make it current."""
        self._target = _load_target_file(path, fmt=fmt)
        self._invalidate()
        return self._target

    def set_target(self, target: BioSequence) -> None:
        """Set the target from an in-memory sequence (e.g. an API upload)."""
        self._target = target
        self._invalidate()

    def load_panel(
        self, path: str | Path, name: str | None = None, fmt: str | None = None
    ) -> SequencePanel:
        """Load the related-sequence panel from a file and make it current."""
        self._panel = _load_panel_file(path, name=name, fmt=fmt)
        self._invalidate()
        return self._panel

    def set_panel(self, panel: SequencePanel | None) -> None:
        """Set (or clear, with None) the panel from an in-memory object."""
        self._panel = panel
        self._invalidate()

    def clear_panel(self) -> None:
        """Drop the panel, returning the session to single-target design."""
        self._panel = None
        self._invalidate()

    def _invalidate(self) -> None:
        # The cached alignment and the last result both belong to the old inputs.
        self._variation = None
        self._last_design = None

    # -- accessors --------------------------------------------------------- #

    @property
    def target(self) -> BioSequence | None:
        return self._target

    @property
    def panel(self) -> SequencePanel | None:
        return self._panel

    @property
    def has_panel(self) -> bool:
        return self._panel is not None and len(self._panel) > 0

    @property
    def last_design(self) -> DesignResult | None:
        return self._last_design

    def _require_target(self) -> BioSequence:
        if self._target is None:
            raise HushError("No target loaded; load a target sequence first.")
        return self._target

    def design_input(self) -> DesignInput:
        """The current inputs as a DesignInput hand-off."""
        return DesignInput(target=self._require_target(), panel=self._panel)

    # -- operations -------------------------------------------------------- #

    def inspect(self) -> InspectResult:
        """Preview the loaded inputs (target role, panel, readiness)."""
        return InspectResult(design_input=self.design_input())

    def _panel_variation(self, mrna: BioSequence) -> PanelVariation:
        # Align the panel to the target once; window evaluation (which depends
        # on the design params) is cheap and re-run per design call.
        if self._variation is None:
            self._variation = PanelVariation.from_alignment(
                align_star(mrna, self._panel)
            )
        return self._variation

    def design(self, params: DesignParams | None = None) -> DesignResult:
        """Enumerate candidates for the target and filter them by the params.

        With a panel loaded, the allele gate (aligned once, cached) is ANDed with
        the siRNA filters. Raises HushError for a protein target or an unknown
        allele id, and ValueError for an incoherent request or bad criteria.
        """
        params = params or DesignParams()
        design_input = self.design_input()
        mrna = design_input.require_nucleotide_target()
        candidates = list(
            enumerate_candidates(mrna, length=params.length, overhang=params.overhang)
        )

        # Allele gate (MSA / stage 3). Windows and candidates share the target's
        # coordinates one-for-one, so a candidate is looked up by start position.
        targeting = None
        window_by_pos = None
        if design_input.has_panel:
            targeting = resolve_targeting(params.targeting, design_input.panel)
            window_criteria = WindowCriteria(
                targeting=targeting,
                length=params.length,
                min_spare_mismatches=params.min_spare_mismatches,
            )
            window_by_pos = {
                ev.position: ev
                for ev in evaluate_windows(
                    self._panel_variation(mrna), window_criteria
                )
            }

        # Off-target scanning is on by default but needs a reference index; with
        # none configured it is skipped (a flag the caller can warn about).
        offtarget_searcher = None
        self_accession = None
        offtarget_active = params.use_offtarget_filter and params.offtarget_db is not None
        offtarget_skipped = params.use_offtarget_filter and params.offtarget_db is None
        if offtarget_active:
            self_accession = params.offtarget_self_accession or _detect_accession(mrna)
            backend = BowtieSearcher(
                params.offtarget_db,
                max_mismatch=params.offtarget_max_mismatch,
                bowtie_path=params.bowtie_path,
            )
            reads = [c.sense.seq for c in candidates]
            reads += [c.antisense.seq for c in candidates]
            offtarget_searcher = backend.precompute(reads)

        criteria = SelectionCriteria(
            use_ui_tei=params.use_ui_tei,
            use_reynolds=params.use_reynolds,
            use_amarzguioui=params.use_amarzguioui,
            rule_mode=params.rule_mode,
            reynolds_threshold=params.reynolds_threshold,
            amarzguioui_threshold=params.amarzguioui_threshold,
            use_seed_filter=params.use_seed_filter,
            seed_tm_threshold=params.seed_tm_threshold,
            use_offtarget_filter=offtarget_active,
            offtarget_min_mismatch=params.offtarget_min_mismatch,
            offtarget_searcher=offtarget_searcher,
            offtarget_self_accession=self_accession,
            use_gc_filter=params.use_gc_filter,
            gc_min=params.gc_min,
            gc_max=params.gc_max,
            use_hairpin_filter=params.use_hairpin_filter,
            hairpin_tm_max=params.hairpin_tm_max,
            use_homopolymer_filter=params.use_homopolymer_filter,
            max_homopolymer_run=params.max_homopolymer_run,
            use_asymmetry_filter=params.use_asymmetry_filter,
            asymmetry_window=params.asymmetry_window,
            asymmetry_temp_c=params.asymmetry_temp_c,
            thermo_params=params.thermo_params,
        )
        evaluations = list(evaluate_all(candidates, criteria))

        def _allele_ok(evaluation) -> bool:
            return (
                window_by_pos is None
                or window_by_pos[evaluation.candidate.position].satisfies
            )

        selected = [e for e in evaluations if e.selected and _allele_ok(e)]
        selected = sort_selection(selected, params.sort, window_by_pos, targeting)

        result = DesignResult(
            target=mrna,
            length=params.length,
            overhang=params.overhang,
            criteria=criteria,
            evaluations=evaluations,
            selected=selected,
            targeting=targeting,
            window_by_pos=window_by_pos,
            min_spare_mismatches=params.min_spare_mismatches,
            offtarget_skipped=offtarget_skipped,
        )
        self._last_design = result
        return result

    def convert(
        self,
        options: ShrnaOptions | None = None,
        candidates: list[SiRNACandidate] | None = None,
    ) -> ConvertResult:
        """Graft siRNA candidates onto a stem-loop scaffold as shRNA constructs.

        With no candidates given, the selected candidates of the most recent
        design run are used, so a session can design then convert without a TSV
        round-trip. The 19-bp core is preserved unchanged.
        """
        options = options or ShrnaOptions()
        if candidates is None:
            if self._last_design is None:
                raise HushError(
                    "Nothing to convert; run design first, or pass candidates."
                )
            candidates = [e.candidate for e in self._last_design.selected]
        constructs = [convert_candidate(c, options) for c in candidates]
        return ConvertResult(constructs=constructs, options=options)
