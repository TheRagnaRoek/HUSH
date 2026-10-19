"""HushSession: a loaded target and panel, with design and convert operations.

The panel alignment is cached across design calls. Methods return the results
of service.results and raise HushError, or ValueError for an incoherent
request; they print nothing. See docs/architecture.md.
"""

from __future__ import annotations

import re
from pathlib import Path

from core.errors import HushError, NothingToConvertError
from core.inputs import load_panel as _load_panel_file
from core.inputs import load_sequences as _load_sequences_file
from core.inputs import load_target as _load_target_file
from core.offtarget import BowtieSearcher, SelfTranscripts, resolve_self_transcripts
from core.reference import default_index_prefix, index_exists
from core.selection import SelectionCriteria, evaluate_all
from core.sequence import BioSequence, DesignInput, SequencePanel
from core.sirna import SiRNACandidate, enumerate_candidates
from msa.align_star import align_star
from msa.targeting import AlleleTargeting
from msa.variation import PanelVariation
from msa.windows import WindowCriteria, evaluate_windows
from shrna.convert import ShrnaConstruct, ShrnaOptions, convert_candidate

from service.params import AlleleTargetingSpec, DesignParams
from service.results import ConvertResult, DesignResult, InspectResult

# RefSeq accession, e.g. NM_005018.3 / NR_ / XM_ / XR_, with optional version.
_REFSEQ_ACCESSION = re.compile(r"\b([NX][MR]_\d+(?:\.\d+)?)\b")


def _detect_accession(seq: BioSequence) -> str | None:
    """RefSeq accession from the target's id or description, or None (e.g. IMGT headers)."""
    match = _REFSEQ_ACCESSION.search(f"{seq.id} {seq.description}")
    return match.group(1) if match else None


def resolve_targeting(
    spec: AlleleTargetingSpec | None, panel: SequencePanel
) -> AlleleTargeting:
    """Resolve an intent against a panel; None means pan-allele.

    Raises ValueError for mutually exclusive intents or spare without hit, and
    TargetingError or PanelError for ids not in the panel.
    """
    spec = spec or AlleleTargetingSpec()
    intents = [
        name
        for name, chosen in (
            ("pan_allele", spec.pan_allele),
            ("allele_specific", spec.allele_specific is not None),
            ("target_specific", spec.target_specific),
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
    if spec.target_specific:
        return AlleleTargeting.target_specific(panel)
    if spec.hit is not None:
        if spec.spare is not None:
            return AlleleTargeting(hit=spec.hit, spare=spec.spare)
        return AlleleTargeting.subset(spec.hit, panel)
    return AlleleTargeting.pan_allele(panel)


def sort_selection(selected, sort, window_by_pos, targeting):
    """Order the selected candidates.

    "auto" ranks by discrimination when the request spares something, else by
    position. "position" keeps target order. "discrimination" raises ValueError
    when nothing is spared. Ties break by position.
    """
    discriminates = (
        window_by_pos is not None
        and targeting is not None
        and targeting.discriminating
    )
    if sort == "auto":
        sort = "discrimination" if discriminates else "position"
    if sort == "position":
        return selected
    if sort != "discrimination":
        raise ValueError(
            f"unknown sort {sort!r}; choose 'auto', 'position' or 'discrimination'."
        )
    if not discriminates:
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


def apply_convert_filters(
    constructs: list[ShrnaConstruct], options: ShrnaOptions
) -> ConvertResult:
    """Withhold constructs per the options and package the result.

    Shared by the session and the CLI.
    """
    dropped: list[ShrnaConstruct] = []
    if options.drop_premature_terminators:
        kept = [c for c in constructs if c.premature_terminator_free]
        dropped = [c for c in constructs if not c.premature_terminator_free]
        constructs = kept
    return ConvertResult(
        constructs=constructs, options=options, dropped=dropped
    )


class HushSession:
    """A loaded target and optional panel.

    The alignment is computed once per (target, panel) and reused by design();
    changing either input invalidates it.
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

    def load_target(self, path: str | Path, fmt: str | None = None) -> BioSequence:
        """Load the target from a file."""
        self._target = _load_target_file(path, fmt=fmt)
        self._invalidate()
        return self._target

    def set_target(self, target: BioSequence) -> None:
        """Set the target from an in-memory sequence."""
        self._target = target
        self._invalidate()

    def load_panel(
        self, path: str | Path, name: str | None = None, fmt: str | None = None
    ) -> SequencePanel:
        """Load the panel from a file."""
        self._panel = _load_panel_file(path, name=name, fmt=fmt)
        self._invalidate()
        return self._panel

    def set_panel(self, panel: SequencePanel | None) -> None:
        """Set the panel from an in-memory object, or clear it with None."""
        self._panel = panel
        self._invalidate()

    def clear_panel(self) -> None:
        """Drop the panel."""
        self._panel = None
        self._invalidate()

    def _invalidate(self) -> None:
        # The cached alignment and last result belong to the old inputs.
        self._variation = None
        self._last_design = None

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
        """The current inputs as a DesignInput."""
        return DesignInput(target=self._require_target(), panel=self._panel)

    def inspect(self) -> InspectResult:
        """Preview the loaded inputs."""
        return InspectResult(design_input=self.design_input())

    def _panel_variation(self, mrna: BioSequence) -> PanelVariation:
        # Align once; window evaluation depends on the params and is cheap.
        if self._variation is None:
            self._variation = PanelVariation.from_alignment(
                align_star(mrna, self._panel)
            )
        return self._variation

    def _offtarget_self(
        self, mrna: BioSequence, params: DesignParams, offtarget_db
    ) -> SelfTranscripts:
        """The design's own transcripts: declared accessions, the target's header
        accession, and (if enabled) those resolved by homology to the target,
        panel and self FASTA. See docs/design-rules.md (What the design owns).
        """
        declared = set(params.offtarget_self_accessions)
        detected = _detect_accession(mrna)
        if detected:
            declared.add(detected)

        resolved: frozenset[str] = frozenset()
        if params.resolve_offtarget_self:
            sequences = [mrna.seq]
            if self._panel is not None:
                sequences += [member.seq for member in self._panel]
            if params.offtarget_self_fasta is not None:
                sequences += [
                    seq.seq
                    for seq in _load_sequences_file(params.offtarget_self_fasta)
                ]
            resolved = resolve_self_transcripts(
                BowtieSearcher(
                    offtarget_db,
                    max_mismatch=0,
                    bowtie_path=params.bowtie_path,
                ),
                sequences,
                length=params.length,
                min_coverage=params.offtarget_self_coverage,
            )
        return SelfTranscripts.of(
            declared=declared,
            resolved=resolved,
            min_coverage=(
                params.offtarget_self_coverage
                if params.resolve_offtarget_self
                else None
            ),
        )

    def design(self, params: DesignParams | None = None) -> DesignResult:
        """Enumerate candidates and apply the gates.

        With a panel, the allele gate is ANDed with the siRNA gates. Raises
        HushError for a protein target or unknown allele id, ValueError for an
        incoherent request.
        """
        params = params or DesignParams()
        design_input = self.design_input()
        mrna = design_input.require_nucleotide_target()
        candidates = list(
            enumerate_candidates(mrna, length=params.length, overhang=params.overhang)
        )

        # Allele gate: windows and candidates share target coordinates, so a
        # candidate is looked up by start position.
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

        # Explicit --offtarget-db wins, else the built reference. With neither,
        # the run is unreferenced and offtarget_skipped is set.
        offtarget_searcher = None
        self_transcripts = None
        offtarget_db = params.offtarget_db
        if offtarget_db is None and index_exists(default_index_prefix()):
            offtarget_db = default_index_prefix()
        offtarget_active = params.use_offtarget_filter and offtarget_db is not None
        offtarget_skipped = params.use_offtarget_filter and offtarget_db is None
        if offtarget_active:
            self_transcripts = self._offtarget_self(mrna, params, offtarget_db)
            backend = BowtieSearcher(
                offtarget_db,
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
            offtarget_self_transcripts=self_transcripts,
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
        """Convert candidates to shRNA constructs.

        Defaults to the selected candidates of the last design, in that order.
        Withheld constructs are returned on the result.
        """
        options = options or ShrnaOptions()
        if candidates is None:
            if self._last_design is None:
                raise NothingToConvertError(
                    "Nothing to convert; run design first, or pass candidates."
                )
            candidates = [e.candidate for e in self._last_design.selected]
        constructs = [convert_candidate(c, options) for c in candidates]
        return apply_convert_filters(constructs, options)
