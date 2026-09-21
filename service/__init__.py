"""HUSH service layer: the orchestration seam between the library and its frontends.

The pure pipeline lives in core/msa/shrna; the CLI and the future HTTP API are
thin renderers. This package is what sits between them - a HushSession that holds
a loaded target and panel, runs the pipeline, and returns structured results
(never printed text, never JSON). Both frontends drive HUSH through it, so the
pipeline is wired together in exactly one place.

* params  - DesignParams / AlleleTargetingSpec: parsed, framework-agnostic requests.
* results - InspectResult / DesignResult / ConvertResult: lossless computed state.
* session - HushSession: load inputs, inspect / design / convert.
"""

from service.params import AlleleTargetingSpec, DesignParams
from service.results import ConvertResult, DesignResult, InspectResult
from service.session import HushSession, resolve_targeting, sort_selection

__all__ = [
    "AlleleTargetingSpec",
    "DesignParams",
    "InspectResult",
    "DesignResult",
    "ConvertResult",
    "HushSession",
    "resolve_targeting",
    "sort_selection",
]
