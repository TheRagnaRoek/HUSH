"""HUSH service layer: HushSession, parameter objects and results.

The CLI and the API both run the pipeline through this package. See
docs/architecture.md.
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
