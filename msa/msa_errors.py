"""Exception types for the MSA / allele module (pipeline stage 3).

Kept local to the module so the stage owns its own failure vocabulary, but every
type derives from core.errors.HushError: callers that catch HushError to handle
"anything HUSH raised" keep working, and the stage stays catchable as a unit via
MsaError.

The dependency runs msa -> core, the same direction as the rest of the module.
"""

from __future__ import annotations

from core.errors import HushError


class MsaError(HushError):
    """Base class for every error raised by the MSA / allele stage."""


class PanelError(MsaError):
    """A sequence panel is unsuitable for allele analysis.

    e.g. members without an id, duplicate ids, mixed molecule types, or too few
    sequences to compare.
    """


class TargetingError(MsaError):
    """An allele targeting request is malformed or does not fit its panel.

    e.g. an empty target set, an allele required to be both hit and spared, or a
    requested allele that is not in the panel.
    """


class AlignmentError(MsaError):
    """A panel could not be aligned to the target, or the alignment is unusable.

    Reserved for the alignment unit (msa.alignment).
    """
