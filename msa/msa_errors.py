"""Exception types for the allele stage; all derive from HushError."""

from __future__ import annotations

from core.errors import HushError


class MsaError(HushError):
    """Base class for allele stage errors."""


class PanelError(MsaError):
    """A panel is unsuitable: missing or duplicate ids, mixed types, too few members."""


class TargetingError(MsaError):
    """A targeting request is malformed or does not fit its panel."""


class AlignmentError(MsaError):
    """A panel could not be aligned to the target, or the alignment is unusable."""
