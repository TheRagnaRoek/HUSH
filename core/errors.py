"""Exception types for HUSH."""

from __future__ import annotations


class HushError(Exception):
    """Base class for all HUSH errors."""


class SequenceValidationError(HushError):
    """A sequence contains residues invalid for its molecule type."""


class MolTypeError(HushError):
    """An operation is invalid for the sequence's molecule type."""


class InputParseError(HushError):
    """An input file could not be parsed into sequences."""


class NothingToConvertError(HushError):
    """A conversion was requested before any design exists to convert."""


class OffTargetError(HushError):
    """The transcriptome off-target search could not be run."""


class ReferenceBuildError(HushError):
    """The transcriptome reference index could not be built."""
