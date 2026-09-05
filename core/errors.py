"""Shared exception types for HUSH."""

from __future__ import annotations


class HushError(Exception):
    """Base class for all HUSH-specific errors."""


class SequenceValidationError(HushError):
    """A sequence contains residues not valid for its molecule type."""


class MolTypeError(HushError):
    """An operation was requested that is invalid for a sequence's molecule type.

    e.g. asking for the RNA form of a protein sequence, or running siRNA
    design against a target that is not a nucleotide sequence.
    """


class InputParseError(HushError):
    """A FASTA/PDB (or other) input file could not be parsed into sequences."""


class OffTargetError(HushError):
    """A transcriptome off-target search could not be run.

    e.g. the bowtie binary is not on PATH, the reference index is missing, or
    the aligner exited with an error.
    """
