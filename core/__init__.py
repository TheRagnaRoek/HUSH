"""HUSH core: shared domain models and the siRNA design engine.

This package holds the foundational pieces the rest of the pipeline builds on:

* core.sequence - the sequence domain model (BioSequence, SequencePanel,
  DesignInput) used everywhere downstream.
* core.inputs   - the Input Layer: FASTA/PDB parsing into that model.

Later stages (siRNA design, MSA/allele logic, shRNA conversion) consume these
types so that every module speaks the same vocabulary.
"""

from core.sequence import BioSequence, DesignInput, MolType, SequencePanel

__all__ = ["BioSequence", "DesignInput", "MolType", "SequencePanel"]
