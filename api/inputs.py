"""Parse request text into the sequence model; failures raise HushError."""

from __future__ import annotations

from core.errors import InputParseError
from core.inputs import read_fasta_text
from core.sequence import BioSequence, SequencePanel


def parse_target(text: str, seq_id: str = "target") -> BioSequence:
    """Target from FASTA text (first record) or a bare sequence."""
    text = text.strip()
    if not text:
        raise InputParseError("The target sequence is empty.")
    if text.startswith(">"):
        return read_fasta_text(text)[0]
    return BioSequence(text, id=seq_id)


def parse_panel(text: str, name: str = "panel") -> SequencePanel:
    """Panel from FASTA text; members are addressed by id."""
    if not text or not text.strip():
        raise InputParseError("The panel is empty.")
    return SequencePanel(name=name, sequences=read_fasta_text(text))
