"""Turn request text into the HUSH domain model.

The browser client sends the target and panel as text in the JSON body, not as
files. A target may be a full FASTA record or a bare sequence (the client can
paste either); a panel must be FASTA, since allele targeting addresses its
members by id. Parsing failures raise HushError, which the app maps to 400.
"""

from __future__ import annotations

from core.errors import InputParseError
from core.inputs import read_fasta_text
from core.sequence import BioSequence, SequencePanel


def parse_target(text: str, seq_id: str = "target") -> BioSequence:
    """Parse a target from FASTA text or a bare sequence string.

    A leading '>' is taken as FASTA (the first record is the target, mirroring
    the CLI's load_target); anything else is treated as a raw sequence whose
    molecule type is auto-detected.
    """
    text = text.strip()
    if not text:
        raise InputParseError("The target sequence is empty.")
    if text.startswith(">"):
        return read_fasta_text(text)[0]
    return BioSequence(text, id=seq_id)


def parse_panel(text: str, name: str = "panel") -> SequencePanel:
    """Parse a related-sequence panel from FASTA text (members need ids)."""
    if not text or not text.strip():
        raise InputParseError("The panel is empty.")
    return SequencePanel(name=name, sequences=read_fasta_text(text))
