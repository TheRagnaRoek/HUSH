"""Stage 1 input: read FASTA and PDB/mmCIF files into the sequence model.

Protein sequences from structures are labelled as protein; the nucleotide
requirement is enforced by DesignInput.require_nucleotide_target.
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path

from Bio import SeqIO

from core.errors import InputParseError
from core.sequence import BioSequence, DesignInput, SequencePanel

_FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn", ".fas", ".mpfa"}
_PDB_SUFFIXES = {".pdb", ".ent"}
_CIF_SUFFIXES = {".cif", ".mmcif"}


def detect_format(path: Path) -> str:
    """Input format ("fasta", "pdb" or "cif") from the file suffix."""
    suffix = path.suffix.lower()
    if suffix in _FASTA_SUFFIXES:
        return "fasta"
    if suffix in _PDB_SUFFIXES:
        return "pdb"
    if suffix in _CIF_SUFFIXES:
        return "cif"
    raise InputParseError(
        f"Cannot infer input format from {path.name!r}. Supported extensions: "
        f"FASTA ({', '.join(sorted(_FASTA_SUFFIXES))}), "
        f"PDB ({', '.join(sorted(_PDB_SUFFIXES))}), "
        f"mmCIF ({', '.join(sorted(_CIF_SUFFIXES))}). "
        "Pass an explicit format to override."
    )


def _record_to_bioseq(record, source: Path, source_format: str, **extra) -> BioSequence:
    """Convert a Biopython SeqRecord to a BioSequence."""
    seq_id = record.id or ""
    # Biopython's description starts with the id; drop it.
    description = record.description or ""
    if description.startswith(seq_id):
        description = description[len(seq_id):].strip()
    metadata = {"source": str(source), "source_format": source_format, **extra}
    return BioSequence(
        seq=str(record.seq),
        id=seq_id,
        description=description,
        metadata=metadata,
    )


def read_fasta(path: str | Path) -> list[BioSequence]:
    """Read all records of a FASTA file."""
    path = Path(path)
    if not path.is_file():
        raise InputParseError(f"FASTA file not found: {path}")
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except Exception as exc:  # pragma: no cover - Biopython raises many types
        raise InputParseError(f"Failed to parse FASTA {path}: {exc}") from exc
    if not records:
        raise InputParseError(f"No sequences found in FASTA file: {path}")
    return [_record_to_bioseq(r, path, "fasta") for r in records]


def read_fasta_text(text: str, source: str = "<uploaded>") -> list[BioSequence]:
    """Parse FASTA records from a string (used by the API)."""
    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
    except Exception as exc:  # pragma: no cover - Biopython raises many types
        raise InputParseError(f"Failed to parse FASTA text: {exc}") from exc
    if not records:
        raise InputParseError("No sequences found in the FASTA text.")
    return [_record_to_bioseq(r, Path(source), "fasta") for r in records]


def read_pdb(path: str | Path, fmt: str | None = None) -> list[BioSequence]:
    """One sequence per chain of a PDB/mmCIF file.

    Uses SEQRES records, falling back to ATOM residues when there are none.
    """
    path = Path(path)
    if not path.is_file():
        raise InputParseError(f"Structure file not found: {path}")

    fmt = fmt or detect_format(path)
    seqres_format = "cif-seqres" if fmt == "cif" else "pdb-seqres"
    atom_format = "cif-atom" if fmt == "cif" else "pdb-atom"

    records = _try_structure_parse(path, seqres_format)
    if not records:
        # No SEQRES: use the modelled (ATOM) residues.
        records = _try_structure_parse(path, atom_format)
    if not records:
        raise InputParseError(
            f"No sequences could be extracted from structure file: {path}"
        )

    seqs = []
    for record in records:
        chain = record.annotations.get("chain", "")
        seqs.append(
            _record_to_bioseq(record, path, fmt, chain=chain)
        )
    return seqs


def _try_structure_parse(path: Path, biopython_format: str) -> list:
    try:
        return list(SeqIO.parse(str(path), biopython_format))
    except Exception:
        return []


def load_sequences(path: str | Path, fmt: str | None = None) -> list[BioSequence]:
    """Load all sequences from a file; the format defaults to the suffix."""
    path = Path(path)
    fmt = fmt or detect_format(path)
    if fmt == "fasta":
        return read_fasta(path)
    if fmt in ("pdb", "cif"):
        return read_pdb(path, fmt=fmt)
    raise InputParseError(f"Unsupported input format: {fmt!r}")


def load_target(path: str | Path, fmt: str | None = None) -> BioSequence:
    """Load the first sequence of a file as the target.

    With several records, the count is stored in metadata["selected_from"].
    """
    seqs = load_sequences(path, fmt=fmt)
    target = seqs[0]
    if len(seqs) > 1:
        object_metadata = dict(target.metadata)
        object_metadata["selected_from"] = len(seqs)
        # BioSequence is frozen: rebuild with the new metadata.
        from dataclasses import replace

        target = replace(target, metadata=object_metadata)
    return target


def load_panel(
    path: str | Path, name: str | None = None, fmt: str | None = None
) -> SequencePanel:
    """Load all sequences of a file as a panel."""
    seqs = load_sequences(path, fmt=fmt)
    panel_name = name or Path(path).stem
    return SequencePanel(name=panel_name, sequences=seqs)


def build_design_input(
    target_path: str | Path,
    panel_path: str | Path | None = None,
    *,
    target_fmt: str | None = None,
    panel_fmt: str | None = None,
    panel_name: str | None = None,
) -> DesignInput:
    """Load a target and optional panel into a DesignInput."""
    target = load_target(target_path, fmt=target_fmt)
    panel = (
        load_panel(panel_path, name=panel_name, fmt=panel_fmt)
        if panel_path is not None
        else None
    )
    return DesignInput(target=target, panel=panel)
