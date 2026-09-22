"""The Input Layer (pipeline stage 1).

Parses target and panel sequences from FASTA or PDB files into the HUSH domain
model (core.sequence). Biopython does the low-level file parsing; this module
adapts its records into validated BioSequence objects and assembles the
DesignInput hand-off.

FASTA is the primary path (a target mRNA/CDS, or a multi-record allele panel).
PDB is supported for convenience, but note: a PDB structure of a protein yields
an amino-acid sequence, which cannot be silenced directly - siRNA/shRNA design
needs the coding transcript. Such sequences are parsed and labelled faithfully
as protein; the requirement is enforced downstream (see
DesignInput.require_nucleotide_target).
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path

from Bio import SeqIO

from core.errors import InputParseError
from core.sequence import BioSequence, DesignInput, SequencePanel

# Map file extensions to a logical input format.
_FASTA_SUFFIXES = {".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn", ".fas", ".mpfa"}
_PDB_SUFFIXES = {".pdb", ".ent"}
_CIF_SUFFIXES = {".cif", ".mmcif"}


def detect_format(path: Path) -> str:
    """Infer the input format ("fasta", "pdb", "cif") from the suffix."""
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
    """Adapt a Biopython SeqRecord into a validated BioSequence."""
    seq_id = record.id or ""
    # Biopython's description includes the id as its first token; drop it.
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
    """Read every record from a FASTA file into a list of BioSequence."""
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
    """Parse FASTA records from an in-memory string (e.g. an HTTP upload).

    The string counterpart to read_fasta: same record adaptation, no file. Used
    by the API, where the target/panel arrive as request text rather than paths.
    """
    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
    except Exception as exc:  # pragma: no cover - Biopython raises many types
        raise InputParseError(f"Failed to parse FASTA text: {exc}") from exc
    if not records:
        raise InputParseError("No sequences found in the FASTA text.")
    return [_record_to_bioseq(r, Path(source), "fasta") for r in records]


def read_pdb(path: str | Path, fmt: str | None = None) -> list[BioSequence]:
    """Extract per-chain sequences from a PDB/mmCIF structure file.

    Prefers SEQRES records (the deposited sequence); falls back to the sequence
    observed in ATOM records when no SEQRES is present. Each chain becomes one
    BioSequence; molecule type is auto-detected, so nucleic-acid chains are
    labelled RNA/DNA and protein chains protein.
    """
    path = Path(path)
    if not path.is_file():
        raise InputParseError(f"Structure file not found: {path}")

    fmt = fmt or detect_format(path)
    seqres_format = "cif-seqres" if fmt == "cif" else "pdb-seqres"
    atom_format = "cif-atom" if fmt == "cif" else "pdb-atom"

    records = _try_structure_parse(path, seqres_format)
    if not records:
        # No SEQRES header - fall back to the modelled (ATOM) residues.
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
    """Load all sequences from path, dispatching on format (auto by suffix)."""
    path = Path(path)
    fmt = fmt or detect_format(path)
    if fmt == "fasta":
        return read_fasta(path)
    if fmt in ("pdb", "cif"):
        return read_pdb(path, fmt=fmt)
    raise InputParseError(f"Unsupported input format: {fmt!r}")


def load_target(path: str | Path, fmt: str | None = None) -> BioSequence:
    """Load a single target sequence.

    If the file holds multiple records (e.g. a multi-chain PDB or multi-FASTA),
    the first is returned and its selection is recorded in metadata under
    'selected_from' so the choice is auditable. Use load_panel when you mean to
    load a set.
    """
    seqs = load_sequences(path, fmt=fmt)
    target = seqs[0]
    if len(seqs) > 1:
        object_metadata = dict(target.metadata)
        object_metadata["selected_from"] = len(seqs)
        # BioSequence is frozen; rebuild with the annotated metadata.
        from dataclasses import replace

        target = replace(target, metadata=object_metadata)
    return target


def load_panel(
    path: str | Path, name: str | None = None, fmt: str | None = None
) -> SequencePanel:
    """Load a related sequence set (e.g. an HLA allele family) as a panel."""
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
    """Assemble the DesignInput handed to the rest of the pipeline."""
    target = load_target(target_path, fmt=target_fmt)
    panel = (
        load_panel(panel_path, name=panel_name, fmt=panel_fmt)
        if panel_path is not None
        else None
    )
    return DesignInput(target=target, panel=panel)
