"""The HUSH sequence domain model.

Every stage of the pipeline speaks in terms of the types defined here, so this
module is intentionally dependency-light (pure standard library, no Biopython)
and heavily validated.

Core types:

* MolType       - DNA / RNA / protein / unknown.
* BioSequence   - an immutable, validated biological sequence with the
  nucleotide operations (transcription, reverse-complement, GC content, ...)
  the siRNA/shRNA stages rely on.
* SequencePanel - a named set of related sequences (e.g. an HLA allele family)
  fed to the MSA/allele module.
* DesignInput   - the bundle handed from the Input Layer to the rest of the
  pipeline: a target plus an optional panel.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum

from core.errors import MolTypeError, SequenceValidationError

# --------------------------------------------------------------------------- #
# Alphabets
# --------------------------------------------------------------------------- #

DNA_BASES = frozenset("ACGT")
RNA_BASES = frozenset("ACGU")
# IUPAC nucleotide ambiguity codes (shared between DNA and RNA).
IUPAC_NUC_AMBIGUITY = frozenset("RYSWKMBDHVN")
PROTEIN_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
# Extended protein symbols: ambiguous (B,Z,J,X), selenocysteine (U),
# pyrrolysine (O), and stop (*).
PROTEIN_EXTRA = frozenset("XBZJUO*")
GAP_CHARS = frozenset("-.")

# IUPAC complement (uppercase, DNA convention: T rather than U).
_DNA_COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
    "-": "-", ".": ".",
}


class MolType(str, Enum):
    """The molecule type of a BioSequence."""

    DNA = "dna"
    RNA = "rna"
    PROTEIN = "protein"
    UNKNOWN = "unknown"

    @property
    def is_nucleotide(self) -> bool:
        return self in (MolType.DNA, MolType.RNA)


def _clean(seq: str) -> str:
    """Uppercase and strip all whitespace from a raw sequence string."""
    return "".join(seq.split()).upper()


def detect_moltype(seq: str) -> MolType:
    """Best-effort guess of a sequence's molecule type from its residues.

    Heuristic: if 90% or more of the (gapless) residues are valid nucleotide
    letters and the sequence does not mix 'T' and 'U', it is treated as nucleic
    acid - RNA when a 'U' is present, otherwise DNA. Everything else is protein.

    Short, ambiguous sequences (e.g. 'ACGA', valid as both DNA and a peptide)
    resolve to DNA under this rule; pass an explicit moltype to BioSequence when
    the type is known.
    """
    s = "".join(c for c in _clean(seq) if c not in GAP_CHARS)
    if not s:
        return MolType.UNKNOWN

    nucleotide_letters = DNA_BASES | RNA_BASES | IUPAC_NUC_AMBIGUITY
    nuc_fraction = sum(c in nucleotide_letters for c in s) / len(s)
    has_t = "T" in s
    has_u = "U" in s

    if nuc_fraction >= 0.90 and not (has_t and has_u):
        return MolType.RNA if has_u else MolType.DNA
    return MolType.PROTEIN


def _allowed_alphabet(moltype: MolType) -> frozenset[str] | None:
    """Return the set of legal residues for moltype, or None to skip checks."""
    if moltype is MolType.DNA:
        return DNA_BASES | IUPAC_NUC_AMBIGUITY | GAP_CHARS
    if moltype is MolType.RNA:
        return RNA_BASES | IUPAC_NUC_AMBIGUITY | GAP_CHARS
    if moltype is MolType.PROTEIN:
        return PROTEIN_AA | PROTEIN_EXTRA | GAP_CHARS
    return None  # UNKNOWN - no validation


@dataclass(frozen=True)
class BioSequence:
    """An immutable, validated biological sequence.

    The 'seq' string is normalised to uppercase with whitespace stripped on
    construction, and validated against the alphabet for its 'moltype'. When
    'moltype' is left as MolType.UNKNOWN it is auto-detected via detect_moltype.

    Parameters:
        seq:         Raw residue string.
        id:          Short identifier (e.g. a FASTA accession or 'HLA-A*02:01').
        description: Free-text description (e.g. the remainder of a FASTA
                     header line).
        moltype:     Molecule type; auto-detected when UNKNOWN.
        metadata:    Provenance and any extra key/value data (source file,
                     chain id, ...).
    """

    seq: str
    id: str = ""
    description: str = ""
    moltype: MolType = MolType.UNKNOWN
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        cleaned = _clean(self.seq)
        object.__setattr__(self, "seq", cleaned)

        moltype = self.moltype
        if moltype is MolType.UNKNOWN:
            moltype = detect_moltype(cleaned)
            object.__setattr__(self, "moltype", moltype)

        self._validate(cleaned, moltype)

    # -- validation -------------------------------------------------------- #

    @staticmethod
    def _validate(seq: str, moltype: MolType) -> None:
        allowed = _allowed_alphabet(moltype)
        if allowed is None:
            return
        bad = sorted({c for c in seq if c not in allowed})
        if bad:
            raise SequenceValidationError(
                f"Invalid {moltype.value} residue(s) {bad} in sequence "
                f"{seq!r:.40}"
            )

    # -- basic accessors --------------------------------------------------- #

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def is_nucleotide(self) -> bool:
        return self.moltype.is_nucleotide

    @property
    def is_protein(self) -> bool:
        return self.moltype is MolType.PROTEIN

    @property
    def label(self) -> str:
        """A human-friendly one-line label: id plus description if present."""
        return f"{self.id} {self.description}".strip() or "<unnamed>"

    def subseq(self, start: int, end: int) -> "BioSequence":
        """Return the half-open slice [start, end) as a new sequence.

        The new id is suffixed with ':start-end' (1-based, inclusive) so the
        provenance of a design window stays legible downstream.
        """
        if start < 0 or end > len(self.seq) or start >= end:
            raise ValueError(
                f"Invalid slice [{start}, {end}) for length {len(self.seq)}"
            )
        return replace(
            self,
            seq=self.seq[start:end],
            id=f"{self.id}:{start + 1}-{end}" if self.id else self.id,
        )

    # -- nucleotide operations -------------------------------------------- #

    def _require_nucleotide(self, op: str) -> None:
        if not self.is_nucleotide:
            raise MolTypeError(
                f"{op} requires a nucleotide sequence, but "
                f"{self.label!r} is {self.moltype.value}."
            )

    @property
    def gc_content(self) -> float:
        """Fraction of G/C among unambiguous bases (0.0 for empty/none)."""
        self._require_nucleotide("gc_content")
        gc = sum(c in "GC" for c in self.seq)
        acgt = sum(c in "ACGTU" for c in self.seq)
        return gc / acgt if acgt else 0.0

    def transcribe(self) -> "BioSequence":
        """DNA to RNA (coding-strand transcription: 'T' becomes 'U')."""
        if self.moltype is MolType.RNA:
            return self
        self._require_nucleotide("transcribe")
        return replace(self, seq=self.seq.replace("T", "U"), moltype=MolType.RNA)

    def reverse_transcribe(self) -> "BioSequence":
        """RNA to DNA ('U' becomes 'T')."""
        if self.moltype is MolType.DNA:
            return self
        self._require_nucleotide("reverse_transcribe")
        return replace(self, seq=self.seq.replace("U", "T"), moltype=MolType.DNA)

    def as_rna(self) -> "BioSequence":
        """Return this sequence as RNA (identity for RNA, transcribe for DNA).

        Convenience for the siRNA stage, which operates on the target mRNA.
        """
        return self.transcribe()

    def complement(self) -> "BioSequence":
        """IUPAC complement (not reversed). Works on DNA or RNA."""
        self._require_nucleotide("complement")
        as_dna = self.seq.replace("U", "T")
        comp = "".join(_DNA_COMPLEMENT[c] for c in as_dna)
        if self.moltype is MolType.RNA:
            comp = comp.replace("T", "U")
        return replace(self, seq=comp)

    def reverse_complement(self) -> "BioSequence":
        """IUPAC reverse-complement. Works on DNA or RNA."""
        rc = self.complement()
        return replace(rc, seq=rc.seq[::-1])

    def __str__(self) -> str:
        return self.seq

    def __repr__(self) -> str:
        preview = self.seq if len(self.seq) <= 30 else f"{self.seq[:27]}..."
        return (
            f"BioSequence(id={self.id!r}, moltype={self.moltype.value}, "
            f"len={len(self.seq)}, seq={preview!r})"
        )


@dataclass
class SequencePanel:
    """A named collection of related sequences (e.g. an HLA allele family).

    Fed to the MSA/allele module, which aligns the members and derives
    conserved vs. polymorphic windows for broad or allele-specific targeting.
    """

    name: str
    sequences: list[BioSequence] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)

    def __getitem__(self, index: int) -> BioSequence:
        return self.sequences[index]

    @property
    def ids(self) -> list[str]:
        return [s.id for s in self.sequences]

    @property
    def moltypes(self) -> set[MolType]:
        return {s.moltype for s in self.sequences}

    @property
    def is_uniform(self) -> bool:
        """True when every member shares one molecule type (required for MSA)."""
        return len(self.moltypes) <= 1

    def get(self, seq_id: str) -> BioSequence | None:
        """Return the first member whose id equals seq_id, or None."""
        return next((s for s in self.sequences if s.id == seq_id), None)


@dataclass
class DesignInput:
    """The Input Layer's output: a target sequence and an optional panel.

    This is the hand-off object consumed by every downstream stage. The panel is
    optional - allele-aware targeting is only engaged when the user supplies a
    related sequence set.
    """

    target: BioSequence
    panel: SequencePanel | None = None

    @property
    def has_panel(self) -> bool:
        return self.panel is not None and len(self.panel) > 0

    def require_nucleotide_target(self) -> BioSequence:
        """Return the target, guaranteeing it is a nucleotide sequence.

        siRNA/shRNA design needs the target mRNA/CDS. A protein target (e.g.
        parsed from a PDB structure) cannot be silenced directly and must first
        be resolved to a coding sequence - this raises a clear error in that
        case rather than letting a protein leak into the design stage.
        """
        if not self.target.is_nucleotide:
            raise MolTypeError(
                f"Target {self.target.label!r} is {self.target.moltype.value}; "
                "siRNA/shRNA design requires a nucleotide (mRNA/CDS) target. "
                "Provide the coding sequence, or resolve the protein to its "
                "transcript first."
            )
        return self.target
