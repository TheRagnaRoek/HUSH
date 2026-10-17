"""Sequence model: MolType, BioSequence, SequencePanel and DesignInput.

Standard library only, so it stays easy to test.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum

from core.errors import MolTypeError, SequenceValidationError

DNA_BASES = frozenset("ACGT")
RNA_BASES = frozenset("ACGU")
IUPAC_NUC_AMBIGUITY = frozenset("RYSWKMBDHVN")
PROTEIN_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
# Ambiguous (B, Z, J, X), selenocysteine (U), pyrrolysine (O), stop (*).
PROTEIN_EXTRA = frozenset("XBZJUO*")
GAP_CHARS = frozenset("-.")

# Complement table in DNA letters; RNA converts T back to U.
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
    """Uppercase a sequence string and strip all whitespace."""
    return "".join(seq.split()).upper()


def detect_moltype(seq: str) -> MolType:
    """Guess the molecule type from the residues.

    Nucleotide when at least 90% of the gapless residues are nucleotide letters
    and T and U do not both occur (RNA if U is present, else DNA); otherwise
    protein. Short strings valid as both, such as 'ACGA', resolve to DNA.
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
    """Legal residues for the molecule type, or None for no validation."""
    if moltype is MolType.DNA:
        return DNA_BASES | IUPAC_NUC_AMBIGUITY | GAP_CHARS
    if moltype is MolType.RNA:
        return RNA_BASES | IUPAC_NUC_AMBIGUITY | GAP_CHARS
    if moltype is MolType.PROTEIN:
        return PROTEIN_AA | PROTEIN_EXTRA | GAP_CHARS
    return None  # UNKNOWN - no validation


@dataclass(frozen=True)
class BioSequence:
    """An immutable sequence, uppercased and validated against its moltype.

    A moltype of UNKNOWN is detected with detect_moltype. 'metadata' holds
    provenance such as the source file or chain id.
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
        """Id and description on one line."""
        return f"{self.id} {self.description}".strip() or "<unnamed>"

    def subseq(self, start: int, end: int) -> "BioSequence":
        """Half-open slice [start, end); the id gains a 1-based ':start-end' suffix."""
        if start < 0 or end > len(self.seq) or start >= end:
            raise ValueError(
                f"Invalid slice [{start}, {end}) for length {len(self.seq)}"
            )
        return replace(
            self,
            seq=self.seq[start:end],
            id=f"{self.id}:{start + 1}-{end}" if self.id else self.id,
        )

    def _require_nucleotide(self, op: str) -> None:
        if not self.is_nucleotide:
            raise MolTypeError(
                f"{op} requires a nucleotide sequence, but "
                f"{self.label!r} is {self.moltype.value}."
            )

    @property
    def gc_content(self) -> float:
        """Fraction of G/C among A/C/G/T/U bases; 0.0 when there are none."""
        self._require_nucleotide("gc_content")
        gc = sum(c in "GC" for c in self.seq)
        acgt = sum(c in "ACGTU" for c in self.seq)
        return gc / acgt if acgt else 0.0

    def transcribe(self) -> "BioSequence":
        """DNA to RNA (T becomes U); RNA is returned unchanged."""
        if self.moltype is MolType.RNA:
            return self
        self._require_nucleotide("transcribe")
        return replace(self, seq=self.seq.replace("T", "U"), moltype=MolType.RNA)

    def reverse_transcribe(self) -> "BioSequence":
        """RNA to DNA (U becomes T); DNA is returned unchanged."""
        if self.moltype is MolType.DNA:
            return self
        self._require_nucleotide("reverse_transcribe")
        return replace(self, seq=self.seq.replace("U", "T"), moltype=MolType.DNA)

    def as_rna(self) -> "BioSequence":
        """Alias for transcribe()."""
        return self.transcribe()

    def complement(self) -> "BioSequence":
        """IUPAC complement, not reversed."""
        self._require_nucleotide("complement")
        as_dna = self.seq.replace("U", "T")
        comp = "".join(_DNA_COMPLEMENT[c] for c in as_dna)
        if self.moltype is MolType.RNA:
            comp = comp.replace("T", "U")
        return replace(self, seq=comp)

    def reverse_complement(self) -> "BioSequence":
        """IUPAC reverse complement."""
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
    """A named set of related sequences, e.g. an allele family."""

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
        """True when all members share one molecule type."""
        return len(self.moltypes) <= 1

    def get(self, seq_id: str) -> BioSequence | None:
        """First member with this id, or None."""
        return next((s for s in self.sequences if s.id == seq_id), None)


@dataclass
class DesignInput:
    """Stage 1 output: a target and an optional panel.

    Allele targeting engages only when a panel is present.
    """

    target: BioSequence
    panel: SequencePanel | None = None

    @property
    def has_panel(self) -> bool:
        return self.panel is not None and len(self.panel) > 0

    def require_nucleotide_target(self) -> BioSequence:
        """The target, or MolTypeError if it is not nucleotide (e.g. from a PDB)."""
        if not self.target.is_nucleotide:
            raise MolTypeError(
                f"Target {self.target.label!r} is {self.target.moltype.value}; "
                "siRNA/shRNA design requires a nucleotide (mRNA/CDS) target. "
                "Provide the coding sequence, or resolve the protein to its "
                "transcript first."
            )
        return self.target
