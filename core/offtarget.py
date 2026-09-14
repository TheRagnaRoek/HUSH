"""Transcriptome homology off-target filter.

An siRNA guide silences any transcript its sequence is complementary to, not only
its intended target. siDirect BLASTs each guide against the human transcriptome
and drops candidates whose nearest off-target is a near-perfect match; this module
is HUSH's open, offline equivalent of that search. For each candidate it looks for
near-complementary sites in a reference transcriptome, reports the off-target
profile (nearest-match mismatch count and hit counts by mismatch level, on both
strands, with the intended transcript excluded), and exposes an optional gate that
rejects a candidate whose nearest off-target sits below a mismatch threshold.

This is distinct from core.seed, which is a *seed-mediated* off-target proxy (a
sequence-intrinsic Tm). Here we search for direct full-length homology.

Two searchers sit behind one protocol, sharing this module's result model, the
on-target exclusion and the gate:

* BruteForceSearcher - a pure-Python Hamming scan. Correct by construction, no
  external tool; suited to small references (a curated panel such as an HLA allele
  family, and the tests). It is the oracle the bowtie backend is checked against.
* BowtieSearcher - shells out to bowtie1 in `-v` mode (end-to-end alignment with at
  most N mismatches, all hits reported). bowtie1 is the open, offline counterpart of
  the GGGenome k-mismatch search siDirect relies on; it scales to the whole
  transcriptome. Only the binary is external - bowtie is reached via stdlib
  subprocess, so importing this module never needs it.

Both real engines are ungapped, Hamming-distance and forward-orientation, so their
results agree exactly. Design rationale (why bowtie1 rather than BLAST, the RefSeq
pinning, strand semantics and the allele-aware caveat) is in
docs/offtarget_design_note.md.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Protocol, runtime_checkable

from core.errors import OffTargetError
from core.sirna import SiRNACandidate

# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #

# Reject a candidate if its nearest off-target has fewer than this many mismatches
# (i.e. a 0/1-mismatch off-target exists). This matches siDirect on the PDCD1
# benchmark: siDirect keeps candidates whose nearest off-target is >= 2 mismatches
# and drops only the near-perfect (0/1) ones (see tools/compare_sidirect.py).
DEFAULT_OFFTARGET_MIN_MISMATCH = 2

# How many mismatches the search enumerates. bowtie1 -v mode allows at most 3, so
# this is capped there; 3 also gives siDirect-style counts at every level 0..3 and
# lets the reported nearest-mismatch match siDirect's column (which reaches 3).
DEFAULT_OFFTARGET_MAX_MISMATCH = 3
BOWTIE_MAX_V = 3


def _to_dna(seq: str) -> str:
    """Uppercase and spell a read as DNA (U -> T) for matching against a database.

    The reference transcriptome and bowtie both work in DNA letters; candidate
    strands are RNA. Matching is on the letters, so both sides must agree.
    """
    return seq.upper().replace("U", "T")


def _base_accession(accession: str) -> str:
    """An accession without its version suffix, lowercased (NM_1.3 -> nm_1)."""
    return accession.split(".", 1)[0].strip().lower()


def _same_accession(a: str, b: str) -> bool:
    """True when two accessions name the same transcript, version aside."""
    return bool(a) and bool(b) and _base_accession(a) == _base_accession(b)


# --------------------------------------------------------------------------- #
# Hit + profile model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class OffTargetHit:
    """One near-complementary site of a read in the reference transcriptome.

    Fields:
        accession:  the transcript the read aligned to (reference name).
        position:   0-based offset of the match on that transcript.
        mismatches: number of mismatched positions over the full read length.
    """

    accession: str
    position: int
    mismatches: int


@dataclass(frozen=True)
class OffTargetProfile:
    """The off-target picture for one candidate, both strands, self excluded.

    A candidate is loaded into RISC as either strand, so both are searched. The
    'guide' side is the near-matches of the guide's target site (its off-targets
    when the guide loads); the 'passenger' side is the same for the passenger.

    counts[i] is the number of retained off-target hits with exactly i mismatches,
    for i in 0..max_mismatch. min_mismatch is the nearest hit, or None if the read
    has no off-target within max_mismatch.

    A candidate passes only when BOTH strands' nearest off-target is at least
    'threshold' mismatches away (or absent).
    """

    guide_min_mismatch: int | None
    guide_counts: tuple[int, ...]
    passenger_min_mismatch: int | None
    passenger_counts: tuple[int, ...]
    max_mismatch: int
    threshold: int

    @property
    def min_mismatch(self) -> int | None:
        """The nearer of the two strands' nearest off-targets (None if neither)."""
        mins = [
            m
            for m in (self.guide_min_mismatch, self.passenger_min_mismatch)
            if m is not None
        ]
        return min(mins) if mins else None

    @property
    def guide_passes(self) -> bool:
        return (
            self.guide_min_mismatch is None
            or self.guide_min_mismatch >= self.threshold
        )

    @property
    def passenger_passes(self) -> bool:
        return (
            self.passenger_min_mismatch is None
            or self.passenger_min_mismatch >= self.threshold
        )

    @property
    def passes(self) -> bool:
        """True when the nearest off-target on either strand clears the threshold."""
        m = self.min_mismatch
        return m is None or m >= self.threshold

    def __repr__(self) -> str:
        return (
            f"OffTargetProfile(passes={self.passes}, "
            f"min_mismatch={self.min_mismatch}, threshold={self.threshold}, "
            f"guide_counts={self.guide_counts}, "
            f"passenger_counts={self.passenger_counts})"
        )


def _min_mismatch(hits: list[OffTargetHit]) -> int | None:
    return min((h.mismatches for h in hits), default=None)


def _counts(hits: list[OffTargetHit], max_mismatch: int) -> tuple[int, ...]:
    """Number of hits at each mismatch level 0..max_mismatch."""
    buckets = [0] * (max_mismatch + 1)
    for hit in hits:
        if 0 <= hit.mismatches <= max_mismatch:
            buckets[hit.mismatches] += 1
    return tuple(buckets)


# --------------------------------------------------------------------------- #
# Searcher protocol + backends
# --------------------------------------------------------------------------- #


@runtime_checkable
class OffTargetSearcher(Protocol):
    """Finds near-matches of a read in a reference, forward strand, up to N mismatches.

    'read' is a DNA-spelled string (see _to_dna). Implementations return every hit
    with at most 'max_mismatch' mismatches over the full read length, ungapped, on
    the reference's forward strand only (transcripts are single-stranded mRNA). The
    on-target transcript is not excluded here - that is the module's job, not the
    searcher's.
    """

    max_mismatch: int

    def search(self, read: str) -> list[OffTargetHit]: ...


class BruteForceSearcher:
    """Exhaustive pure-Python Hamming search over an in-memory reference.

    O(references x length) per read, so it is meant for small references - a
    curated allele panel, or the tests - not the whole transcriptome. Use
    BowtieSearcher for that. Serves as the reference implementation the bowtie
    backend is checked against.
    """

    def __init__(
        self,
        references: Iterable[tuple[str, str]],
        max_mismatch: int = DEFAULT_OFFTARGET_MAX_MISMATCH,
    ) -> None:
        if max_mismatch < 0:
            raise ValueError(f"max_mismatch must be >= 0, got {max_mismatch}.")
        self.max_mismatch = max_mismatch
        self._refs: list[tuple[str, str]] = [
            (accession, _to_dna(seq)) for accession, seq in references
        ]

    @classmethod
    def from_fasta(
        cls,
        path: str | Path,
        max_mismatch: int = DEFAULT_OFFTARGET_MAX_MISMATCH,
    ) -> "BruteForceSearcher":
        """Build a searcher from a FASTA reference (each record's id is its accession)."""
        # Imported lazily so this module stays importable without touching the
        # Biopython-backed input layer.
        from core.inputs import load_sequences

        refs = [(s.id, s.seq) for s in load_sequences(path)]
        return cls(refs, max_mismatch=max_mismatch)

    def search(self, read: str) -> list[OffTargetHit]:
        read = _to_dna(read)
        k = len(read)
        cap = self.max_mismatch
        hits: list[OffTargetHit] = []
        for accession, ref in self._refs:
            for start in range(len(ref) - k + 1):
                mm = _hamming_leq(read, ref, start, cap)
                if mm is not None:
                    hits.append(OffTargetHit(accession, start, mm))
        return hits


def _hamming_leq(read: str, ref: str, start: int, cap: int) -> int | None:
    """Mismatch count of 'read' against ref[start:], or None if it exceeds 'cap'."""
    mm = 0
    for i, base in enumerate(read):
        if base != ref[start + i]:
            mm += 1
            if mm > cap:
                return None
    return mm


class PrecomputedSearcher:
    """A searcher backed by an already-computed {read: hits} map (dict lookup).

    BowtieSearcher.precompute() returns one of these after a single batched bowtie
    run, so the per-candidate selection loop pays no per-read process cost.
    """

    def __init__(
        self, hits_by_read: dict[str, list[OffTargetHit]], max_mismatch: int
    ) -> None:
        self.max_mismatch = max_mismatch
        self._hits = {_to_dna(read): hits for read, hits in hits_by_read.items()}

    def search(self, read: str) -> list[OffTargetHit]:
        return self._hits.get(_to_dna(read), [])


class BowtieSearcher:
    """Off-target search via bowtie1 in `-v` mode over a prebuilt index.

    Runs `bowtie -v N -a --norc`: end-to-end alignment allowing at most N
    mismatches (bowtie1 caps N at 3), reporting all hits on the forward strand.
    build the index with tools/build_offtarget_db.py.

    search() runs one query; precompute() runs a whole batch in a single bowtie
    call and returns a PrecomputedSearcher - the path the CLI uses.
    """

    def __init__(
        self,
        index_prefix: str | Path,
        max_mismatch: int = DEFAULT_OFFTARGET_MAX_MISMATCH,
        bowtie_path: str = "bowtie",
    ) -> None:
        if not 0 <= max_mismatch <= BOWTIE_MAX_V:
            raise ValueError(
                f"bowtie -v mode allows 0-{BOWTIE_MAX_V} mismatches, got {max_mismatch}."
            )
        self.max_mismatch = max_mismatch
        self.index_prefix = str(index_prefix)
        self.bowtie_path = bowtie_path
        if shutil.which(bowtie_path) is None:
            raise OffTargetError(
                f"bowtie executable {bowtie_path!r} not found on PATH. Install "
                "bowtie1 (e.g. via bioconda or your package manager)."
            )
        if not Path(f"{self.index_prefix}.1.ebwt").is_file():
            raise OffTargetError(
                f"No bowtie index at {self.index_prefix!r} (expected "
                f"{self.index_prefix}.1.ebwt). Build one with "
                "tools/build_offtarget_db.py."
            )

    def search(self, read: str) -> list[OffTargetHit]:
        return self.precompute([read])._hits.get(_to_dna(read), [])

    def precompute(self, reads: Iterable[str]) -> PrecomputedSearcher:
        """Align every read in one bowtie run; return a lookup-backed searcher."""
        unique = sorted({_to_dna(r) for r in reads})
        hits_by_read: dict[str, list[OffTargetHit]] = {r: [] for r in unique}
        if not unique:
            return PrecomputedSearcher(hits_by_read, self.max_mismatch)

        with tempfile.TemporaryDirectory() as tmp:
            reads_fa = Path(tmp) / "reads.fa"
            # The read sequence doubles as its FASTA name so hits regroup by read.
            reads_fa.write_text(
                "".join(f">{r}\n{r}\n" for r in unique), encoding="ascii"
            )
            stdout = self._run_bowtie(reads_fa)

        for line in stdout.splitlines():
            name, hit = _parse_bowtie_line(line)
            if name in hits_by_read:
                hits_by_read[name].append(hit)
        return PrecomputedSearcher(hits_by_read, self.max_mismatch)

    def _run_bowtie(self, reads_fa: Path) -> str:
        cmd = [
            self.bowtie_path,
            "-v", str(self.max_mismatch),
            "-a",              # report all alignments, not just the best
            "--norc",          # forward strand only (transcripts are sense mRNA)
            "--quiet",
            "-f", self.index_prefix, str(reads_fa),
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise OffTargetError(f"Could not run bowtie: {exc}") from exc
        # A non-zero exit with no output is a real failure; bowtie can also exit
        # non-zero merely because no read aligned, in which case stdout is empty
        # and there is nothing to parse anyway.
        if proc.returncode != 0 and not proc.stdout:
            detail = proc.stderr.strip() or f"exit code {proc.returncode}"
            raise OffTargetError(f"bowtie failed: {detail}")
        return proc.stdout


def _parse_bowtie_line(line: str) -> tuple[str, OffTargetHit]:
    """Parse one default-format bowtie hit into (read_name, OffTargetHit).

    Columns: name, strand, reference, 0-based offset, read, quals, reserved,
    mismatch descriptors. The mismatch count is the number of comma-separated
    descriptors in the final field (empty means a perfect match).
    """
    fields = line.rstrip("\n").split("\t")
    name = fields[0]
    accession = fields[2]
    position = int(fields[3])
    descriptors = fields[7] if len(fields) > 7 else ""
    mismatches = descriptors.count(",") + 1 if descriptors else 0
    return name, OffTargetHit(accession, position, mismatches)


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def _exclude_self(
    hits: list[OffTargetHit], self_accession: str | None
) -> list[OffTargetHit]:
    """Drop hits to the candidate's own source transcript (siDirect's behaviour)."""
    if not self_accession:
        return hits
    return [h for h in hits if not _same_accession(h.accession, self_accession)]


def evaluate_offtarget(
    candidate: SiRNACandidate,
    searcher: OffTargetSearcher,
    threshold: int = DEFAULT_OFFTARGET_MIN_MISMATCH,
    self_accession: str | None = None,
) -> OffTargetProfile:
    """Off-target profile of a candidate against a reference, both strands, self out.

    The guide silences transcripts complementary to it, i.e. those carrying its
    target site (the sense strand); the passenger, if mis-loaded, silences those
    carrying the antisense strand. So the guide side searches the sense read and
    the passenger side the antisense read.
    """
    max_mm = searcher.max_mismatch
    guide_hits = _exclude_self(searcher.search(candidate.sense.seq), self_accession)
    passenger_hits = _exclude_self(
        searcher.search(candidate.antisense.seq), self_accession
    )
    return OffTargetProfile(
        guide_min_mismatch=_min_mismatch(guide_hits),
        guide_counts=_counts(guide_hits, max_mm),
        passenger_min_mismatch=_min_mismatch(passenger_hits),
        passenger_counts=_counts(passenger_hits, max_mm),
        max_mismatch=max_mm,
        threshold=threshold,
    )
