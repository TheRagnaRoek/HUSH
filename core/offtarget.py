"""Transcriptome homology off-target filter.

Searches a reference transcriptome for near-complementary sites of each
candidate on both strands, reports the off-target profile and gates on the
nearest hit. Design transcripts are excluded as SelfTranscripts. Rationale,
the owned-set sources and the coverage default: docs/design-rules.md
(Transcriptome off-target).

Two searchers share one protocol: BruteForceSearcher (pure Python, for small
references and as the test oracle) and BowtieSearcher (bowtie1 in -v mode,
reached through subprocess so importing this module never needs it). Both are
ungapped, Hamming-distance and forward-orientation.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Protocol, runtime_checkable

from core.errors import OffTargetError
from core.sirna import SiRNACandidate

# Reject when the nearest off-target has fewer mismatches than this.
DEFAULT_OFFTARGET_MIN_MISMATCH = 2

# Mismatch depth searched; bowtie1 -v mode allows at most 3.
DEFAULT_OFFTARGET_MAX_MISMATCH = 3
BOWTIE_MAX_V = 3

# Fraction of a design sequence's k-mers a transcript must carry to count as self.
DEFAULT_SELF_MIN_COVERAGE = 0.10


def _to_dna(seq: str) -> str:
    """Uppercase and read U as T; references and bowtie use DNA letters."""
    return seq.upper().replace("U", "T")


def _base_accession(accession: str) -> str:
    """Accession without version suffix, lowercased (NM_1.3 -> nm_1)."""
    return accession.split(".", 1)[0].strip().lower()


def _accession_set(values: Iterable[str] | str | None) -> frozenset[str]:
    """Normalise accessions to base form; a bare string counts as one."""
    if values is None:
        return frozenset()
    if isinstance(values, str):
        values = [values]
    return frozenset(_base_accession(v) for v in values if v and v.strip())


@dataclass(frozen=True)
class SelfTranscripts:
    """Reference transcripts the design owns, excluded from off-target hits.

    'declared' holds accessions named by the caller, 'resolved' those found by
    resolve_self_transcripts. Both are version-insensitive.
    """

    declared: frozenset[str] = frozenset()
    resolved: frozenset[str] = frozenset()
    min_coverage: float | None = None

    @staticmethod
    def of(
        declared: Iterable[str] | str | None = None,
        resolved: Iterable[str] | str | None = None,
        min_coverage: float | None = None,
    ) -> "SelfTranscripts":
        """Build from raw accessions."""
        return SelfTranscripts(
            declared=_accession_set(declared),
            resolved=_accession_set(resolved),
            min_coverage=min_coverage,
        )

    @property
    def accessions(self) -> frozenset[str]:
        return self.declared | self.resolved

    def __len__(self) -> int:
        return len(self.accessions)

    def __contains__(self, accession: str) -> bool:
        return _base_accession(accession) in self.accessions

    def __repr__(self) -> str:
        return (
            f"SelfTranscripts(n={len(self)}, declared={len(self.declared)}, "
            f"resolved={len(self.resolved)}, min_coverage={self.min_coverage})"
        )


@dataclass(frozen=True)
class OffTargetHit:
    """One near-match: transcript accession, 0-based position and mismatch count."""

    accession: str
    position: int
    mismatches: int


@dataclass(frozen=True)
class OffTargetProfile:
    """Off-target hits of one candidate, per strand, with self transcripts removed.

    counts[i] is the number of hits with exactly i mismatches (i up to
    max_mismatch); min_mismatch is the nearest hit, or None for no hit. The
    candidate passes when each strand's nearest hit is at least 'threshold'
    mismatches away or absent.
    """

    guide_min_mismatch: int | None
    guide_counts: tuple[int, ...]
    passenger_min_mismatch: int | None
    passenger_counts: tuple[int, ...]
    max_mismatch: int
    threshold: int

    @property
    def min_mismatch(self) -> int | None:
        """The nearer of the two strands' nearest hits, or None."""
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
        """True when the nearest hit on both strands clears the threshold."""
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


@runtime_checkable
class OffTargetSearcher(Protocol):
    """Finds every ungapped forward-strand match of a read within max_mismatch.

    The read is DNA-spelled (see _to_dna). Self exclusion is not done here.
    """

    max_mismatch: int

    def search(self, read: str) -> list[OffTargetHit]: ...


class BruteForceSearcher:
    """Exhaustive Hamming search over an in-memory reference.

    O(references x length) per read: for small references and tests, and the
    oracle for BowtieSearcher.
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
        """Build from a FASTA file; each record id is an accession."""
        # Lazy import keeps this module free of the Biopython-backed input layer.
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
    """Mismatches of read against ref[start:], or None when above cap."""
    mm = 0
    for i, base in enumerate(read):
        if base != ref[start + i]:
            mm += 1
            if mm > cap:
                return None
    return mm


class PrecomputedSearcher:
    """Searcher over a precomputed {read: hits} map."""

    def __init__(
        self, hits_by_read: dict[str, list[OffTargetHit]], max_mismatch: int
    ) -> None:
        self.max_mismatch = max_mismatch
        self._hits = {_to_dna(read): hits for read, hits in hits_by_read.items()}

    def search(self, read: str) -> list[OffTargetHit]:
        return self._hits.get(_to_dna(read), [])


class BowtieSearcher:
    """Off-target search with bowtie1: -v N -a --norc over a prebuilt index.

    Build the index with 'hush reference build' (core.reference). precompute()
    aligns a batch of reads in one bowtie run.
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
                f"bowtie executable {bowtie_path!r} not found on PATH. The "
                "HUSH container ships bowtie1; on a host install, install "
                "bowtie1 from your package manager or bioconda."
            )
        if not Path(f"{self.index_prefix}.1.ebwt").is_file():
            raise OffTargetError(
                f"No bowtie index at {self.index_prefix!r} (expected "
                f"{self.index_prefix}.1.ebwt). Build one with "
                "'hush reference build'."
            )

    def search(self, read: str) -> list[OffTargetHit]:
        return self.precompute([read])._hits.get(_to_dna(read), [])

    def precompute(self, reads: Iterable[str]) -> PrecomputedSearcher:
        """Align all reads in one bowtie run."""
        unique = sorted({_to_dna(r) for r in reads})
        hits_by_read: dict[str, list[OffTargetHit]] = {r: [] for r in unique}
        if not unique:
            return PrecomputedSearcher(hits_by_read, self.max_mismatch)

        with tempfile.TemporaryDirectory() as tmp:
            reads_fa = Path(tmp) / "reads.fa"
            # The sequence doubles as the FASTA name so hits regroup by read.
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
            "-a",              # all alignments, not just the best
            "--norc",          # forward strand only
            "--quiet",
            "-f", self.index_prefix, str(reads_fa),
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError as exc:
            raise OffTargetError(f"Could not run bowtie: {exc}") from exc
        # bowtie may exit non-zero when no read aligns; that leaves stdout empty.
        if proc.returncode != 0 and not proc.stdout:
            detail = proc.stderr.strip() or f"exit code {proc.returncode}"
            raise OffTargetError(f"bowtie failed: {detail}")
        return proc.stdout


def _parse_bowtie_line(line: str) -> tuple[str, OffTargetHit]:
    """Parse one default-format bowtie line into (read name, hit).

    The mismatch count is the number of comma-separated descriptors in the
    last column.
    """
    fields = line.rstrip("\n").split("\t")
    name = fields[0]
    accession = fields[2]
    position = int(fields[3])
    descriptors = fields[7] if len(fields) > 7 else ""
    mismatches = descriptors.count(",") + 1 if descriptors else 0
    return name, OffTargetHit(accession, position, mismatches)


def _tiles(sequence: str, length: int) -> set[str]:
    """Distinct k-mers of the given length, DNA-spelled."""
    seq = _to_dna(sequence)
    return {seq[i : i + length] for i in range(len(seq) - length + 1)}


def resolve_self_transcripts(
    searcher: OffTargetSearcher,
    sequences: Iterable[str],
    length: int,
    min_coverage: float = DEFAULT_SELF_MIN_COVERAGE,
) -> frozenset[str]:
    """Transcripts that carry at least min_coverage of any input's k-mers exactly.

    Each input is tiled into k-mers of 'length' and judged on its own, so a
    long target does not dilute a short panel member.
    """
    if not 0 < min_coverage <= 1:
        raise ValueError(f"min_coverage must be in (0, 1], got {min_coverage}.")
    if length < 1:
        raise ValueError(f"length must be >= 1, got {length}.")

    tiled = [tiles for tiles in (_tiles(s, length) for s in sequences) if tiles]
    if not tiled:
        return frozenset()
    pooled = sorted(set().union(*tiled))

    # One batched pass when the backend supports it, else per-read search.
    lookup = (
        searcher.precompute(pooled) if hasattr(searcher, "precompute") else searcher
    )
    exact = {
        tile: {
            _base_accession(h.accession)
            for h in lookup.search(tile)
            if h.mismatches == 0
        }
        for tile in pooled
    }

    resolved: set[str] = set()
    for tiles in tiled:
        covered: Counter[str] = Counter()
        for tile in tiles:
            covered.update(exact[tile])
        floor = min_coverage * len(tiles)
        resolved.update(acc for acc, n in covered.items() if n >= floor)
    return frozenset(resolved)


def _as_self(
    spec: "SelfTranscripts | Iterable[str] | str | None",
) -> SelfTranscripts | None:
    """Accept a SelfTranscripts, an accession or an iterable of accessions."""
    if spec is None or isinstance(spec, SelfTranscripts):
        return spec
    return SelfTranscripts.of(declared=spec)


def _exclude_self(
    hits: list[OffTargetHit], self_transcripts: SelfTranscripts | None
) -> list[OffTargetHit]:
    """Drop hits to the design's own transcripts."""
    if not self_transcripts:
        return hits
    return [h for h in hits if h.accession not in self_transcripts]


def evaluate_offtarget(
    candidate: SiRNACandidate,
    searcher: OffTargetSearcher,
    threshold: int = DEFAULT_OFFTARGET_MIN_MISMATCH,
    self_transcripts: "SelfTranscripts | Iterable[str] | str | None" = None,
) -> OffTargetProfile:
    """Off-target profile of a candidate, both strands, self transcripts excluded.

    Guide side: sense read; passenger side: antisense read. See
    docs/design-rules.md (Transcriptome off-target, Both strands are searched).
    """
    max_mm = searcher.max_mismatch
    own = _as_self(self_transcripts)
    guide_hits = _exclude_self(searcher.search(candidate.sense.seq), own)
    passenger_hits = _exclude_self(searcher.search(candidate.antisense.seq), own)
    return OffTargetProfile(
        guide_min_mismatch=_min_mismatch(guide_hits),
        guide_counts=_counts(guide_hits, max_mm),
        passenger_min_mismatch=_min_mismatch(passenger_hits),
        passenger_counts=_counts(passenger_hits, max_mm),
        max_mismatch=max_mm,
        threshold=threshold,
    )
