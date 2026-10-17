"""Loop registry for stem-loop shRNAs: sequence and source of each loop.

Selection rationale and references: docs/shrna-conversion.md (Loops).
"""

from __future__ import annotations

from dataclasses import dataclass

from core.sequence import BioSequence, MolType


@dataclass(frozen=True)
class LoopScaffold:
    """A named loop, stored as RNA 5'->3' with its citation."""

    name: str
    sequence: BioSequence
    citation: str

    @property
    def seq(self) -> str:
        """The loop nucleotides, 5'->3'."""
        return self.sequence.seq

    def __len__(self) -> int:
        return len(self.sequence)


def _loop(name: str, seq: str, citation: str) -> LoopScaffold:
    return LoopScaffold(
        name=name,
        sequence=BioSequence(seq, id=f"loop:{name}", moltype=MolType.RNA),
        citation=citation,
    )


LOOPS: dict[str, LoopScaffold] = {
    "jensen": _loop(
        "jensen", "UGUGCUU", "Jensen et al. 2012, PLoS ONE 7(8):e43095"
    ),
    "psuper": _loop(
        "psuper", "UUCAAGAGA", "Brummelkamp et al. 2002, Science 296:550 (pSUPER)"
    ),
    "trc": _loop(
        "trc", "CUCGAG", "Moffat et al. 2006, Cell 124:1283 (pLKO.1 / TRC)"
    ),
}

DEFAULT_LOOP = "jensen"


def get_loop(name: str) -> LoopScaffold:
    """The loop for a registry name; ValueError if unknown."""
    key = name.strip().lower()
    if key not in LOOPS:
        choices = ", ".join(sorted(LOOPS))
        raise ValueError(f"unknown loop {name!r}; choose from {choices}.")
    return LOOPS[key]
