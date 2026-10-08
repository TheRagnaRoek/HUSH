"""The stem-loop scaffold loop registry.

A stem-loop shRNA closes its duplex with a short unpaired loop. Which loop is
used measurably affects knockdown, so HUSH offers a small set of established,
literature-backed loops rather than one hard-coded choice. Each entry records
its sequence (RNA, 5'->3') and the source it comes from.

The default is Jensen et al. (2012)'s library-selected UGUGCUU: it was chosen
specifically because it supports robust knockdown with the least
stem-sequence dependence, which is the safest property for an automated tool
that does not know the stem in advance.

References:
* Jensen, S.M.R., Schmitz, A., Pedersen, F.S., Kjems, J., Bramsen, J.B. (2012).
  Functional selection of shRNA loops from randomized retroviral libraries. PLoS
  ONE 7(8), e43095.
* Brummelkamp, T.R., Bernards, R. & Agami, R. (2002). A system for stable
  expression of short interfering RNAs in mammalian cells (pSUPER). Science 296,
  550-553.
* Moffat, J., Grueneberg, D.A., Yang, X., Kim, S.Y., Kloepfer, A.M.,
  Hinkle, G. et al. (2006). A lentiviral RNAi library for human and mouse genes
  applied to an arrayed viral high-content screen (pLKO.1 / TRC). Cell 124(6),
  1283-1298. doi:10.1016/j.cell.2006.01.040
"""

from __future__ import annotations

from dataclasses import dataclass

from core.sequence import BioSequence, MolType


@dataclass(frozen=True)
class LoopScaffold:
    """One named stem-loop loop, stored as RNA 5'->3'.

    Attributes:
        name:     short registry key (e.g. 'jensen').
        sequence: the loop as a validated RNA BioSequence, 5'->3'.
        citation: the source the loop is taken from.
    """

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


# The registry, keyed by short name. Sequences are RNA (the transcribed loop).
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

# Least stem-dependent of the set - the safe default for an automated converter.
DEFAULT_LOOP = "jensen"


def get_loop(name: str) -> LoopScaffold:
    """Return the LoopScaffold for a registry name, or raise ValueError."""
    key = name.strip().lower()
    if key not in LOOPS:
        choices = ", ".join(sorted(LOOPS))
        raise ValueError(f"unknown loop {name!r}; choose from {choices}.")
    return LOOPS[key]
