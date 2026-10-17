"""Seed-duplex Tm off-target filter (Naito et al. 2009).

The 21.5 C threshold is valid only against the paper's exact constants, so they
are pinned. Seed definition, equation and validation: docs/design-rules.md
(Seed-duplex off-target Tm).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from core.sequence import BioSequence
from core.sirna import SiRNACandidate
from core.thermo import (
    NN_PARAMS_FREIER_1986,
    ZERO_CELSIUS_K,
    as_wc_rna,
    stacking_terms,
)

SEED_START = 2
SEED_END = 8
SEED_LENGTH = SEED_END - SEED_START + 1  # 7 nt


def seed_region(strand: BioSequence | str) -> str:
    """Seed heptamer (positions 2-8, 5'->3') of a strand, as RNA.

    Ambiguity codes outside the seed are ignored. Raises ValueError when the
    strand is shorter than 8 nt, SequenceValidationError for an ambiguous seed.
    """
    if isinstance(strand, BioSequence):
        text = strand.as_rna().seq
    else:
        text = "".join(str(strand).split()).upper()
    if len(text) < SEED_END:
        raise ValueError(
            f"A seed spans positions {SEED_START}-{SEED_END}, so the strand needs "
            f"at least {SEED_END} nt; got {len(text)}."
        )
    return as_wc_rna(text[SEED_START - 1 : SEED_END])


# Constants of Naito et al. 2009, taken verbatim. Do not tune.
SIDIRECT_HELIX_INIT_DS = -10.8

# The paper's rounded gas constant, not core.thermo's R_CAL.
SIDIRECT_R = 1.987

SIDIRECT_TOTAL_STRAND_CONC_M = 100e-6
SIDIRECT_SODIUM_M = 0.1
SALT_CORRECTION_COEFFICIENT = 16.6

# A seed passes when its Tm is strictly below this.
SEED_TM_THRESHOLD_C = 21.5


def seed_duplex_tm(
    seq: BioSequence | str,
    total_strand_conc: float = SIDIRECT_TOTAL_STRAND_CONC_M,
    sodium: float = SIDIRECT_SODIUM_M,
) -> float:
    """Tm in degrees C of a seed duplex (Naito et al. 2009).

        Tm = 1000*dH / (S_init + dS + R*ln(Ct/4)) - 273.15 + 16.6*log10[Na+]

    dH and dS are the Freier stacking sums of the strand. The concentration and
    sodium arguments are for sensitivity checks only. Raises ValueError for a
    strand under 2 nt or a non-positive concentration.
    """
    if total_strand_conc <= 0:
        raise ValueError(
            f"Total strand concentration must be positive, got {total_strand_conc}."
        )
    if sodium <= 0:
        raise ValueError(f"Sodium concentration must be positive, got {sodium}.")

    rna = as_wc_rna(seq)
    if len(rna) < 2:
        raise ValueError(
            f"A seed duplex needs at least 2 nt to have a nearest-neighbor "
            f"step; got {len(rna)}."
        )
    dh, ds = stacking_terms(rna, params=NN_PARAMS_FREIER_1986)

    denominator = (
        SIDIRECT_HELIX_INIT_DS
        + ds
        + SIDIRECT_R * math.log(total_strand_conc / 4.0)
    )
    kelvin = (1000.0 * dh) / denominator
    return (
        kelvin
        - ZERO_CELSIUS_K
        + SALT_CORRECTION_COEFFICIENT * math.log10(sodium)
    )


@dataclass(frozen=True)
class SeedTmResult:
    """Seed Tm of both strands of one candidate; passes only if both are below threshold."""

    guide_seed: str
    guide_tm: float
    passenger_seed: str
    passenger_tm: float
    threshold: float

    @property
    def max_tm(self) -> float:
        """The higher of the two strands' Tm."""
        return max(self.guide_tm, self.passenger_tm)

    @property
    def guide_passes(self) -> bool:
        return self.guide_tm < self.threshold

    @property
    def passenger_passes(self) -> bool:
        return self.passenger_tm < self.threshold

    @property
    def passes(self) -> bool:
        """True when both strands are below the threshold."""
        return self.max_tm < self.threshold

    def __repr__(self) -> str:
        return (
            f"SeedTmResult(passes={self.passes}, "
            f"guide={self.guide_seed} {self.guide_tm:.1f}C, "
            f"passenger={self.passenger_seed} {self.passenger_tm:.1f}C, "
            f"threshold={self.threshold:.1f}C)"
        )


def seed_tm(
    candidate: SiRNACandidate, threshold: float = SEED_TM_THRESHOLD_C
) -> SeedTmResult:
    """Seed Tm of both strands of a candidate. Raises ValueError for strands under 8 nt."""
    guide_seed = seed_region(candidate.antisense)
    passenger_seed = seed_region(candidate.sense)
    return SeedTmResult(
        guide_seed=guide_seed,
        guide_tm=seed_duplex_tm(guide_seed),
        passenger_seed=passenger_seed,
        passenger_tm=seed_duplex_tm(passenger_seed),
        threshold=threshold,
    )
