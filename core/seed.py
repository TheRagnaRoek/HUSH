"""Seed-duplex Tm off-target filter (siDirect 2.0).

An siRNA also behaves like a microRNA: a match to its seed (guide positions 2-8)
in a transcript's 3' UTR is enough to silence it. Off-target strength tracks the
seed-target duplex stability (Ui-Tei et al. 2008), so siDirect filters out seeds
that pair too stably - Tm < 21.5 C, on BOTH strands (either can end up loaded).
This module reproduces that filter exactly.

The Tm here is deliberately NOT core.thermo's: it uses the Freier (1986)
parameters and siDirect's own conditions (100 uM strands, 100 mM Na+), because
the 21.5 C threshold is only meaningful against that exact combination. The
filter is a benchmark, believed over-conservative; treat it as fixed.

Primary source: Naito et al. 2009 (siDirect 2.0). Design rationale, the exact
formula, the validation against siDirect's published figures, and the
over-conservatism critique are in docs/design-rules.md; the dated design
notes behind it are on the archive/design-notes branch. Read those before
changing any constant here.
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

# --------------------------------------------------------------------------- #
# Seed region
# --------------------------------------------------------------------------- #

# The seed is positions 2-8 (1-based, inclusive) from the 5' end of the strand.
# Position 1 is excluded because it is buried in the Mid-domain binding pocket
# of Argonaute and does not pair with the target.
SEED_START = 2
SEED_END = 8
SEED_LENGTH = SEED_END - SEED_START + 1  # 7 nt


def seed_region(strand: BioSequence | str) -> str:
    """Return the seed heptamer (positions 2-8, 5'->3') of a strand, as RNA.

    Only the extracted heptamer is validated as unambiguous RNA. An ambiguity
    code elsewhere in the strand is irrelevant to this filter, which never looks
    outside the seed, so it is not grounds for rejecting the candidate here.

    Raises:
        ValueError: if the strand is too short to contain a full seed.
        SequenceValidationError: if the seed itself is not unambiguous RNA.
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


# --------------------------------------------------------------------------- #
# siDirect Tm constants
# --------------------------------------------------------------------------- #
#
# Every constant below is taken verbatim from the siDirect 2.0 equation. They
# are named rather than inlined so a future recalibration is a visible edit, but
# they should not be tuned casually: the 21.5 C threshold is only meaningful in
# combination with exactly these values.

# Helix initiation, applied to the entropy sum. Freier et al. (1986).
SIDIRECT_HELIX_INIT_DS = -10.8

# The paper's rounded gas constant, cal/(mol*K). Deliberately not core.thermo's
# more precise 1.9872 - see the module docstring.
SIDIRECT_R = 1.987

# Total strand concentration, mol/L (100 uM in the paper).
SIDIRECT_TOTAL_STRAND_CONC_M = 100e-6

# Sodium concentration the correction is evaluated at, mol/L (100 mM).
SIDIRECT_SODIUM_M = 0.1

# Coefficient of the log10[Na+] salt correction term.
SALT_CORRECTION_COEFFICIENT = 16.6

# Seeds melting below this are treated as off-target-free (Ui-Tei et al. 2008).
# The comparison is strict: the papers say "below 21.5 C".
SEED_TM_THRESHOLD_C = 21.5


def seed_duplex_tm(
    seq: BioSequence | str,
    total_strand_conc: float = SIDIRECT_TOTAL_STRAND_CONC_M,
    sodium: float = SIDIRECT_SODIUM_M,
) -> float:
    """Melting temperature (degrees C) of a seed duplex, per siDirect 2.0.

    'seq' is one strand of the duplex, 5'->3' - normally a seed heptamer from
    seed_region(). The target strand is its perfect complement, so only one
    strand is needed.

    The defaults reproduce siDirect exactly. The concentration and salt
    arguments are exposed for sensitivity checks, not for routine use: move them
    and the 21.5 C threshold no longer means what the papers measured.

    Raises:
        SequenceValidationError: if the sequence is not unambiguous RNA/DNA.
        ValueError: if the sequence is shorter than 2 nt, or a concentration is
            not positive.
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


# --------------------------------------------------------------------------- #
# Result
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class SeedTmResult:
    """Seed-duplex stability of both strands of one candidate.

    A candidate passes only when BOTH strands' seeds melt below the threshold:
    either strand can end up RISC-loaded, and whichever does brings its own seed
    with it.

    Fields:
        guide_seed:     the guide (antisense) heptamer, positions 2-8.
        guide_tm:       its seed-duplex Tm in degrees C.
        passenger_seed: the passenger (sense) heptamer, positions 2-8.
        passenger_tm:   its seed-duplex Tm in degrees C.
        threshold:      the cut-off applied, in degrees C.
    """

    guide_seed: str
    guide_tm: float
    passenger_seed: str
    passenger_tm: float
    threshold: float

    @property
    def max_tm(self) -> float:
        """The worse (higher, more off-target-prone) of the two strands' Tm."""
        return max(self.guide_tm, self.passenger_tm)

    @property
    def guide_passes(self) -> bool:
        return self.guide_tm < self.threshold

    @property
    def passenger_passes(self) -> bool:
        return self.passenger_tm < self.threshold

    @property
    def passes(self) -> bool:
        """True when both strands melt below the threshold."""
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
    """Evaluate the seed-duplex Tm filter for both strands of a candidate.

    Args:
        candidate: the siRNA candidate to filter.
        threshold: Tm cut-off in degrees C. The default 21.5 is siDirect's;
            the paper also reports 15 and 10 as progressively stricter settings
            that trade candidate yield for specificity.

    Raises:
        ValueError: if either strand is too short to contain a seed (< 8 nt).
    """
    guide_seed = seed_region(candidate.antisense)
    passenger_seed = seed_region(candidate.sense)
    return SeedTmResult(
        guide_seed=guide_seed,
        guide_tm=seed_duplex_tm(guide_seed),
        passenger_seed=passenger_seed,
        passenger_tm=seed_duplex_tm(passenger_seed),
        threshold=threshold,
    )
