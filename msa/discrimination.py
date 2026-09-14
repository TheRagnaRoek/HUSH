"""How well a single mismatch discriminates between two alleles (stage 3).

msa.windows answers whether a window can tell alleles apart at all. This module
answers how *well* - by placing each mismatch against published measurements of
which mismatches, at which guide positions, block silencing of the non-target
allele.

Two levels of evidence, used in that order:

* Huang et al. (2009) profiled every mismatch at every position and published a
  lookup keyed by the exact base change (Table 1). Being listed there is the
  most specific statement available: it is measured for this base pair at this
  position, not inferred from a general rule.
* Schwarz et al. (2006) established the general principles - purine:purine
  mismatches discriminate best, and mismatches 3' of the seed discriminate
  robustly while seed mismatches are overwhelmed by high siRNA concentration.
  These order the candidates Huang's table already endorses.

Above both sits the number of mismatches. Both papers address the case of a
*single* nucleotide difference, which is where placement matters. Two or more
differences disrupt pairing with the spared allele regardless of where they sit,
and Ohnishi et al. (2008) showed that adding a second mismatch - even inside the
seed - markedly improves allele discrimination.

Nothing here is a weighted score. Ranking is a tuple of published facts compared
lexicographically (see discrimination_key): each position in the tuple is one
finding, and their order is a stated priority rather than a set of weights that
would have to be calibrated.

Scope and deferrals:

* The *order* of the tuple's entries is only partly evidenced. Huang before
  Schwarz is well founded (specific beats general); region before position is
  founded on a difference in kind rather than degree. Position before class, and
  the count threshold above everything, are judgements - see the comments on
  DISCRIMINATION_ORDER_NOTE.
* Huang's table and Schwarz's regional finding disagree in places: A>C lists
  guide position 8, which lies inside the seed. Disagreements are reported via
  MismatchAssessment.sources_agree rather than resolved here.
* Schwarz lists P5 among the consistently discriminating positions although it
  lies inside the seed; the authors attribute strong seed discrimination in that
  scaffold to an unusually stable seed. The published set is reproduced intact
  and region is checked first, which places P5 in the seed tier.
* Within purine:purine there is further structure - G:G and A:G were more
  selective than A:A and G:A at P10, and A:G alone held up at 20 nM - which is
  not modelled.
* Guide position 1 is held in Argonaute's MID pocket and does not pair with the
  target, so a "mismatch" there is not one. It is reported as GuideRegion.ANCHOR
  and ranks with the seed rather than below it; no measurement separates them.
* IUPAC ambiguity codes have no row in Huang's table and no purine/pyrimidine
  identity. They yield None rather than a guess.

References:

* Huang, H. et al. (2009). Profiling of mismatch discrimination in RNAi enabled
  rational design of allele-specific siRNAs. Nucleic Acids Research 37(22),
  7560-7569. doi:10.1093/nar/gkp835
* Schwarz, D.S. et al. (2006). Designing siRNA that distinguish between genes
  that differ by a single nucleotide. PLoS Genetics 2(9), e140.
  doi:10.1371/journal.pgen.0020140
* Ohnishi, Y. et al. (2008). Enhancement of allele discrimination by
  introduction of nucleotide mismatches into siRNA in allele-specific gene
  silencing by RNAi. PLoS ONE 3(5), e2248. doi:10.1371/journal.pone.0002248
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence

from core.seed import SEED_END, SEED_START
from msa.msa_errors import MsaError

# --------------------------------------------------------------------------- #
# Alphabet
# --------------------------------------------------------------------------- #

PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CTU")

# Single-base complement, DNA convention. A local table rather than a trip
# through BioSequence: this is one character, and the surrounding code works in
# plain letters.
_COMPLEMENT = {"A": "T", "T": "A", "U": "A", "G": "C", "C": "G"}


def _to_dna(base: str) -> str:
    """Uppercase, with U folded to T, so RNA and DNA share one lookup key."""
    base = base.upper()
    return "T" if base == "U" else base


def _is_unambiguous(base: str) -> bool:
    return _to_dna(base) in {"A", "C", "G", "T"}


# --------------------------------------------------------------------------- #
# Level 1: Huang et al. (2009), Table 1
# --------------------------------------------------------------------------- #

# "Recommended positions for allele-specific siRNA design", reproduced verbatim.
#
# Read a row as: the spared allele carries the first base, the targeted allele
# the second. The siRNA guide therefore carries the complement of the second,
# and mismatches the first. The positions are the guide positions (1-based,
# 5'->3') at which that mismatch was least tolerated.
HUANG_RECOMMENDED_POSITIONS: dict[tuple[str, str], frozenset[int]] = {
    ("A", "C"): frozenset({8, 9, 10, 11, 16}),
    ("A", "G"): frozenset({10, 11, 12, 16}),
    ("A", "T"): frozenset({9, 10, 14, 15, 16}),
    ("C", "A"): frozenset({11, 12, 13, 14, 16}),
    ("C", "G"): frozenset({9, 10, 12, 13, 16}),
    ("C", "T"): frozenset({9, 10, 11, 12, 14}),
    ("G", "A"): frozenset({10, 13, 15, 16}),
    ("G", "C"): frozenset({10, 11, 12, 13}),
    ("G", "T"): frozenset({9, 10, 11, 12, 14}),
    ("T", "A"): frozenset({9, 11, 12, 15}),
    ("T", "C"): frozenset({10, 11, 12, 13}),
    ("T", "G"): frozenset({10, 11, 12, 13, 16}),
}


# --------------------------------------------------------------------------- #
# Level 2: Schwarz et al. (2006)
# --------------------------------------------------------------------------- #

# "P5, P9, P10, P12, P13, and P16 mismatches consistently discriminated between
#  two alleles, as tested with three different experimental strategies."
CONSISTENT_ACROSS_ASSAYS = frozenset({5, 9, 10, 12, 13, 16})

# "most apparent for mismatches at positions 10 and 16, but also P9, P12, P13,
#  and P14"
MOST_APPARENT = frozenset({10, 16})
ALSO_ROBUST = frozenset({9, 12, 13, 14})

# Derived for ranking: the two singled out, and everything else the paper names.
STRONGEST_POSITIONS = MOST_APPARENT
CONSISTENT_POSITIONS = (CONSISTENT_ACROSS_ASSAYS | ALSO_ROBUST) - STRONGEST_POSITIONS


# --------------------------------------------------------------------------- #
# Mismatch count
# --------------------------------------------------------------------------- #

# From how many differences a window counts as robust irrespective of placement.
# Two, because Ohnishi et al. (2008) found that a second mismatch - including one
# placed inside the seed - markedly improves allele discrimination, and because
# the placement literature addresses the single-difference case by construction.
MULTI_MISMATCH = 2


# --------------------------------------------------------------------------- #
# Classification
# --------------------------------------------------------------------------- #


class MismatchClass(str, Enum):
    """The pair formed by the guide residue and the spared allele's residue.

    Named as the papers report them, siRNA nucleotide first. Purine:purine
    disrupts the A-form helix most and discriminates best; the other three are
    accommodated in the helix and discriminate only moderately.
    """

    PURINE_PURINE = "purine:purine"
    PURINE_PYRIMIDINE = "purine:pyrimidine"
    PYRIMIDINE_PURINE = "pyrimidine:purine"
    PYRIMIDINE_PYRIMIDINE = "pyrimidine:pyrimidine"


class GuideRegion(str, Enum):
    """Where on the guide a position falls.

    ANCHOR:    position 1, held in Argonaute's MID pocket and not paired.
    SEED:      positions SEED_START..SEED_END, responsible for target binding.
               Mismatches here are overwhelmed by high siRNA concentration.
    POST_SEED: 3' of the seed, where pairing is read for catalysis rather than
               binding. Mismatches here discriminate robustly.
    """

    ANCHOR = "anchor"
    SEED = "seed"
    POST_SEED = "post_seed"


def guide_region(guide_position: int) -> GuideRegion:
    """The region a 1-based guide position falls in."""
    if guide_position < 1:
        raise MsaError(f"Guide positions are 1-based; got {guide_position}.")
    if guide_position == 1:
        return GuideRegion.ANCHOR
    if SEED_START <= guide_position <= SEED_END:
        return GuideRegion.SEED
    return GuideRegion.POST_SEED


def classify_mismatch(guide_base: str, other_base: str) -> MismatchClass | None:
    """The mismatch class of one guide/target pair, or None if either is ambiguous."""
    if not (_is_unambiguous(guide_base) and _is_unambiguous(other_base)):
        return None
    guide_purine = _to_dna(guide_base) in PURINES
    other_purine = _to_dna(other_base) in PURINES
    if guide_purine and other_purine:
        return MismatchClass.PURINE_PURINE
    if guide_purine:
        return MismatchClass.PURINE_PYRIMIDINE
    if other_purine:
        return MismatchClass.PYRIMIDINE_PURINE
    return MismatchClass.PYRIMIDINE_PYRIMIDINE


def huang_recommends(
    spared_base: str, target_base: str, guide_position: int
) -> bool | None:
    """Whether Huang's Table 1 lists this position for this base change.

    None when either base is ambiguous or the two are equal - neither is a row in
    the table, and neither should be guessed at.
    """
    if not (_is_unambiguous(spared_base) and _is_unambiguous(target_base)):
        return None
    key = (_to_dna(spared_base), _to_dna(target_base))
    positions = HUANG_RECOMMENDED_POSITIONS.get(key)
    if positions is None:  # identical bases: not a mismatch, not a table row
        return None
    return guide_position in positions


# --------------------------------------------------------------------------- #
# One mismatch, assessed
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class MismatchAssessment:
    """One difference between the target and a spared allele, placed in evidence.

    Everything here is measured or looked up; the ranking is derived from it by
    discrimination_key rather than stored, so the ordering can be changed without
    touching the data.

    Parameters:
        window_position: 0-based position inside the candidate window.
        guide_position:  1-based position on the guide strand, 5'->3', as both
                         papers count. The guide is the reverse complement of the
                         window, so the two run in opposite directions.
        target_base:     The base the siRNA is designed against.
        spared_base:     The base the allele to be spared carries.
        guide_base:      The guide residue, the complement of target_base.
        huang_recommended: Level 1. None when the bases are ambiguous.
        mismatch_class:  Level 2, class. None when the bases are ambiguous.
        region:          Level 2, region.
    """

    window_position: int
    guide_position: int
    target_base: str
    spared_base: str
    guide_base: str
    huang_recommended: bool | None
    mismatch_class: MismatchClass | None
    region: GuideRegion

    @property
    def mutation_type(self) -> str:
        """The Huang Table 1 row this pair belongs to, e.g. 'A>C'."""
        return f"{_to_dna(self.spared_base)}>{_to_dna(self.target_base)}"

    @property
    def sources_agree(self) -> bool:
        """Whether Huang's table and Schwarz's regional finding point the same way.

        False for a position Huang lists that lies inside the seed (A>C position
        8), or one it omits that lies 3' of the seed. Reported rather than
        resolved: a disagreement between two studies is information, not an error.
        """
        if self.huang_recommended is None:
            return False
        return self.huang_recommended == (self.region is GuideRegion.POST_SEED)

    def __repr__(self) -> str:
        cls = "?" if self.mismatch_class is None else self.mismatch_class.value
        return (
            f"MismatchAssessment(P{self.guide_position}, {self.mutation_type}, "
            f"{cls}, {self.region.value}, huang={self.huang_recommended})"
        )


def guide_position_1based(window_position: int, window_length: int) -> int:
    """The guide position a window position pairs with, counted from 1.

    Both source papers number guide positions from 1 at the 5' end (P1..P19), so
    this module does too. msa.windows.guide_position expresses the same flip
    0-based for internal use; the two differ by exactly one and a test holds them
    in step.

    The guide is the reverse complement of the window, so the two run in opposite
    directions: window position 0 pairs with the guide's 3' end, and the guide's
    5' end pairs with the last window position.

    Raises:
        MsaError: When window_position lies outside the window.
    """
    if not 0 <= window_position < window_length:
        raise MsaError(
            f"Window position {window_position} outside 0..{window_length - 1}"
        )
    return window_length - window_position


def assess_mismatch(
    window_position: int, window_length: int, target_base: str, spared_base: str
) -> MismatchAssessment:
    """Place one window-relative difference in the published evidence.

    Args:
        window_position: 0-based position inside the window.
        window_length:   Length of the window, needed for the coordinate flip.
        target_base:     The base at that position on the targeted allele.
        spared_base:     The base the spared allele carries there.

    Raises:
        MsaError: When window_position lies outside the window.
    """
    position = guide_position_1based(window_position, window_length)
    guide_base = _COMPLEMENT.get(_to_dna(target_base), target_base.upper())
    return MismatchAssessment(
        window_position=window_position,
        guide_position=position,
        target_base=target_base.upper(),
        spared_base=spared_base.upper(),
        guide_base=guide_base,
        huang_recommended=huang_recommends(spared_base, target_base, position),
        mismatch_class=classify_mismatch(guide_base, spared_base),
        region=guide_region(position),
    )


# --------------------------------------------------------------------------- #
# Ranking
# --------------------------------------------------------------------------- #

DISCRIMINATION_ORDER_NOTE = """\
The entries of discrimination_key, in order, and what each rests on:

1. count >= MULTI_MISMATCH - Ohnishi et al. 2008; also, the placement literature
   addresses the single-difference case by construction, so placement matters
   less once there are several. JUDGEMENT: that this outranks a single
   well-placed mismatch is not measured.
2. huang_recommended - measured for this exact base pair at this exact position.
   Specific beats general, so it sits above Schwarz's rules.
3. region is POST_SEED - Schwarz: a seed mismatch blocks binding and is
   overwhelmed by concentration; one 3' of the seed blocks catalysis and is not.
   A difference in kind, hence high placement.
4. position in STRONGEST_POSITIONS - Schwarz: "particularly P10 and P16"; P16
   gave an IC50 shift of at least 20-fold in 7 of 10 pairs.
5. position in CONSISTENT_POSITIONS - Schwarz's remaining named positions.
6. class is PURINE_PURINE - Schwarz's general principle. JUDGEMENT: placed below
   position because the paper's own design advice is positional, and because
   Huang's key is already base-pair specific. The cross-axis comparison - a
   purine:purine at P12 against a milder class at P16 - is not measured.
7. count - a tiebreaker within an otherwise equal rank, not a weight.

Entries 1 and 6 are the two that a reader may reasonably want reordered. Doing
so means moving a line in discrimination_key; nothing else changes.
"""


def _single_key(assessment: MismatchAssessment) -> tuple[bool, bool, bool, bool, bool]:
    """The published facts about one mismatch, ordered. Higher sorts better."""
    return (
        assessment.huang_recommended is True,
        assessment.region is GuideRegion.POST_SEED,
        assessment.guide_position in STRONGEST_POSITIONS,
        assessment.guide_position in CONSISTENT_POSITIONS,
        assessment.mismatch_class is MismatchClass.PURINE_PURINE,
    )


def best_mismatch(
    assessments: Sequence[MismatchAssessment],
) -> MismatchAssessment | None:
    """The most discriminating of several mismatches against one allele."""
    if not assessments:
        return None
    return max(assessments, key=_single_key)


def discrimination_key(
    assessments: Sequence[MismatchAssessment],
) -> tuple[bool, bool, bool, bool, bool, bool, int]:
    """Ranking key for all mismatches against one spared allele.

    Compared lexicographically; higher sorts better. See
    DISCRIMINATION_ORDER_NOTE for what each entry rests on.

    An empty sequence - no difference at all, so no discrimination - sorts below
    everything else.
    """
    count = len(assessments)
    if count == 0:
        return (False, False, False, False, False, False, 0)
    return (count >= MULTI_MISMATCH, *_single_key(best_mismatch(assessments)), count)
