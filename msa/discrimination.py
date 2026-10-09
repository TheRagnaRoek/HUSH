"""How well a mismatch's *placement* discriminates between two alleles (stage 3).

msa.windows answers whether a window can tell alleles apart at all. This module
answers how well the difference is placed - by putting each mismatch against
published measurements of which mismatches, at which guide positions, block
silencing of the non-target allele.

Three sources, read together:

* Huang et al. (2009) profiled every mismatch at every position and published a
  lookup keyed by the exact base change (Table 1). Being listed there is the
  most specific statement available: it is measured for this base pair at this
  position, not inferred from a general rule.
* Schwarz et al. (2006) established the general principles - purine:purine
  mismatches discriminate best, and mismatches 3' of the seed discriminate
  robustly while seed mismatches only block binding and are overwhelmed by high
  siRNA concentration. They name P10, P13 and P16 as the positions that held up
  across their assays.
* Du et al. (2005) mutated the *target site* rather than the siRNA, at all 57
  single-nucleotide permutations, and reported a low-tolerance region: the
  positions where a mismatch abolishes most silencing. Converted to guide
  coordinates it is P9-P15.

The three converge on a low-tolerance window of P9-P16, which is the backbone of
the ranking. Nothing here is a weighted score: the key is a tuple of published
facts compared lexicographically (see discrimination_key), each entry one
finding, their order a stated priority rather than weights to be calibrated.

What this module deliberately does NOT do:

* It does not rank by the *number* of differences. All three sources measure
  single mismatches, so none of them says whether several poorly placed
  differences outweigh one well-placed difference. Any threshold would assert
  something the literature does not support. The count is reported next to the
  ranking (MULTI_MISMATCH, WindowEvaluation.min_spare_mismatches) and left to
  the user; see docs/msa_ranking_note.md for the full argument.
* It does not claim a highly ranked candidate will be silent on the spared
  allele. Du found 73% of singly mismatched targets still significantly
  silenced: a difference is necessary for discrimination, not sufficient.
  Ranking never removes a candidate - filtering is a separate, earlier stage.

Scope and deferrals:

* The *order* of the tuple's entries is only partly evidenced. Huang first is
  well founded (specific beats general). Window before single positions, and
  position before class, are judgements - see DISCRIMINATION_ORDER_NOTE.
* Huang's table and Schwarz's regional finding disagree in places: A>C lists
  guide position 8, which lies inside the seed. Disagreements are reported via
  MismatchAssessment.sources_agree rather than resolved here.
* The sources disagree about P16: Schwarz's best-replicated position falls in
  Du's *medium*-tolerance region. It is kept, on the strength of Schwarz's
  ten-pair generality test and of Huang listing it in 7 of 12 rows.
* Schwarz also saw strong discrimination at P3, P4 and P5, but from a single
  siRNA whose seed pairing they describe as the most thermodynamically stable of
  their set, and whose strand-loading analysis shows it loaded predominantly its
  *sense* strand. Those seed positions are not used as a criterion.
* Within purine:purine there is further structure - G:G and A:G were more
  selective than A:A and G:A at P10, and A:G alone held up at 20 nM - which is
  not modelled. Nor is Du's finding that A:C, C:A and the G:U wobble are the
  *best* tolerated, i.e. the weakest discriminators; they currently rank only as
  "not purine:purine".
* Guide position 1 is held in Argonaute's MID pocket and does not pair with the
  target, so a "mismatch" there is not one. It is reported as GuideRegion.ANCHOR
  and ranks with the seed rather than below it; no measurement separates them.
* IUPAC ambiguity codes have no row in Huang's table and no purine/pyrimidine
  identity. They yield None rather than a guess.

References:

* Huang, H., Qiao, R., Zhao, D., Zhang, T., Li, Y., Yi, F. et al. (2009).
  Profiling of mismatch discrimination in RNAi enabled rational design of
  allele-specific siRNAs. Nucleic Acids Research 37(22), 7560-7569.
  doi:10.1093/nar/gkp835
* Schwarz, D.S., Ding, H., Kennington, L., Moore, J.T., Schelter, J., Burchard,
  J. et al. (2006). Designing siRNA that distinguish between genes that differ
  by a single nucleotide. PLoS Genetics 2(9), e140.
  doi:10.1371/journal.pgen.0020140
* Ohnishi, Y., Tamura, Y., Yoshida, M., Tokunaga, K., Hohjoh, H. (2008).
  Enhancement of allele discrimination by introduction of nucleotide mismatches
  into siRNA in allele-specific gene silencing by RNAi. PLoS ONE 3(5), e2248.
  doi:10.1371/journal.pone.0002248
  NO LONGER USED. Formerly cited for "a second mismatch improves
  discrimination", which put mismatch count at the top of the key. The second
  mismatch there is engineered into the siRNA and mismatches both alleles, and
  the natural SNP is central, so it does not support ranking a multi-mismatch
  window above a well-placed single one. See docs/msa_ranking_note.md section 6.
* Du, Q., Thonberg, H., Wang, J., Wahlestedt, C., Liang, Z. (2005). A systematic
  analysis of the silencing effects of an active siRNA at all single-nucleotide
  mismatched target sites. Nucleic Acids Research 33(5), 1671-1677.
  doi:10.1093/nar/gki312
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

# Schwarz report positions in five places. Counting how often each one appears
# is what separates the merely named from the replicated:
#
#   Fig 1B   Drosophila lysate, G:G mismatch        P5 P9 P10 P13 P14 P15 P16
#   Fig 3    HEK293, discriminating in both systems P5 P10 P12 P13 P16
#   Fig 3C   G:U wobble, 3' of the seed             P8 P11 P13 P14 P16
#   Fig 6A   purine:purine scan, P10 scaffold       P4 P7 P9 P10 P11 P13 P16
#   Fig 6B   purine:purine scan, P4 scaffold        P3 P4 P5 P9 P10 P11 P12 P13 P16
#
# P13 and P16 appear in all five, P10 in four. P16 additionally carries the
# paper's own generality test: ten independent siRNA-mRNA pairs against two
# targets, a higher IC50 for the mismatched combination in every pair and at
# least 20-fold in seven of them.
STRONGEST_POSITIONS = frozenset({10, 13, 16})

# The remaining positions Schwarz names, excluding the seed positions discussed
# in the module docstring. Reported (api.schemas.DiscriminationResult,
# app.cli._position_tier) but NOT a ranking tier: the window below already
# covers them, and the paper gives no basis for a tier between the two.
CONSISTENT_POSITIONS = frozenset({9, 12, 14})


# --------------------------------------------------------------------------- #
# Level 3: Du et al. (2005), and the window all three agree on
# --------------------------------------------------------------------------- #

# Du's low-tolerance region - "The low tolerance region comprises position 5-11
# of the target site. Mismatches at these positions abolish most of the
# suppression activities." Those are *target-site* positions; the conversion to
# guide coordinates is verified against the sequences the paper prints:
#
#   target site 5'-CTTATTGGAGAGAGCACGA-3'   (siCD46, 19 nt)
#   guide       5'-TCGTGCTCTCTCCAATAAG-3'   = exact reverse complement
#
# so guide = 20 - target, and target 5-11 is guide P9-P15.
DU_LOW_TOLERANCE = frozenset(range(9, 16))

# The window used for ranking. P9-P15 is Du's low-tolerance region; P16 is added
# because Schwarz replicate it in every assay and Huang list it in 7 of 12 rows,
# even though Du place it one tier lower. 54 of the 55 position entries in
# Huang's table fall inside P9-P16.
LOW_TOLERANCE_WINDOW = DU_LOW_TOLERANCE | frozenset({16})


# --------------------------------------------------------------------------- #
# Mismatch count - reported, never ranked
# --------------------------------------------------------------------------- #

# siDirect treats two mismatches as enough to rule out an off-target, so a window
# with two or more differences is flagged for comparison with that convention.
#
# This is NOT a ranking criterion and must not become one without a source. All
# three papers measure single mismatches; none of them weighs several poorly
# placed differences against one well-placed difference. See the module
# docstring and docs/msa_ranking_note.md.
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

1. huang_recommended - measured for this exact base pair at this exact position.
   The only criterion conditioned on the actual substitution rather than on
   position alone, so specific beats general and it sits first.
2. position in LOW_TOLERANCE_WINDOW (P9-P16) - the one window all three sources
   support: Du's low-tolerance region is P9-P15, Schwarz name P9, P10, P12, P13,
   P14 and P16, and 54 of the 55 entries in Huang's table fall inside it.
3. position in STRONGEST_POSITIONS (P10, P13, P16) - Schwarz, refining within
   the window: P13 and P16 discriminated in all five of their assays and P10 in
   four, and P16 carries their ten-pair generality test.
4. region is POST_SEED - Schwarz: a seed mismatch blocks binding and is
   overwhelmed by concentration; one 3' of the seed blocks catalysis and is not.
   A difference in kind, but ranked below the window because Du place P16-P19 in
   the medium- and high-tolerance regions, so "anywhere 3' of the seed" carries
   less than the window does.
5. class is PURINE_PURINE - Schwarz test this as a class across their assays;
   Du do not use the term but report the same pattern per base pair, with A:C,
   C:A and the G:U wobble best tolerated and A:A and A:G least. JUDGEMENT:
   placed below position because both papers' design advice is positional. The
   cross-axis comparison - a purine:purine at P12 against a milder class at P16
   - is not measured.

The number of mismatches is deliberately absent; see MULTI_MISMATCH.

Entry 5 is the one a reader may reasonably want reordered. Doing so means moving
a line in _single_key; nothing else changes.
"""


def _single_key(assessment: MismatchAssessment) -> tuple[bool, bool, bool, bool, bool]:
    """The published facts about one mismatch, ordered. Higher sorts better."""
    return (
        assessment.huang_recommended is True,
        assessment.guide_position in LOW_TOLERANCE_WINDOW,
        assessment.guide_position in STRONGEST_POSITIONS,
        assessment.region is GuideRegion.POST_SEED,
        assessment.mismatch_class is MismatchClass.PURINE_PURINE,
    )


def best_mismatch(
    assessments: Sequence[MismatchAssessment],
) -> MismatchAssessment | None:
    """The most discriminating of several mismatches against one allele."""
    if not assessments:
        return None
    return max(assessments, key=_single_key)


def ordered_by_placement(
    assessments: Sequence[MismatchAssessment],
) -> tuple[MismatchAssessment, ...]:
    """All mismatches against one allele, best-placed first.

    The ranking turns on the best-placed difference alone, so every other
    difference would otherwise be invisible. Renderers use this to show the whole
    set and let the user weigh count against placement themselves - the trade-off
    this module refuses to make for them.
    """
    return tuple(sorted(assessments, key=_single_key, reverse=True))


def discrimination_key(
    assessments: Sequence[MismatchAssessment],
) -> tuple[bool, bool, bool, bool, bool]:
    """Ranking key for all mismatches against one spared allele.

    Compared lexicographically; higher sorts better. See
    DISCRIMINATION_ORDER_NOTE for what each entry rests on.

    Ranks the *best-placed* difference only. The number of differences is not
    part of the key and is reported separately - no source measures several
    poorly placed differences against one well-placed difference.

    An empty sequence - no difference at all, so no discrimination - sorts below
    everything else.
    """
    if not assessments:
        return (False, False, False, False, False)
    return _single_key(best_mismatch(assessments))
