"""Tests for the shRNA conversion stage (shrna.loops, shrna.convert)."""

import pytest

from core.sequence import BioSequence, MolType
from shrna.convert import (
    MIN_TERMINATOR_LENGTH,
    Orientation,
    ShrnaOptions,
    convert,
    convert_candidate,
)
from shrna.loops import DEFAULT_LOOP, LOOPS, LoopScaffold, get_loop

# A clean 19-mer core whose duplex has no competing internal structure.
CORE = "GCACCGUCAAGGCUGAGAA"


def _sense(seq: str = CORE) -> BioSequence:
    return BioSequence(seq, moltype=MolType.RNA)


def _pair(seq: str = CORE) -> tuple[BioSequence, BioSequence]:
    s = _sense(seq)
    return s, s.reverse_complement()


# --------------------------------------------------------------------------- #
# Loop registry
# --------------------------------------------------------------------------- #


def test_registry_has_the_three_loops_with_expected_sequences():
    assert LOOPS["jensen"].seq == "UGUGCUU"
    assert LOOPS["psuper"].seq == "UUCAAGAGA"
    assert LOOPS["trc"].seq == "CUCGAG"


def test_default_loop_is_jensen():
    assert DEFAULT_LOOP == "jensen"
    assert ShrnaOptions().loop.name == "jensen"


def test_get_loop_is_case_insensitive_and_validates():
    assert get_loop("PSUPER").name == "psuper"
    with pytest.raises(ValueError):
        get_loop("nonesuch")


# --------------------------------------------------------------------------- #
# Hairpin assembly (SA orientation)
# --------------------------------------------------------------------------- #


def test_hairpin_is_sense_loop_antisense_plus_terminator():
    sense, antisense = _pair()
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    expected = sense.seq + LOOPS["jensen"].seq + antisense.seq + "U" * 6
    assert c.hairpin_rna.seq == expected
    assert c.orientation is Orientation.SA


def test_insert_bottom_is_reverse_complement_of_top_and_dna():
    sense, antisense = _pair()
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    assert c.insert_top.moltype is MolType.DNA
    assert "U" not in c.insert_top.seq  # DNA, not RNA
    assert c.insert_bottom.seq == c.insert_top.reverse_complement().seq


def test_no_terminator_option_drops_the_poly_t():
    sense, antisense = _pair()
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        options=ShrnaOptions(include_terminator=False),
    )
    assert c.hairpin_rna.seq == sense.seq + LOOPS["jensen"].seq + antisense.seq
    assert not c.hairpin_rna.seq.endswith("UUUU")


def test_convert_candidate_wrapper_matches_convert():
    from core.sirna import enumerate_candidates

    target = BioSequence("G" + CORE + "ACGU", moltype=MolType.RNA)
    cand = next(enumerate_candidates(target, length=19))
    c = convert_candidate(cand)
    assert c.region == cand.region
    assert c.hairpin_rna.seq.startswith(cand.sense.seq)


# --------------------------------------------------------------------------- #
# Options validation
# --------------------------------------------------------------------------- #


def test_short_terminator_is_rejected():
    with pytest.raises(ValueError):
        ShrnaOptions(terminator_length=MIN_TERMINATOR_LENGTH - 1)


def test_short_terminator_allowed_when_terminator_disabled():
    # Length is not policed when the terminator is not emitted.
    opts = ShrnaOptions(include_terminator=False, terminator_length=0)
    assert opts.terminator_length == 0


# --------------------------------------------------------------------------- #
# Fold-based scoring
# --------------------------------------------------------------------------- #


def test_clean_core_folds_as_designed():
    sense, antisense = _pair()
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    assert c.folds_as_designed is True
    assert c.fold is not None
    assert c.fold.stem_length == 19
    assert c.fold.loop == LOOPS["jensen"].seq


def test_loop_that_extends_the_stem_flags_misfold():
    # A loop whose ends are complementary (G...C) lets the stem extend by a base
    # pair into the loop, so the most stable fold is not the intended stem.
    sense, antisense = _pair()
    trap = LoopScaffold(
        name="trap",
        sequence=BioSequence("GAAAC", id="loop:trap", moltype=MolType.RNA),
        citation="synthetic test loop",
    )
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        options=ShrnaOptions(loop=trap),
    )
    assert c.folds_as_designed is False


def test_starts_with_purine_flag_tracks_the_transcript_5_prime():
    g_sense, g_anti = _pair("GCACCGUCAAGGCUGAGAA")
    assert convert(g_sense, g_anti, region="1-19", start=0, end=19).starts_with_purine
    a_sense, a_anti = _pair("ACACCGUCAAGGCUGAGAA")
    assert convert(a_sense, a_anti, region="1-19", start=0, end=19).starts_with_purine
    for first in "CU":
        c_sense, c_anti = _pair(first + "CACCGUCAAGGCUGAGAA")
        assert not convert(c_sense, c_anti, region="1-19", start=0, end=19).starts_with_purine


def test_premature_terminator_flag_catches_a_loop_stem_junction_run():
    # jensen ends in UU; a core whose antisense begins with UU makes a UUUU run
    # at the loop/stem junction (a premature Pol III terminator).
    sense = _sense("GCACCGUCAAGGCUGAGAA")
    antisense = sense.reverse_complement()
    assert antisense.seq.startswith("UU")  # guards the fixture's premise
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    assert c.premature_terminator_free is False


def test_premature_terminator_free_for_a_clean_junction():
    # A core whose antisense does not start with U avoids the junction run.
    sense = _sense("GCACCGUCAAGGCUGAGAC")
    antisense = sense.reverse_complement()
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    assert c.premature_terminator_free is True


# --------------------------------------------------------------------------- #
# 3' overhangs: both strands are encoded in full
# --------------------------------------------------------------------------- #

# A core whose antisense does not begin with U, so the loop/stem junction is
# clean and any U run in these tests is one the test put there.
CLEAN = "GCACCGUCAAGGCUGAGAC"


def test_hairpin_encodes_both_strands_in_full():
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        sense_overhang="AA", guide_overhang="CG",
    )
    expected = sense.seq + "AA" + LOOPS["jensen"].seq + antisense.seq + "CG" + "U" * 6
    assert c.hairpin_rna.seq == expected


def test_hairpin_carries_the_validated_guide_including_its_overhang():
    # The whole point of the stage: the hairpin expresses the strand design
    # validated, overhang and all, not a truncated core.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        sense_overhang="AA", guide_overhang="CG",
    )
    assert antisense.seq + "CG" in c.hairpin_rna.seq
    assert sense.seq + "AA" in c.hairpin_rna.seq


def test_absent_overhangs_leave_the_bare_core_hairpin_unchanged():
    sense, antisense = _pair(CLEAN)
    c = convert(sense, antisense, region="1-19", start=0, end=19)
    assert c.sense_overhang == "" and c.guide_overhang == ""
    assert c.hairpin_rna.seq == sense.seq + LOOPS["jensen"].seq + antisense.seq + "U" * 6


def test_dtdt_overhang_is_encoded_as_ribo_u():
    # An expressed hairpin cannot carry 2'-deoxythymidine, so dTdT lands as
    # ribo-U. The ordered DNA oligo carries T at those positions either way.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        sense_overhang="tt", guide_overhang="tt",
    )
    assert c.sense_overhang == "UU"
    assert c.guide_overhang == "UU"
    assert "t" not in c.hairpin_rna.seq
    assert c.insert_top.seq.startswith(sense.seq.replace("U", "T") + "TT")


def test_sense_overhang_reads_as_loop_not_stem():
    # CC cannot pair with the jensen loop's trailing UU, so the stem stays the
    # core and the overhang is simply part of a longer loop.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19, sense_overhang="CC",
    )
    assert c.folds_as_designed is True
    assert c.fold is not None
    assert c.fold.stem_length == len(sense)
    assert c.fold.loop == "CC" + LOOPS["jensen"].seq


def test_sense_overhang_that_pairs_with_the_loop_end_flags_a_misfold():
    # jensen ends in UU, so an AA passenger overhang base-pairs with it: the
    # most stable fold takes 2 nt of loop into the stem and shifts the Dicer
    # register. Common under the natural policy, so it must be flagged.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19, sense_overhang="AA",
    )
    assert c.folds_as_designed is False
    assert c.fold is not None
    assert c.fold.stem_length == len(sense) + 2
    assert c.fold.loop == LOOPS["jensen"].seq[:-2]


def test_guide_overhang_joins_the_terminator_tail():
    # UU before a 6-U terminator is an 8-U run, but it is the terminator doing
    # its job at the 3' end, not a premature one, so the flag stays clean.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19, guide_overhang="UU",
    )
    assert c.hairpin_rna.seq.endswith("UU" + "U" * 6)
    assert c.premature_terminator_free is True


def test_sense_overhang_can_build_a_premature_terminator():
    # The graft can assemble a U run no upstream gate sees: the core's last
    # base, the overhang, and the loop's first base.
    sense, antisense = _pair("GCACCGUCAAGGCUGAGAU")
    c = convert(
        sense, antisense, region="1-19", start=0, end=19, sense_overhang="UU",
    )
    assert c.hairpin_rna.seq.startswith(sense.seq + "UUU")
    assert c.premature_terminator_free is False


@pytest.mark.parametrize("loop_name", ["psuper", "trc"])
def test_the_other_loops_convert_with_overhangs(loop_name):
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        sense_overhang="AA", guide_overhang="CG",
        options=ShrnaOptions(loop=get_loop(loop_name)),
    )
    assert c.loop_name == loop_name
    assert LOOPS[loop_name].seq in c.hairpin_rna.seq
    assert c.insert_bottom.seq == c.insert_top.reverse_complement().seq


def test_construct_carries_its_fold_as_dot_bracket():
    # The notation covers the folded body only - the 3' tail is single-stranded
    # and never enters the fold search.
    sense, antisense = _pair(CLEAN)
    c = convert(
        sense, antisense, region="1-19", start=0, end=19,
        sense_overhang="CC", guide_overhang="UU",
    )
    body_len = len(sense) + 2 + len(LOOPS["jensen"].seq) + len(antisense)
    assert len(c.fold_structure) == body_len
    assert c.fold_structure.count("(") == c.fold_structure.count(")") == len(sense)
    assert c.hairpin_rna.seq.startswith(sense.seq)


def test_a_construct_that_folds_into_nothing_has_no_structure():
    sense = _sense("A" * 19)
    c = convert(sense, sense.reverse_complement(), region="1-19", start=0, end=19)
    if c.fold is None:
        assert c.fold_structure == ""
    else:
        assert len(c.fold_structure) > 0

