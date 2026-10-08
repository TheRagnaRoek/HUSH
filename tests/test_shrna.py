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
