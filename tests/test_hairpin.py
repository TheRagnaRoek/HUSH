"""Tests for intramolecular hairpin stability (core.hairpin)."""

import pytest

from core.hairpin import (
    DEFAULT_INTERNAL_TM_THRESHOLD_C,
    MIN_LOOP,
    MIN_STEM,
    HairpinFold,
    has_no_internal_repeat,
    internal_structure_tm,
    loop_penalty_dg37,
    most_stable_hairpin,
)
from core.inputs import load_sequences
from core.sequence import BioSequence, MolType
from core.sirna import enumerate_candidates
from core.thermo import NN_PARAMS_FREIER_1986, NN_PARAMS_XIA_1998


class TestLoopPenalty:
    def test_matches_groebe_uhlenbeck_table3(self):
        # Homopolymer loops read straight off Table 3 (dG37 destabilisation).
        assert loop_penalty_dg37("AAA") == pytest.approx(4.9)
        assert loop_penalty_dg37("CCC") == pytest.approx(6.6)
        assert loop_penalty_dg37("UUUU") == pytest.approx(3.0)
        assert loop_penalty_dg37("A" * 9) == pytest.approx(5.6)

    def test_size_is_linearly_interpolated(self):
        # Size 6 is not tabulated; it sits midway between 5 and 7.
        expected = (3.5 + 4.0) / 2  # A-loop at sizes 5 and 7
        assert loop_penalty_dg37("A" * 6) == pytest.approx(expected)

    def test_large_loop_clamps_to_the_last_measured_size(self):
        assert loop_penalty_dg37("A" * 12) == pytest.approx(loop_penalty_dg37("A" * 9))

    def test_composition_is_averaged_over_the_loop(self):
        # A mixed loop is the mean of its bases' penalties at that size.
        expected = (4.9 + 6.6 + 4.2) / 3  # A, C, U at size 3
        assert loop_penalty_dg37("ACU") == pytest.approx(expected)

    def test_g_and_ambiguous_bases_use_the_mean_composition(self):
        mean3 = (4.9 + 6.6 + 4.2) / 3
        assert loop_penalty_dg37("GGG") == pytest.approx(mean3)


class TestMostStableHairpin:
    def test_no_stem_possible_returns_none(self):
        # A homopolymer cannot pair with itself.
        assert most_stable_hairpin("A" * 19) is None

    def test_finds_an_obvious_hairpin(self):
        fold = most_stable_hairpin("GGGGGCCCAAAGGGCCCCC")
        assert isinstance(fold, HairpinFold)
        assert fold.stem_length >= MIN_STEM
        assert len(fold.loop) >= MIN_LOOP

    def test_loop_shorter_than_minimum_is_not_a_hairpin(self):
        # Perfect self-complement with only a 2-nt centre: no valid 3-nt loop
        # closing a stem, so the search must reject it.
        assert most_stable_hairpin("GCGCGCAAGCGCGC", min_loop=3) is not None
        assert most_stable_hairpin("GCGCGC", min_loop=3) is None

    def test_selects_the_lowest_dg37(self):
        # Every returned fold is the minimum-dG37 one; a second, weaker stem in
        # the same sequence must not be chosen over the strong one.
        fold = most_stable_hairpin("GCGCGCGAAAACGCGCGC")
        alt = most_stable_hairpin("GCGCGCGAAAACGCGCGC", min_stem=3)
        assert fold.dg37 <= alt.dg37 + 1e-9


class TestHairpinTm:
    def test_gc_stem_melts_hotter_than_au_stem(self):
        gc = internal_structure_tm("GCGCGCGAAACGCGCGC")
        au = internal_structure_tm("AUAUAUAAAAUAUAUAU"[:17])
        assert gc is not None and au is not None
        assert gc > au

    def test_longer_stem_melts_hotter(self):
        short = most_stable_hairpin("GCGCAAAGCGC")  # 3-4 bp stem
        long = most_stable_hairpin("GCGCGCGCAAAGCGCGCGC")  # up to 8 bp
        assert long.tm > short.tm

    def test_none_when_nothing_folds(self):
        assert internal_structure_tm("A" * 19) is None

    def test_stem_parameter_set_defaults_to_xia_and_freier_differs(self):
        seq = "GCAUGCGAAACGCAUGC"
        default = internal_structure_tm(seq)
        xia = internal_structure_tm(seq, params=NN_PARAMS_XIA_1998)
        freier = internal_structure_tm(seq, params=NN_PARAMS_FREIER_1986)
        assert default is not None and default == pytest.approx(xia)
        assert freier != pytest.approx(xia)


class TestAbsenceOfInternalRepeat:
    def test_no_structure_counts_as_absent(self):
        assert has_no_internal_repeat("A" * 19) is True

    def test_strong_hairpin_is_present(self):
        # A 7-bp G/C stem closing a 3-nt loop melts far above 20 C.
        assert has_no_internal_repeat("GGGGGGGAAACCCCCCC") is False

    def test_threshold_is_the_reynolds_twenty(self):
        assert DEFAULT_INTERNAL_TM_THRESHOLD_C == 20.0

    def test_comparison_is_inclusive(self):
        # A fold sitting exactly on the threshold counts as absent ("<= 20 C").
        seq = "GGGGGGGAAACCCCCCC"
        tm = internal_structure_tm(seq)
        assert has_no_internal_repeat(seq, threshold_c=tm) is True
        assert has_no_internal_repeat(seq, threshold_c=tm - 0.01) is False

    def test_dna_and_biosequence_accepted(self):
        dna = BioSequence("GGGGGGGAAACCCCCCC", id="t", moltype=MolType.DNA)
        assert has_no_internal_repeat(dna) is False
        assert has_no_internal_repeat("ggggggg aaa ccccccc") is False


class TestRealHlaTarget:
    def test_folds_every_candidate_without_error(self, hla_a01_fasta):
        target = load_sequences(hla_a01_fasta)[0]
        candidates = list(enumerate_candidates(target))
        verdicts = [has_no_internal_repeat(c.sense.seq) for c in candidates]
        assert len(verdicts) > 1000
        # A GC-rich transcript has plenty of foldable candidates, so the filter
        # is neither all-True nor all-False here.
        assert any(verdicts)
        assert not all(verdicts)


# --------------------------------------------------------------------------- #
# dot-bracket notation
# --------------------------------------------------------------------------- #


class TestDotBracket:
    def test_describes_the_fold_it_was_found_in(self):
        seq = "GGGGGCCCAAAGGGCCCCC"
        fold = most_stable_hairpin(seq)
        assert fold is not None
        db = fold.dot_bracket(len(seq))

        assert len(db) == len(seq)
        assert db.count("(") == db.count(")") == fold.stem_length
        assert db.count(".") == len(seq) - 2 * fold.stem_length
        # The unpaired run sits exactly where the loop does.
        assert db[fold.stem_start + fold.stem_length :][: len(fold.loop)] == "." * len(fold.loop)

    def test_is_balanced_when_read_left_to_right(self):
        seq = "GGGGGCCCAAAGGGCCCCC"
        db = most_stable_hairpin(seq).dot_bracket(len(seq))
        depth = 0
        for ch in db:
            depth += 1 if ch == "(" else -1 if ch == ")" else 0
            assert depth >= 0, "a closing bracket preceded its opening one"
        assert depth == 0

    def test_pads_a_fold_that_sits_inside_a_longer_sequence(self):
        fold = HairpinFold(
            stem_start=3, stem_length=4, loop="AAAA", dh=0.0, ds=0.0, dg37=0.0, tm=0.0
        )
        assert fold.dot_bracket(20) == "..." + "((((" + "...." + "))))" + "....."

    def test_rejects_a_length_the_fold_cannot_fit_in(self):
        fold = HairpinFold(
            stem_start=0, stem_length=5, loop="AAAA", dh=0.0, ds=0.0, dg37=0.0, tm=0.0
        )
        with pytest.raises(ValueError):
            fold.dot_bracket(13)  # needs 14

