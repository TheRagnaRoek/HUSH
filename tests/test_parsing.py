"""Tests for the shared flag-string parsing (app.parsing).

The CLI and the shell both hand user-typed strings to these five functions, so
the failure messages here are what a user actually reads. Each parser is
checked on its accepted spellings, its whitespace and case handling, and the
exact wording of its error.
"""

import pytest

from app.parsing import (
    parse_mode,
    parse_overhang,
    parse_params,
    parse_rules,
    split_ids,
)
from core.selection import RuleMode
from core.sirna import OverhangPolicy
from core.thermo import ThermoParams


class TestParseRules:
    @pytest.mark.parametrize("spec", ["all", "ALL", "  all  ", ""])
    def test_all_and_empty_enable_everything(self, spec):
        # An empty --rules means "do not narrow the default", not "no rules".
        assert parse_rules(spec) == (True, True, True)

    @pytest.mark.parametrize("spec", ["none", "NONE", " none "])
    def test_none_disables_everything(self, spec):
        assert parse_rules(spec) == (False, False, False)

    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("ui-tei", (True, False, False)),
            ("uitei", (True, False, False)),
            ("ui_tei", (True, False, False)),
            ("reynolds", (False, True, False)),
            ("amarzguioui", (False, False, True)),
            ("amarz", (False, False, True)),
        ],
    )
    def test_every_alias_selects_its_rule(self, spec, expected):
        assert parse_rules(spec) == expected

    def test_a_comma_list_selects_several(self):
        assert parse_rules("ui-tei,amarz") == (True, False, True)

    def test_spacing_and_case_are_ignored(self):
        assert parse_rules("  UI-Tei ,  REYNOLDS ") == (True, True, False)

    def test_repeating_a_rule_is_harmless(self):
        # Flags are set, not counted, so a doubled token means the same thing.
        assert parse_rules("reynolds,reynolds") == (False, True, False)

    @pytest.mark.parametrize(
        "spec,expected",
        [("ui-tei,", (True, False, False)), (",reynolds", (False, True, False))],
    )
    def test_trailing_and_leading_commas_are_skipped(self, spec, expected):
        # Common when a list is edited by hand; they must not fail.
        assert parse_rules(spec) == expected

    @pytest.mark.parametrize("spec", [",", ",,"])
    def test_a_list_of_only_commas_enables_nothing(self, spec):
        # Pinned rather than endorsed: an empty string enables everything while
        # a string of separators enables nothing, so a stray comma silently
        # turns all rules off. Change this deliberately, not by accident.
        assert parse_rules(spec) == (False, False, False)

    def test_an_unknown_rule_names_itself_and_the_choices(self):
        with pytest.raises(ValueError) as excinfo:
            parse_rules("ui-tei,bogus")
        message = str(excinfo.value)
        assert "'bogus'" in message
        for choice in ("ui-tei", "reynolds", "amarzguioui", "all", "none"):
            assert choice in message

    def test_all_and_none_are_not_list_tokens(self):
        # They are whole-string spellings; inside a list they are unknown rules.
        with pytest.raises(ValueError):
            parse_rules("all,reynolds")


class TestParseMode:
    @pytest.mark.parametrize(
        "spec,expected",
        [("all", RuleMode.ALL), ("any", RuleMode.ANY), ("  ANY  ", RuleMode.ANY)],
    )
    def test_accepted_spellings(self, spec, expected):
        assert parse_mode(spec) is expected

    def test_an_unknown_mode_names_itself_and_the_choices(self):
        with pytest.raises(ValueError) as excinfo:
            parse_mode("either")
        message = str(excinfo.value)
        assert "'either'" in message
        assert "'all'" in message and "'any'" in message

    def test_the_enum_error_does_not_leak(self):
        # RuleMode() raises its own ValueError; the message must be ours.
        with pytest.raises(ValueError) as excinfo:
            parse_mode("")
        assert "is not a valid RuleMode" not in str(excinfo.value)


class TestParseParams:
    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("xia", ThermoParams.XIA_1998),
            ("xia1998", ThermoParams.XIA_1998),
            ("freier", ThermoParams.FREIER_1986),
            ("freier1986", ThermoParams.FREIER_1986),
            ("  FREIER  ", ThermoParams.FREIER_1986),
        ],
    )
    def test_accepted_spellings(self, spec, expected):
        assert parse_params(spec) is expected

    def test_an_unknown_set_names_itself_and_the_choices(self):
        with pytest.raises(ValueError) as excinfo:
            parse_params("santalucia")
        message = str(excinfo.value)
        assert "'santalucia'" in message
        assert "'xia'" in message and "'freier'" in message


class TestParseOverhang:
    @pytest.mark.parametrize(
        "spec,expected",
        [
            ("natural", OverhangPolicy.NATURAL),
            ("mrna", OverhangPolicy.NATURAL),
            ("dtdt", OverhangPolicy.DTDT),
            ("tt", OverhangPolicy.DTDT),
            ("uu", OverhangPolicy.UU),
            ("  UU  ", OverhangPolicy.UU),
        ],
    )
    def test_accepted_spellings(self, spec, expected):
        assert parse_overhang(spec) is expected

    def test_an_unknown_policy_names_itself_and_the_choices(self):
        with pytest.raises(ValueError) as excinfo:
            parse_overhang("blunt")
        message = str(excinfo.value)
        assert "'blunt'" in message
        for choice in ("'natural'", "'dtdt'", "'uu'"):
            assert choice in message

    def test_every_policy_is_reachable_by_at_least_one_spelling(self):
        # A new OverhangPolicy without an alias would be unreachable from the
        # CLI, which is the kind of gap this guards against.
        reachable = {parse_overhang(s) for s in ("natural", "dtdt", "uu")}
        assert reachable == set(OverhangPolicy)


class TestSplitIds:
    def test_a_plain_list(self):
        assert split_ids("HLA-A*01:01,HLA-A*02:01") == ["HLA-A*01:01", "HLA-A*02:01"]

    def test_surrounding_space_is_stripped(self):
        assert split_ids(" HLA-A*01:01 , HLA-A*02:01 ") == [
            "HLA-A*01:01",
            "HLA-A*02:01",
        ]

    @pytest.mark.parametrize("spec", ["", "   ", ",", ",,", " , , "])
    def test_blank_input_gives_an_empty_list(self, spec):
        assert split_ids(spec) == []

    def test_blank_tokens_are_dropped_not_kept_as_empty_ids(self):
        # An empty id would later fail an allele lookup with a confusing
        # message, so it is dropped here.
        assert split_ids("A,,B,") == ["A", "B"]

    def test_order_is_preserved_and_duplicates_are_kept(self):
        # Deduplication belongs to AlleleTargeting, not to the parser.
        assert split_ids("B,A,B") == ["B", "A", "B"]

    def test_case_is_not_folded(self):
        # Allele ids are case-sensitive identifiers, unlike the option words.
        assert split_ids("hla-a*01:01") == ["hla-a*01:01"]
