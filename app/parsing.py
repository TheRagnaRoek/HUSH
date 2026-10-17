"""Flag-string parsing shared by the CLI and the shell.

Raises ValueError with a user-facing message; the CLI turns it into a usage
error, the shell prints it.
"""

from __future__ import annotations

from core.selection import RuleMode
from core.sirna import OverhangPolicy
from core.thermo import ThermoParams

_RULE_ALIASES = {
    "ui-tei": "ui_tei",
    "uitei": "ui_tei",
    "ui_tei": "ui_tei",
    "reynolds": "reynolds",
    "amarzguioui": "amarzguioui",
    "amarz": "amarzguioui",
}

_PARAMS_ALIASES = {
    "xia": ThermoParams.XIA_1998,
    "xia1998": ThermoParams.XIA_1998,
    "freier": ThermoParams.FREIER_1986,
    "freier1986": ThermoParams.FREIER_1986,
}

_OVERHANG_ALIASES = {
    "natural": OverhangPolicy.NATURAL,
    "mrna": OverhangPolicy.NATURAL,
    "dtdt": OverhangPolicy.DTDT,
    "tt": OverhangPolicy.DTDT,
    "uu": OverhangPolicy.UU,
}


def parse_rules(spec: str) -> tuple[bool, bool, bool]:
    """(use_ui_tei, use_reynolds, use_amarzguioui) from a comma list, 'all' or 'none'."""
    text = spec.strip().lower()
    if text in ("", "all"):
        return True, True, True
    if text == "none":
        return False, False, False
    flags = {"ui_tei": False, "reynolds": False, "amarzguioui": False}
    for token in (t.strip() for t in text.split(",")):
        if not token:
            continue
        key = _RULE_ALIASES.get(token)
        if key is None:
            raise ValueError(
                f"unknown rule {token!r}; choose from ui-tei, reynolds, "
                "amarzguioui (comma-separated), or 'all'/'none'."
            )
        flags[key] = True
    return flags["ui_tei"], flags["reynolds"], flags["amarzguioui"]


def parse_mode(spec: str) -> RuleMode:
    """Rule combination mode: 'all' or 'any'."""
    try:
        return RuleMode(spec.strip().lower())
    except ValueError:
        raise ValueError(f"unknown mode {spec!r}; choose 'all' or 'any'.")


def parse_params(spec: str) -> ThermoParams:
    """Nearest-neighbor parameter set: 'xia' or 'freier'."""
    key = _PARAMS_ALIASES.get(spec.strip().lower())
    if key is None:
        raise ValueError(
            f"unknown parameter set {spec!r}; choose 'xia' or 'freier'."
        )
    return key


def parse_overhang(spec: str) -> OverhangPolicy:
    """3' overhang policy: 'natural', 'dtdt' or 'uu'."""
    key = _OVERHANG_ALIASES.get(spec.strip().lower())
    if key is None:
        raise ValueError(
            f"unknown overhang {spec!r}; choose 'natural', 'dtdt', or 'uu'."
        )
    return key


def split_ids(spec: str) -> list[str]:
    """Split a comma-separated id list, dropping blanks."""
    return [token.strip() for token in spec.split(",") if token.strip()]
