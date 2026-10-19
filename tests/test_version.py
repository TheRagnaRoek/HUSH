"""The version is declared in three files and read at runtime; they must agree."""

import json
import re
import tomllib
from pathlib import Path

from typer.testing import CliRunner

from app.cli import app
from core.version import __version__

ROOT = Path(__file__).resolve().parent.parent
runner = CliRunner()


def _pyproject_version() -> str:
    with open(ROOT / "pyproject.toml", "rb") as handle:
        return tomllib.load(handle)["project"]["version"]


def test_citation_version_matches_pyproject():
    text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"^version:\s*(\S+)\s*$", text, re.MULTILINE)
    assert match, "CITATION.cff has no top-level version"
    assert match.group(1) == _pyproject_version()


def test_web_package_version_matches_pyproject():
    package = json.loads((ROOT / "web" / "package.json").read_text(encoding="utf-8"))
    assert package["version"] == _pyproject_version()


def test_cli_reports_the_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == f"hush {__version__}"
