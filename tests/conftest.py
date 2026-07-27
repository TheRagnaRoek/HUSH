"""Shared pytest fixtures: paths to example and test data."""

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES_DIR = REPO_ROOT / "data" / "examples"
TEST_DATA_DIR = Path(__file__).resolve().parent / "data"


@pytest.fixture
def examples_dir() -> Path:
    return EXAMPLES_DIR


@pytest.fixture
def test_data_dir() -> Path:
    return TEST_DATA_DIR


@pytest.fixture
def target_fasta() -> Path:
    return EXAMPLES_DIR / "target_example.fasta"


@pytest.fixture
def panel_fasta() -> Path:
    return EXAMPLES_DIR / "hla_panel_example.fasta"


@pytest.fixture
def hla_a01_fasta() -> Path:
    """A real HLA-A*01:01 coding sequence (IPD-IMGT/HLA HLA00001)."""
    return EXAMPLES_DIR / "hla_a01.fasta"


@pytest.fixture(scope="session")
def amarzguioui_table3() -> list[dict]:
    """Table 3 of Amarzguioui & Prydz (2004): 46 siRNAs scored by the authors."""
    path = TEST_DATA_DIR / "amarzguioui_table3.tsv"
    rows = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("name\t"):
                continue
            name, sense, functional, dt3, score = line.split("\t")
            rows.append(
                {
                    "name": name,
                    "sense": sense,
                    "functional": functional,
                    "dt3": int(dt3),
                    "score": int(score),
                }
            )
    return rows


@pytest.fixture
def mini_pdb() -> Path:
    return TEST_DATA_DIR / "mini_protein.pdb"
