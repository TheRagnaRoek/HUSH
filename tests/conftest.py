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
def mini_pdb() -> Path:
    return TEST_DATA_DIR / "mini_protein.pdb"
