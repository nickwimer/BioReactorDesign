"""Test setup: import path and optional local inputs.

Run from anywhere with::

    python -m pytest scripts/ma_reactions/kinetics_fit/tests

Tests that need the surrogate files or the Exp 17 extract are skipped when
those are not available (see local.example.yaml).
"""

import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))

from fit_kinetics import local_config  # noqa: E402
from kinetics0d import Surrogate  # noqa: E402


def pytest_addoption(parser):
    parser.addoption("--surrogate-dir", default=None)
    parser.addoption("--exp17-dir", default=None)


@pytest.fixture(scope="session")
def surrogate_dir(request):
    path = request.config.getoption("--surrogate-dir") or local_config().get(
        "surrogate_dir"
    )
    if not path or not (Path(path) / "tri_points.txt").exists():
        pytest.skip("surrogate files not available")
    return Path(path)


@pytest.fixture(scope="session")
def sur(surrogate_dir):
    return Surrogate(surrogate_dir)


@pytest.fixture(scope="session")
def exp17_dir(request):
    path = request.config.getoption("--exp17-dir") or local_config().get(
        "exp17_dir"
    )
    if not path or not (Path(path) / "samples.csv").exists():
        pytest.skip("Exp 17 extract not available")
    return Path(path)


@pytest.fixture(scope="session")
def expected():
    """Expected values for the unpublished inputs, from local.yaml.

    They are kept out of the repository with the inputs themselves. Tests
    that compare against them are skipped when the block is missing.
    """
    return local_config().get("expected") or {}
