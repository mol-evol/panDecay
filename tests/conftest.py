#!/usr/bin/env python3
"""
Pytest configuration and shared fixtures for panDecay tests.
"""

import pytest
from pathlib import Path
import tempfile
import shutil


@pytest.fixture
def fixtures_dir():
    """Return path to fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def dna_alignment_small(fixtures_dir):
    """Return path to small DNA alignment fixture."""
    path = fixtures_dir / "alignments" / "dna_small.fas"
    if not path.exists():
        pytest.skip(f"Fixture not found: {path}")
    return path


@pytest.fixture
def dna_alignment_nexus(fixtures_dir):
    """Return path to NEXUS DNA alignment fixture."""
    path = fixtures_dir / "alignments" / "dna_medium.nex"
    if not path.exists():
        pytest.skip(f"Fixture not found: {path}")
    return path


@pytest.fixture
def protein_alignment_small(fixtures_dir):
    """Return path to small protein alignment fixture."""
    path = fixtures_dir / "alignments" / "protein_small.fas"
    if not path.exists():
        pytest.skip(f"Fixture not found: {path}")
    return path


@pytest.fixture
def simple_tree(fixtures_dir):
    """Return path to simple tree fixture."""
    path = fixtures_dir / "trees" / "simple_tree.nwk"
    if not path.exists():
        pytest.skip(f"Fixture not found: {path}")
    return path


@pytest.fixture
def test_config(fixtures_dir):
    """Return path to test configuration file."""
    path = fixtures_dir / "configs" / "test_config.txt"
    if not path.exists():
        pytest.skip(f"Fixture not found: {path}")
    return path


@pytest.fixture
def temp_workspace():
    """Create a temporary workspace for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def temp_output_dir(temp_workspace):
    """Create a temporary output directory."""
    output_dir = temp_workspace / "output"
    output_dir.mkdir()
    return output_dir


def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "requires_paup: tests that require PAUP* to be installed"
    )
    config.addinivalue_line(
        "markers", "requires_mrbayes: tests that require MrBayes to be installed"
    )
    config.addinivalue_line(
        "markers", "slow: tests that take significant time to run"
    )
    config.addinivalue_line(
        "markers", "integration: integration tests that test full workflows"
    )


def pytest_collection_modifyitems(config, items):
    """Automatically mark integration tests."""
    for item in items:
        # Mark integration tests
        if "integration" in str(item.fspath):
            item.add_marker(pytest.mark.integration)

        # Mark slow tests based on name patterns
        if "workflow" in item.name.lower() or "full" in item.name.lower():
            item.add_marker(pytest.mark.slow)
