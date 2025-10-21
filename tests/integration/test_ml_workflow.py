#!/usr/bin/env python3
"""
Integration tests for ML workflow.

These tests run end-to-end ML analysis workflows.
Tests requiring PAUP* will be skipped if PAUP* is not available.
"""

import pytest
from pathlib import Path
import tempfile
import shutil
import subprocess

from pandecay.core.analysis_engine import panDecayIndices


def is_paup_available():
    """Check if PAUP* is available in PATH."""
    try:
        result = subprocess.run(
            ['paup', '-v'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


requires_paup = pytest.mark.skipif(
    not is_paup_available(),
    reason="PAUP* not available"
)


class TestMLWorkflowWithoutPaup:
    """Tests that can run without PAUP* (initialization only)."""

    def test_ml_workflow_initialization(self):
        """Test ML workflow can be initialized."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6,
                analysis_mode="ml",
                threads=1
            )

            # Check initialization
            assert instance.alignment is not None
            assert len(instance.alignment) == 5  # 5 taxa
            assert instance.alignment.get_alignment_length() == 32
            assert instance.do_ml is True
            assert instance.nexus_file_path.exists()

    def test_ml_workflow_with_gamma(self):
        """Test ML workflow initialization with gamma rate variation."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6,
                gamma=True,
                analysis_mode="ml",
                threads=1
            )

            assert instance.gamma_arg is True
            assert "+G" in instance.model_str

    def test_ml_workflow_protein(self):
        """Test ML workflow initialization with protein data."""
        alignment_file = Path("tests/fixtures/alignments/protein_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="protein",
                protein="wag",
                analysis_mode="ml",
                threads=1
            )

            assert instance.protein_arg == "wag"
            assert "WAG" in instance.model_str

    def test_ml_workflow_nexus_input(self):
        """Test ML workflow with NEXUS input."""
        alignment_file = Path("tests/fixtures/alignments/dna_medium.nex")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="nexus",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6,
                analysis_mode="ml",
                threads=1
            )

            assert instance.alignment is not None
            assert len(instance.alignment) == 6  # 6 taxa


@requires_paup
class TestMLWorkflowWithPaup:
    """
    Tests that require PAUP* to be installed.

    These tests are skipped if PAUP* is not available in PATH.
    """

    def test_build_ml_tree(self):
        """Test building ML tree with PAUP*."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6,
                analysis_mode="ml",
                threads=1,
                keep_files=True
            )

            try:
                instance.build_ml_tree()

                # Check that tree was built
                assert instance.ml_tree is not None
                assert instance.ml_likelihood is not None

                # ML likelihood should be negative
                assert instance.ml_likelihood < 0

            except Exception as e:
                pytest.fail(f"ML tree building failed: {e}")

    def test_ml_tree_with_gamma(self):
        """Test building ML tree with gamma rate variation."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6,
                gamma=True,
                analysis_mode="ml",
                threads=1,
                keep_files=True
            )

            try:
                instance.build_ml_tree()

                assert instance.ml_tree is not None
                assert instance.ml_likelihood is not None

            except Exception as e:
                pytest.fail(f"ML tree with gamma failed: {e}")

    def test_ml_tree_different_nst(self):
        """Test building ML tree with different nst values."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        for nst_value in [1, 2, 6]:
            with tempfile.TemporaryDirectory() as tmpdir:
                instance = panDecayIndices(
                    alignment_file=alignment_file,
                    alignment_format="fasta",
                    temp_dir=Path(tmpdir),
                    data_type="dna",
                    nst=nst_value,
                    analysis_mode="ml",
                    threads=1,
                    keep_files=True
                )

                try:
                    instance.build_ml_tree()

                    assert instance.ml_tree is not None
                    assert instance.ml_likelihood is not None

                except Exception as e:
                    pytest.fail(f"ML tree with nst={nst_value} failed: {e}")


class TestMLWorkflowErrorHandling:
    """Test error handling in ML workflow."""

    def test_invalid_alignment_format(self):
        """Test handling of invalid alignment format."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with pytest.raises(Exception):  # BioPython will raise an exception
            with tempfile.TemporaryDirectory() as tmpdir:
                panDecayIndices(
                    alignment_file=alignment_file,
                    alignment_format="invalid_format",
                    temp_dir=Path(tmpdir),
                    data_type="dna",
                    nst=6
                )

    def test_empty_temp_directory_handling(self):
        """Test that temporary directory is created if it doesn't exist."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            temp_path = Path(tmpdir) / "new_temp_dir"

            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=temp_path,
                data_type="dna",
                nst=6
            )

            assert temp_path.exists()
            assert temp_path.is_dir()
