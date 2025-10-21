#!/usr/bin/env python3
"""
Unit tests for panDecay analysis engine core functionality.

These tests focus on testable components that don't require external tools.
"""

import pytest
from pathlib import Path
import tempfile
import shutil

from pandecay.core.analysis_engine import panDecayIndices, AnalysisEngineError


class TestPanDecayIndicesInit:
    """Test panDecayIndices initialization."""

    def test_init_with_valid_alignment(self):
        """Test initialization with valid alignment file."""
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
                keep_files=True
            )

            assert instance.alignment_file == alignment_file
            assert instance.alignment_format == "fasta"
            assert instance.data_type == "dna"
            assert instance.nst_arg == 6
            assert instance.threads >= 1

    def test_init_with_nonexistent_file(self):
        """Test initialization with non-existent alignment file."""
        with pytest.raises(AnalysisEngineError, match="Alignment file not found"):
            panDecayIndices(
                alignment_file="/nonexistent/file.fas",
                alignment_format="fasta",
                data_type="dna",
                nst=6
            )

    def test_init_creates_temp_directory(self):
        """Test that initialization creates temporary directory."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        instance = panDecayIndices(
            alignment_file=alignment_file,
            alignment_format="fasta",
            data_type="dna",
            nst=6,
            keep_files=True
        )

        assert instance.temp_path.exists()
        assert instance.temp_path.is_dir()

        # Cleanup
        if instance.temp_path.exists():
            shutil.rmtree(instance.temp_path)

    def test_init_thread_auto(self):
        """Test thread count with 'auto' setting."""
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
                threads="auto"
            )

            # Auto should select a reasonable number
            assert instance.threads >= 1
            assert isinstance(instance.threads, int)

    def test_init_thread_custom(self):
        """Test thread count with custom value."""
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
                threads=4
            )

            assert instance.threads == 4

    def test_init_invalid_thread_count(self):
        """Test thread count validation."""
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
                threads=0  # Invalid
            )

            # Should default to 1
            assert instance.threads == 1

    def test_init_data_type_validation(self):
        """Test data type normalization."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="DNA",  # Uppercase
                nst=6
            )

            # Should be normalized to lowercase
            assert instance.data_type == "dna"

    def test_init_creates_nexus_file(self):
        """Test that NEXUS file is created."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6
            )

            assert instance.nexus_file_path.exists()
            assert instance.nexus_file_path.suffix == ".nex"

            # Check NEXUS file has content
            assert instance.nexus_file_path.stat().st_size > 0


class TestModelConversion:
    """Test model string conversion to PAUP* commands."""

    def test_build_model_string_dna_nst1(self):
        """Test model string for JC-type DNA model."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=1
            )

            assert "JC-type" in instance.model_str

    def test_build_model_string_dna_nst2(self):
        """Test model string for HKY-type DNA model."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=2
            )

            assert "HKY-type" in instance.model_str

    def test_build_model_string_dna_nst6(self):
        """Test model string for GTR-type DNA model."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6
            )

            assert "GTR-type" in instance.model_str

    def test_build_model_string_with_gamma(self):
        """Test model string with gamma rate variation."""
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
                gamma=True
            )

            assert "+G" in instance.model_str

    def test_build_model_string_with_invariable(self):
        """Test model string with invariable sites."""
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
                invariable=True
            )

            assert "+I" in instance.model_str

    def test_build_model_string_protein(self):
        """Test model string for protein data."""
        alignment_file = Path("tests/fixtures/alignments/protein_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="protein",
                protein="wag"
            )

            assert "WAG" in instance.model_str


class TestAnalysisMode:
    """Test analysis mode parsing."""

    def test_ml_only(self):
        """Test ML-only analysis mode."""
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
                analysis_mode="ml"
            )

            assert instance.do_ml is True
            assert instance.do_bayesian is False
            assert instance.do_parsimony is False

    def test_bayesian_only(self):
        """Test Bayesian-only analysis mode."""
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
                analysis_mode="bayesian"
            )

            assert instance.do_ml is False
            assert instance.do_bayesian is True
            assert instance.do_parsimony is False

    def test_parsimony_only(self):
        """Test parsimony-only analysis mode."""
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
                analysis_mode="parsimony"
            )

            assert instance.do_ml is False
            assert instance.do_bayesian is False
            assert instance.do_parsimony is True

    def test_all_analyses(self):
        """Test all analysis modes enabled."""
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
                analysis_mode="all"
            )

            assert instance.do_ml is True
            assert instance.do_bayesian is True
            assert instance.do_parsimony is True

    def test_combined_ml_parsimony(self):
        """Test combined ML and parsimony mode."""
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
                analysis_mode="ml+parsimony"
            )

            assert instance.do_ml is True
            assert instance.do_bayesian is False
            assert instance.do_parsimony is True


class TestNexusConversion:
    """Test alignment to NEXUS conversion."""

    def test_nexus_conversion_dna(self):
        """Test NEXUS conversion for DNA alignment."""
        alignment_file = Path("tests/fixtures/alignments/dna_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6
            )

            # Check NEXUS file was created
            assert instance.nexus_file_path.exists()

            # Read and check content
            content = instance.nexus_file_path.read_text()
            assert "#NEXUS" in content
            assert "DATATYPE=DNA" in content
            assert "NTAX=" in content
            assert "NCHAR=" in content

    def test_nexus_conversion_protein(self):
        """Test NEXUS conversion for protein alignment."""
        alignment_file = Path("tests/fixtures/alignments/protein_small.fas")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="fasta",
                temp_dir=Path(tmpdir),
                data_type="protein",
                protein="wag"
            )

            content = instance.nexus_file_path.read_text()
            assert "#NEXUS" in content
            assert "DATATYPE=PROTEIN" in content

    def test_nexus_already_nexus(self):
        """Test handling of already-NEXUS input."""
        alignment_file = Path("tests/fixtures/alignments/dna_medium.nex")

        if not alignment_file.exists():
            pytest.skip("Test fixture not found")

        with tempfile.TemporaryDirectory() as tmpdir:
            instance = panDecayIndices(
                alignment_file=alignment_file,
                alignment_format="nexus",
                temp_dir=Path(tmpdir),
                data_type="dna",
                nst=6
            )

            # Should still create a NEXUS file (copy)
            assert instance.nexus_file_path.exists()

            content = instance.nexus_file_path.read_text()
            assert "#NEXUS" in content
