#!/usr/bin/env python3
"""
Unit tests for panDecay utility functions.
"""

import pytest
from pathlib import Path
import tempfile
import argparse

from pandecay.core.utils import (
    get_display_path,
    format_taxon_for_paup,
    format_support_symbol,
    build_effective_model_display,
    detect_format_from_extension,
    detect_sequence_type,
    validate_file_path,
    truncate_string
)


class TestDisplayPath:
    """Test get_display_path function."""

    def test_relative_path(self):
        """Test relative path display."""
        path = Path("tests/fixtures/alignments/dna_small.fas")
        result = get_display_path(path)
        assert "dna_small.fas" in result

    def test_absolute_path(self):
        """Test absolute path display."""
        path = Path("/home/user/panDecay/tests")
        result = get_display_path(path)
        assert isinstance(result, str)

    def test_string_input(self):
        """Test string path input."""
        result = get_display_path("tests/fixtures")
        assert isinstance(result, str)


class TestFormatTaxonForPaup:
    """Test taxon name formatting for PAUP*."""

    def test_simple_name(self):
        """Test simple taxon name without special characters."""
        result = format_taxon_for_paup("Taxon1")
        assert result == "Taxon1"

    def test_name_with_spaces(self):
        """Test taxon name with spaces gets quoted."""
        result = format_taxon_for_paup("Taxon 1")
        assert result.startswith("'")
        assert result.endswith("'")
        assert "Taxon 1" in result or "Taxon_1" in result

    def test_name_with_special_chars(self):
        """Test taxon name with NEXUS special characters."""
        result = format_taxon_for_paup("Taxon(1)")
        assert result.startswith("'")
        assert result.endswith("'")

    def test_name_with_colon(self):
        """Test taxon name with colon."""
        result = format_taxon_for_paup("Taxon:1")
        assert result.startswith("'")
        assert result.endswith("'")

    def test_name_with_single_quote(self):
        """Test taxon name with single quote gets replaced."""
        result = format_taxon_for_paup("Taxon'1")
        assert "'" in result
        # Single quote should be replaced with underscore
        assert "Taxon_1" in result


class TestFormatSupportSymbol:
    """Test support value formatting."""

    def test_high_pvalue(self):
        """Test high p-value (non-significant)."""
        result = format_support_symbol(0.99)
        assert "0.99" in result
        assert "*" not in result  # Not significant

    def test_low_pvalue(self):
        """Test low p-value (significant)."""
        result = format_support_symbol(0.01)
        assert "*" in result  # Significant

    def test_very_low_pvalue(self):
        """Test very low p-value."""
        result = format_support_symbol(0.0001)
        assert "<0.001" in result
        assert "*" in result

    def test_none_value(self):
        """Test None value."""
        result = format_support_symbol(None)
        assert result == "ns"

    def test_boundary_values(self):
        """Test boundary p-values."""
        assert "*" in format_support_symbol(0.049)
        assert "*" not in format_support_symbol(0.051)


class TestBuildEffectiveModelDisplay:
    """Test model display string building."""

    def test_dna_nst1(self):
        """Test JC-type DNA model display."""
        args = argparse.Namespace(
            data_type="dna",
            nst=1,
            protein=None,
            gamma=False,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "JC-type" in result

    def test_dna_nst2(self):
        """Test HKY-type DNA model display."""
        args = argparse.Namespace(
            data_type="dna",
            nst=2,
            protein=None,
            gamma=False,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "HKY-type" in result

    def test_dna_nst6(self):
        """Test GTR-type DNA model display."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma=False,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "GTR-type" in result

    def test_dna_with_gamma(self):
        """Test DNA model with gamma rate variation."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma=True,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "+G" in result

    def test_dna_with_invariable(self):
        """Test DNA model with invariable sites."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma=False,
            invariable=True,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "+I" in result

    def test_dna_with_gamma_and_invariable(self):
        """Test DNA model with both gamma and invariable."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma=True,
            invariable=True,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "+G" in result
        assert "+I" in result

    def test_protein_wag(self):
        """Test protein model display."""
        args = argparse.Namespace(
            data_type="protein",
            nst=None,
            protein="wag",
            gamma=False,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "WAG" in result

    def test_protein_with_gamma(self):
        """Test protein model with gamma."""
        args = argparse.Namespace(
            data_type="protein",
            nst=None,
            protein="lg",
            gamma=True,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "LG" in result
        assert "+G" in result

    def test_discrete_model(self):
        """Test discrete data model display."""
        args = argparse.Namespace(
            data_type="discrete",
            nst=None,
            protein=None,
            gamma=False,
            invariable=False,
            gamma_shape=None,
            prop_invar=None,
            base_freq=None
        )
        result = build_effective_model_display(args)
        assert "Mk" in result

    def test_model_with_fixed_params(self):
        """Test model with fixed parameters."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma=True,
            invariable=True,
            gamma_shape=0.5,
            prop_invar=0.2,
            base_freq="equal"
        )
        result = build_effective_model_display(args)
        assert "shape=0.5" in result
        assert "pinvar=0.2" in result
        assert "basefreq=equal" in result


class TestDetectFormatFromExtension:
    """Test alignment format detection from file extension."""

    def test_fasta_extensions(self):
        """Test various FASTA extensions."""
        assert detect_format_from_extension(Path("test.fas")) == "fasta"
        assert detect_format_from_extension(Path("test.fasta")) == "fasta"
        assert detect_format_from_extension(Path("test.fa")) == "fasta"
        assert detect_format_from_extension(Path("test.fna")) == "fasta"

    def test_nexus_extensions(self):
        """Test NEXUS extensions."""
        assert detect_format_from_extension(Path("test.nex")) == "nexus"
        assert detect_format_from_extension(Path("test.nexus")) == "nexus"
        assert detect_format_from_extension(Path("test.nxs")) == "nexus"

    def test_phylip_extensions(self):
        """Test PHYLIP extensions."""
        assert detect_format_from_extension(Path("test.phy")) == "phylip"
        assert detect_format_from_extension(Path("test.phylip")) == "phylip"

    def test_clustal_extensions(self):
        """Test Clustal extensions."""
        assert detect_format_from_extension(Path("test.aln")) == "clustal"
        assert detect_format_from_extension(Path("test.clustal")) == "clustal"

    def test_case_insensitive(self):
        """Test case-insensitive extension matching."""
        assert detect_format_from_extension(Path("test.FAS")) == "fasta"
        assert detect_format_from_extension(Path("test.NEX")) == "nexus"


class TestDetectSequenceType:
    """Test DNA vs protein sequence detection."""

    def test_detect_dna(self):
        """Test DNA sequence detection."""
        dna_file = Path("tests/fixtures/alignments/dna_small.fas")
        if dna_file.exists():
            result = detect_sequence_type(dna_file, "fasta")
            assert result in ["dna", None]  # None if detection fails

    def test_detect_protein(self):
        """Test protein sequence detection."""
        protein_file = Path("tests/fixtures/alignments/protein_small.fas")
        if protein_file.exists():
            result = detect_sequence_type(protein_file, "fasta")
            assert result in ["protein", None]


class TestValidateFilePath:
    """Test file path validation."""

    def test_existing_file(self):
        """Test validation of existing file."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write("test")
            temp_path = f.name

        try:
            result = validate_file_path(temp_path, must_exist=True)
            assert isinstance(result, Path)
            assert result.exists()
        finally:
            Path(temp_path).unlink()

    def test_nonexistent_file_must_exist(self):
        """Test validation of non-existent file when must_exist=True."""
        with pytest.raises(FileNotFoundError):
            validate_file_path("/nonexistent/path/file.txt", must_exist=True)

    def test_nonexistent_file_optional(self):
        """Test validation of non-existent file when must_exist=False."""
        result = validate_file_path("/nonexistent/path/file.txt", must_exist=False)
        assert isinstance(result, Path)


class TestTruncateString:
    """Test string truncation."""

    def test_short_string(self):
        """Test string shorter than max length."""
        result = truncate_string("Hello", max_length=10)
        assert result == "Hello"

    def test_long_string(self):
        """Test string longer than max length."""
        result = truncate_string("Hello World!", max_length=8)
        assert len(result) <= 8
        assert "..." in result

    def test_custom_suffix(self):
        """Test custom truncation suffix."""
        result = truncate_string("Hello World!", max_length=10, suffix=">>")
        assert ">>" in result
        assert len(result) <= 10

    def test_exact_length(self):
        """Test string exactly at max length."""
        result = truncate_string("12345", max_length=5)
        assert result == "12345"
