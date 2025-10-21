#!/usr/bin/env python3
"""
Unit tests for panDecay CLI argument parsing and validation.
"""

import pytest
import argparse
from pathlib import Path

from pandecay.cli import (
    setup_argument_parser,
    validate_model_flags,
    normalize_data_type,
    setup_output_structure
)


class TestNormalizeDataType:
    """Test data type normalization."""

    def test_dna_variations(self):
        """Test various DNA type inputs."""
        assert normalize_data_type("dna") == "dna"
        assert normalize_data_type("DNA") == "dna"
        assert normalize_data_type("nucleotide") == "dna"
        assert normalize_data_type("nt") == "dna"

    def test_protein_variations(self):
        """Test various protein type inputs."""
        assert normalize_data_type("protein") == "protein"
        assert normalize_data_type("Protein") == "protein"
        assert normalize_data_type("PROTEIN") == "protein"
        assert normalize_data_type("amino") == "protein"
        assert normalize_data_type("aa") == "protein"

    def test_discrete_variations(self):
        """Test various discrete type inputs."""
        assert normalize_data_type("discrete") == "discrete"
        assert normalize_data_type("Discrete") == "discrete"
        assert normalize_data_type("DISCRETE") == "discrete"
        assert normalize_data_type("morphology") == "discrete"
        assert normalize_data_type("morph") == "discrete"

    def test_none_input(self):
        """Test None input returns default."""
        result = normalize_data_type(None)
        assert result == "dna"  # Default fallback

    def test_empty_string(self):
        """Test empty string returns default."""
        result = normalize_data_type("")
        assert result == "dna"


class TestValidateModelFlags:
    """Test model flag validation."""

    def test_dna_without_nst(self):
        """Test DNA data without nst flag raises error."""
        args = argparse.Namespace(
            data_type="dna",
            nst=None,
            protein=None
        )

        with pytest.raises(ValueError, match="--nst is REQUIRED for DNA data"):
            validate_model_flags(args)

    def test_dna_with_protein_flag(self):
        """Test DNA data with protein flag raises error."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein="wag"
        )

        with pytest.raises(ValueError, match="--protein flag cannot be used with DNA data"):
            validate_model_flags(args)

    def test_protein_without_protein_flag(self):
        """Test protein data without protein flag raises error."""
        args = argparse.Namespace(
            data_type="protein",
            nst=None,
            protein=None
        )

        with pytest.raises(ValueError, match="--protein is REQUIRED for protein data"):
            validate_model_flags(args)

    def test_protein_with_nst_flag(self):
        """Test protein data with nst flag raises error."""
        args = argparse.Namespace(
            data_type="protein",
            nst=6,
            protein="wag"
        )

        with pytest.raises(ValueError, match="--nst flag cannot be used with protein data"):
            validate_model_flags(args)

    def test_discrete_with_nst(self):
        """Test discrete data with nst flag raises error."""
        args = argparse.Namespace(
            data_type="discrete",
            nst=1,
            protein=None
        )

        with pytest.raises(ValueError, match="--nst not needed for discrete data"):
            validate_model_flags(args)

    def test_discrete_with_protein(self):
        """Test discrete data with protein flag raises error."""
        args = argparse.Namespace(
            data_type="discrete",
            nst=None,
            protein="wag"
        )

        with pytest.raises(ValueError, match="--protein not needed for discrete data"):
            validate_model_flags(args)

    def test_valid_dna_nst1(self):
        """Test valid DNA with nst=1."""
        args = argparse.Namespace(
            data_type="dna",
            nst=1,
            protein=None,
            gamma_shape=None,
            prop_invar=None
        )

        # Should not raise
        validate_model_flags(args)

    def test_valid_dna_nst2(self):
        """Test valid DNA with nst=2."""
        args = argparse.Namespace(
            data_type="dna",
            nst=2,
            protein=None,
            gamma_shape=None,
            prop_invar=None
        )

        # Should not raise
        validate_model_flags(args)

    def test_valid_dna_nst6(self):
        """Test valid DNA with nst=6."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma_shape=None,
            prop_invar=None
        )

        # Should not raise
        validate_model_flags(args)

    def test_valid_protein_wag(self):
        """Test valid protein with WAG model."""
        args = argparse.Namespace(
            data_type="protein",
            nst=None,
            protein="wag",
            gamma_shape=None,
            prop_invar=None
        )

        # Should not raise
        validate_model_flags(args)

    def test_valid_discrete(self):
        """Test valid discrete data."""
        args = argparse.Namespace(
            data_type="discrete",
            nst=None,
            protein=None,
            gamma_shape=None,
            prop_invar=None
        )

        # Should not raise
        validate_model_flags(args)

    def test_negative_gamma_shape(self):
        """Test negative gamma shape raises error."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma_shape=-0.5,
            prop_invar=None
        )

        with pytest.raises(ValueError, match="Gamma shape must be positive"):
            validate_model_flags(args)

    def test_invalid_prop_invar(self):
        """Test invalid proportion invariable raises error."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma_shape=None,
            prop_invar=1.5
        )

        with pytest.raises(ValueError, match="Proportion invariable must be between 0 and 1"):
            validate_model_flags(args)

    def test_prop_invar_negative(self):
        """Test negative proportion invariable raises error."""
        args = argparse.Namespace(
            data_type="dna",
            nst=6,
            protein=None,
            gamma_shape=None,
            prop_invar=-0.1
        )

        with pytest.raises(ValueError, match="Proportion invariable must be between 0 and 1"):
            validate_model_flags(args)

    def test_invalid_data_type(self):
        """Test invalid data type raises error."""
        args = argparse.Namespace(
            data_type="invalid_type",
            nst=None,
            protein=None,
            gamma_shape=None,
            prop_invar=None
        )

        with pytest.raises(ValueError, match="Invalid data type"):
            validate_model_flags(args)


class TestSetupArgumentParser:
    """Test argument parser setup."""

    def test_parser_creation(self):
        """Test parser is created successfully."""
        parser = setup_argument_parser()
        assert isinstance(parser, argparse.ArgumentParser)

    def test_parser_has_alignment_argument(self):
        """Test parser has alignment argument."""
        parser = setup_argument_parser()
        # Test that we can parse an alignment argument
        args = parser.parse_args(['test.fas', '--nst', '6'])
        assert args.alignment == 'test.fas'

    def test_parser_nst_choices(self):
        """Test nst argument accepts valid choices."""
        parser = setup_argument_parser()

        # Valid choices
        args = parser.parse_args(['test.fas', '--nst', '1'])
        assert args.nst == 1

        args = parser.parse_args(['test.fas', '--nst', '2'])
        assert args.nst == 2

        args = parser.parse_args(['test.fas', '--nst', '6'])
        assert args.nst == 6

    def test_parser_protein_choices(self):
        """Test protein argument accepts valid choices."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--protein', 'wag'])
        assert args.protein == 'wag'

        args = parser.parse_args(['test.fas', '--protein', 'lg'])
        assert args.protein == 'lg'

        args = parser.parse_args(['test.fas', '--protein', 'jtt'])
        assert args.protein == 'jtt'

    def test_parser_analysis_mode_choices(self):
        """Test analysis mode choices."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--nst', '6', '--analysis', 'ml'])
        assert args.analysis == 'ml'

        args = parser.parse_args(['test.fas', '--nst', '6', '--analysis', 'bayesian'])
        assert args.analysis == 'bayesian'

        args = parser.parse_args(['test.fas', '--nst', '6', '--analysis', 'all'])
        assert args.analysis == 'all'

    def test_parser_gamma_flag(self):
        """Test gamma flag."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--nst', '6', '--gamma'])
        assert args.gamma is True

        args = parser.parse_args(['test.fas', '--nst', '6'])
        assert args.gamma is False

    def test_parser_invariable_flag(self):
        """Test invariable flag."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--nst', '6', '--invariable'])
        assert args.invariable is True

        args = parser.parse_args(['test.fas', '--nst', '6'])
        assert args.invariable is False

    def test_parser_threads_default(self):
        """Test threads default value."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--nst', '6'])
        assert args.threads == 'auto'

    def test_parser_threads_custom(self):
        """Test custom threads value."""
        parser = setup_argument_parser()

        args = parser.parse_args(['test.fas', '--nst', '6', '--threads', '8'])
        assert args.threads == '8'


class TestSetupOutputStructure:
    """Test output directory structure setup."""

    def test_output_structure_default(self):
        """Test default output structure creation."""
        args = argparse.Namespace(
            project_name=None,
            output_dir=None
        )
        alignment_file = Path("test_alignment.fas")

        result = setup_output_structure(args, alignment_file)

        assert 'basename' in result
        assert result['basename'] == 'test_alignment'
        assert 'results_dir' in result
        assert 'trees_dir' in result
        assert 'site_analysis_dir' in result
        assert 'supplementary_dir' in result
        assert 'viz_dir' in result

    def test_output_structure_custom_project_name(self):
        """Test output structure with custom project name."""
        args = argparse.Namespace(
            project_name="my_project",
            output_dir=None
        )
        alignment_file = Path("test_alignment.fas")

        result = setup_output_structure(args, alignment_file)

        assert result['basename'] == 'my_project'
        assert 'my_project' in str(result['summary_txt'])

    def test_output_structure_file_paths(self):
        """Test all expected file paths are created."""
        args = argparse.Namespace(
            project_name=None,
            output_dir=None
        )
        alignment_file = Path("test.fas")

        result = setup_output_structure(args, alignment_file)

        # Check main files
        assert 'summary_txt' in result
        assert 'report_md' in result
        assert 'data_csv' in result

        # Check tree files
        assert 'tree_basic' in result
        assert 'tree_comprehensive' in result
        assert 'tree_pd' in result
        assert 'tree_ld' in result
        assert 'tree_bd' in result

        # Check site analysis files
        assert 'site_likelihoods' in result
        assert 'supporting_sites' in result
        assert 'conflicting_sites' in result

        # Check all are Path objects
        for key, value in result.items():
            if key != 'basename':
                assert isinstance(value, Path)
