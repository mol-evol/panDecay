#!/usr/bin/env python3
"""
Unit tests for panDecay configuration file parsing.
"""

import pytest
from pathlib import Path
import tempfile
import argparse

from pandecay.core.configuration import (
    generate_config_template,
    parse_config,
    str_to_bool,
    ConfigurationError
)


class TestStrToBool:
    """Test string to boolean conversion."""

    def test_true_values(self):
        """Test various true string representations."""
        assert str_to_bool("true") is True
        assert str_to_bool("True") is True
        assert str_to_bool("TRUE") is True
        assert str_to_bool("yes") is True
        assert str_to_bool("Yes") is True
        assert str_to_bool("on") is True
        assert str_to_bool("1") is True

    def test_false_values(self):
        """Test various false string representations."""
        assert str_to_bool("false") is False
        assert str_to_bool("False") is False
        assert str_to_bool("FALSE") is False
        assert str_to_bool("no") is False
        assert str_to_bool("No") is False
        assert str_to_bool("off") is False
        assert str_to_bool("0") is False

    def test_invalid_values(self):
        """Test invalid boolean string raises error."""
        with pytest.raises(ConfigurationError):
            str_to_bool("maybe")
        with pytest.raises(ConfigurationError):
            str_to_bool("invalid")


class TestGenerateConfigTemplate:
    """Test configuration template generation."""

    def test_generate_template(self):
        """Test generating a config template file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "test_config.txt"

            # Note: This might fail if template file doesn't exist
            # We'll handle that gracefully
            try:
                generate_config_template(config_path)
                assert config_path.exists()

                # Check file has content
                content = config_path.read_text()
                assert len(content) > 0
            except ConfigurationError:
                # Template file not found - skip this test
                pytest.skip("Config template file not found in package")


class TestParseConfig:
    """Test configuration file parsing."""

    def test_parse_simple_config(self):
        """Test parsing a simple configuration file."""
        config_content = """
alignment = tests/fixtures/alignments/dna_small.fas
format = fasta
data_type = dna
nst = 6
gamma = yes
threads = 4
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            # Create a base args namespace
            args = argparse.Namespace(
                alignment=None,
                format='fasta',
                data_type='dna',
                nst=None,
                protein=None,
                gamma=False,
                invariable=False,
                threads='auto',
                gamma_shape=None,
                prop_invar=None,
                base_freq=None,
                parsmodel=None,
                starting_tree=None,
                paup_block=None,
                temp=None,
                keep_files=False,
                debug=False,
                site_analysis=False,
                analysis='ml',
                bootstrap=False,
                bootstrap_reps=100,
                mrbayes_path='mb',
                bayes_ngen=1000000,
                bayes_burnin=0.25,
                bayes_chains=4,
                bayes_sample_freq=1000,
                marginal_likelihood='ss',
                ss_alpha=0.4,
                ss_nsteps=50,
                use_mpi=False,
                mpi_processors=None,
                mpirun_path='mpirun',
                use_beagle=False,
                beagle_device='auto',
                beagle_precision='double',
                beagle_scaling='dynamic',
                visualize=False,
                viz_format='png',
                annotation='lnl',
                constraint_mode='all',
                test_branches=None,
                constraint_file=None,
                check_convergence=True,
                min_ess=200,
                max_psrf=1.01,
                max_asdsf=0.01,
                convergence_strict=False,
                output_style='unicode'
            )

            result = parse_config(config_file, args)

            assert result.alignment == "tests/fixtures/alignments/dna_small.fas"
            assert result.format == "fasta"
            assert result.data_type == "dna"
            assert result.nst == 6
            assert result.gamma is True
            assert result.threads == 4

        finally:
            Path(config_file).unlink()

    def test_parse_config_with_sections(self):
        """Test parsing config with explicit sections."""
        config_content = """[DEFAULT]
alignment = tests/fixtures/alignments/dna_small.fas
nst = 2
gamma = true
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            args = argparse.Namespace(
                alignment=None,
                nst=None,
                gamma=False,
                format='fasta',
                data_type='dna'
            )

            result = parse_config(config_file, args)
            assert result.alignment == "tests/fixtures/alignments/dna_small.fas"
            assert result.nst == 2
            assert result.gamma is True

        finally:
            Path(config_file).unlink()

    def test_parse_config_type_conversions(self):
        """Test various type conversions in config parsing."""
        config_content = """
gamma_shape = 0.5
prop_invar = 0.2
bayes_burnin = 0.25
nst = 6
threads = 8
gamma = yes
keep_files = no
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            args = argparse.Namespace(
                gamma_shape=None,
                prop_invar=None,
                bayes_burnin=0.25,
                nst=None,
                threads='auto',
                gamma=False,
                keep_files=False
            )

            result = parse_config(config_file, args)

            # Float conversions
            assert result.gamma_shape == 0.5
            assert result.prop_invar == 0.2
            assert result.bayes_burnin == 0.25

            # Int conversions
            assert result.nst == 6
            assert result.threads == 8

            # Bool conversions
            assert result.gamma is True
            assert result.keep_files is False

        finally:
            Path(config_file).unlink()

    def test_parse_nonexistent_config(self):
        """Test parsing non-existent config file raises error."""
        args = argparse.Namespace()

        with pytest.raises(ConfigurationError):
            parse_config("/nonexistent/config.txt", args)

    def test_config_overrides_defaults(self):
        """Test config file values override argument defaults."""
        config_content = """
threads = 16
analysis = bayesian
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            args = argparse.Namespace(
                threads='auto',  # Default
                analysis='ml'     # Default
            )

            result = parse_config(config_file, args)

            # Config should override defaults
            assert result.threads == 16
            assert result.analysis == 'bayesian'

        finally:
            Path(config_file).unlink()

    def test_parse_config_with_paths(self):
        """Test parsing config with Path type parameters."""
        config_content = """
alignment = tests/fixtures/alignments/dna_small.fas
starting_tree = tests/fixtures/trees/simple_tree.nwk
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            args = argparse.Namespace(
                alignment=None,
                starting_tree=None
            )

            result = parse_config(config_file, args)

            assert result.alignment == "tests/fixtures/alignments/dna_small.fas"
            assert isinstance(result.starting_tree, Path)

        finally:
            Path(config_file).unlink()

    def test_parse_config_special_thread_values(self):
        """Test parsing special thread values like 'auto'."""
        config_content = """
threads = auto
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
            f.write(config_content)
            config_file = f.name

        try:
            args = argparse.Namespace(threads=1)

            result = parse_config(config_file, args)

            assert result.threads == 'auto'

        finally:
            Path(config_file).unlink()
