#!/usr/bin/env python3
"""
Configuration loader for panDecay supporting YAML, TOML, and legacy INI formats.

This module provides utilities to load and validate configuration files
using the Pydantic models defined in config_models.py.
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, Any, Union, Optional
import configparser

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    
try:
    import toml
    HAS_TOML = True
except ImportError:
    HAS_TOML = False

from .config_models import (
    PanDecayConfig, InputOutputConfig, AnalysisConfig, ModelConfig,
    ComputationalConfig, BayesianConfig, VisualizationConfig, ConstraintConfig
)

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Exception raised for configuration-related errors."""
    pass


def detect_config_format(config_path: Path) -> str:
    """Detect configuration file format based on extension and content."""
    
    suffix = config_path.suffix.lower()
    
    # Prioritize file extension for known formats
    if suffix in ['.yaml', '.yml']:
        return 'yaml'
    elif suffix in ['.toml']:
        return 'toml'
    elif suffix in ['.ini', '.cfg', '.config', '.conf']:
        return 'ini'
    
    # Try to detect by content only if extension is unknown
    try:
        with open(config_path, 'r') as f:
            content = f.read().strip()
        
        # Check for YAML indicators
        if content.startswith(('---', '%YAML')) or ':\n' in content or ': ' in content:
            return 'yaml'
        
        # Check for TOML indicators
        if '[' in content and ']' in content and '=' in content:
            # Could be TOML or INI, check for TOML-specific syntax
            if '[[' in content or content.count('"') > content.count("'"):
                return 'toml'
            else:
                return 'ini'
                
        # Default to INI for legacy compatibility
        return 'ini'
        
    except Exception:
        # If we can't read the file, default to INI
        return 'ini'


def load_yaml_config(config_path: Path) -> Dict[str, Any]:
    """Load YAML configuration file."""
    
    if not HAS_YAML:
        raise ConfigurationError(
            "PyYAML is required for YAML configuration files. "
            "Install with: pip install pyyaml"
        )
    
    try:
        with open(config_path, 'r') as f:
            data = yaml.safe_load(f)
        
        if data is None:
            data = {}
            
        return data
        
    except yaml.YAMLError as e:
        error_msg = f"Error parsing YAML configuration in {config_path}:"
        error_msg += f"\n  → {e}"
        error_msg += f"\n  → Make sure the file uses proper YAML syntax (check indentation, colons, etc.)"
        error_msg += f"\n  → If this is an INI file, use .ini or .conf extension for proper detection"
        raise ConfigurationError(error_msg)
    except Exception as e:
        raise ConfigurationError(f"Error loading YAML configuration from {config_path}: {e}")


def load_toml_config(config_path: Path) -> Dict[str, Any]:
    """Load TOML configuration file."""
    
    if not HAS_TOML:
        raise ConfigurationError(
            "toml is required for TOML configuration files. "
            "Install with: pip install toml"
        )
    
    try:
        with open(config_path, 'r') as f:
            data = toml.load(f)
        
        return data
        
    except toml.TomlDecodeError as e:
        raise ConfigurationError(f"Error parsing TOML configuration: {e}")
    except Exception as e:
        raise ConfigurationError(f"Error loading TOML configuration: {e}")


def load_ini_config(config_path: Path) -> Dict[str, Any]:
    """Load legacy INI configuration file and convert to new format."""
    
    config = configparser.ConfigParser()
    
    try:
        config.read(config_path)
    except configparser.Error as e:
        error_msg = f"Error parsing INI configuration in {config_path}:"
        error_msg += f"\n  → {e}"
        if "no section headers" in str(e).lower():
            error_msg += f"\n  → INI files require section headers like [DEFAULT]"
            error_msg += f"\n  → Try regenerating the config with: --generate-ini-config"
        elif "duplicate section" in str(e).lower():
            error_msg += f"\n  → Duplicate section names found in INI file"
        error_msg += f"\n  → For better validation, consider using YAML format instead"
        raise ConfigurationError(error_msg)
    except Exception as e:
        raise ConfigurationError(f"Error loading INI configuration from {config_path}: {e}")
    
    # Convert INI format to new nested dictionary format
    data = {}
    
    # Helper function to convert string values
    def convert_value(value: str) -> Union[str, bool, int, float]:
        """Convert string values to appropriate types."""
        value = value.strip()
        
        # Boolean conversion
        if value.lower() in ('true', 'yes', 'on', '1'):
            return True
        elif value.lower() in ('false', 'no', 'off', '0'):
            return False
        
        # Numeric conversion
        try:
            if '.' in value:
                return float(value)
            else:
                return int(value)
        except ValueError:
            pass
        
        # String value
        return value
    
    # Map INI sections and keys to new configuration structure
    
    # Input/Output settings (default section)
    input_output = {}
    main_section = config['DEFAULT'] if 'DEFAULT' in config else config
    
    for key, value in main_section.items():
        if key in ['alignment', 'alignment_file']:
            input_output['alignment_file'] = value
        elif key == 'format':
            input_output['alignment_format'] = value
        elif key == 'data_type':
            input_output['data_type'] = value
        elif key == 'output':
            input_output['output_prefix'] = value.replace('.txt', '')
        elif key == 'tree':
            input_output['tree_prefix'] = value
        elif key == 'temp':
            input_output['temp_directory'] = value
        elif key == 'keep_files':
            input_output['keep_files'] = convert_value(value)
        elif key == 'debug':
            input_output['debug'] = convert_value(value)
    
    if input_output:
        data['input_output'] = input_output
    
    # Analysis settings
    analysis = {}
    if 'analysis' in main_section:
        analysis_value = main_section['analysis']
        if '+' in analysis_value:
            analysis['analysis_types'] = analysis_value.split('+')
        elif analysis_value == 'all':
            analysis['analysis_types'] = ['ml', 'bayesian', 'parsimony']
        else:
            analysis['analysis_types'] = [analysis_value]
    
    for key, value in main_section.items():
        if key == 'bootstrap':
            analysis['bootstrap'] = convert_value(value)
        elif key == 'bootstrap_reps':
            analysis['bootstrap_reps'] = convert_value(value)
        elif key == 'site_analysis':
            analysis['site_analysis'] = convert_value(value)
    
    if analysis:
        data['analysis'] = analysis
    
    # Model settings
    model = {}
    for key, value in main_section.items():
        if key == 'model':
            # Determine data type and set appropriate model
            data_type = input_output.get('data_type', 'dna')
            if data_type == 'dna':
                model['dna_model'] = value
            elif data_type == 'protein':
                model['protein_model'] = value
            elif data_type == 'discrete':
                model['discrete_model'] = value
        elif key == 'gamma':
            model['gamma'] = convert_value(value)
        elif key == 'gamma_shape':
            model['gamma_shape'] = convert_value(value)
        elif key == 'invariant':
            model['invariant'] = convert_value(value)
        elif key == 'prop_invar':
            model['prop_invar'] = convert_value(value)
        elif key == 'base_freq':
            model['base_freq'] = value
        elif key == 'nst':
            model['nst'] = convert_value(value)
    
    if model:
        data['model'] = model
    
    # Computational settings
    computational = {}
    for key, value in main_section.items():
        if key == 'threads':
            computational['threads'] = convert_value(value) if value.isdigit() else value
        elif key == 'paup':
            computational['paup_path'] = value
        elif key == 'starting_tree':
            computational['starting_tree'] = value
        elif key == 'paup_block':
            computational['paup_block'] = value
    
    if computational:
        data['computational'] = computational
    
    # Bayesian settings
    bayesian = {}
    for key, value in main_section.items():
        if key == 'bayesian_software':
            bayesian['software'] = value
        elif key == 'mrbayes_path':
            bayesian['mrbayes_path'] = value
        elif key == 'bayes_ngen':
            bayesian['ngen'] = convert_value(value)
        elif key == 'bayes_burnin':
            bayesian['burnin'] = convert_value(value)
        elif key == 'bayes_chains':
            bayesian['chains'] = convert_value(value)
        elif key == 'bayes_sample_freq':
            bayesian['sample_freq'] = convert_value(value)
        elif key == 'marginal_likelihood':
            bayesian['marginal_likelihood'] = value
        elif key == 'ss_alpha':
            bayesian['ss_alpha'] = convert_value(value)
        elif key == 'ss_nsteps':
            bayesian['ss_nsteps'] = convert_value(value)
        elif key == 'check_convergence':
            bayesian['check_convergence'] = convert_value(value)
        elif key == 'min_ess':
            bayesian['min_ess'] = convert_value(value)
        elif key == 'max_psrf':
            bayesian['max_psrf'] = convert_value(value)
        elif key == 'max_asdsf':
            bayesian['max_asdsf'] = convert_value(value)
        elif key == 'convergence_strict':
            bayesian['convergence_strict'] = convert_value(value)
        elif key == 'use_mpi':
            bayesian['use_mpi'] = convert_value(value)
        elif key == 'mpi_processors':
            bayesian['mpi_processors'] = convert_value(value)
        elif key == 'mpirun_path':
            bayesian['mpirun_path'] = value
        elif key == 'use_beagle':
            bayesian['use_beagle'] = convert_value(value)
        elif key == 'beagle_device':
            bayesian['beagle_device'] = value
        elif key == 'beagle_precision':
            bayesian['beagle_precision'] = value
        elif key == 'beagle_scaling':
            bayesian['beagle_scaling'] = value
    
    if bayesian:
        data['bayesian'] = bayesian
    
    # Visualization settings
    visualization = {}
    for key, value in main_section.items():
        if key == 'visualize':
            visualization['enable'] = convert_value(value)
        elif key == 'viz_format':
            # Map old format to new format
            if value in ['png', 'pdf', 'svg']:
                visualization['format'] = 'static'
            else:
                visualization['format'] = value
        elif key == 'annotation':
            visualization['annotation'] = value
    
    if visualization:
        data['visualization'] = visualization
    
    # Constraint settings
    constraints = {}
    for key, value in main_section.items():
        if key == 'constraint_mode':
            constraints['mode'] = value
        elif key == 'test_branches':
            if value:
                constraints['test_branches'] = [b.strip() for b in value.split(';')]
        elif key == 'constraint_file':
            constraints['constraint_file'] = value
    
    # Handle constraints section
    if 'constraints' in config:
        custom_constraints = {}
        for clade_name, taxa_list in config['constraints'].items():
            custom_constraints[clade_name] = [t.strip() for t in taxa_list.split(',')]
        if custom_constraints:
            constraints['custom_constraints'] = custom_constraints
    
    if constraints:
        data['constraints'] = constraints
    
    return data


def load_configuration(config_path: Union[str, Path]) -> PanDecayConfig:
    """Load and validate configuration from file."""
    
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise ConfigurationError(f"Configuration file not found: {config_path}")
    
    # Detect format and load data
    format_type = detect_config_format(config_path)
    logger.info(f"Loading {format_type.upper()} configuration from: {config_path}")
    
    if format_type == 'yaml':
        data = load_yaml_config(config_path)
    elif format_type == 'toml':
        data = load_toml_config(config_path)
    elif format_type == 'ini':
        data = load_ini_config(config_path)
        logger.warning(
            "INI configuration format is deprecated. "
            "Consider migrating to YAML or TOML format for better features."
        )
    else:
        raise ConfigurationError(f"Unsupported configuration format: {format_type}")
    
    # Validate and create configuration object
    try:
        config = PanDecayConfig(**data)
        logger.info("Configuration loaded and validated successfully")
        return config
        
    except Exception as e:
        # Provide more helpful error messages based on the type of error
        error_msg = f"Configuration validation failed for {config_path} ({format_type.upper()} format)"
        
        if "alignment_file" in str(e):
            error_msg += f"\n  → The specified alignment file was not found. Please check the file path."
        elif "analysis_types" in str(e):
            error_msg += f"\n  → Invalid analysis type specified. Valid options: ml, bayesian, parsimony"
        elif "normalization" in str(e):
            error_msg += f"\n  → Invalid normalization level. Valid options: none, basic, full"
        elif "threads" in str(e):
            error_msg += f"\n  → Invalid thread count. Use 'auto', 'all', or a positive integer"
        else:
            error_msg += f"\n  → {e}"
            
        error_msg += f"\n  → See the configuration guide for help: Use --generate-yaml-config for examples"
        
        raise ConfigurationError(error_msg)


def create_example_yaml_config(output_path: Path) -> None:
    """Create an example YAML configuration file."""
    
    example_config = {
        'input_output': {
            'alignment_file': 'my_sequences.fasta',
            'alignment_format': 'fasta',
            'data_type': 'dna',
            'output_prefix': 'pan_decay_analysis',
            'tree_prefix': 'annotated_tree',
            'keep_files': False,
            'debug': False
        },
        'analysis': {
            'analysis_types': ['ml', 'bayesian'],
            'bootstrap': False,
            'bootstrap_reps': 100,
            'site_analysis': False
        },
        'model': {
            'dna_model': 'GTR',
            'gamma': True,
            'invariant': False,
            'base_freq': 'estimate'
        },
        'computational': {
            'threads': 'auto',
            'paup_path': 'paup',
            'ml_timeout': 3600,
            'constraint_timeout': 600
        },
        'bayesian': {
            'software': 'mrbayes',
            'mrbayes_path': 'mb',
            'ngen': 1000000,
            'burnin': 0.25,
            'chains': 4,
            'sample_freq': 1000,
            'marginal_likelihood': 'ss',
            'check_convergence': True,
            'min_ess': 200,
            'max_psrf': 1.01,
            'max_asdsf': 0.01
        },
        'visualization': {
            'enable': True,
            'format': 'both',
            'static': {
                'dpi': 300,
                'formats': ['png', 'pdf'],
                'style': 'publication'
            },
            'interactive': {
                'theme': 'plotly_white',
                'export_html': True,
                'include_controls': True
            }
        },
        'constraints': {
            'mode': 'all'
        }
    }
    
    with open(output_path, 'w') as f:
        yaml.dump(example_config, f, default_flow_style=False, sort_keys=False, indent=2)
    
    logger.info(f"Example YAML configuration created: {output_path}")


def create_example_toml_config(output_path: Path) -> None:
    """Create an example TOML configuration file."""
    
    if not HAS_TOML:
        raise ConfigurationError(
            "toml is required for TOML configuration files. "
            "Install with: pip install toml"
        )
    
    example_config = {
        'input_output': {
            'alignment_file': 'my_sequences.fasta',
            'alignment_format': 'fasta',
            'data_type': 'dna',
            'output_prefix': 'pan_decay_analysis',
            'tree_prefix': 'annotated_tree',
            'keep_files': False,
            'debug': False
        },
        'analysis': {
            'analysis_types': ['ml', 'bayesian'],
            'bootstrap': False,
            'bootstrap_reps': 100,
            'site_analysis': False
        },
        'model': {
            'dna_model': 'GTR',
            'gamma': True,
            'invariant': False,
            'base_freq': 'estimate'
        },
        'computational': {
            'threads': 'auto',
            'paup_path': 'paup',
            'ml_timeout': 3600,
            'constraint_timeout': 600
        },
        'bayesian': {
            'software': 'mrbayes',
            'mrbayes_path': 'mb',
            'ngen': 1000000,
            'burnin': 0.25,
            'chains': 4,
            'sample_freq': 1000,
            'marginal_likelihood': 'ss',
            'check_convergence': True,
            'min_ess': 200,
            'max_psrf': 1.01,
            'max_asdsf': 0.01
        },
        'visualization': {
            'enable': True,
            'format': 'both',
            'static': {
                'dpi': 300,
                'formats': ['png', 'pdf'],
                'style': 'publication'
            },
            'interactive': {
                'theme': 'plotly_white',
                'export_html': True,
                'include_controls': True
            }
        },
        'constraints': {
            'mode': 'all'
        }
    }
    
    with open(output_path, 'w') as f:
        toml.dump(example_config, f)
    
    logger.info(f"Example TOML configuration created: {output_path}")


# Command-line utility functions
def main():
    """Command-line interface for configuration utilities."""
    import argparse
    
    parser = argparse.ArgumentParser(description="panDecay configuration utilities")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate configuration file')
    validate_parser.add_argument('config_file', help='Configuration file to validate')
    
    # Example command
    example_parser = subparsers.add_parser('example', help='Create example YAML configuration')
    example_parser.add_argument('output_file', help='Output YAML configuration file')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        try:
            config = load_configuration(args.config_file)
            print(f"✓ Configuration file {args.config_file} is valid")
            print(f"  Analysis types: {config.analysis.analysis_types}")
            print(f"  Data type: {config.input_output.data_type}")
            print(f"  Model: {config.get_model_for_data_type(config.input_output.data_type)}")
        except ConfigurationError as e:
            print(f"✗ Configuration validation failed: {e}")
            sys.exit(1)
    
    elif args.command == 'example':
        try:
            create_example_yaml_config(Path(args.output_file))
            print(f"✓ Example configuration created: {args.output_file}")
        except Exception as e:
            print(f"✗ Failed to create example: {e}")
            sys.exit(1)
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()