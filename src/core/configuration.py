#!/usr/bin/env python3
"""
Configuration file handling for panDecay.

Extracted from the original monolithic system to preserve
exact configuration parsing behavior.
"""

import sys
import logging
import configparser
import argparse
from pathlib import Path
from typing import Union

logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    """Custom exception for configuration-related errors."""
    pass


def generate_config_template(filepath: Union[str, Path]) -> None:
    """Generate a template configuration file with all options and comments."""
    # Get the path to the template file relative to this module
    template_path = Path(__file__).parent / "config_template.txt"
    
    try:
        # Read the template from file
        with open(template_path, 'r') as template_file:
            template = template_file.read()
    except FileNotFoundError:
        msg = f"Configuration template file not found: {template_path}"
        logger.error(msg)
        raise ConfigurationError(msg)
    except Exception as e:
        msg = f"Failed to read configuration template: {e}"
        logger.error(msg)
        raise ConfigurationError(msg)
    
    try:
        with open(filepath, 'w') as f:
            f.write(template)
        logger.info(f"Template configuration file generated: {filepath}")
        logger.info("Edit this file with your parameters and run:")
        logger.info(f"  python panDecay.py --config {filepath}")
    except Exception as e:
        msg = f"Failed to generate config template: {e}"
        logger.error(msg)
        raise ConfigurationError(msg)


def str_to_bool(value: str) -> bool:
    """Convert string to boolean."""
    if value.lower() in ('true', 'yes', 'on', '1'):
        return True
    elif value.lower() in ('false', 'no', 'off', '0'):
        return False
    else:
        raise ConfigurationError(f"Cannot convert '{value}' to boolean")


def parse_config(config_file: Union[str, Path], args: argparse.Namespace) -> argparse.Namespace:
    """Parse configuration file and update args namespace with values."""
    config = configparser.ConfigParser(allow_no_value=True)
    
    try:
        # Read config file, prepending a [DEFAULT] section if none exists
        with open(config_file, 'r') as f:
            config_content = f.read()
        
        # Check if file starts with a section header
        if not config_content.strip().startswith('['):
            # Prepend DEFAULT section for files without explicit sections
            config_content = '[DEFAULT]\n' + config_content
        
        config.read_string(config_content)
        
        # Process main section (unnamed section at top of file)
        if config.has_section('DEFAULT') or len(config.sections()) > 0 or len(config.defaults()) > 0:
            # Get all items from DEFAULT section and any named sections
            items = dict(config.defaults())
            
            # Map config file parameters to argparse arguments
            param_map = {
                'alignment': 'alignment',
                'format': 'format',
                'model': 'model',
                'gamma': 'gamma',
                'invariable': 'invariable',
                'paup': 'paup',
                'output': 'output',
                'tree': 'tree',
                'data_type': 'data_type',
                'gamma_shape': 'gamma_shape',
                'prop_invar': 'prop_invar',
                'base_freq': 'base_freq',
                'rates': 'rates',
                'protein_model': 'protein_model',
                'nst': 'nst',
                'parsmodel': 'parsmodel',
                'threads': 'threads',
                'starting_tree': 'starting_tree',
                'paup_block': 'paup_block',
                'temp': 'temp',
                'keep_files': 'keep_files',
                'debug': 'debug',
                'site_analysis': 'site_analysis',
                'analysis': 'analysis',
                'bootstrap': 'bootstrap',
                'bootstrap_reps': 'bootstrap_reps',
                'bayesian_software': 'bayesian_software',
                'mrbayes_path': 'mrbayes_path',
                'bayes_model': 'bayes_model',
                'bayes_ngen': 'bayes_ngen',
                'bayes_burnin': 'bayes_burnin',
                'bayes_chains': 'bayes_chains',
                'bayes_sample_freq': 'bayes_sample_freq',
                'marginal_likelihood': 'marginal_likelihood',
                'ss_alpha': 'ss_alpha',
                'ss_nsteps': 'ss_nsteps',
                'use_mpi': 'use_mpi',
                'mpi_processors': 'mpi_processors',
                'mpirun_path': 'mpirun_path',
                'use_beagle': 'use_beagle',
                'beagle_device': 'beagle_device',
                'beagle_precision': 'beagle_precision',
                'beagle_scaling': 'beagle_scaling',
                'visualize': 'visualize',
                'viz_format': 'viz_format',
                'annotation': 'annotation',
                'constraint_mode': 'constraint_mode',
                'test_branches': 'test_branches',
                'constraint_file': 'constraint_file',
                'check_convergence': 'check_convergence',
                'min_ess': 'min_ess',
                'max_psrf': 'max_psrf',
                'max_asdsf': 'max_asdsf',
                'convergence_strict': 'convergence_strict',
                'mrbayes_parse_timeout': 'mrbayes_parse_timeout',
                'output_style': 'output_style',
            }
            
            # Process each parameter
            for config_param, arg_param in param_map.items():
                if config_param in items:
                    value = items[config_param]
                    
                    # Skip if command line already provided this argument
                    if hasattr(args, arg_param):
                        current_value = getattr(args, arg_param)
                        # For alignment, only skip if it was provided on command line (not None)
                        if arg_param == 'alignment' and current_value is not None:
                            continue  # Command line alignment takes precedence
                        
                    # Convert types as needed
                    if arg_param in ['gamma', 'invariable', 'keep_files', 'debug', 
                                     'site_analysis', 'bootstrap', 'use_mpi', 'use_beagle',
                                     'visualize', 'check_convergence', 'convergence_strict']:
                        value = str_to_bool(value)
                    elif arg_param in ['gamma_shape', 'prop_invar', 'bayes_burnin', 'ss_alpha', 
                                       'max_psrf', 'max_asdsf', 'mrbayes_parse_timeout']:
                        value = float(value)
                    elif arg_param in ['nst', 'threads', 'bootstrap_reps', 'bayes_ngen', 
                                       'bayes_chains', 'bayes_sample_freq', 'ss_nsteps', 
                                       'mpi_processors', 'min_ess']:
                        if value != 'auto' and value != 'all':  # Special thread values
                            value = int(value)
                    elif arg_param in ['starting_tree', 'paup_block', 'constraint_file']:
                        if value:
                            value = Path(value)
                    
                    # Set the argument value
                    setattr(args, arg_param, value)
                    
            # Handle constraints section
            if config.has_section('constraints'):
                constraints_dict = dict(config.items('constraints'))
                if constraints_dict:
                    # Convert to the format expected by the main code
                    args.config_constraints = constraints_dict
                    logger.info(f"Loaded {len(constraints_dict)} constraint definitions from config file")
                    
        logger.info(f"Configuration loaded from {config_file}")
        
    except Exception as e:
        msg = f"Error parsing config file {config_file}: {e}"
        logger.error(msg)
        raise ConfigurationError(msg)
    
    return args


def read_paup_block(paup_block_file_path: Path) -> str:
    """Read PAUP* block commands from file."""
    try:
        with open(paup_block_file_path, 'r') as f:
            content = f.read().strip()
        return content
    except Exception as e:
        msg = f"Failed to read PAUP* block file {paup_block_file_path}: {e}"
        logger.error(msg)
        raise ConfigurationError(msg)