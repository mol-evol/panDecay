#!/usr/bin/env python3
"""
Configuration file handling for panDecay.

Extracted from the original monolithic system to preserve
exact configuration parsing behavior.
"""

import sys
import logging
import configparser
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_config_template(filepath):
    """Generate a template configuration file with all options and comments."""
    template = '''# panDecay Configuration File
# Lines starting with # are comments and will be ignored
# This file allows you to specify all parameters for a panDecay analysis
# Format: parameter = value (spaces around = are optional)

# ==============================================================================
# BASIC INPUT/OUTPUT SETTINGS
# ==============================================================================

# Input alignment file (required)
alignment = my_sequences.fasta

# Alignment format (default: fasta)
# Options: fasta, phylip, nexus, clustal, stockholm
format = fasta

# Data type (default: dna)
# Options: dna, protein, discrete
data_type = dna

# Output file prefix (default: pan_decay_indices)
# This will generate files like: myanalysis.txt, myanalysis_tree.nwk, etc.
output = pan_decay_indices.txt

# Base name for annotated tree files (default: annotated_tree)
tree = annotated_tree

# Keep temporary files after analysis (default: false)
# Options: true, false
keep_files = false

# Enable debug mode with detailed logging (default: false)
debug = false

# Custom temporary directory (optional)
# If not specified, uses system temp directory
# temp = /path/to/temp

# ==============================================================================
# ANALYSIS MODE SETTINGS
# ==============================================================================

# Analysis type (default: ml)
# Options: ml, bayesian, parsimony, ml+parsimony, bayesian+parsimony, ml+bayesian, all
# Note: 'all' runs all three analysis types (ML + Bayesian + Parsimony)
# Use '+' to combine multiple analyses (e.g., ml+parsimony)
analysis = ml

# Perform bootstrap analysis (default: false)
bootstrap = false

# Number of bootstrap replicates (default: 100)
bootstrap_reps = 100

# Perform site-specific analysis (default: false)
site_analysis = false

# Generate visualizations (default: false)
# Requires matplotlib and seaborn
visualize = false

# Visualization format (default: png)
# Options: png, pdf, svg
viz_format = png

# Annotation type for visualization (default: lnl)
# Options: au, lnl
annotation = lnl


# ==============================================================================
# MODEL SETTINGS (for ML and Bayesian analyses)
# ==============================================================================

# Base substitution model
# DNA options: JC, K80, HKY, GTR
# Protein options: JTT, WAG, LG, Dayhoff
# Discrete options: Mk
model = GTR

# Add gamma rate heterogeneity (default: false)
gamma = false

# Fixed gamma shape parameter (optional)
# If not specified, will be estimated
# gamma_shape = 0.5

# Add proportion of invariable sites (default: false)
invariable = false

# Fixed proportion of invariable sites (optional)
# If not specified, will be estimated
# prop_invar = 0.2

# Base frequencies (optional)
# Options: equal, estimate, empirical
# base_freq = empirical

# Number of substitution types for DNA (optional)
# Options: 1, 2, 6
# nst = 6

# Protein-specific model (overrides base model for protein data)
# protein_model = WAG

# ==============================================================================
# COMPUTATIONAL SETTINGS
# ==============================================================================

# Number of threads for PAUP* (default: auto)
# Options: auto (uses cores-2), all, or specific number (e.g., 8)
threads = auto

# Path to PAUP* executable (default: paup)
# Specify full path if not in system PATH
paup = paup

# Starting tree file in Newick format (optional)
# starting_tree = my_starting_tree.nwk

# Custom PAUP* commands file (optional)
# This overrides most model settings above
# paup_block = custom_paup_commands.txt

# ==============================================================================
# BAYESIAN-SPECIFIC SETTINGS
# ==============================================================================

# Bayesian software (required if analysis includes bayesian)
# Options: mrbayes
bayesian_software = mrbayes

# Path to MrBayes executable (default: mb)
mrbayes_path = mb


# Model for Bayesian analysis (optional)
# If not specified, uses same as ML model
# bayes_model = GTR+G+I

# Number of MCMC generations (default: 1000000)
bayes_ngen = 1000000

# Burnin fraction 0-1 (default: 0.25)
bayes_burnin = 0.25

# Number of MCMC chains (default: 4)
bayes_chains = 4

# Sample frequency (default: 1000)
bayes_sample_freq = 1000

# Marginal likelihood estimation method (default: ss)
# Options: ss (stepping-stone), ps (path sampling), hm (harmonic mean)
marginal_likelihood = ss

# Stepping-stone sampling parameters
ss_alpha = 0.4
ss_nsteps = 50

# ==============================================================================
# CONVERGENCE CHECKING (MrBayes only)
# ==============================================================================

# Check MCMC convergence diagnostics (default: true)
check_convergence = true

# Minimum ESS (Effective Sample Size) threshold (default: 200)
# Values below this indicate poor mixing
min_ess = 200

# Maximum PSRF (Potential Scale Reduction Factor) threshold (default: 1.01)
# Values above this indicate lack of convergence between chains
max_psrf = 1.01

# Maximum ASDSF (Average Standard Deviation of Split Frequencies) (default: 0.01)
# Values above this indicate lack of convergence between independent runs
max_asdsf = 0.01

# Fail analysis if convergence criteria not met (default: false)
# If false, warnings are issued but analysis continues
convergence_strict = false

# ==============================================================================
# PARALLEL PROCESSING (MrBayes only)
# ==============================================================================

# Use MPI version of MrBayes (default: false)
use_mpi = false

# Number of MPI processors (default: number of chains)
# mpi_processors = 8

# Path to mpirun executable (default: mpirun)
mpirun_path = mpirun

# Use BEAGLE library for acceleration (default: false)
use_beagle = false

# BEAGLE device (default: auto)
# Options: cpu, gpu, auto
beagle_device = auto

# BEAGLE precision (default: double)
# Options: single, double
beagle_precision = double

# BEAGLE scaling (default: dynamic)
# Options: none, dynamic, always
beagle_scaling = dynamic

# ==============================================================================
# CONSTRAINT SETTINGS
# ==============================================================================

# Constraint mode (default: all)
# Options: all (test all branches), specific (test listed branches), exclude
constraint_mode = all

# Specify branches to test (optional)
# Format: semicolon-separated clades, each clade is comma-separated taxa
# Example: Homo_sapiens,Pan_troglodytes;Mus_musculus,Rattus_norvegicus
# test_branches = 

# Constraint file (optional)
# File containing constraint definitions (one per line)
# constraint_file = constraints.txt

# ==============================================================================
# CONSTRAINT DEFINITIONS
# ==============================================================================
# Define specific clades to test when constraint_mode = specific
# Format: clade_name = comma-separated list of taxa
# Note: These are only used when constraint_mode = specific

[constraints]
# Example constraints (uncomment and modify as needed):
# clade1 = Homo_sapiens,Pan_troglodytes,Gorilla_gorilla
# clade2 = Rattus_norvegicus,Mus_musculus
# clade3 = Gallus_gallus,Taeniopygia_guttata

# ==============================================================================
# ADVANCED SETTINGS
# ==============================================================================

# Parsimony model for discrete data (default: true for discrete)
# parsmodel = true

# Site rate variation (overrides gamma setting)
# Options: equal, gamma
# rates = gamma

# ==============================================================================
# NOTES
# ==============================================================================
# 1. Command-line arguments override settings in this file
# 2. Paths can be absolute or relative to the working directory
# 3. Boolean values can be: true/false, yes/no, on/off, 1/0
# 4. Comments can appear on their own line or after values
# 5. Section headers like [constraints] are optional for organization
# 6. Empty lines are ignored
'''
    
    try:
        with open(filepath, 'w') as f:
            f.write(template)
        logger.info(f"Template configuration file generated: {filepath}")
        logger.info("Edit this file with your parameters and run:")
        logger.info(f"  python panDecay.py --config {filepath}")
    except Exception as e:
        logger.error(f"Failed to generate config template: {e}")
        sys.exit(1)


def str_to_bool(value):
    """Convert string to boolean."""
    if value.lower() in ('true', 'yes', 'on', '1'):
        return True
    elif value.lower() in ('false', 'no', 'off', '0'):
        return False
    else:
        raise ValueError(f"Cannot convert '{value}' to boolean")


def parse_config(config_file, args):
    """Parse configuration file and update args namespace with values."""
    config = configparser.ConfigParser()
    
    try:
        config.read(config_file)
        
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
            }
            
            # Process each parameter
            for config_param, arg_param in param_map.items():
                if config_param in items:
                    value = items[config_param]
                    
                    # Skip if command line already provided this argument
                    if hasattr(args, arg_param):
                        current_value = getattr(args, arg_param)
                        # Check if it's a default value or explicitly set
                        if arg_param == 'alignment':
                            continue  # Special case - required positional arg
                        
                    # Convert types as needed
                    if arg_param in ['gamma', 'invariable', 'keep_files', 'debug', 
                                     'site_analysis', 'bootstrap', 'use_mpi', 'use_beagle',
                                     'visualize', 'check_convergence', 'convergence_strict']:
                        value = str_to_bool(value)
                    elif arg_param in ['gamma_shape', 'prop_invar', 'bayes_burnin', 'ss_alpha', 
                                       'max_psrf', 'max_asdsf']:
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
        logger.error(f"Error parsing config file {config_file}: {e}")
        sys.exit(1)
    
    return args


def read_paup_block(paup_block_file_path: Path):
    """Read PAUP* block commands from file."""
    try:
        with open(paup_block_file_path, 'r') as f:
            content = f.read().strip()
        return content
    except Exception as e:
        logger.error(f"Failed to read PAUP* block file {paup_block_file_path}: {e}")
        sys.exit(1)