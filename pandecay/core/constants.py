#!/usr/bin/env python3
"""
Constants and filename definitions for panDecay.

All constant values extracted from the original monolithic system
to ensure exact preservation of behavior.
"""

# Version
VERSION = "1.1.0"

# --- Constants for Filenames ---
NEXUS_ALIGNMENT_FN = "alignment.nex"
ML_TREE_FN = "ml_tree.tre"
ML_SCORE_FN = "ml_score.txt"
ML_SEARCH_NEX_FN = "ml_search.nex"
ML_LOG_FN = "paup_ml.log"

AU_TEST_NEX_FN = "au_test.nex"
AU_TEST_SCORE_FN = "au_test_results.txt"
AU_LOG_FN = "paup_au.log"

# Additional filename patterns used throughout the codebase
BOOTSTRAP_TREES_FN = "bootstrap_trees.tre"
BOOTSTRAP_SEARCH_NEX_FN = "bootstrap_search.nex"
BOOTSTRAP_LOG_FN = "paup_bootstrap.log"

PARSIMONY_TREES_FN = "parsimony_trees.tre"
PARSIMONY_SCORE_FN = "parsimony_score.txt"
PARSIMONY_SEARCH_NEX_FN = "parsimony_search.nex"
PARSIMONY_LOG_FN = "paup_parsimony.log"

# Constraint analysis filenames
CONSTRAINT_TREE_PREFIX = "constraint_tree_"
CONSTRAINT_SCORE_PREFIX = "constraint_score_"
CONSTRAINT_SEARCH_PREFIX = "constraint_search_"
CONSTRAINT_LOG_PREFIX = "paup_constraint_"

# Site analysis filenames
SITE_ANALYSIS_PREFIX = "site_analysis_"
SITE_LNL_PREFIX = "site_lnl_"

# MrBayes filenames
MRBAYES_LOG_FN = "mrbayes.log"

# Output filename patterns
ANNOTATED_TREE_AU = "_au.nwk"
ANNOTATED_TREE_DELTA_LNL = "_delta_lnl.nwk"
ANNOTATED_TREE_COMBINED = "_combined.nwk"
ANNOTATED_TREE_BOOTSTRAP = "_bootstrap.nwk"
ANNOTATED_TREE_COMPREHENSIVE = "_comprehensive.nwk"

# Default values for analysis parameters
DEFAULT_THREADS = "auto"
DEFAULT_MODEL = "GTR"
DEFAULT_DATA_TYPE = "dna"
DEFAULT_ALIGNMENT_FORMAT = "fasta"
DEFAULT_ANALYSIS_MODE = "ml"
DEFAULT_OUTPUT_FILE = "pan_decay_indices.txt"
DEFAULT_TREE_BASE = "annotated_tree"

# MrBayes defaults
DEFAULT_BAYES_NGEN = 1000000
DEFAULT_BAYES_BURNIN = 0.25
DEFAULT_BAYES_CHAINS = 4
DEFAULT_BAYES_SAMPLE_FREQ = 1000
DEFAULT_MARGINAL_LIKELIHOOD = "ss"
DEFAULT_SS_ALPHA = 0.4
DEFAULT_SS_NSTEPS = 50

# Bootstrap defaults
DEFAULT_BOOTSTRAP_REPS = 100

# Convergence defaults
DEFAULT_MIN_ESS = 200
DEFAULT_MAX_PSRF = 1.01
DEFAULT_MAX_ASDSF = 0.01
DEFAULT_MRBAYES_PARSE_TIMEOUT = 30.0

# Visualization defaults
DEFAULT_VIZ_FORMAT = "png"
DEFAULT_ANNOTATION = "lnl"
DEFAULT_OUTPUT_STYLE = "unicode"

# Constraint defaults
DEFAULT_CONSTRAINT_MODE = "all"

# Tool paths
DEFAULT_PAUP_PATH = "paup"
DEFAULT_MRBAYES_PATH = "mb"
DEFAULT_MPIRUN_PATH = "mpirun"

# BEAGLE defaults
DEFAULT_BEAGLE_DEVICE = "auto"
DEFAULT_BEAGLE_PRECISION = "double"
DEFAULT_BEAGLE_SCALING = "dynamic"