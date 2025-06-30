# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

panDecay is a Python command-line tool for calculating phylogenetic decay indices across multiple analysis frameworks. It can compute parsimony-based decay indices (traditional Bremer support), Maximum Likelihood (ML)-based decay indices, and Bayesian decay indices using MrBayes, assessing the robustness of clades in a phylogenetic tree through multiple statistical frameworks.

## Required Dependencies

- Python 3.8+
- PAUP* command-line executable (required for parsimony and ML analyses)
- MrBayes command-line executable (required for Bayesian analysis)
- Python packages:
  - biopython>=1.79
  - numpy>=1.20.0
  - matplotlib>=3.5.0 (optional, for visualizations)
  - seaborn>=0.11.0 (optional, for visualizations)

## Common Commands

### Installing Dependencies

```bash
pip install -r requirements.txt
```

### Running panDecay

Basic command:
```bash
python panDecay.py <alignment_file> --model <model_name> [options...]
```

Example with parsimony analysis (traditional Bremer support):
```bash
python panDecay.py alignment.fas --analysis parsimony
```

Example with ML analysis (default):
```bash
python panDecay.py alignment.fas --model GTR --gamma --invariable --data-type dna
```

Example with combined ML and parsimony:
```bash
python panDecay.py alignment.fas --analysis ml+parsimony --model GTR --gamma
```

Example with Bayesian analysis only:
```bash
python panDecay.py alignment.fas --analysis bayesian --bayesian-software mrbayes --bayes-model GTR --gamma
```

Example with all three analysis types:
```bash
python panDecay.py alignment.fas --analysis all --model GTR --gamma --bayesian-software mrbayes
```

Example with protein data:
```bash
python panDecay.py proteins.phy --format phylip --data-type protein --protein-model WAG --gamma
```

Example with site-specific analysis:
```bash
python panDecay.py alignment.fas --model GTR --gamma --site-analysis --visualize
```

Example with configuration file:
```bash
# Generate template configuration
python panDecay.py --generate-config my_config.ini
# Edit my_config.ini, then run:
python panDecay.py --config my_config.ini
```

Example testing specific branches:
```bash
# Test only specified clades
python panDecay.py alignment.fas --constraint-mode specific --test-branches "Homo,Pan,Gorilla;Mus,Rattus"
# Test all except specified clades
python panDecay.py alignment.fas --constraint-mode exclude --test-branches "outgroup1,outgroup2"
```

## Code Architecture

panDecay is organized around a main Python class `panDecayIndices` in panDecay.py, which handles:

1. **Initialization and Setup**: 
   - Parsing arguments and initializing parameters
   - Setting up temporary directories
   - Converting alignment to NEXUS format for PAUP*

2. **Tree Building and Analysis**:
   - Building the initial tree (parsimony or ML) with PAUP*
   - For each internal branch:
     - Creating a constraint to force non-monophyly
     - Searching for the best tree under this constraint
     - Calculating parsimony steps difference or likelihood differences

3. **Results Processing**:
   - Running the Approximately Unbiased (AU) test
   - Optional site-specific likelihood analysis
   - Generating annotated trees
   - Creating results files and visualizations

4. **PAUP* and MrBayes Interaction**:
   - Generating PAUP* command files for parsimony and ML analyses
   - Generating MrBayes command blocks for Bayesian analyses
   - Executing PAUP*/MrBayes as subprocesses
   - Parsing output logs and score files

The code uses subprocess calls to execute PAUP* and MrBayes, and depends on proper installation of these executables.

## Important Features

- Three analysis frameworks: parsimony (traditional Bremer support), maximum likelihood, and Bayesian
- Supports DNA, protein, and binary (0/1) discrete morphological data
- Can use various evolutionary models (GTR, HKY, JTT, WAG, Mk, etc.)
- Supports gamma-distributed rate heterogeneity (+G) and proportion of invariable sites (+I)
- Can analyze which specific sites support or conflict with each branch
- Generates annotated trees and comprehensive results files
- Provides optional visualizations of support values
- Can perform combined analyses with flexible combinations (e.g., ML+Parsimony, ML+Bayesian, Bayesian+Parsimony, or all three types)
- Configuration file support (INI format) for reproducible analyses
- Selective branch testing: test specific clades, exclude clades, or test all branches
- Constraint definitions can be provided via command line, config file, or external files
- Template configuration file generation with comprehensive documentation
- MCMC convergence checking for Bayesian analyses (ESS, PSRF, ASDSF metrics)

## Bayesian Analysis Considerations

### Bayes Decay vs Bayes Factor
- **Bayes Decay (BD)** is the primary metric: marginal log-likelihood difference (unconstrained - constrained)
- **Bayes Factor (BF)** = e^BD, but can become astronomically large
- BD is more interpretable: 0-1 (weak), 1-3 (positive), 3-5 (strong), >5 (very strong)
- BF values are capped at 10^6 for display to avoid numerical issues

### Model Dimension Effects
When comparing constrained vs unconstrained trees, be aware that:
- Unconstrained trees have more degrees of freedom (more possible topologies)
- This can bias Bayes factors toward unconstrained models
- Very large BD values (>10) may partially reflect this dimensional difference
- Interpret extreme values with caution

### Prior Sensitivity
Bayes factors are sensitive to:
- Prior distributions on tree topologies
- Prior distributions on branch lengths
- Model specification choices

### Convergence Diagnostics
panDecay automatically checks:
- **ESS (Effective Sample Size)**: Should be ≥200 for all parameters
- **PSRF (Potential Scale Reduction Factor)**: Should be ≤1.01
- **ASDSF (Average Standard Deviation of Split Frequencies)**: Should be <0.01

Use `--check-convergence` (default: on) and adjust thresholds with:
- `--min-ess 200`
- `--max-psrf 1.01`
- `--max-asdsf 0.01`
- `--convergence-strict` to fail if criteria not met