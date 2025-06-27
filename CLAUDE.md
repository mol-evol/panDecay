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

Example with DNA data and GTR+G+I model:
```bash
python panDecay.py alignment.fas --model GTR --gamma --invariable --data-type dna
```

Example with protein data:
```bash
python panDecay.py proteins.phy --format phylip --data-type protein --protein-model WAG --gamma
```

Example with site-specific analysis:
```bash
python panDecay.py alignment.fas --model GTR --gamma --site-analysis --visualize
```

Example with Bayesian analysis:
```bash
python panDecay.py alignment.fas --analysis bayesian --bayesian-software mrbayes --bayes-model GTR --gamma
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
- Can perform combined analyses (e.g., ML+Bayesian, or all three types)