# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MLDecay is a Python command-line tool for calculating Maximum Likelihood (ML)-based phylogenetic decay indices, also known as ML-Bremer support or ML branch support. It leverages the phylogenetic software PAUP* to assess the robustness of clades in a phylogenetic tree.

## Required Dependencies

- Python 3.8+
- PAUP* command-line executable (must be installed separately)
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

### Running MLDecay

Basic command:
```bash
python MLDecay.py <alignment_file> --model <model_name> [options...]
```

Example with DNA data and GTR+G+I model:
```bash
python MLDecay.py alignment.fas --model GTR --gamma --invariable --data-type dna
```

Example with protein data:
```bash
python MLDecay.py proteins.phy --format phylip --data-type protein --protein-model WAG --gamma
```

Example with site-specific analysis:
```bash
python MLDecay.py alignment.fas --model GTR --gamma --site-analysis --visualize
```

## Code Architecture

MLDecay is organized around a main Python class `MLDecayIndices` in MLDecay.py, which handles:

1. **Initialization and Setup**: 
   - Parsing arguments and initializing parameters
   - Setting up temporary directories
   - Converting alignment to NEXUS format for PAUP*

2. **Tree Building and Analysis**:
   - Building the initial ML tree with PAUP*
   - For each internal branch:
     - Creating a constraint to force non-monophyly
     - Searching for the best tree under this constraint
     - Calculating likelihood differences

3. **Results Processing**:
   - Running the Approximately Unbiased (AU) test
   - Optional site-specific likelihood analysis
   - Generating annotated trees
   - Creating results files and visualizations

4. **PAUP* Interaction**:
   - Generating PAUP* command files
   - Executing PAUP* as a subprocess
   - Parsing output logs and score files

The code uses subprocess calls to execute PAUP* and depends on proper installation of the PAUP* executable.

## Important Features

- Supports DNA, protein, and binary (0/1) discrete morphological data
- Can use various evolutionary models (GTR, HKY, JTT, WAG, Mk, etc.)
- Supports gamma-distributed rate heterogeneity (+G) and proportion of invariable sites (+I)
- Can analyze which specific sites support or conflict with each branch
- Generates annotated trees and comprehensive results files
- Provides optional visualizations of support values