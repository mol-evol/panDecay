#!/usr/bin/env python3
"""
panDecay - Phylogenetic Decay Indices Calculator

A comprehensive implementation of phylogenetic decay indices (Bremer support) 
using Maximum Likelihood, Bayesian, and Parsimony approaches.

This is the main entry point for the panDecay analysis system.
"""

import sys
from pathlib import Path

# Add the src directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from src.main import main

if __name__ == "__main__":
    main()