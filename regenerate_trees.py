#!/usr/bin/env python3
"""
Quick script to regenerate annotated trees with Bayesian results included.
This reads the existing ML decay indices file and regenerates the annotated trees.
"""

import sys
import logging
from pathlib import Path

# Add the MLDecay directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

from MLDecay import MLDecayIndices

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Find the most recent debug run directory
    debug_runs_dir = Path("debug_runs")
    if not debug_runs_dir.exists():
        logger.error("No debug_runs directory found")
        return
    
    # Use the specific run directory from the analysis
    latest_run = debug_runs_dir / "mldecay_20250627_124554"
    if not latest_run.exists():
        logger.error(f"Run directory not found: {latest_run}")
        return
    
    logger.info(f"Using run directory: {latest_run}")
    
    # Create a minimal MLDecayIndices instance
    # Use the actual alignment file from the run
    alignment_file = latest_run / "alignment.nex"
    if not alignment_file.exists():
        # Try the original file
        alignment_file = Path("Primate10.nex")
        if not alignment_file.exists():
            logger.error("Could not find alignment file")
            return
    
    decay_calc = MLDecayIndices(
        alignment_file=str(alignment_file),
        alignment_format="nexus" if alignment_file.suffix == ".nex" else "fasta",
        paup_path="paup",
        debug=True,
        keep_files=True,
        analysis_mode="both"  # Important for including Bayesian results
    )
    
    # Set the temp path to the existing run directory
    decay_calc.temp_path = latest_run
    
    # Load the ML tree
    ml_tree_file = latest_run / "ml_tree.tre"
    if not ml_tree_file.exists():
        logger.error(f"ML tree file not found: {ml_tree_file}")
        return
    
    decay_calc.ml_tree_file = ml_tree_file
    
    # Read the decay indices from the results file
    results_file = Path("ml_decay_indices.txt")
    if not results_file.exists():
        logger.error("ml_decay_indices.txt not found")
        return
    
    # Parse the results file to reconstruct decay_indices
    decay_indices = {}
    with open(results_file, 'r') as f:
        lines = f.readlines()
        
    # Find the data section
    in_data = False
    for line in lines:
        line = line.strip()
        if line.startswith("Clade_ID"):
            in_data = True
            continue
        if in_data and line and not line.startswith("-"):
            parts = line.split('\t')
            if len(parts) >= 9:
                clade_id = parts[0]
                decay_indices[clade_id] = {
                    'num_taxa': int(parts[1]),
                    'constrained_lnl': float(parts[2]),
                    'lnl_diff': float(parts[3]),
                    'AU_pvalue': float(parts[4]),
                    'significant_au': parts[5] == 'Yes',
                    'bayes_decay': float(parts[6]) if parts[6] != 'N/A' else None,
                    'bayes_factor': float(parts[7]) if parts[7] != 'N/A' else None,
                    'taxa': parts[8].split(',')
                }
    
    decay_calc.decay_indices = decay_indices
    
    # Analysis mode already set to 'both' in initialization
    
    # Regenerate annotated trees
    output_dir = Path(".")
    tree_files = decay_calc.annotate_trees(output_dir, "annotated_tree")
    
    if tree_files:
        logger.info(f"Successfully regenerated {len(tree_files)} annotated trees with Bayesian results:")
        for tree_type, path in tree_files.items():
            logger.info(f"  - {tree_type}: {path}")
    else:
        logger.error("Failed to regenerate annotated trees")

if __name__ == "__main__":
    main()