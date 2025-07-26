#!/usr/bin/env python3
"""
Test MrBayes and keep files for examination.
"""

from src.external_tools.external_tool_runner import ExternalToolRunner
from pathlib import Path
import shutil

def test_mrbayes_keep_files():
    """Test MrBayes and keep output files."""
    
    # Create temporary directory
    temp_dir = Path("test_mb_keep")
    temp_dir.mkdir(exist_ok=True)
    
    # Copy test data
    shutil.copy("test_simple.nex", temp_dir)
    
    # Create runner
    runner = ExternalToolRunner(temp_dir, debug=True)
    
    # Test analysis with stepping stone
    model_settings = {'nst': 1, 'rates': 'equal'}
    mcmc_settings = {
        'ngen': 200, 
        'nchains': 1, 
        'samplefreq': 100, 
        'printfreq': 100,
        'marginal_likelihood': 'ss',  # Enable stepping stone
        'use_ss': True,
        'ss_alpha': 0.3,
        'ss_nsteps': 5
    }
    
    print("Running MrBayes analysis with stepping stone...")
    result = runner.run_mrbayes_analysis(
        alignment_file=Path("test_simple.nex"),
        analysis_id="test_ss",
        model_settings=model_settings,
        mcmc_settings=mcmc_settings,
        timeout=120
    )
    
    print(f"Success: {result['success']}")
    if 'marginal_likelihood' in result:
        print(f"Marginal likelihood: {result['marginal_likelihood']}")
    
    # Check output files
    print("\nOutput files:")
    output_files = list(temp_dir.glob("*"))
    for f in sorted(output_files):
        if f.is_file():
            print(f"  {f.name}: {f.stat().st_size} bytes")
    
    print(f"\nFiles kept in: {temp_dir}")
    return result.get('success', False)

if __name__ == "__main__":
    success = test_mrbayes_keep_files()
    print(f"\nTest {'PASSED' if success else 'FAILED'}")