#!/usr/bin/env python3
"""
Test the fixed MrBayes execution with filtering.
"""

from src.external_tools.external_tool_runner import ExternalToolRunner
from pathlib import Path
import shutil

def test_fixed_mrbayes():
    """Test MrBayes execution with the fixed system."""
    
    # Create temporary directory
    temp_dir = Path("test_fixed_mb_temp")
    temp_dir.mkdir(exist_ok=True)
    
    try:
        # Copy test data
        shutil.copy("test_simple.nex", temp_dir)
        
        # Create runner
        runner = ExternalToolRunner(temp_dir, debug=True)
        
        # Test filtering
        filtered_file = runner._create_filtered_alignment_file(
            Path("test_simple.nex"), "test"
        )
        
        print(f"Created filtered file: {filtered_file}")
        
        # Show first few lines of filtered content
        content = filtered_file.read_text()
        lines = content.split('\n')
        print("\nFiltered content (first 15 lines):")
        for i, line in enumerate(lines[:15], 1):
            print(f"{i:2d}: {line}")
        
        # Test analysis
        model_settings = {'nst': 1, 'rates': 'equal'}
        mcmc_settings = {'ngen': 100, 'nchains': 1, 'samplefreq': 100, 'printfreq': 100}
        
        print("\nRunning MrBayes analysis...")
        result = runner.run_mrbayes_analysis(
            alignment_file=Path("test_simple.nex"),
            analysis_id="test",
            model_settings=model_settings,
            mcmc_settings=mcmc_settings,
            timeout=60
        )
        
        print(f"Success: {result['success']}")
        if 'marginal_likelihood' in result:
            print(f"Marginal likelihood: {result['marginal_likelihood']}")
        
        if not result['success']:
            print(f"Error: {result.get('error', 'Unknown error')}")
        
        # Check output files
        print("\nOutput files:")
        output_files = list(temp_dir.glob("*"))
        for f in output_files:
            if f.is_file():
                print(f"  {f.name}: {f.stat().st_size} bytes")
        
        # Look for MrBayes output patterns
        run_files = list(temp_dir.glob("*.run*.p")) + list(temp_dir.glob("*.run*.t"))
        con_files = list(temp_dir.glob("*.con.tre"))
        lstat_files = list(temp_dir.glob("*.lstat"))
        
        print(f"\nMrBayes output files:")
        print(f"  Run files: {len(run_files)}")
        print(f"  Consensus files: {len(con_files)}")
        print(f"  Lstat files: {len(lstat_files)}")
        
        return len(con_files) > 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        # Clean up
        shutil.rmtree(temp_dir, ignore_errors=True)

if __name__ == "__main__":
    success = test_fixed_mrbayes()
    print(f"\nTest {'PASSED' if success else 'FAILED'}")