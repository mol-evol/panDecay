#!/usr/bin/env python3
"""
Comprehensive end-to-end test of the complete panDecay pipeline.
"""

import sys
import subprocess
import time
from pathlib import Path

def test_complete_pipeline():
    """Test the complete panDecay pipeline with all functionality."""
    
    print("=== TESTING COMPLETE panDecay PIPELINE ===")
    print()
    
    # Test 1: Simple analysis with all modes
    print("TEST 1: Complete analysis with all modes (ml, bayesian, parsimony)")
    print("-" * 60)
    
    start_time = time.time()
    
    cmd = [
        "python3", "src/workflow_main.py",
        "test_simple.nex",
        "--analysis", "all",
        "--nst", "1",
        "--data-type", "dna",
        "--debug",
        "--keep-files",
        "--visualize",
        "--site-analysis"
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        
        execution_time = time.time() - start_time
        
        print(f"Execution time: {execution_time:.2f} seconds")
        print(f"Return code: {result.returncode}")
        print()
        
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
            print()
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
            print()
        
        # Check for success indicators
        success_indicators = [
            "ML analysis completed",
            "Bayesian analysis completed", 
            "Parsimony analysis completed",
            "Generated output files",
            "Analysis completed successfully"
        ]
        
        found_indicators = []
        for indicator in success_indicators:
            if indicator.lower() in result.stdout.lower():
                found_indicators.append(indicator)
        
        print(f"Success indicators found: {len(found_indicators)}/{len(success_indicators)}")
        for indicator in found_indicators:
            print(f"  ✓ {indicator}")
        
        missing_indicators = set(success_indicators) - set(found_indicators)
        for indicator in missing_indicators:
            print(f"  ✗ {indicator}")
        
        print()
        
        # Check output files
        print("CHECKING OUTPUT FILES:")
        print("-" * 30)
        
        # Look for output directories
        output_dirs = list(Path.cwd().glob("*panDecay*"))
        if output_dirs:
            output_dir = max(output_dirs, key=lambda p: p.stat().st_mtime)
            print(f"Latest output directory: {output_dir}")
            
            # Check for expected files
            expected_files = [
                "results/pan_decay_indices.txt",
                "reports/*.md",
                "trees/*.nwk", 
                "visualizations/*.png"
            ]
            
            for pattern in expected_files:
                files = list(output_dir.glob(pattern))
                print(f"  {pattern}: {len(files)} files")
                for f in files[:3]:  # Show first 3 files
                    print(f"    - {f.name}")
        else:
            print("  No output directories found")
        
        print()
        
        # Overall assessment
        if result.returncode == 0 and len(found_indicators) >= 3:
            print("🎉 TEST PASSED: Complete pipeline executed successfully!")
            return True
        else:
            print("❌ TEST FAILED: Pipeline did not complete successfully")
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ TEST FAILED: Pipeline timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"❌ TEST FAILED: Exception during execution: {e}")
        return False

def test_individual_components():
    """Test individual components to isolate any issues."""
    
    print()
    print("=== TESTING INDIVIDUAL COMPONENTS ===")
    print()
    
    # Test external tools
    print("Testing external tool availability:")
    try:
        # Test PAUP*
        paup_result = subprocess.run(["which", "paup"], capture_output=True, text=True)
        if paup_result.returncode == 0:
            print("  ✓ PAUP* found")
        else:
            print("  ✗ PAUP* not found")
        
        # Test MrBayes
        mb_result = subprocess.run(["which", "mb"], capture_output=True, text=True)
        if mb_result.returncode == 0:
            print("  ✓ MrBayes found")
        else:
            print("  ✗ MrBayes not found")
            
    except Exception as e:
        print(f"  Error checking external tools: {e}")
    
    print()
    
    # Test Python imports
    print("Testing Python imports:")
    test_imports = [
        "src.workflow.phylogenetic_workflow",
        "src.external_tools.external_tool_runner", 
        "src.analysis.ml_analysis",
        "src.analysis.bayesian_analysis",
        "src.visualization.plot_manager",
        "src.utils.output_manager"
    ]
    
    for module in test_imports:
        try:
            __import__(module)
            print(f"  ✓ {module}")
        except ImportError as e:
            print(f"  ✗ {module}: {e}")
    
    print()

if __name__ == "__main__":
    print("panDecay Pipeline Integration Test")
    print("=" * 50)
    print()
    
    # Test individual components first
    test_individual_components()
    
    # Run complete pipeline test
    success = test_complete_pipeline()
    
    print()
    if success:
        print("🎉 ALL TESTS PASSED! panDecay pipeline is working correctly.")
        sys.exit(0)
    else:
        print("❌ TESTS FAILED! Check the output above for issues.")
        sys.exit(1)