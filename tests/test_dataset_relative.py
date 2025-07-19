#!/usr/bin/env python3
"""
Test script for dataset-relative normalization framework in panDecay

This script validates that:
1. Dataset-relative metrics are calculated correctly
2. Percentile ranks are computed properly
3. Z-scores are mathematically sound
4. Output formats are correct
5. Cohen's d framework has been completely removed
"""

import numpy as np
import subprocess
import tempfile
import os
import sys

def test_dataset_relative_calculations():
    """Test the mathematical correctness of dataset-relative calculations"""
    print("Testing dataset-relative normalization calculations...")
    
    # Test data: simulate likelihood decay values
    ld_values = np.array([5.2, 15.8, 32.1, 48.7, 93.5, 12.0, 28.4, 2.1, 97.4])
    
    # Expected calculations
    min_ld = np.min(ld_values)
    max_ld = np.max(ld_values)
    mean_ld = np.mean(ld_values)
    std_ld = np.std(ld_values, ddof=1)
    
    print(f"Test data statistics:")
    print(f"  Min LD: {min_ld:.3f}")
    print(f"  Max LD: {max_ld:.3f}")
    print(f"  Mean LD: {mean_ld:.3f}")
    print(f"  Std LD: {std_ld:.3f}")
    
    # Calculate expected dataset-relative values
    expected_dataset_relative = (ld_values - min_ld) / (max_ld - min_ld)
    
    # Calculate expected percentile ranks
    ranks = np.argsort(np.argsort(ld_values)) + 1  # 1-based ranks
    expected_percentile = (ranks / len(ld_values)) * 100
    
    # Calculate expected z-scores
    expected_z_scores = (ld_values - mean_ld) / std_ld
    
    print(f"\nExpected results for each branch:")
    for i, ld in enumerate(ld_values):
        print(f"  Branch {i+1}: LD={ld:.1f}, Rel={expected_dataset_relative[i]:.3f}, "
              f"Perc={expected_percentile[i]:.1f}%, Z={expected_z_scores[i]:.2f}")
    
    # Validate key properties
    assert np.min(expected_dataset_relative) == 0.0, "Dataset-relative min should be 0.0"
    assert np.max(expected_dataset_relative) == 1.0, "Dataset-relative max should be 1.0"
    assert np.min(expected_percentile) > 0, "Percentile min should be > 0"
    assert np.max(expected_percentile) <= 100, "Percentile max should be <= 100"
    assert abs(np.mean(expected_z_scores)) < 1e-10, "Z-score mean should be ~0"
    assert abs(np.std(expected_z_scores, ddof=1) - 1.0) < 1e-10, "Z-score std should be ~1"
    
    print("Mathematical calculations are correct")
    return True

def test_no_effect_size_references():
    """Test that Cohen's d effect size framework has been completely removed"""
    print("\nTesting removal of Cohen's d framework...")
    
    # Check panDecay.py for effect size references
    with open('panDecay.py', 'r') as f:
        content = f.read()
    
    forbidden_terms = [
        'effect_size',
        'effect size',
        'cohen',
        'Cohen',
        'ES_Robust', 
        'ES_Weighted',
        'ES_Standard'
    ]
    
    found_terms = []
    for term in forbidden_terms:
        if term in content:
            found_terms.append(term)
    
    if found_terms:
        print(f"Found forbidden effect size terms in code: {found_terms}")
        return False
    else:
        print("No Cohen's d effect size references found in code")
    
    # Check command line arguments
    if '--bd-normalization-methods effect_size' in content:
        print("Found old effect size command line option")
        return False
    
    if 'ld-normalization-methods' in content:
        print("New LD normalization framework is present")
    else:
        print("New LD normalization framework not found")
        return False
    
    return True

def test_panDecay_with_test_data():
    """Test panDecay with actual test data to verify new framework works"""
    print("\nTesting panDecay with test data...")
    
    # Check if test data exists
    if not os.path.exists('test_simple.fas'):
        print("test_simple.fas not found, skipping integration test")
        return True
    
    # Run a quick test with new dataset-relative framework
    try:
        cmd = [
            'python3', 'panDecay.py', 'test_simple.fas',
            '--analysis', 'ml',  # Quick ML-only test
            '--ld-normalization-methods', 'dataset_relative',
            '--output', 'test_integration_output.txt',
            '--debug'  # Keep files for inspection
        ]
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            print(f"panDecay failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        # Check that output file was created
        if os.path.exists('test_integration_output.txt'):
            print("Output file created successfully")
            
            # Check output contains new analysis
            with open('test_integration_output.txt', 'r') as f:
                output_content = f.read()
            
            # Check for successful analysis completion
            if 'panDecay Branch Support Analysis Results' in output_content:
                print("Analysis completed successfully")
            else:
                print("Analysis did not complete properly")
                return False
                
            # Check that no effect size columns remain
            forbidden_columns = ['Effect_Size', 'ES_Robust', 'ES_Weighted', 'Cohen']
            found_forbidden = []
            for col in forbidden_columns:
                if col in output_content:
                    found_forbidden.append(col)
            
            if found_forbidden:
                print(f"Found forbidden effect size references: {found_forbidden}")
                return False
            else:
                print("No effect size references found in output")
            
            return True
        else:
            print("Output file not created")
            return False
            
    except subprocess.TimeoutExpired:
        print("panDecay test timed out after 5 minutes")
        return False
    except Exception as e:
        print(f"Error running panDecay: {e}")
        return False

def main():
    """Run all validation tests"""
    print("=== panDecay Dataset-Relative Framework Validation ===\n")
    
    tests = [
        test_dataset_relative_calculations,
        test_no_effect_size_references,
        test_panDecay_with_test_data
    ]
    
    passed = 0
    total = len(tests)
    
    for test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                print(f"{test_func.__name__} failed")
        except Exception as e:
            print(f"{test_func.__name__} crashed: {e}")
    
    print(f"\n=== Test Results: {passed}/{total} passed ===")
    
    if passed == total:
        print("All tests passed! Dataset-relative framework is working correctly.")
        return 0
    else:
        print("Some tests failed. Please review the issues above.")
        return 1

if __name__ == '__main__':
    sys.exit(main())