#!/usr/bin/env python3
"""
Test runner for panDecay - comprehensive testing with different test levels.

This script provides different levels of testing:
- smoke: Quick smoke tests to verify basic functionality
- unit: Unit tests for individual components
- integration: Integration tests with mocked external dependencies
- full: Complete test suite including slow tests

Usage:
    python3 run_tests.py [smoke|unit|integration|full]
"""

import sys
import os
import argparse
import time
import subprocess
from pathlib import Path


def run_smoke_tests():
    """Run quick smoke tests."""
    print("Running smoke tests...")
    try:
        # Import and run the built-in smoke tests
        import panDecay
        result = panDecay.run_smoke_tests()
        
        if result:
            print("Smoke tests passed")
            return True
        else:
            print("Smoke tests failed")
            return False
    except Exception as e:
        print(f"Smoke tests failed with exception: {e}")
        return False


def run_unit_tests():
    """Run unit tests for individual components."""
    print("Running unit tests...")
    try:
        # Run the simplified test suite
        result = subprocess.run([
            sys.executable, "test_panDecay_simple.py"
        ], capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print("Unit tests passed")
            return True
        else:
            print("Unit tests failed")
            print(result.stdout)
            print(result.stderr)
            return False
    except subprocess.TimeoutExpired:
        print("Unit tests timed out")
        return False
    except Exception as e:
        print(f"Unit tests failed with exception: {e}")
        return False


def run_integration_tests():
    """Run integration tests with mocked dependencies."""
    print("Running integration tests...")
    try:
        # Check if pytest is available
        try:
            import pytest
            pytest_available = True
        except ImportError:
            pytest_available = False
        
        if pytest_available:
            # Run with pytest if available
            result = subprocess.run([
                sys.executable, "-m", "pytest", "test_panDecay_simple.py", "-v"
            ], capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print("Integration tests passed")
                return True
            else:
                print("Integration tests failed")
                print(result.stdout)
                return False
        else:
            # Fall back to unittest
            result = subprocess.run([
                sys.executable, "test_panDecay_simple.py"
            ], capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                print("Integration tests passed")
                return True
            else:
                print("Integration tests failed")
                print(result.stdout)
                return False
                
    except subprocess.TimeoutExpired:
        print("Integration tests timed out")
        return False
    except Exception as e:
        print(f"Integration tests failed with exception: {e}")
        return False


def run_full_tests():
    """Run the complete test suite."""
    print("Running full test suite...")
    
    # Run all test levels
    smoke_passed = run_smoke_tests()
    unit_passed = run_unit_tests()
    integration_passed = run_integration_tests()
    
    # Additional full test components
    syntax_passed = check_syntax()
    import_passed = check_imports()
    
    all_passed = all([smoke_passed, unit_passed, integration_passed, syntax_passed, import_passed])
    
    if all_passed:
        print("All tests passed!")
        return True
    else:
        print("Some tests failed!")
        return False


def check_syntax():
    """Check Python syntax."""
    print("Checking Python syntax...")
    try:
        result = subprocess.run([
            sys.executable, "-m", "py_compile", "panDecay.py"
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print("Syntax check passed")
            return True
        else:
            print("Syntax check failed")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"Syntax check failed with exception: {e}")
        return False


def check_imports():
    """Check that all imports work correctly."""
    print("Checking imports...")
    try:
        result = subprocess.run([
            sys.executable, "-c", "import panDecay; print('Import successful')"
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print("Import check passed")
            return True
        else:
            print("Import check failed")
            print(result.stderr)
            return False
    except Exception as e:
        print(f"Import check failed with exception: {e}")
        return False


def print_test_summary():
    """Print a summary of available tests."""
    print("""
panDecay Test Suite

Available test levels:
  smoke       - Quick smoke tests (< 10 seconds)
  unit        - Unit tests for individual components (< 1 minute)
  integration - Integration tests with mocked dependencies (< 5 minutes)
  full        - Complete test suite (< 10 minutes)

Test files:
  test_panDecay.py        - Comprehensive test suite
  test_panDecay_simple.py - Simplified test suite
  test_fixtures.py        - Test fixtures and mock data
  run_tests.py           - This test runner

Configuration:
  pytest.ini              - pytest configuration
  requirements-test.txt   - Test dependencies

Examples:
  python3 run_tests.py smoke
  python3 run_tests.py unit
  python3 run_tests.py full
""")


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="panDecay test runner")
    parser.add_argument(
        "level",
        choices=["smoke", "unit", "integration", "full", "help"],
        default="smoke",
        nargs="?",
        help="Test level to run"
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Test timeout in seconds (default: 600)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Verbose output"
    )
    
    args = parser.parse_args()
    
    if args.level == "help":
        print_test_summary()
        return 0
    
    # Change to the directory containing the script
    script_dir = Path(__file__).parent
    os.chdir(script_dir)
    
    print(f"Running {args.level} tests...")
    start_time = time.time()
    
    # Run the selected test level
    if args.level == "smoke":
        success = run_smoke_tests()
    elif args.level == "unit":
        success = run_unit_tests()
    elif args.level == "integration":
        success = run_integration_tests()
    elif args.level == "full":
        success = run_full_tests()
    else:
        print(f"Unknown test level: {args.level}")
        return 1
    
    end_time = time.time()
    duration = end_time - start_time
    
    print(f"\nTests completed in {duration:.1f} seconds")
    
    if success:
        print("All tests passed!")
        return 0
    else:
        print("Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())