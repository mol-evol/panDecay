#!/usr/bin/env python3
"""
Test runner for panDecay.

This script provides different ways to run the test suite with various
options for coverage, performance testing, and selective test execution.
"""

import sys
import subprocess
import argparse
from pathlib import Path


def run_command(cmd, description=""):
    """Run a command and return the result."""
    if description:
        print(f"\n{'='*60}")
        print(f"Running: {description}")
        print(f"Command: {' '.join(cmd)}")
        print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        if result.stdout:
            print(result.stdout)
        
        if result.stderr:
            print("STDERR:", result.stderr, file=sys.stderr)
        
        return result.returncode
    except FileNotFoundError:
        print(f"Error: Command not found: {cmd[0]}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error running command: {e}", file=sys.stderr)
        return 1


def run_unit_tests(verbose=False, coverage=True):
    """Run unit tests."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend(["--cov=src", "--cov-report=html", "--cov-report=term-missing"])
    
    cmd.extend(["-m", "unit", "tests/"])
    
    return run_command(cmd, "Unit Tests")


def run_integration_tests(verbose=False):
    """Run integration tests."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    cmd.extend(["-m", "integration", "tests/"])
    
    return run_command(cmd, "Integration Tests")


def run_all_tests(verbose=False, coverage=True, parallel=False):
    """Run all tests."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    if coverage:
        cmd.extend([
            "--cov=src", 
            "--cov-report=html:htmlcov", 
            "--cov-report=term-missing",
            "--cov-report=xml",
            "--cov-fail-under=80"
        ])
    
    if parallel:
        cmd.extend(["-n", "auto"])
    
    cmd.append("tests/")
    
    return run_command(cmd, "All Tests")


def run_specific_tests(test_pattern, verbose=False):
    """Run specific tests matching a pattern."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    cmd.extend(["-k", test_pattern, "tests/"])
    
    return run_command(cmd, f"Tests matching: {test_pattern}")


def run_performance_tests(verbose=False):
    """Run performance tests."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    cmd.extend(["-m", "performance", "--benchmark-only", "tests/"])
    
    return run_command(cmd, "Performance Tests")


def run_memory_tests(verbose=False):
    """Run memory-related tests."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    cmd.extend(["-m", "memory", "tests/"])
    
    return run_command(cmd, "Memory Tests")


def run_external_tests(verbose=False):
    """Run tests requiring external tools."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    cmd.extend(["-m", "external", "tests/"])
    
    return run_command(cmd, "External Tool Tests")


def run_code_quality_checks():
    """Run code quality checks."""
    exit_code = 0
    
    # Run flake8
    print("\n" + "="*60)
    print("Running flake8 (code style)")
    print("="*60)
    result = run_command(["flake8", "src/", "tests/", "--max-line-length=100"])
    if result != 0:
        exit_code = result
    
    # Run black check
    print("\n" + "="*60)
    print("Running black (code formatting check)")
    print("="*60)
    result = run_command(["black", "--check", "--diff", "src/", "tests/"])
    if result != 0:
        exit_code = result
    
    # Run isort check
    print("\n" + "="*60)
    print("Running isort (import sorting check)")
    print("="*60)
    result = run_command(["isort", "--check-only", "--diff", "src/", "tests/"])
    if result != 0:
        exit_code = result
    
    # Run mypy
    print("\n" + "="*60)
    print("Running mypy (type checking)")
    print("="*60)
    result = run_command(["mypy", "src/"])
    if result != 0:
        exit_code = result
    
    return exit_code


def generate_coverage_report():
    """Generate detailed coverage report."""
    cmd = ["python", "-m", "pytest", 
           "--cov=src", 
           "--cov-report=html:htmlcov",
           "--cov-report=xml:coverage.xml",
           "--cov-report=term-missing",
           "tests/"]
    
    result = run_command(cmd, "Coverage Report Generation")
    
    if result == 0:
        print("\n" + "="*60)
        print("Coverage report generated!")
        print("HTML report: ./htmlcov/index.html")
        print("XML report: ./coverage.xml")
        print("="*60)
    
    return result


def check_dependencies():
    """Check if all test dependencies are installed."""
    required_packages = [
        "pytest",
        "pytest-cov", 
        "pytest-mock",
        "pytest-xdist",
        "coverage",
        "flake8",
        "black",
        "isort",
        "mypy"
    ]
    
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package.replace("-", "_"))
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print("Missing test dependencies:", file=sys.stderr)
        for package in missing_packages:
            print(f"  - {package}", file=sys.stderr)
        print("\nInstall with: pip install -r dev/requirements-test.txt", file=sys.stderr)
        return False
    
    return True


def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(description="Run panDecay tests")
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    parser.add_argument("--integration", action="store_true", help="Run integration tests only")
    parser.add_argument("--performance", action="store_true", help="Run performance tests only")
    parser.add_argument("--memory", action="store_true", help="Run memory tests only")
    parser.add_argument("--external", action="store_true", help="Run external tool tests only")
    parser.add_argument("--quality", action="store_true", help="Run code quality checks")
    parser.add_argument("--coverage", action="store_true", help="Generate coverage report")
    parser.add_argument("--pattern", "-k", help="Run tests matching pattern")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--no-coverage", action="store_true", help="Disable coverage reporting")
    parser.add_argument("--parallel", "-n", action="store_true", help="Run tests in parallel")
    parser.add_argument("--check-deps", action="store_true", help="Check test dependencies")
    
    args = parser.parse_args()
    
    # Check dependencies if requested
    if args.check_deps:
        if check_dependencies():
            print("All test dependencies are installed.")
            return 0
        else:
            return 1
    
    # Ensure we're in the right directory
    project_root = Path(__file__).parent
    if not (project_root / "src").exists():
        print("Error: Must run from project root directory", file=sys.stderr)
        return 1
    
    coverage = not args.no_coverage
    exit_code = 0
    
    try:
        if args.quality:
            exit_code = run_code_quality_checks()
        elif args.coverage:
            exit_code = generate_coverage_report()
        elif args.unit:
            exit_code = run_unit_tests(args.verbose, coverage)
        elif args.integration:
            exit_code = run_integration_tests(args.verbose)
        elif args.performance:
            exit_code = run_performance_tests(args.verbose)
        elif args.memory:
            exit_code = run_memory_tests(args.verbose)
        elif args.external:
            exit_code = run_external_tests(args.verbose)
        elif args.pattern:
            exit_code = run_specific_tests(args.pattern, args.verbose)
        else:
            # Run all tests by default
            exit_code = run_all_tests(args.verbose, coverage, args.parallel)
    
    except KeyboardInterrupt:
        print("\nTest run interrupted by user")
        return 130
    
    if exit_code == 0:
        print("\n" + "="*60)
        print("🎉 All tests passed!")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("❌ Some tests failed or code quality issues found")
        print("="*60)
    
    return exit_code


if __name__ == "__main__":
    sys.exit(main())