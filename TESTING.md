# Testing Guide for panDecay

## Test Suite Summary

This document provides an overview of the comprehensive test suite for panDecay.

### Test Statistics

- **Total Tests**: 111 tests (107 passing, 4 skipped when external tools unavailable)
- **Unit Tests**: 97 tests covering utilities, CLI, configuration, and core engine
- **Integration Tests**: 14 tests covering full ML workflows
- **Test Coverage**: Core utilities, argument parsing, configuration, model validation

### Quick Start

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run all tests
pytest

# Run with coverage report
pytest --cov=pandecay --cov-report=html

# Run only unit tests
pytest tests/unit/

# Run only integration tests
pytest tests/integration/
```

## Test Organization

```
tests/
├── unit/                       # Fast unit tests (no external dependencies)
│   ├── test_utils.py          # 35 tests for utility functions
│   ├── test_configuration.py  # 10 tests for config parsing
│   ├── test_cli.py            # 30 tests for CLI argument validation
│   └── test_analysis_engine.py # 22 tests for core engine initialization
├── integration/               # End-to-end workflow tests
│   └── test_ml_workflow.py    # 14 tests for ML analysis workflows
├── fixtures/                  # Test data
│   ├── alignments/           # Sample DNA, protein alignments
│   ├── trees/                # Sample tree files
│   └── configs/              # Sample configuration files
├── conftest.py               # Shared pytest fixtures
└── README.md                 # Detailed testing documentation
```

## What's Tested

### ✅ Fully Covered

1. **Utility Functions** (`test_utils.py`)
   - Path handling and display
   - Taxon name formatting for PAUP*
   - Support value formatting
   - Model display string building
   - File format auto-detection
   - Sequence type detection (DNA vs protein)
   - File path validation
   - String truncation

2. **Configuration Parsing** (`test_configuration.py`)
   - String to boolean conversion
   - Config file parsing with type conversion
   - Default value handling
   - Path parameter handling
   - Config override behavior

3. **CLI Argument Validation** (`test_cli.py`)
   - Data type normalization
   - Model flag validation (--nst, --protein)
   - DNA/protein/discrete model requirements
   - Parameter range validation (gamma_shape, prop_invar)
   - Argument parser setup
   - Output structure creation

4. **Analysis Engine Initialization** (`test_analysis_engine.py`)
   - Instance creation with various parameters
   - Thread count validation
   - Data type normalization
   - NEXUS file conversion
   - Model string building
   - Analysis mode parsing

5. **Integration Workflows** (`test_ml_workflow.py`)
   - ML workflow initialization
   - Gamma rate variation
   - Protein data handling
   - NEXUS input handling
   - Error handling
   - *Full ML tree building (requires PAUP*)*

### ⏭️ Skipped Tests (Require External Tools)

- `test_build_ml_tree`: Requires PAUP* in PATH
- `test_ml_tree_with_gamma`: Requires PAUP* in PATH
- `test_ml_tree_different_nst`: Requires PAUP* in PATH
- `test_generate_template`: Requires config template file

These tests are automatically skipped when dependencies aren't available.

## Continuous Integration

GitHub Actions workflow (`.github/workflows/tests.yml`) runs automatically on:
- Pushes to `main`, `develop`, or `claude/*` branches
- Pull requests to `main` or `develop`

The CI pipeline:
1. Tests on Python 3.8, 3.9, 3.10, 3.11, 3.12
2. Runs unit tests with coverage reporting
3. Runs integration tests (skips tests requiring PAUP*/MrBayes)
4. Checks code formatting with black
5. Lints code with flake8
6. Uploads coverage reports to Codecov

## Test Development

### Adding New Tests

When adding new functionality:

1. **Write tests first** (TDD) or alongside your code
2. **Use appropriate test type**:
   - Unit test for isolated functions
   - Integration test for workflows
3. **Add fixtures** if you need sample data
4. **Mark appropriately**: Use `@pytest.mark.requires_paup` for tests needing external tools

Example:
```python
import pytest

def test_my_new_feature():
    """Test description."""
    result = my_function("input")
    assert result == "expected"

@pytest.mark.requires_paup
def test_with_paup():
    """Test requiring PAUP*."""
    # Test code here
```

### Running Specific Tests

```bash
# Run specific test file
pytest tests/unit/test_utils.py

# Run specific test class
pytest tests/unit/test_utils.py::TestFormatTaxonForPaup

# Run specific test function
pytest tests/unit/test_utils.py::TestFormatTaxonForPaup::test_simple_name

# Run tests matching pattern
pytest -k "format_taxon"

# Run with verbose output
pytest -v

# Run with extra verbosity
pytest -vv

# Show local variables on failure
pytest -l
```

## Test Fixtures

Shared fixtures are defined in `tests/conftest.py`:

- `fixtures_dir`: Path to fixtures directory
- `dna_alignment_small`: Small DNA alignment (5 taxa, 32 sites)
- `dna_alignment_nexus`: NEXUS DNA alignment (6 taxa, 60 sites)
- `protein_alignment_small`: Small protein alignment (5 taxa, 20 sites)
- `simple_tree`: Simple Newick tree
- `test_config`: Test configuration file
- `temp_workspace`: Temporary directory for test outputs
- `temp_output_dir`: Temporary output directory

## Coverage Goals

Current coverage by module:

- **pandecay.core.utils**: ~85% (high coverage of utility functions)
- **pandecay.core.configuration**: ~80% (config parsing well-tested)
- **pandecay.cli**: ~75% (CLI validation thoroughly tested)
- **pandecay.core.analysis_engine**: ~30% (initialization tested, analysis workflows require external tools)

Target: 80% overall coverage for testable code (excluding PAUP*/MrBayes integration).

## Future Test Additions

Recommended additions:

1. **Bayesian workflow tests** (requires MrBayes)
2. **Parsimony workflow tests** (requires PAUP*)
3. **Bootstrap analysis tests** (requires PAUP*)
4. **AU test workflow** (requires PAUP*)
5. **Site analysis tests**
6. **Tree annotation tests**
7. **Output file generation tests**
8. **Constraint testing workflow**

## Troubleshooting

### Common Issues

**Tests not found:**
```bash
# Ensure you're in the repo root
cd /path/to/panDecay
pytest
```

**Import errors:**
```bash
# Install package in development mode
pip install -e .

# Or install just dependencies
pip install numpy biopython pytest
```

**Fixtures not found:**
```bash
# Check fixtures exist
ls tests/fixtures/alignments/
```

**CI failing:**
- Check that all imports are available without external dependencies
- Ensure tests marked with `@pytest.mark.requires_paup` for PAUP*-dependent tests
- Verify tests work in a fresh environment

## Questions?

- See `tests/README.md` for detailed testing documentation
- Check test files for examples of testing patterns
- Open an issue on GitHub for testing-related questions
