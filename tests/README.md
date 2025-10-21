# panDecay Test Suite

This directory contains comprehensive tests for the panDecay phylogenetic analysis package.

## Test Structure

```
tests/
├── unit/                       # Unit tests for individual components
│   ├── test_utils.py          # Utility function tests
│   ├── test_configuration.py  # Configuration file parsing tests
│   ├── test_cli.py            # CLI argument parsing tests
│   └── test_analysis_engine.py # Core engine tests (no external tools)
├── integration/               # Integration tests for full workflows
│   └── test_ml_workflow.py    # ML analysis workflow tests
├── fixtures/                  # Test data and expected outputs
│   ├── alignments/           # Sample alignment files
│   ├── trees/                # Sample tree files
│   ├── configs/              # Sample configuration files
│   └── expected_outputs/     # Expected results for validation
└── conftest.py               # Pytest configuration and fixtures

## Running Tests

### Run all tests
```bash
pytest
```

### Run only unit tests
```bash
pytest tests/unit/
```

### Run only integration tests
```bash
pytest tests/integration/
```

### Run tests with verbose output
```bash
pytest -v
```

### Run tests with coverage report
```bash
pytest --cov=pandecay --cov-report=html
```

### Run specific test file
```bash
pytest tests/unit/test_utils.py
```

### Run specific test class or function
```bash
pytest tests/unit/test_utils.py::TestFormatTaxonForPaup
pytest tests/unit/test_utils.py::TestFormatTaxonForPaup::test_simple_name
```

## Test Categories

Tests are organized using pytest markers:

- **`unit`**: Unit tests for individual components
- **`integration`**: Integration tests for full workflows
- **`requires_paup`**: Tests that require PAUP* (skipped if not available)
- **`requires_mrbayes`**: Tests that require MrBayes (skipped if not available)
- **`slow`**: Tests that take significant time to run

### Running tests by marker

```bash
# Run only unit tests
pytest -m unit

# Run only integration tests
pytest -m integration

# Skip slow tests
pytest -m "not slow"

# Run only tests that require PAUP*
pytest -m requires_paup

# Run all tests except those requiring external tools
pytest -m "not requires_paup and not requires_mrbayes"
```

## External Tool Requirements

Some integration tests require external phylogenetic tools:

- **PAUP***: Required for ML and parsimony analysis tests
  - Tests requiring PAUP* will be automatically skipped if not available
  - To enable these tests, ensure PAUP* is in your PATH

- **MrBayes**: Required for Bayesian analysis tests
  - Tests requiring MrBayes will be automatically skipped if not available
  - To enable these tests, ensure MrBayes is in your PATH

## Test Fixtures

The `fixtures/` directory contains sample data files for testing:

### Alignments
- `dna_small.fas`: Small DNA alignment (5 taxa, 32 sites)
- `protein_small.fas`: Small protein alignment (5 taxa, 20 sites)
- `dna_medium.nex`: Medium DNA alignment in NEXUS format (6 taxa, 60 sites)
- `dna_phylip.phy`: DNA alignment in PHYLIP format

### Trees
- `simple_tree.nwk`: Simple Newick tree for testing tree operations

### Configurations
- `test_config.txt`: Sample configuration file

## Writing New Tests

### Unit Test Template

```python
import pytest
from pandecay.core.utils import my_function

class TestMyFunction:
    """Test my_function behavior."""

    def test_basic_case(self):
        """Test basic functionality."""
        result = my_function("input")
        assert result == "expected_output"

    def test_edge_case(self):
        """Test edge case."""
        with pytest.raises(ValueError):
            my_function(None)
```

### Integration Test Template

```python
import pytest
from pathlib import Path
from pandecay.core.analysis_engine import panDecayIndices

@pytest.mark.integration
@pytest.mark.requires_paup
def test_ml_workflow(dna_alignment_small, temp_workspace):
    """Test complete ML workflow."""
    instance = panDecayIndices(
        alignment_file=dna_alignment_small,
        alignment_format="fasta",
        temp_dir=temp_workspace,
        data_type="dna",
        nst=6
    )

    instance.build_ml_tree()

    assert instance.ml_tree is not None
    assert instance.ml_likelihood is not None
```

## Test Coverage

To generate a coverage report:

```bash
pytest --cov=pandecay --cov-report=html --cov-report=term
```

This will create an HTML coverage report in `htmlcov/index.html`.

### Current Coverage Goals

- **Unit tests**: Aim for >80% coverage of utility functions, configuration parsing, and CLI
- **Integration tests**: Cover main workflows (ML, Bayesian, parsimony)
- **Critical paths**: 100% coverage of model validation and data type detection

## Continuous Integration

Tests are designed to work in CI/CD environments:

- Tests requiring external tools (PAUP*, MrBayes) are automatically skipped if tools aren't available
- Fast unit tests run on every commit
- Slower integration tests can be run on a schedule or before releases

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.11
      - name: Install dependencies
        run: |
          pip install -e ".[dev]"
      - name: Run unit tests
        run: pytest tests/unit/ -v
      - name: Run integration tests (without external tools)
        run: pytest tests/integration/ -v -m "not requires_paup and not requires_mrbayes"
```

## Troubleshooting

### Tests fail with "fixture not found"
- Ensure you're running pytest from the repository root directory
- Check that test fixtures exist in `tests/fixtures/`

### Import errors
- Install the package in development mode: `pip install -e .`
- Ensure all dependencies are installed: `pip install -e ".[dev]"`

### Tests requiring PAUP*/MrBayes are not running
- These tests are skipped by default if tools aren't available
- Install PAUP* or MrBayes and ensure they're in your PATH
- Use `-v` flag to see skip reasons: `pytest -v`

## Contributing Tests

When contributing new code:

1. **Write tests first** (TDD approach) or alongside your code
2. **Ensure all tests pass** before submitting a PR
3. **Add fixtures** if your tests need new sample data
4. **Document test purpose** with clear docstrings
5. **Use appropriate markers** for test categorization
6. **Check coverage** to ensure new code is tested

## Questions?

For questions about tests, please open an issue on GitHub or consult the main README.md.
