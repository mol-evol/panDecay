# panDecay Developer Guide

This guide provides comprehensive information for developers who want to understand, extend, or contribute to the panDecay codebase.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Project Structure](#project-structure)  
3. [Core Modules](#core-modules)
4. [Key Classes and Components](#key-classes-and-components)
5. [Data Flow and Processing Pipeline](#data-flow-and-processing-pipeline)
6. [Extension Points](#extension-points)
7. [Development Setup](#development-setup)
8. [Testing Framework](#testing-framework)
9. [Coding Standards](#coding-standards)
10. [Contributing Guidelines](#contributing-guidelines)

## Architecture Overview

panDecay uses a **modular architecture** designed for extensibility, maintainability, and performance. The application has evolved from a single-file script to a multi-module system that supports:

- **Multiple analysis frameworks** (ML, Bayesian, Parsimony)
- **Async processing** for performance optimization
- **Dual visualization systems** (static + interactive)
- **Flexible configuration** (YAML, TOML, INI)
- **Container deployment** with Docker

### Design Principles

1. **Separation of Concerns**: Each module handles a specific aspect (config, async processing, visualization)
2. **Backward Compatibility**: All existing workflows continue to work
3. **Extensibility**: Clear interfaces for adding new analysis methods and external tools
4. **Performance**: Async processing and optimized algorithms for large datasets
5. **User Experience**: Progressive enhancement with sensible defaults

## Project Structure

```
panDecay/
├── src/                          # Core application modules
│   ├── panDecay.py               # Main application logic (7,967 lines)
│   ├── config_loader.py          # Configuration parsing (593 lines)
│   ├── config_models.py          # Pydantic configuration models (414 lines)
│   ├── async_constraint_processor.py  # Async processing (395 lines)
│   └── dual_visualization.py     # Visualization systems (715 lines)
├── tests/                        # Test suite
│   ├── data/                     # Test datasets
│   ├── test_*.py                 # Unit and integration tests
│   └── test_fixtures.py          # Test utilities and mocks
├── examples/                     # Example data and configurations
│   ├── data/                     # Sample alignments and trees
│   └── configs/                  # Example configuration files
├── docs/                         # Documentation
│   ├── DEVELOPER_GUIDE.md        # This file
│   ├── API_REFERENCE.md          # API documentation
│   ├── CONFIGURATION_GUIDE.md    # Configuration guide
│   ├── DOCKER_GUIDE.md           # Docker deployment
│   ├── MIGRATION_GUIDE.md        # Feature migration guide
│   └── style_guide.md            # Code style guide
├── docker/                       # Container configuration
│   ├── Dockerfile                # Production container
│   ├── docker-compose.yml        # Multi-service deployment
│   └── docker-deploy.sh          # Deployment scripts
├── dev/                          # Development tools
│   ├── pytest.ini               # Test configuration
│   ├── requirements-dev.txt      # Development dependencies
│   └── run_tests.py              # Test runner
├── panDecay.py                   # CLI entry point
└── README.md                     # User documentation
```

## Core Modules

### 1. `src/panDecay.py` - Main Application Logic

**Primary Classes:**
- `panDecayIndices`: Core analysis engine
- `ExternalToolRunner`: Interface to PAUP* and MrBayes
- `TreeManager`: Tree manipulation and constraint generation
- `OutputManager`: Result formatting and file generation
- `DatasetNormalizer`: Cross-study normalization calculations

**Key Features:**
- Multi-framework phylogenetic analysis (ML, Bayesian, Parsimony)
- Constraint tree generation and analysis
- External software integration (PAUP*, MrBayes)
- Results aggregation and statistical testing

### 2. `src/config_loader.py` - Configuration Management

**Primary Functions:**
- `load_configuration()`: Universal config loader for YAML/TOML/INI
- `detect_config_format()`: Automatic format detection
- `validate_config_file()`: Configuration validation
- `ConfigurationError`: Custom exception handling

**Features:**
- Multi-format support (YAML, TOML, INI)
- Automatic format detection via file extensions and content
- Comprehensive error handling with user-friendly messages
- Pydantic integration for schema validation

### 3. `src/config_models.py` - Configuration Schema

**Primary Classes:**
- `PanDecayConfig`: Root configuration model
- `AnalysisConfig`: Analysis-specific settings
- `ComputationalConfig`: Performance and resource settings
- `VisualizationConfig`: Plotting and output configuration
- `OutputConfig`: File output settings

**Features:**
- Pydantic-based validation with detailed error messages
- Type safety and automatic conversion
- Comprehensive default values
- Documentation strings for all fields

### 4. `src/async_constraint_processor.py` - Async Processing

**Primary Classes:**
- `AsyncConstraintProcessor`: Main async orchestrator
- `ConstraintTask`: Individual constraint analysis task
- `ConstraintResult`: Analysis result container

**Features:**
- ThreadPoolExecutor-based parallelism for I/O-bound operations
- Configurable worker limits and timeouts
- Progress tracking and error handling
- Graceful degradation to sequential processing

### 5. `src/dual_visualization.py` - Visualization Systems

**Primary Classes:**
- `PlotManager`: Matplotlib-based visualization system
- `AlignmentVisualizer`: Site-specific alignment visualization
- `FileTracker`: Organized output directory management

**Features:**
- Static publication-ready plots (PNG format)
- Organized timestamp-based directory structure
- Site-specific alignment visualization
- Progress tracking with clean console output

## Key Classes and Components

### panDecayIndices Class

The core analysis engine that orchestrates all phylogenetic decay calculations.

**Key Methods:**
```python
def run_analysis()                    # Main entry point
def _ml_analysis()                    # Maximum likelihood analysis  
def _bayesian_analysis()              # Bayesian analysis with MrBayes
def _parsimony_analysis()             # Traditional Bremer support
def _generate_constraint_trees()      # Constraint generation
def _site_specific_analysis()         # Optional per-site analysis
def _effect_size_analysis()           # Normalization calculations
```

**Integration Points:**
- External tool execution via `ExternalToolRunner`
- Async processing via `AsyncConstraintProcessor`
- Tree operations via `TreeManager`
- Output generation via `OutputManager`

### ExternalToolRunner Class

Handles execution of external phylogenetic software with robust error handling.

**Supported Software:**
- **PAUP\***: ML searches, AU tests, parsimony analysis
- **MrBayes**: Bayesian analysis, marginal likelihood estimation

**Features:**
- Subprocess management with timeouts
- Output parsing and error detection
- Platform-specific executable detection
- Debug mode with file retention

### TreeManager Class

Manages phylogenetic tree operations and constraint generation.

**Key Functionality:**
- Tree parsing and format conversion
- Constraint generation for monophyly testing
- Tree annotation with support values
- Newick format validation and cleaning

## Data Flow and Processing Pipeline

### 1. Initialization Phase
```
Input Files → Format Detection → Validation → Temporary Directory Setup
```

### 2. Configuration Phase
```
CLI Args → Config File Loading → Schema Validation → Parameter Merge
```

### 3. Analysis Phase
```
Alignment Processing → Tree Search → Constraint Generation → Analysis Execution
```

### 4. Processing Phase (Async-enabled)
```
Constraint Tasks → Thread Pool → Parallel Execution → Result Aggregation
```

### 5. Statistical Testing Phase
```
Result Collection → AU Test / Marginal Likelihood → P-value Calculation
```

### 6. Output Phase
```
Result Formatting → File Generation → Visualization → Report Creation
```

### 7. Cleanup Phase
```
Temporary File Management → Debug File Retention → Resource Cleanup
```

## Extension Points

### Adding New Analysis Methods

1. **Extend panDecayIndices class:**
```python
def _new_method_analysis(self):
    """Implement new analysis method."""
    # Method-specific logic
    return results
```

2. **Update run_analysis() method:**
```python
if "new_method" in self.analysis_type:
    self._new_method_analysis()
```

3. **Add configuration support:**
```python
# In config_models.py
class AnalysisConfig(BaseModel):
    new_method_param: float = Field(default=1.0, description="New method parameter")
```

### Adding External Software Support

1. **Extend ExternalToolRunner:**
```python
def run_new_software(self, input_file, output_file):
    """Execute new phylogenetic software."""
    # Implementation
```

2. **Add software detection:**
```python
def find_new_software_executable(self):
    """Locate new software executable."""
    # Platform-specific detection
```

3. **Update configuration schema:**
```python
class ExternalToolsConfig(BaseModel):
    new_software_path: Optional[str] = None
    new_software_timeout: int = 3600
```

### Adding Visualization Types

1. **Extend DualVisualizationSystem:**
```python
def create_new_plot_type(self, data, config):
    """Create new visualization type."""
    # Static version (matplotlib)
    # Static matplotlib visualization
```

2. **Update configuration:**
```python
class VisualizationConfig(BaseModel):
    new_plot_enabled: bool = True
    new_plot_params: Dict[str, Any] = Field(default_factory=dict)
```

## Development Setup

### 1. Environment Setup
```bash
# Clone repository
git clone https://github.com/mol-evol/panDecay.git
cd panDecay

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # Linux/macOS
# or
venv\Scripts\activate     # Windows

# Install development dependencies
pip install -r dev/requirements-dev.txt
pip install -r requirements.txt
```

### 2. External Software (Optional)
```bash
# Install PAUP* (required for ML/parsimony)
# Download from: https://phylosolutions.com/paup-test/

# Install MrBayes (required for Bayesian)
# Download from: https://nbisweden.github.io/MrBayes/
```

### 3. Development Tools
```bash
# Run tests
python dev/run_tests.py

# Run with pytest directly
pytest tests/ -v

# Run specific test category
pytest tests/ -m "unit"

# Code formatting (if using black)
black src/ tests/

# Type checking (if using mypy)
mypy src/
```

## Testing Framework

### Test Structure
```
tests/
├── test_panDecay.py          # Core functionality tests
├── test_panDecay_simple.py   # Simple integration tests  
├── test_dataset_relative.py  # Normalization tests
├── test_effect_sizes.py      # Effect size calculation tests
├── test_fixtures.py          # Test utilities and mocks
└── data/                     # Test datasets
```

### Test Categories

**Unit Tests:**
- Individual method testing
- Mock external dependencies
- Configuration validation
- Utility function testing

**Integration Tests:**
- End-to-end workflow testing
- External software integration (when available)
- File I/O operations
- Error handling scenarios

**Mock Tests:**
- External software simulation
- Network-independent testing
- Deterministic results
- Fast execution

### Running Tests
```bash
# All tests
pytest tests/

# Unit tests only
pytest tests/ -m "unit"

# Integration tests only  
pytest tests/ -m "integration"

# With coverage
pytest tests/ --cov=src --cov-report=html

# Specific test file
pytest tests/test_panDecay.py -v
```

## Coding Standards

### Code Style
- **Line Length**: 100 characters (120 for tables/comments)
- **Indentation**: 4 spaces
- **Naming**: snake_case for functions/variables, PascalCase for classes
- **Documentation**: Comprehensive docstrings for all public methods

### Import Organization
```python
# Standard library imports
import os
import sys
from pathlib import Path

# Third-party imports
import numpy as np
import pandas as pd
from pydantic import BaseModel

# Local imports  
from .config_models import PanDecayConfig
from .async_constraint_processor import AsyncConstraintProcessor
```

### Error Handling
```python
# Use specific exception types
try:
    result = risky_operation()
except FileNotFoundError as e:
    logger.error(f"Input file not found: {e}")
    raise ConfigurationError(f"Cannot locate input file: {e}")
except subprocess.TimeoutExpired:
    logger.warning("Operation timed out, retrying...")
    # Implement retry logic
```

### Logging
```python
# Use module-level logger
logger = logging.getLogger(__name__)

# Appropriate log levels
logger.debug("Detailed debugging information")
logger.info("General information about progress")  
logger.warning("Something unexpected happened")
logger.error("Error occurred but execution continues")
logger.critical("Serious error, execution cannot continue")
```

## Contributing Guidelines

### 1. Development Workflow
1. **Fork** the repository
2. **Create** feature branch: `git checkout -b feature/new-analysis-method`
3. **Implement** changes with tests
4. **Test** thoroughly: `pytest tests/`
5. **Document** changes and new features
6. **Submit** pull request with clear description

### 2. Pull Request Requirements
- **All tests pass**: Unit and integration tests
- **Code coverage**: Maintain or improve coverage
- **Documentation**: Update relevant docs
- **Backward compatibility**: Existing workflows must continue working
- **Performance**: No significant performance regressions

### 3. Code Review Process
- **Automated checks**: CI/CD pipeline validation
- **Peer review**: Code quality and design review
- **Integration testing**: Full workflow validation
- **Documentation review**: Clarity and completeness

### 4. Release Process
1. **Version bump**: Update version in relevant files
2. **Changelog**: Document new features and fixes
3. **Testing**: Comprehensive test suite execution
4. **Documentation**: Update user and developer docs
5. **Container rebuild**: Update Docker images
6. **Tag release**: Git tag with semantic versioning

---

## Getting Help

- **Issues**: Report bugs and request features on GitHub Issues
- **Discussions**: Ask questions in GitHub Discussions
- **Documentation**: Check existing docs in the `docs/` directory
- **Examples**: Review examples in the `examples/` directory

For specific implementation questions, refer to the [API Reference](API_REFERENCE.md) or examine the test suite for usage examples.