# panDecay API Reference (v1.1.0)

This document provides API documentation for the panDecay v1.1.0 system, featuring a comprehensive monolithic architecture. Updated August 1, 2025.

## Architecture Overview

The panDecay system uses a proven monolithic architecture:

```
src/
├── main.py                       # Main entry point and argument parsing
└── core/                        # Core functionality
    ├── analysis_engine.py       # Main panDecayIndices class (4,963 lines)
    ├── configuration.py         # Configuration management
    ├── constants.py            # System constants and defaults
    └── utils.py                # Utility functions and progress tracking
```

### Design Principles

- **Comprehensive Integration**: All analysis types (ML, Bayesian, Parsimony) in one cohesive system
- **Resource Management**: Proper cleanup and temporary file handling
- **Error Handling**: Robust exception handling throughout
- **Professional Presentation**: Clean UI with citations and organized output
- **External Tool Integration**: Seamless PAUP* and MrBayes integration

## Core Components

### panDecayIndices Class

The main analysis class in `src/core/analysis_engine.py` provides all functionality:

```python
from src.core.analysis_engine import panDecayIndices

# Initialize analysis
decay_calc = panDecayIndices(
    alignment_file="alignment.fas",
    alignment_format="fasta",
    model="GTR+G",
    data_type="dna"
)

# Run analysis
decay_calc.build_ml_tree()
decay_calc.calculate_decay_indices()
```

### Key Methods

#### Tree Building
- `build_ml_tree()`: Construct optimal ML tree using PAUP*
- `run_bootstrap_analysis()`: Perform bootstrap analysis
- `run_bayesian_analysis()`: Run MrBayes analysis

#### Decay Calculation
- `calculate_decay_indices()`: Calculate decay indices for all branches
- `calculate_au_test()`: Perform AU test for ML support
- `calculate_marginal_likelihood()`: Bayesian marginal likelihood comparison

#### Output Generation
- `write_formatted_results()`: Generate main results table
- `generate_detailed_report()`: Create comprehensive markdown report
- `annotate_trees()`: Create annotated tree files
- `export_results_csv()`: Export data in CSV format

### Configuration Management

Constants and configuration handled in `src/core/constants.py`:

```python
from src.core.constants import VERSION, DEFAULT_MODEL, DEFAULT_THREADS

# Version information
VERSION = "1.1.0"

# Default parameters
DEFAULT_MODEL = "GTR"
DEFAULT_THREADS = "auto"
DEFAULT_DATA_TYPE = "dna"
```

### Utility Functions

Progress tracking and utilities in `src/core/utils.py`:

```python
from src.core.utils import print_runtime_parameters, ProgressIndicator

# Display runtime parameters
print_runtime_parameters(args, model_display)

# Progress tracking
progress = ProgressIndicator(total_steps=10)
progress.update(1, "Building ML tree...")
```

## Usage Examples

### Basic Analysis

```python
from src.core.analysis_engine import panDecayIndices

# DNA analysis with GTR+G model
decay_calc = panDecayIndices(
    alignment_file="dna_alignment.fas",
    model="GTR+G",
    data_type="dna"
)

decay_calc.build_ml_tree()
decay_calc.calculate_decay_indices()
```

### Multi-Analysis

```python
# Combined ML and Bayesian analysis
decay_calc = panDecayIndices(
    alignment_file="alignment.fas",
    model="GTR+G",
    analysis_mode="ml+bayesian",
    bayes_ngen=1000000
)

decay_calc.build_ml_tree()
decay_calc.run_bayesian_analysis()
decay_calc.calculate_decay_indices()
```

### Protein Analysis

```python
# Protein data with WAG model
decay_calc = panDecayIndices(
    alignment_file="proteins.phy",
    alignment_format="phylip",
    data_type="protein",
    protein_model="WAG"
)

decay_calc.build_ml_tree()
decay_calc.calculate_decay_indices()
```

## Error Handling

The system includes comprehensive error handling:

```python
from src.core.analysis_engine import AnalysisEngineError, ExternalToolError

try:
    decay_calc = panDecayIndices("alignment.fas")
    decay_calc.build_ml_tree()
except ExternalToolError as e:
    print(f"External tool error: {e}")
except AnalysisEngineError as e:
    print(f"Analysis error: {e}")
```

## External Tool Integration

### PAUP* Integration
- Automatic NEXUS file generation
- ML tree search and AU testing
- Result parsing and validation

### MrBayes Integration
- Bayesian analysis with MCMC
- MPI support for parallel processing
- Convergence diagnostics and validation

## Output Structure

Results are organized in timestamped directories:

```
{basename}_pandecay_results/
├── {basename}_summary.txt          # Main results table
├── {basename}_report.md            # Detailed markdown report
├── {basename}_data.csv             # CSV export
├── trees/                          # Annotated tree files
├── site_analysis/                  # Site-specific analysis (optional)
├── supplementary/                  # Configuration and logs
└── visualizations/                 # Plots (optional)
```

## API Compatibility

This monolithic architecture provides:
- **Stability**: Proven, extensively tested system
- **Integration**: All analysis types work together seamlessly
- **Reliability**: Robust error handling and resource management
- **Performance**: Optimized for memory usage and processing speed
- **Professional Output**: Clean presentation with proper citations

For detailed usage examples and command-line interface documentation, see the main README.md file.