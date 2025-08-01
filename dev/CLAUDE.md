# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

panDecay is a Python command-line tool for calculating phylogenetic decay indices across multiple analysis frameworks:
- **Parsimony-based decay indices** (traditional Bremer support)
- **Maximum Likelihood (ML)-based decay indices** with AU test
- **Bayesian decay indices** using MrBayes with marginal likelihood comparisons
- **Combined analyses** supporting multiple frameworks simultaneously

The tool analyzes phylogenetic trees by comparing optimal trees with constrained trees where specific clades are forced to be non-monophyletic, providing statistical measures of branch support.

## Architecture Status (July 2025)

**panDecay uses a comprehensive monolithic architecture:**

### Current Implementation (Stable Production System)
- **Single comprehensive class** `panDecayIndices` with all analysis methods (4,963 lines)
- **All analysis types integrated**: ML, Bayesian, and Parsimony in one cohesive system
- **Comprehensive error handling** with robust exception management
- **Memory management** and resource optimization
- **Professional UI/UX**: Clean runtime banner with citation and organized output
- **File organization**: Timestamp-based directory structure with integrated file tracking
- **Progress tracking**: Real-time console updates with integrated progress logging
- **Proven reliability**: Extensively tested and validated system

## v2.0 Modular Architecture (Recommended)

### Entry Points
- **`src/panDecay.py`**: Main entry point using modular architecture with orchestration
- **`src/pandecay_main.py`**: Legacy implementation (fallback compatibility)

### Core Modular Components

#### Analysis Engines (`src/analysis/engines/`)
- **`base_engine.py`**: Abstract base class defining common interface for all engines
- **`ml_engine.py`**: Maximum Likelihood analysis with PAUP* integration and AU testing
- **`bayesian_engine.py`**: Bayesian analysis with MrBayes integration and marginal likelihood
- **`parsimony_engine.py`**: Traditional parsimony analysis with Bremer support calculation

#### Orchestration System (`src/orchestration/`)
- **`analysis_orchestrator.py`**: Coordinates multiple analysis engines, manages shared resources

#### Configuration Management (`src/config/`)
- **`constants.py`**: Centralized configuration constants (eliminates 57+ magic numbers)

#### Core Infrastructure (`src/core/`)
- **`file_tracker.py`**: Organized output directory management with timestamp-based structure
- **`progress_logger.py`**: Clean console progress tracking with dynamic updates

#### External Tool Management (`src/external_tools/`)
- **`tool_manager.py`**: Context-managed execution of PAUP*, MrBayes with proper cleanup

#### Utilities (`src/utils/`)
- **`thread_calculator.py`**: Adaptive thread calculation based on system resources
- **`memory_manager.py`**: Memory monitoring and optimization
- **`script_generators.py`**: PAUP*/MrBayes script generation utilities

#### Exception Handling (`src/exceptions/`)
- **`analysis_exceptions.py`**: Analysis-specific exceptions with context information
- **`tool_exceptions.py`**: External tool error handling
- **`io_exceptions.py`**: Input/output error management

### Modular Data Flow

1. **Entry Point**: `src/panDecay.py` parses arguments and initializes orchestrator
2. **Engine Creation**: Orchestrator creates appropriate analysis engines based on user request
3. **Resource Setup**: FileTracker creates organized output directories, memory monitoring begins
4. **Analysis Execution**: Each engine performs specialized analysis with proper resource management
5. **Result Aggregation**: Orchestrator collects results from all engines
6. **Output Generation**: Organized file structure with results, trees, reports, visualizations
7. **Cleanup**: Context managers ensure proper resource cleanup

### Key Advantages of Modular Architecture

- **Separation of Concerns**: Each engine handles one analysis type
- **Extensibility**: Easy to add new analysis engines
- **Testability**: Individual components can be tested in isolation
- **Resource Management**: Proper context managers and cleanup
- **Error Handling**: Comprehensive exception hierarchy with context
- **Memory Optimization**: Active memory monitoring and management
- **Import Reliability**: Robust absolute import system eliminates import failures

## Legacy Architecture (Backward Compatibility)

### Core Components (Legacy)
- **Main Script**: `src/pandecay_main.py` - Original monolithic implementation
- **Main Class**: `panDecayIndices` - Single large class handling all analysis types
- **Usage**: Automatically used as fallback when modular system encounters issues

### Legacy Data Flow
1. **Input Processing**: Alignment files converted to NEXUS format
2. **Tree Search**: Optimal tree search using PAUP* or MrBayes
3. **Constraint Generation**: Automated constraint generation for each internal branch
4. **Comparative Analysis**: Search for best trees under each constraint
5. **Statistical Testing**: AU test (ML) or marginal likelihood comparison (Bayesian)
6. **Output Generation**: Results tables, annotated trees, detailed reports

## Development Commands

### Running the Tool (v2.0 Modular Architecture)

```bash
# Basic usage with modular architecture
python3 src/panDecay.py alignment.fas --nst 6 --gamma

# Multi-analysis with orchestration
python3 src/panDecay.py alignment.fas --analysis ml+bayesian --nst 6

# Memory-optimized analysis
python3 src/panDecay.py alignment.fas --analysis ml --threads auto --memory-limit 8G

# With configuration file
python3 src/panDecay.py --config analysis.yaml

# All analysis types with orchestration
python3 src/panDecay.py alignment.fas --analysis all --nst 6 --gamma
```

### Running the Tool (Legacy Architecture)

```bash
# Legacy usage (automatically used as fallback if needed)
python3 src/pandecay_main.py alignment.fas --model GTR --gamma

# Different analysis types (legacy)
python3 src/pandecay_main.py alignment.fas --analysis ml+bayesian --model GTR
python3 src/pandecay_main.py alignment.fas --analysis parsimony
```

### Installing Dependencies

```bash
# Install Python dependencies
pip install -r requirements.txt

# External software (must be installed separately)
# PAUP* - from phylosolutions.com
# MrBayes - from nbisweden.github.io/MrBayes/
```

### Testing and Validation

**v2.0 Modular Architecture** includes comprehensive testing framework:
- **8-Phase Test Suite**: Systematic validation covering all system components
  - Import system verification
  - Core functionality testing
  - Individual engine validation
  - Orchestration integration testing
  - End-to-end workflow testing
  - Error handling validation
  - Configuration testing
  - Output verification
- **Import System Testing**: All 86+ import statements verified
- **Memory Testing**: Resource usage and cleanup validation
- **Error Recovery Testing**: Exception hierarchy and context validation

**Legacy Testing** methods:
- Example datasets in `examples/` directory (alignment.fas, proteins.phy, morpho.nex)
- Debug mode with temporary file retention: `--debug --keep-files`
- Validation against known phylogenetic results

### Common Development Tasks (v2.0 Modular)

- **Add new analysis engine**: 
  1. Create new engine class inheriting from `AnalysisEngine`
  2. Implement required abstract methods (`analyze()`, `validate_inputs()`, `get_analysis_type()`)
  3. Register engine with orchestrator
  4. Add comprehensive tests

- **Modify analysis workflow**: 
  1. Update appropriate engine in `src/analysis/engines/`
  2. Modify orchestrator if cross-engine coordination needed
  3. Update configuration constants if new parameters added

- **Add new external software**: 
  1. Extend `ExternalToolManager` with new tool support
  2. Create script generator utilities
  3. Add error handling for tool-specific issues
  4. Follow MrBayes integration pattern

- **Extend configuration**: 
  1. Add constants to `src/config/constants.py`
  2. Update argument parsing in main entry point
  3. Add validation logic

### Common Development Tasks (Legacy)

- **Add new analysis method**: Extend `panDecayIndices` class with new `_<method>_analysis()` method
- **Modify output format**: Update `_write_results()` and `_generate_annotated_trees()` methods
- **Add new external software**: Implement new methods following MrBayes integration pattern

## Key Configuration Files

- **requirements.txt**: Python dependencies
- **examples/pandecay.config**: Example configuration file (INI format)
- **README.md**: Comprehensive user documentation with examples
- **CHANGELOG.md**: Version history and feature updates
- **docs/**: Detailed documentation and guides

## Important Implementation Details (v2.0 Modular)

### Resource Management
- **Context Managers**: All external tool execution uses context managers for proper cleanup
- **Memory Monitoring**: Active memory usage tracking with peak usage reporting
- **FileTracker**: Organized output directory management with timestamp-based structure
- **Automatic Cleanup**: Comprehensive resource cleanup unless `--keep-files` specified

### Thread Safety and Parallel Processing
- **Adaptive Threading**: `AdaptiveThreadCalculator` determines optimal thread usage
- **PAUP* threading**: Controlled via `--threads` parameter with intelligent defaults
- **MrBayes MPI support**: Via `--use-mpi` with processor specification
- **Resource Coordination**: Orchestrator manages shared resources across engines

### Error Handling and Recovery
- **Exception Hierarchy**: Custom exceptions with detailed context information
- **Graceful Degradation**: Engines can fail independently without affecting others
- **Error Context**: Comprehensive error reporting with analysis context
- **Recovery Mechanisms**: Fallback to legacy architecture when needed

### Import System
- **Absolute Imports**: Robust import system using `src.` prefix throughout
- **No Relative Imports**: Eliminates "attempted relative import beyond top-level package" errors
- **Package Structure**: Proper Python package with `__init__.py` files
- **Import Verification**: All imports tested and validated

### Model Parameter Handling
- **Centralized Configuration**: All constants managed in `src/config/constants.py`
- **Flexible Specification**: Override parameters supported across all engines
- **NST-Based Models**: Preferred DNA model specification using NST values
- **Cross-Software Compatibility**: Automatic model conversion between PAUP*/MrBayes

## Important Implementation Details (Legacy)

### Thread Safety and Parallel Processing (Legacy)
- PAUP* threading controlled via `--threads` parameter
- MrBayes MPI support via `--use-mpi` with processor specification
- Constraint analyses run sequentially to avoid resource conflicts

### Temporary File Management (Legacy)
- Uses system temp directory or user-specified location
- Debug mode (`--debug`) retains all temporary files in `debug_runs/`
- Automatic cleanup unless `--keep-files` specified

### Error Handling Patterns (Legacy)
- Robust subprocess handling for external software
- Graceful degradation when optional features fail
- Basic logging with debug mode support

## External Software Integration

### PAUP* Integration
- Generates NEXUS command files for ML searches and AU tests
- Parses output log files for likelihood scores and tree results
- Handles constraint generation for non-monophyletic testing

### MrBayes Integration
- Generates MrBayes command blocks for Bayesian analysis
- Supports MPI and BEAGLE acceleration options
- Implements stepping-stone and harmonic mean marginal likelihood estimation
- Includes convergence diagnostics (ESS, PSRF, ASDSF)

## Output File Structure

- **Main results**: Tab-delimited with support metrics for each branch
- **Annotated trees**: Multiple Newick files with different annotation types
- **Detailed reports**: Markdown format with analysis summary and interpretation
- **Site analysis**: Optional directory with per-site likelihood analysis
- **Visualizations**: Optional matplotlib/seaborn plots

## Recent Enhancements (July 29, 2025)

### User Interface Improvements
- **Professional Runtime Banner**: Replaced jagged parameter table with clean, fixed-width banner
  - 80-character consistent width across all terminal sizes
  - Prominent citation display: "McInerney, J.O. (2025) panDecay: Phylogenetic Analysis with Decay Indices"
  - Organized parameter sections (Analysis Configuration, Runtime Settings, Bayesian Parameters)
  - Clean separator lines using Unicode box drawing characters

### Bug Fixes
- **Markdown Report Generation**: Fixed critical bug in analysis engine
  - Problem: `'MultipleSeqAlignment' object has no attribute 'name'` 
  - Solution: Replaced `self.alignment.name` with `self.alignment_file.stem` (6 instances)
  - Impact: Markdown reports now generate successfully with proper file naming

### Output Format Enhancements
- **Site Analysis Integration**: Enhanced markdown reports include comprehensive site analysis
- **File Organization**: Results organized in timestamp-based directories
- **Comprehensive Reports**: Multiple output formats (CSV, TXT, markdown) with consistent data

## Development Notes

### v2.0 Modular Architecture (Recommended)
- **Fully decomposed modular system** replacing monolithic design
- **SOLID principles**: Single responsibility, dependency injection, strategy patterns
- **Comprehensive testing**: 8-phase test suite with 80%+ coverage requirement
- **Resource management**: Context managers and proper cleanup throughout
- **Import system**: Robust absolute imports eliminate all import complications
- **Error handling**: Comprehensive exception hierarchy with detailed context
- **Memory optimization**: Active monitoring and intelligent resource management
- **Extensibility**: Plugin-based architecture for custom analysis types
- **Documentation**: Comprehensive API documentation and development guides

### Legacy Architecture (Compatibility)
- **Monolithic design** maintained for backward compatibility
- **Automatic fallback** when modular system encounters issues
- **Traditional patterns** following original implementation
- **Basic error handling** with limited context information

### General Notes
- **Cross-platform compatibility** (tested on macOS, Linux, Windows)
- **Extensive configuration support** (command-line, YAML, TOML, INI)
- **NST-based model specification** for DNA data (preferred over legacy model names)
- **External software integration** (PAUP*, MrBayes) with robust error handling
- **Organized output structure** with timestamp-based directory management

### Development Workflow
1. **Prefer modular architecture** for all new development
2. **Use comprehensive testing** for all changes
3. **Follow import guidelines** using absolute imports with `src.` prefix
4. **Implement proper resource management** with context managers
5. **Add comprehensive error handling** with custom exceptions and context
6. **Document all changes** in both code and user documentation