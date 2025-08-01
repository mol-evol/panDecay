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

## Current Architecture (Monolithic System)

### Entry Points
- **`panDecay.py`**: Main entry point (root directory)
- **`src/main.py`**: Argument parsing and workflow coordination
- **`src/core/analysis_engine.py`**: Core analysis implementation (4,963 lines)

### Core Components

#### Core Analysis (`src/core/`)
- **`analysis_engine.py`**: Main `panDecayIndices` class with all analysis methods (4,963 lines)
- **`configuration.py`**: Configuration management and validation
- **`constants.py`**: Centralized configuration constants
- **`utils.py`**: Utility functions and progress tracking

### Monolithic Data Flow

1. **Entry Point**: `panDecay.py` imports and calls `src/main.py`
2. **Argument Processing**: `src/main.py` parses arguments and validates configuration
3. **Analysis Setup**: Creates `panDecayIndices` instance with all parameters
4. **Analysis Execution**: Single comprehensive class handles all analysis types
5. **Output Generation**: Organized file structure with results, trees, reports, visualizations
6. **Cleanup**: Comprehensive cleanup of temporary files

### Key Advantages of Monolithic Architecture

- **Simplicity**: Single comprehensive analysis class
- **Reliability**: Proven, extensively tested system
- **Integration**: All analysis types work together seamlessly
- **Resource Management**: Centralized resource management and cleanup
- **Error Handling**: Comprehensive exception handling throughout
- **Memory Management**: Efficient memory usage and monitoring

## Implementation Details

### Core Workflow
1. **Input Processing**: Alignment files converted to NEXUS format
2. **Tree Search**: Optimal tree search using PAUP* or MrBayes
3. **Constraint Generation**: Automated constraint generation for each internal branch
4. **Comparative Analysis**: Search for best trees under each constraint
5. **Statistical Testing**: AU test (ML) or marginal likelihood comparison (Bayesian)
6. **Output Generation**: Results tables, annotated trees, detailed reports

## Development Commands

### Running the Tool

```bash
# Basic usage
python3 panDecay.py alignment.fas --nst 6 --gamma

# Multi-analysis 
python3 panDecay.py alignment.fas --analysis ml+bayesian --nst 6

# With configuration file
python3 panDecay.py --config analysis.yaml

# All analysis types
python3 panDecay.py alignment.fas --analysis all --nst 6 --gamma

# Different data types
python3 panDecay.py proteins.phy --data-type protein --protein-model WAG
python3 panDecay.py morpho.nex --data-type discrete --discrete-model Mk
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

**Testing Framework** includes comprehensive validation:
- **Multi-Phase Test Suite**: Systematic validation covering core functionality
  - Import system verification
  - Core analysis engine testing
  - End-to-end workflow testing
  - Error handling validation
  - Configuration testing
  - Output verification
- **Memory Testing**: Resource usage and cleanup validation
- **Integration Testing**: External tool integration validation

**Testing Methods**:
- Example datasets in `examples/` directory (alignment.fas, proteins.phy, morpho.nex)
- Debug mode with temporary file retention: `--debug --keep-files`
- Validation against known phylogenetic results

### Common Development Tasks

- **Add new analysis method**: Extend `panDecayIndices` class with new analysis method
- **Modify output format**: Update result writing and tree generation methods
- **Add new external software**: Implement new integration following PAUP*/MrBayes pattern
- **Extend configuration**: 
  1. Add constants to `src/core/constants.py`
  2. Update argument parsing in `src/main.py`
  3. Add validation logic
- **Modify analysis workflow**: Update methods in `src/core/analysis_engine.py`

## Key Configuration Files

- **requirements.txt**: Python dependencies
- **examples/pandecay.config**: Example configuration file (INI format)
- **README.md**: Comprehensive user documentation with examples
- **CHANGELOG.md**: Version history and feature updates
- **docs/**: Detailed documentation and guides

## Important Implementation Details

### Resource Management
- **Memory Monitoring**: Active memory usage tracking with peak usage reporting
- **File Organization**: Organized output directory management with timestamp-based structure
- **Automatic Cleanup**: Comprehensive resource cleanup unless `--keep-files` specified

### Thread Safety and Parallel Processing
- **PAUP* threading**: Controlled via `--threads` parameter with intelligent defaults
- **MrBayes MPI support**: Via `--use-mpi` with processor specification
- **Constraint analyses**: Run sequentially to avoid resource conflicts

### Error Handling and Recovery
- **Exception Handling**: Custom exceptions with detailed context information
- **Graceful Degradation**: Robust error handling for external tool failures
- **Error Context**: Comprehensive error reporting with analysis context
- **Logging**: Enhanced logging with debug mode support

### Import System
- **Simple Imports**: Clear import structure: panDecay.py → src/main.py → src/core/analysis_engine.py
- **Package Structure**: Proper Python package with `__init__.py` files
- **Import Reliability**: All imports tested and validated

### Model Parameter Handling
- **Centralized Configuration**: All constants managed in `src/core/constants.py`
- **Flexible Specification**: Override parameters supported for all analysis types
- **NST-Based Models**: Preferred DNA model specification using NST values
- **Cross-Software Compatibility**: Automatic model conversion between PAUP*/MrBayes

### Temporary File Management
- Uses system temp directory or user-specified location
- Debug mode (`--debug`) retains all temporary files
- Automatic cleanup unless `--keep-files` specified

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

### Current Architecture (Monolithic System)
- **Comprehensive monolithic design** with all functionality integrated
- **Proven reliability**: Extensively tested and validated system
- **Resource management**: Proper cleanup and file management throughout
- **Import system**: Simple, reliable import structure
- **Error handling**: Enhanced exception handling with detailed context
- **Memory optimization**: Active monitoring and intelligent resource management
- **Professional presentation**: Clean UI with proper citations and organized output
- **Documentation**: Comprehensive user documentation and development guides

### General Notes
- **Cross-platform compatibility** (tested on macOS, Linux, Windows)
- **Extensive configuration support** (command-line, YAML, TOML, INI)
- **NST-based model specification** for DNA data (preferred over legacy model names)
- **External software integration** (PAUP*, MrBayes) with robust error handling
- **Organized output structure** with timestamp-based directory management

### Development Workflow
1. **Work within monolithic structure** maintaining integration of all components
2. **Use comprehensive testing** for all changes
3. **Follow simple import structure** with clear entry point flow
4. **Implement proper resource management** with cleanup and file handling
5. **Add comprehensive error handling** with custom exceptions and context
6. **Document all changes** in both code and user documentation