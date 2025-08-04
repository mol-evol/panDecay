# Changelog

All notable changes to panDecay will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-08-04

### Added
- **Professional Python Package**: Complete conversion to installable pip package
- **PyPI Support**: Ready for publication on Python Package Index
- **Modern Packaging**: Uses pyproject.toml and PEP 518 standards
- **Entry Point**: Global `pandecay` command after pip installation
- **Backward Compatibility**: Deprecation wrapper for existing `panDecay.py` users

### Changed
- **Package Structure**: Renamed `src/` → `pandecay/` (standard Python naming)
- **CLI Module**: Renamed `main.py` → `cli.py` for clarity 
- **Import System**: Updated to absolute imports for better packaging
- **Installation Method**: Now supports `pip install pandecay`
- **Command Usage**: New command `pandecay` instead of `python3 panDecay.py`
- **Documentation**: Updated README with pip installation instructions

### Technical
- **Build System**: setuptools with pyproject.toml configuration
- **Dependencies**: Managed through requirements.txt and project metadata
- **Distribution**: Source distribution (sdist) and wheel formats
- **Manifest**: MANIFEST.in for proper file inclusion/exclusion
- **Entry Points**: Console script configuration for global command

### Migration
- **Old Usage**: `python3 panDecay.py alignment.fas --model GTR` (still works with warning)
- **New Usage**: `pandecay alignment.fas --model GTR` (recommended)
- **Installation**: `pip install pandecay` instead of manual clone and setup
- **No Breaking Changes**: All CLI arguments and functionality identical

### Package Info
- **Name**: `pandecay` (lowercase, PyPI standard)
- **Version**: 1.1.0
- **Python**: Requires Python ≥3.8
- **License**: MIT
- **Author**: James McInerney

## [1.0.0] - Previous Releases

### Features
- ML-based phylogenetic decay indices
- Bayesian decay analysis with MrBayes integration  
- Traditional parsimony Bremer support
- Approximately Unbiased (AU) test implementation
- Site-specific likelihood analysis
- Bootstrap support analysis
- Multiple analysis modes (ML, Bayesian, Parsimony, combined)
- Comprehensive output formats (text, trees, markdown reports)
- Visualization support (matplotlib/seaborn)
- Configuration file support
- Advanced MrBayes features (MPI, BEAGLE, convergence checking)
- DNA, protein, and discrete morphological data support