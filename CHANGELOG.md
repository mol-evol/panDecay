# Changelog

All notable changes to panDecay will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1] - 2025-08-04

### Breaking Changes
- **Complete Model Specification Overhaul**: Removed all deprecated model flags for a clean, unambiguous system
- **REMOVED FLAGS** (no backward compatibility):
  - `--model` flag completely removed
  - `--protein-model` flag completely removed  
  - `--bayes-model` flag completely removed
- **NEW REQUIRED FLAGS**:
  - `--nst {1,2,6}` now REQUIRED for DNA data (1=JC-like, 2=HKY-like, 6=GTR-like)
  - `--protein {jtt,wag,lg,dayhoff,mtrev,cprev,blosum62,hivb,hivw}` now REQUIRED for protein data
  - Discrete data requires no model flags (uses Mk model automatically)

### Added
- **Foolproof Model Specification**: Zero ambiguity - exactly one way to specify each model component
- **Case-Insensitive Data Types**: Handles `dna`/`DNA`, `protein`/`Protein`, `discrete`/`Discrete` variants
- **Strict Validation**: Clear error messages prevent all invalid flag combinations
- **Clean Model Display**: Shows "GTR-type+G", "WAG+G", "Mk" instead of confusing model names

### Changed
- **Bayesian Analysis**: Now uses same model specification as ML analysis (no separate bayes-model needed)
- **Configuration Files**: Updated parameter mapping to use new flag system
- **Help Output**: Reorganized with clear "Model Specification (REQUIRED)" section
- **Error Messages**: More helpful validation errors for missing required flags

### Migration Guide
**Old → New Examples:**
- `--model GTR --gamma` → `--nst 6 --gamma`
- `--model HKY --gamma` → `--nst 2 --gamma`  
- `--model JC` → `--nst 1`
- `--protein-model WAG --gamma` → `--protein wag --gamma`
- `--data-type discrete --model Mk` → `--data-type discrete` (no model flag needed)

### Documentation
- **Updated README**: All examples converted to new flag system
- **Updated Config Template**: Reflects new model specification requirements
- **Updated Help Text**: Clear descriptions of new required flags

## [1.1.0] - 2025-08-04

### Added
- **Professional Python Package**: Complete conversion to installable pip package
- **PyPI Support**: Ready for publication on Python Package Index
- **Modern Packaging**: Uses pyproject.toml and PEP 518 standards
- **Entry Point**: Global `pandecay` command after pip installation

### Changed
- **Package Structure**: Renamed `src/` → `pandecay/` (standard Python naming)
- **CLI Module**: Renamed `main.py` → `cli.py` for clarity 
- **Import System**: Updated to absolute imports for better packaging
- **Installation Method**: Now supports `pip install pandecay`
- **Command Usage**: New command `pandecay` replaces previous script usage
- **Documentation**: Updated README with pip installation instructions

### Technical
- **Build System**: setuptools with pyproject.toml configuration
- **Dependencies**: Managed through requirements.txt and project metadata
- **Distribution**: Source distribution (sdist) and wheel formats
- **Manifest**: MANIFEST.in for proper file inclusion/exclusion
- **Entry Points**: Console script configuration for global command

### Migration
- **New Usage**: `pandecay alignment.fas --model GTR`
- **Installation**: `pip install pandecay` instead of manual clone and setup
- **Breaking Change**: Previous script-based usage no longer supported

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