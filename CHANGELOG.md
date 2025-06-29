# Changelog

All notable changes to panDecay will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-06-29

### Added
- Configuration file support (INI format) for reproducible analyses
  - `--config` option to read parameters from configuration file
  - `--generate-config` option to create a fully-commented template
  - Command-line arguments override configuration file values
- Selective branch testing functionality
  - `--constraint-mode` option with three modes: all, specific, exclude
  - `--test-branches` option to specify branches via taxa lists, branch IDs, or file
  - `--constraint-file` option for external constraint definitions
  - Support for constraints in configuration file `[constraints]` section
- Enhanced documentation with new examples for configuration files and constraint usage

### Changed
- Updated version to 1.1.0
- Modified branch processing logic in ML, Bayesian, and parsimony analyses to respect constraint settings
- Made alignment argument optional when using configuration file

### Fixed
- Improved error handling for missing constraints in "specific" mode

## [1.0.3] - Previous release

### Features
- ML-based phylogenetic decay indices calculation
- Bayesian decay indices with MrBayes support
- Traditional parsimony Bremer support calculation
- Site-specific likelihood analysis
- Bootstrap support integration
- Visualization capabilities
- MPI and BEAGLE support for Bayesian analyses