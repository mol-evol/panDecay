# Changelog

All notable changes to panDecay will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **Interactive tree visualization now works properly**
  - Fixed file path mismatch where HTML referenced `.nwk` files but only `.nwk.cleaned` files existed
  - Updated JavaScript libraries from phylotree@1.0.0-alpha.3 to phylotree@1.4.2 and D3.js to v6.7.0 for better stability
  - Improved error handling with detailed troubleshooting messages when tree files are missing
  - Fixed file cleanup logic to preserve cleaned tree files needed by HTML visualizations
  - Added better offline mode support with informative fallback messages
- **User interface improvements**
  - Consistent progress display format across ML, Bayesian, and parsimony analyses
  - Fixed confusing progress counters (e.g., "Clade 11 of 9") to show "Testing clade X (Y of Z)"
  - Simplified progress box styling with dashed borders for better terminal compatibility
  - Reduced verbose output while maintaining essential information
  - Changed default output filename from `ml_decay_indices.txt` to `pan_decay_indices.txt`
  - Updated all output to use relative paths instead of full absolute paths

### Changed
- Enhanced error messages for interactive tree visualization with specific troubleshooting steps
- Updated documentation with comprehensive troubleshooting guide for interactive visualizations
- Improved progress tracking logic to pre-calculate testable branches for accurate progress reporting

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