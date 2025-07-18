# Changelog

All notable changes to panDecay will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- **BD/site calculation accuracy**
  - Resolved issue where BD/site (per-site normalization) was incorrectly showing as 0.000000 despite valid BD values
  - Added comprehensive diagnostic logging to track BD/site calculation pipeline
  - Fixed calculation pipeline to ensure proper BD/site value storage and reporting
- **Effect Size Robust (MAD-based) scaling**
  - Added missing 1.4826 scaling factor to Median Absolute Deviation calculations
  - Fixed ES Robust values that were orders of magnitude higher than expected
  - MAD-based effect sizes now provide statistically equivalent robust alternatives to standard effect sizes
- **Parsimony analysis reporting**
  - Fixed variable name bug in parsimony report generation that caused crashes
  - Corrected reference from `bayes_decays` to `pars_decays` in parsimony-specific code sections

### Improved
- **Normalization method consistency**
  - All normalization methods now produce values on consistent, interpretable scales
  - Enhanced mathematical foundation alignment with established statistical frameworks
  - Improved error handling and fallback mechanisms for edge cases in effect size calculations

### Removed
- **Configuration migration functionality**
  - Removed `--migrate-config` command-line argument and associated implementation
  - Removed `migrate_ini_to_yaml()` function and migration utilities
  - Since panDecay has never been released, there are no existing users with legacy configurations to migrate
  - All three configuration formats (INI, YAML, TOML) remain fully supported
- **HTML tree visualization functionality**
  - Removed all HTML-based interactive tree visualization code (~600 lines)
  - HTML visualization was causing persistent browser compatibility and CORS issues
  - Users should use dedicated tree visualization software (e.g., FigTree) with the generated Newick files

### Fixed
- **MrBayes compatibility**
  - Fixed "Could not find command 'options'" error by filtering PAUP*-specific commands from NEXUS files
  - MrBayes now properly processes NEXUS files that contain PAUP* options blocks
- **Program hanging issues**
  - Fixed redundant NEXUS-to-NEXUS conversion that was causing >1000% CPU usage
  - Resolved multiple PAUP* processes spawning issue
- **Documentation accuracy**
  - Updated Bayesian analysis default from harmonic mean to stepping-stone sampling
  - Added missing command-line arguments (MPI, BEAGLE, visualization options)
  - Fixed duplicate default values in help text
  - Added annotated tree visualization example

### Changed
- **Bayesian decay interpretation guidance**
  - Removed arbitrary cutoff values (0-2, 2-5, 5-10, >10)
  - Now emphasizes comparative interpretation across branches
  - Recommends evaluating BD values alongside other support metrics
  - Notes that BD values may scale with alignment properties
- **Documentation improvements**
  - Added visual examples including annotated tree and site-specific analysis plots
  - Updated tree format documentation to match actual FigTree output
  - Clarified that strong support is best identified through concordance across multiple metrics

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