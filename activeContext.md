# Active Context - panDecay Project

## Current Status
**Date**: 2025-07-29  
**Branch**: feature/alignment-visualization  
**Status**: Production refinements completed - professional UI enhancements and comprehensive documentation updates

## Project Overview
panDecay is a phylogenetic analysis tool for studying gene family evolution and decay processes. The project has recently undergone comprehensive architectural improvements to create a cleaner, more intuitive interface.

## Recent Major Changes (July 2025)

### 1. NST-First Migration (Phase 1)
- Migrated from traditional model names (GTR, HKY, JC) to NST-based parameters
- NST (Number of Substitution Types): 1=equal rates, 2=ti/tv rates, 6=all rates differ
- Provides clearer control over DNA substitution complexity

### 2. Case-Insensitive Argument Parsing
- Fixed issue where `--data-type DNA` failed due to case sensitivity
- Implemented custom type functions for case-insensitive handling
- Supports uppercase, lowercase, and mixed case variations

### 3. Interface Simplification
- **Removed redundant flags**: Eliminated `--gamma` and `--invariable` flags
- **Enhanced --rates flag**: Now supports all 4 rate model combinations:
  - `equal`: Equal rates across sites
  - `gamma`: Gamma-distributed rate variation
  - `propinv`: Proportion of invariable sites
  - `invgamma`: Both gamma distribution and invariable sites

### 4. Comprehensive Architecture Cleanup
- **File Organization**: 
  - Renamed root `panDecay.py` → `pandecay` (entry point)
  - Renamed `src/panDecay.py` → `src/pandecay_main.py` (main implementation)
  - Eliminated confusing dual naming
- **Clean Model Specification**:
  - DNA: `--nst` parameter (default: 2)
  - Proteins: `--protein-model` (default: WAG)
  - Discrete/Morphological: `--discrete-model` (default: Mk)
- **Removed Complex Fallbacks**: Eliminated confusing model-to-NST inference logic
- **Consolidated Argument Parsing**: All parsing now in single location

### 5. Help Text Improvements
- Removed misleading software-specific references (e.g., "PAUP*'s default")
- Used clear biological language instead of technical jargon
- Eliminated confusing "+G", "+I" notation
- Made all parameter descriptions consistent

## Major Architectural Refactoring (July 24-26, 2025)

### 6. Complete System Decomposition - v2.0 Architecture
Following the interface improvements, a comprehensive architectural overhaul was undertaken to address fundamental design issues identified in the monolithic codebase.

#### Key Problems Addressed
- **Monolithic Design**: 6,293-line `panDecayIndices` class with 111 methods
- **Code Duplication**: Extensive duplication across analysis types
- **Poor Resource Management**: No context managers or proper cleanup
- **Import Complexities**: Relative import failures preventing modular usage
- **Error Handling**: Inadequate exception hierarchy and context

#### Architectural Transformation
- **Modular Analysis Engines**: Decomposed into specialized ML, Bayesian, and Parsimony engines
- **Orchestration System**: `AnalysisOrchestrator` coordinates multiple analysis types
- **Resource Management**: Context managers and proper cleanup throughout
- **Configuration Management**: Centralized constants and configuration validation
- **Memory Management**: Active memory monitoring and optimization
- **File Organization**: Timestamp-based organized output directories

#### Import System Resolution
- **Root Cause Identified**: Relative imports failing when modules run as top-level packages
- **Solution Implemented**: Complete migration to absolute imports with `src.` prefix
- **Result**: Robust, complication-free import system across 27 files and 86+ import statements

#### Testing and Validation
- **Comprehensive Testing**: 8-phase test suite covering all system components
- **Import Verification**: All modules load correctly without fallback to legacy code
- **Integration Testing**: Full end-to-end analysis workflow verified
- **Error Handling**: Exception hierarchy tested with proper context information

#### Performance and Reliability Improvements
- **Memory Optimization**: Active monitoring with peak usage tracking
- **Thread Management**: Adaptive thread calculation based on system resources
- **Error Recovery**: Graceful degradation and detailed error messages
- **Progress Tracking**: Clean console output with dynamic progress updates

## Current Architecture (v2.0 Modular System)

### Entry Points
- **src/panDecay.py**: Main entry point using modular architecture with orchestration
- **src/pandecay_main.py**: Legacy implementation (fallback when needed)

### Core Modular Components
- **src/analysis/engines/**: Specialized analysis engines
  - `base_engine.py`: Abstract base class for all engines
  - `ml_engine.py`: Maximum Likelihood analysis with PAUP*
  - `bayesian_engine.py`: Bayesian analysis with MrBayes
  - `parsimony_engine.py`: Parsimony analysis with Bremer support
- **src/orchestration/**: Analysis coordination and workflow management
  - `analysis_orchestrator.py`: Coordinates multiple analysis engines
- **src/config/**: Configuration and constants management
  - `constants.py`: Centralized configuration constants
- **src/core/**: Core utilities and infrastructure
  - `file_tracker.py`: Organized output directory management
  - `progress_logger.py`: Clean console progress tracking
- **src/external_tools/**: External software integration
  - `tool_manager.py`: Context-managed tool execution
- **src/utils/**: Shared utilities
  - `thread_calculator.py`: Adaptive thread calculation
  - `memory_manager.py`: Memory monitoring and optimization
- **src/exceptions/**: Comprehensive error handling
  - `analysis_exceptions.py`: Analysis-specific exceptions with context

### Model Specification System
```
DNA Data: --nst [1|2|6] --rates [equal|gamma|propinv|invgamma]
Protein Data: --protein-model [WAG|LG|JTT|...] --rates [equal|gamma|propinv|invgamma]
Discrete Data: --discrete-model [Mk|...] --rates [equal|gamma|propinv|invgamma]
```

## Configuration Files
- **examples/configs/example_config.yaml**: Updated with new parameter structure
- All model specifications follow new clean interface

## Testing Status
- **Interface Changes**: Tested with various argument combinations ✓
- **Case-insensitive Parsing**: Verified for all relevant parameters ✓
- **Model Conversion Logic**: Verified for both PAUP* and MrBayes output ✓
- **Modular Architecture**: Comprehensive 8-phase test suite completed ✓
- **Import System**: All 86+ import statements verified working ✓
- **Integration Testing**: End-to-end analysis workflow validated ✓
- **Error Handling**: Exception hierarchy and context tested ✓
- **Resource Management**: Memory monitoring and cleanup verified ✓

## Known Issues
- **Import complexities**: ✅ RESOLVED - Robust absolute import system implemented
- **Monolithic architecture**: ✅ RESOLVED - Fully decomposed into modular system
- **Resource management**: ✅ RESOLVED - Context managers and proper cleanup
- **Error handling**: ✅ RESOLVED - Comprehensive exception hierarchy
- **Code duplication**: ✅ RESOLVED - Shared utilities and focused components
- No critical issues currently identified

## Next Potential Improvements
- **Performance Optimization**: Further memory and processing optimizations
- **Additional Analysis Engines**: Support for specialized phylogenetic methods
- **Enhanced Visualization**: Extended plotting and output options
- **Configuration System**: Advanced YAML/TOML configuration features
- **Plugin Architecture**: Framework for custom analysis extensions

## Latest Enhancements (July 29, 2025)

### 7. User Interface Polish and Professional Presentation
- **Runtime Parameters Banner**: Replaced jagged table with professional fixed-width banner
  - Consistent 80-character width across all terminal sizes
  - Prominent software citation display
  - Organized parameter sections with clean separators
- **Bug Fixes**: Resolved critical markdown report generation error
  - Fixed BioPython alignment attribute access issue
  - Restored proper file naming in output reports
- **Documentation Synchronization**: Comprehensive update of all documentation files
  - README.md enhanced with current functionality descriptions
  - progress.md updated with July 29 development session details
  - dev/CLAUDE.md expanded with recent feature implementations
  - Planning for comprehensive output interpretation guides

### 8. Production Readiness Features
- **Professional Presentation**: Software now displays with academic-quality formatting
- **Citation Integration**: Proper attribution included in all user-facing output
- **Enhanced Reporting**: Markdown reports now generate successfully with comprehensive content
- **File Organization**: Maintained timestamp-based directory structure for organized output

## Development Notes
- **Architectural Milestone**: Successfully transitioned from monolithic to modular design
- **Import System**: Robust solution eliminates all import complications
- **Testing**: Comprehensive validation ensures reliability and maintainability
- **No backward compatibility**: Clean slate approach enables optimal design
- **Code Quality**: Follows SOLID principles with clear separation of concerns
- **Professional Polish**: Academic-quality presentation with proper citations
- **Documentation**: All changes reflected in comprehensive documentation with current functionality