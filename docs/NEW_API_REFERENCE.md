# panDecay API Reference (v2.0 Modular Architecture)

This document provides comprehensive API documentation for the panDecay v2.0 system, featuring a fully decomposed modular architecture with specialized analysis engines, orchestration, and comprehensive utility components. Updated July 29, 2025.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Analysis Engines](#analysis-engines)
3. [Orchestration System](#orchestration-system)
4. [Configuration Management](#configuration-management)
5. [Utility Modules](#utility-modules)
6. [Exception Handling](#exception-handling)
7. [Memory Management](#memory-management)
8. [Usage Examples](#usage-examples)

## Architecture Overview

The refactored panDecay system follows a modular architecture with clear separation of concerns:

```
src/
├── analysis/engines/          # Specialized analysis engines
├── orchestration/            # Analysis coordination
├── config/                   # Configuration management
├── core/                     # Core infrastructure (FileTracker, ProgressLogger)
├── external_tools/           # Tool management
├── utils/                   # Shared utilities
├── exceptions/              # Error handling
├── visualization/           # Plot and visualization management
└── io/                     # Input/output operations
```

### Design Principles

- **Single Responsibility**: Each module has one clear purpose
- **Dependency Injection**: Easy testing and mocking
- **Strategy Pattern**: Pluggable analysis algorithms
- **Resource Management**: Proper cleanup and limits
- **Error Handling**: Comprehensive exception hierarchy

## Analysis Engines

### Base Engine

All analysis engines inherit from the base `AnalysisEngine` class.

#### `AnalysisEngine` (Abstract Base Class)

```python
from src.analysis.engines.base_engine import AnalysisEngine, AnalysisData, AnalysisResult

class AnalysisEngine(ABC):
    def __init__(self, temp_dir: Path, debug: bool = False):
        """Initialize analysis engine."""
        
    @abstractmethod
    def get_analysis_type(self) -> str:
        """Return the analysis type identifier."""
        
    @abstractmethod
    def validate_inputs(self, data: AnalysisData) -> bool:
        """Validate analysis inputs."""
        
    @abstractmethod
    def analyze(self, data: AnalysisData) -> AnalysisResult:
        """Perform the phylogenetic analysis."""
```

#### Data Classes

##### `AnalysisData`

```python
@dataclass
class AnalysisData:
    alignment_file: Path
    tree_file: Optional[Path]
    temp_dir: Path
    model_settings: Dict[str, Any]
    constraints: List[Dict[str, Any]]
    timeouts: AnalysisTimeouts
```

##### `AnalysisResult`

```python
@dataclass
class AnalysisResult:
    analysis_type: str
    success: bool
    branches_tested: int
    decay_indices: Dict[str, Any]
    execution_time: float
    error_message: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
```

### ML Analysis Engine

Maximum Likelihood analysis using PAUP* with AU testing.

```python
from src.analysis.engines.ml_engine import MLAnalysisEngine

engine = MLAnalysisEngine(
    temp_dir=Path("temp/"),
    paup_path="paup",
    debug=False
)

result = engine.analyze(analysis_data)
```

#### Key Methods

- **`_build_ml_tree()`**: Constructs maximum likelihood tree
- **`_test_constraints()`**: Tests topological constraints
- **`_run_au_test()`**: Performs Approximately Unbiased test
- **`_parse_likelihood_score()`**: Parses PAUP* output

### Bayesian Analysis Engine

Bayesian analysis using MrBayes with marginal likelihood comparison.

```python
from src.analysis.engines.bayesian_engine import BayesianAnalysisEngine

engine = BayesianAnalysisEngine(
    temp_dir=Path("temp/"),
    mrbayes_path="mb",
    debug=False
)

result = engine.analyze(analysis_data)
```

#### Key Methods

- **`_run_unconstrained_analysis()`**: Runs unconstrained Bayesian analysis
- **`_run_constrained_analyses()`**: Runs constrained analyses
- **`_parse_marginal_likelihood()`**: Parses marginal likelihood from output
- **`_calculate_decay_indices()`**: Calculates Bayes factors

### Parsimony Analysis Engine

Parsimony analysis using PAUP* with Bremer support calculation.

```python
from src.analysis.engines.parsimony_engine import ParsimonyAnalysisEngine

engine = ParsimonyAnalysisEngine(
    temp_dir=Path("temp/"),
    paup_path="paup",
    debug=False
)

result = engine.analyze(analysis_data)
```

#### Key Methods

- **`_build_parsimony_trees()`**: Constructs most parsimonious trees
- **`_test_constraints()`**: Tests constraint trees
- **`_parse_parsimony_score()`**: Parses parsimony scores

## Orchestration System

The orchestration system coordinates multiple analysis engines and manages the overall workflow.

### AnalysisOrchestrator

```python
from src.orchestration.analysis_orchestrator import AnalysisOrchestrator

orchestrator = AnalysisOrchestrator(
    temp_dir=Path("temp/"),
    debug=False
)

# Register analysis engines
orchestrator.create_engines(
    analysis_modes=["ml", "bayesian"],
    paup_path="/usr/local/bin/paup",
    mrbayes_path="/usr/local/bin/mb"
)

# Run coordinated analysis
result = orchestrator.run_analysis(
    analysis_modes=["ml", "bayesian"],
    data=analysis_data
)
```

#### Key Methods

##### `register_engine(analysis_type: str, engine_class: Type[AnalysisEngine])`
Register a new analysis engine type.

##### `create_engines(analysis_modes: List[str], **kwargs)`
Create engine instances for specified analysis modes.

##### `run_analysis(analysis_modes: List[str], data: AnalysisData) -> OrchestrationResult`
Execute coordinated analysis across multiple engines.

##### `parse_analysis_mode_string(mode_string: str) -> List[str]`
Parse analysis mode strings like "ml+bayesian" or "all".

#### OrchestrationResult

```python
@dataclass
class OrchestrationResult:
    success: bool
    total_execution_time: float
    analysis_results: Dict[str, AnalysisResult]
    error_message: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None
```

## Configuration Management

Centralized configuration system with constants and type validation.

### Constants

```python
from src.config.constants import (
    FileNames, AnalysisThresholds, AnalysisTimeouts, 
    ResourceLimits, AnalysisModes, ValidationLimits
)

# File naming
alignment_file = temp_dir / FileNames.NEXUS_ALIGNMENT
ml_tree = temp_dir / FileNames.ML_TREE

# Analysis parameters
significance = AnalysisThresholds.AU_SIGNIFICANCE  # 0.05
min_ess = AnalysisThresholds.MIN_ESS  # 200

# Timeouts
timeouts = AnalysisTimeouts(
    ml_timeout=3600,
    bayesian_timeout=7200
)

# Resource limits
limits = ResourceLimits(
    max_memory_mb=2000,
    max_processes=8
)
```

### Analysis Modes

```python
# Single analysis types
AnalysisModes.ML           # "ml"
AnalysisModes.BAYESIAN     # "bayesian" 
AnalysisModes.PARSIMONY    # "parsimony"

# Combined analysis types
AnalysisModes.ML_BAYESIAN  # "ml+bayesian"
AnalysisModes.ALL          # "all"
```

## Utility Modules

### Format Detection

```python
from src.utils.format_detectors import (
    ConfigurationFormatDetector, 
    AlignmentFormatMapper,
    AnalysisTypeConfiguration
)

# Configuration format detection
format_detector = ConfigurationFormatDetector()
config_format = format_detector.detect_format(config_content)

# Alignment format detection
alignment_mapper = AlignmentFormatMapper()
alignment_format = alignment_mapper.detect_format("data.fas", content)

# Analysis type configuration
config = AnalysisTypeConfiguration(
    has_ml=True, 
    has_bayesian=True, 
    has_parsimony=False
)
layout = config.get_layout_configuration()
```

### Configuration Conversion

```python
from src.utils.config_converters import ConfigurationValueConverter

converter = ConfigurationValueConverter()

# Type conversion with validation
gamma = converter.convert('gamma', 'true')      # -> True
threads = converter.convert('threads', '4')     # -> 4
norm = converter.convert('normalization', 'basic')  # -> 'basic'
```

### Thread Calculation

```python
from src.utils.thread_calculator import ThreadCalculator, AdaptiveThreadCalculator

# Basic thread calculation
calculator = ThreadCalculator()
optimal_threads = calculator.calculate_optimal_threads("auto")

# Adaptive calculation
adaptive = AdaptiveThreadCalculator()
ml_threads = adaptive.calculate_for_analysis_type('ml', 'auto')
memory_safe_threads = adaptive.calculate_for_data_size(100.0, 'auto')
```

### Memory Management

```python
from src.utils.memory_manager import (
    MemoryMonitor, MemoryOptimizer, 
    ChunkedFileReader, memory_limited_operation
)

# Memory monitoring
monitor = MemoryMonitor()
monitor.start_monitoring(interval_seconds=5.0)
snapshot = monitor.take_snapshot("operation_start")

# Memory optimization
optimizer = MemoryOptimizer()
recommendations = optimizer.optimize_for_large_alignment(150.0)
estimates = optimizer.estimate_analysis_memory_requirements(
    alignment_size_mb=100.0,
    num_sequences=50,
    analysis_types=['ml', 'bayesian']
)

# Memory-limited operations
with memory_limited_operation(100.0, "large_analysis") as monitor:
    # Perform memory-intensive operations
    result = heavy_computation()

# Chunked file reading
reader = ChunkedFileReader(chunk_size=8192)
for chunk in reader.read_chunks(large_file):
    process_chunk(chunk)
```

### Tree Utilities

```python
from src.utils.tree_utils import (
    prepare_starting_tree, validate_tree_file, 
    clean_tree_file, count_taxa_in_tree
)

# Tree preparation
prepared_tree = prepare_starting_tree(
    starting_tree_path, 
    temp_dir, 
    "prepared.tre"
)

# Tree validation
is_valid = validate_tree_file(tree_file, require_branch_lengths=True)

# Tree cleaning
clean_tree_file(tree_file, remove_metadata=True)

# Taxa counting
taxa_count = count_taxa_in_tree(tree_file)
```

### Script Generation

```python
from src.utils.script_generators import PAUPScriptGenerator, MrBayesScriptGenerator

# PAUP* script generation
paup_generator = PAUPScriptGenerator(debug=False)
ml_script = paup_generator.generate_ml_search_script(
    alignment_file="alignment.nex",
    threads=4,
    model_settings={'nst': 6, 'rates': 'gamma'}
)

# MrBayes script generation
mb_generator = MrBayesScriptGenerator(debug=False)
bayes_script = mb_generator.generate_bayesian_analysis_script(
    alignment_file="alignment.nex",
    model_settings={'nst': 6},
    mcmc_settings={'ngen': 1000000, 'nchains': 4}
)
```

## Exception Handling

Comprehensive exception hierarchy for different error types.

### Exception Hierarchy

```
Exception
└── PanDecayError (base)
    ├── AnalysisError
    ├── ValidationError  
    ├── ConfigurationError
    ├── ExternalToolError
    │   ├── ToolNotFoundError
    │   ├── ToolExecutionError
    │   └── ToolTimeoutError
    └── FileOperationError
        ├── AlignmentParsingError
        └── TreeParsingError
```

### Usage Examples

```python
from src.exceptions.analysis_exceptions import AnalysisError, ValidationError
from src.exceptions.tool_exceptions import ExternalToolError, ToolNotFoundError
from src.exceptions.io_exceptions import FileOperationError

# Analysis errors with context
try:
    result = engine.analyze(data)
except AnalysisError as e:
    print(f"Analysis failed: {e}")
    print(f"Analysis type: {e.analysis_type}")
    print(f"Context: {e.context}")

# Tool errors
try:
    tool_manager.execute_script(script_file)
except ToolNotFoundError as e:
    print(f"Tool not found: {e.tool_name}")
    print(f"Search paths: {e.search_paths}")

# Validation errors
try:
    orchestrator.validate_analysis_request(modes, data)
except ValidationError as e:
    print(f"Validation failed: {e}")
    print(f"Field: {e.field}, Value: {e.value}")
```

## External Tool Management

Centralized management of external tools with proper resource handling.

### ExternalToolManager

```python
from src.external_tools.tool_manager import ExternalToolManager

# Context manager for automatic cleanup
with ExternalToolManager(temp_dir, debug=False) as tool_manager:
    # Check tool availability
    tool_manager.check_tool_availability("paup")
    
    # Execute script file
    tool_manager.execute_script_file(
        tool_path="paup",
        script_file="analysis.nex",
        timeout=3600,
        log_file="output.log"
    )
    
    # Execute direct command
    result = tool_manager.execute_command(
        command=["paup", "-n", "script.nex"],
        timeout=1800
    )
```

## Usage Examples

### Basic Single Analysis

```python
from pathlib import Path
from src.orchestration.analysis_orchestrator import AnalysisOrchestrator
from src.config.constants import AnalysisModes, AnalysisTimeouts

# Setup
temp_dir = Path("temp_analysis")
orchestrator = AnalysisOrchestrator(temp_dir, debug=False)

# Create analysis data
analysis_data = orchestrator.create_analysis_data(
    alignment_file=Path("data/alignment.nex"),
    temp_dir=temp_dir,
    model_settings={
        'nst': 6,  # GTR model
        'rates': 'gamma',
        'threads': 4
    },
    constraints=[
        {'id': 'primates', 'taxa': ['Human', 'Chimp', 'Gorilla']},
        {'id': 'rodents', 'taxa': ['Mouse', 'Rat']}
    ]
)

# Create engines
orchestrator.create_engines([AnalysisModes.ML])

# Run analysis
result = orchestrator.run_analysis([AnalysisModes.ML], analysis_data)

if result.success:
    ml_result = result.analysis_results['ml']
    print(f"Analysis completed in {result.total_execution_time:.1f} seconds")
    print(f"Branches tested: {ml_result.branches_tested}")
    for constraint_id, decay_data in ml_result.decay_indices.items():
        print(f"Constraint {constraint_id}: {decay_data}")
```

### Multi-Analysis Workflow

```python
# Run multiple analysis types
modes = orchestrator.parse_analysis_mode_string("ml+bayesian")
orchestrator.create_engines(modes, 
                           paup_path="/usr/local/bin/paup",
                           mrbayes_path="/usr/local/bin/mb")

result = orchestrator.run_analysis(modes, analysis_data)

# Compare results across analysis types
if result.success:
    for analysis_type, analysis_result in result.analysis_results.items():
        print(f"\n{analysis_type.upper()} Results:")
        for constraint_id, decay_data in analysis_result.decay_indices.items():
            print(f"  {constraint_id}: {decay_data}")
```

### Memory-Optimized Large Dataset Analysis

```python
from src.utils.memory_manager import MemoryOptimizer
from src.utils.thread_calculator import AdaptiveThreadCalculator

# Optimize for large dataset
optimizer = MemoryOptimizer()
calculator = AdaptiveThreadCalculator()

alignment_size_mb = 150.0  # 150MB alignment
recommendations = optimizer.optimize_for_large_alignment(alignment_size_mb)

# Adjust configuration based on recommendations
if recommendations['use_chunked_reading']:
    # Enable chunked processing
    pass

optimal_threads = calculator.calculate_for_data_size(alignment_size_mb, 'auto')

# Create memory-aware analysis data
analysis_data = orchestrator.create_analysis_data(
    alignment_file=large_alignment_file,
    temp_dir=temp_dir,
    model_settings={
        'nst': 2,  # Simpler model for large datasets
        'rates': 'equal',  # No gamma for speed
        'threads': optimal_threads
    },
    constraints=constraints
)

# Run with memory monitoring
with memory_limited_operation(500.0, "large_dataset_analysis") as monitor:
    result = orchestrator.run_analysis([AnalysisModes.ML], analysis_data)
    
    memory_summary = monitor.get_memory_summary()
    print(f"Peak memory usage: {memory_summary['peak_memory_mb']:.1f} MB")
```

### Custom Analysis Engine

```python
from src.analysis.engines.base_engine import AnalysisEngine, AnalysisData, AnalysisResult

class CustomAnalysisEngine(AnalysisEngine):
    def get_analysis_type(self) -> str:
        return "custom"
    
    def validate_inputs(self, data: AnalysisData) -> bool:
        return data.alignment_file.exists()
    
    def analyze(self, data: AnalysisData) -> AnalysisResult:
        # Custom analysis implementation
        return AnalysisResult(
            analysis_type="custom",
            success=True,
            branches_tested=len(data.constraints),
            decay_indices={'custom_metric': 0.5},
            execution_time=1.0
        )

# Register and use custom engine
orchestrator.register_engine("custom", CustomAnalysisEngine)
orchestrator.create_engines(["custom"])
result = orchestrator.run_analysis(["custom"], analysis_data)
```

## Configuration Integration

### Type-Safe Configuration

```python
from src.utils.config_converters import ConfigurationValueConverter
from src.utils.format_detectors import AnalysisTypeConfiguration

converter = ConfigurationValueConverter()

# Process configuration with validation
config_dict = {
    'analysis_type': 'ml+bayesian',
    'gamma': 'true',
    'threads': '4',
    'gamma_shape': '0.5',
    'normalization': 'basic'
}

processed_config = {}
for param, value in config_dict.items():
    processed_config[param] = converter.convert(param, value)

# Create analysis configuration
analysis_modes = orchestrator.parse_analysis_mode_string(
    processed_config['analysis_type']
)
analysis_config = AnalysisTypeConfiguration(
    has_ml='ml' in analysis_modes,
    has_bayesian='bayesian' in analysis_modes,
    has_parsimony='parsimony' in analysis_modes
)

# Use configuration
if analysis_config.requires_model_configuration():
    model_settings = {
        'nst': 6,
        'rates': 'gamma' if processed_config['gamma'] else 'equal',
        'threads': processed_config['threads']
    }
```

## Recent Enhancements (July 29, 2025)

### Professional User Interface
- **Runtime Parameters Banner**: Clean, fixed-width banner format with academic citation
- **Progress Tracking**: Enhanced console output with ProgressLogger integration
- **File Organization**: Timestamp-based directory structure via FileTracker

### Enhanced Reporting
- **Markdown Reports**: Comprehensive analysis reports with site analysis integration
- **Output Organization**: Structured output files with consistent formatting
- **Citation Integration**: Proper academic attribution in all user-facing output

### Core Infrastructure
- **FileTracker**: Organized output directory management with timestamp-based structure
- **ProgressLogger**: Clean console progress tracking with dynamic updates
- **Professional Presentation**: Academic-quality software presentation standards

This API reference provides comprehensive documentation for the panDecay v2.0 system, focusing on the modular architecture and enhanced capabilities for extensibility, testing, maintainability, and professional presentation.