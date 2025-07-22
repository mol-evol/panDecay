# panDecay API Reference

This document provides comprehensive API documentation for panDecay's classes, methods, and configuration schema. Use this reference when extending panDecay or integrating it as a library.

## Table of Contents

1. [Core Classes](#core-classes)
2. [Configuration Models](#configuration-models)  
3. [Async Processing](#async-processing)
4. [Visualization Systems](#visualization-systems)
5. [Utility Functions](#utility-functions)
6. [Configuration Schema](#configuration-schema)
7. [Usage Examples](#usage-examples)

## Core Classes

### panDecayIndices

The main analysis engine for phylogenetic decay calculations.

#### Constructor
```python
panDecayIndices(
    alignment_file: Union[str, Path],
    analysis_type: str = "ml",
    model: str = "GTR", 
    gamma: bool = True,
    invariant: bool = False,
    threads: Union[int, str] = "auto",
    **kwargs
)
```

**Parameters:**
- `alignment_file`: Path to sequence alignment (FASTA, NEXUS, or PHYLIP)
- `analysis_type`: Analysis framework ("ml", "bayesian", "parsimony", or combinations)
- `model`: Evolutionary model (GTR, HKY, JC for DNA; AUTO for proteins)
- `gamma`: Enable gamma rate variation
- `invariant`: Enable invariant sites
- `threads`: Number of threads ("auto" or integer)

#### Key Methods

##### `run_analysis() -> Dict[str, Any]`
Main entry point that orchestrates the complete analysis workflow.

**Returns:** Dictionary containing analysis results with decay indices, statistical tests, and metadata.

##### `_ml_analysis() -> Dict[str, Any]`
Performs maximum likelihood decay analysis with AU testing.

**Process:**
1. ML tree search with PAUP*
2. Constraint generation for each clade
3. Constrained tree searches
4. AU test for statistical significance

##### `_bayesian_analysis() -> Dict[str, Any]`
Performs Bayesian decay analysis using MrBayes.

**Process:**
1. Bayesian tree search with MrBayes
2. Marginal likelihood estimation (stepping-stone sampling)
3. Constrained analyses with topology constraints
4. Bayes factor calculations

##### `_parsimony_analysis() -> Dict[str, Any]`
Traditional parsimony-based Bremer support calculation.

**Process:**
1. Most parsimonious tree search
2. Constrained parsimony searches
3. Step difference calculations (Bremer support)

##### `parse_constraints(constraint_string: str = None) -> List[List[str]]`
Parses user-defined constraint specifications.

**Parameters:**
- `constraint_string`: Semicolon-separated constraint definitions

**Returns:** List of taxon groups for constraint testing.

##### `should_test_clade(clade_taxa: List[str], user_constraints: List[List[str]]) -> bool`
Determines whether a clade should be tested based on constraint mode.

**Parameters:**
- `clade_taxa`: Taxon names in the clade
- `user_constraints`: User-defined constraints

**Returns:** Boolean indicating whether to test the clade.

#### Properties

- `alignment_file: Path` - Input alignment file path
- `alignment_format: str` - Detected alignment format
- `temp_path: Path` - Temporary working directory
- `threads: int` - Number of processing threads
- `analysis_type: str` - Analysis framework(s)
- `constraint_mode: str` - Constraint testing mode ("all" or "specific")

### ExternalToolRunner

Interface for executing external phylogenetic software.

#### Constructor
```python
ExternalToolRunner(
    temp_path: Path,
    paup_path: str = "paup",
    mrbayes_path: str = "mb", 
    debug: bool = False
)
```

#### Key Methods

##### `run_paup_command_file(cmd_file: str, log_file: str, timeout: int = None) -> bool`
Executes PAUP* command file with timeout and error handling.

**Parameters:**
- `cmd_file`: PAUP* command file name (in temp directory)
- `log_file`: Output log file name
- `timeout`: Execution timeout in seconds

**Returns:** True if successful, raises RuntimeError on failure.

##### `run_mrbayes(nexus_file: str, output_prefix: str) -> Optional[float]`
Executes MrBayes analysis and extracts marginal likelihood.

**Parameters:**
- `nexus_file`: MrBayes input file name
- `output_prefix`: Output file prefix

**Returns:** Marginal likelihood value or None if failed.

##### `parse_paup_output(log_file: str) -> Dict[str, Any]`
Parses PAUP* output log for likelihood scores and tree information.

**Returns:** Dictionary with parsed results including likelihood scores, trees, and AU test results.

### TreeManager

Handles tree operations and constraint generation.

#### Constructor
```python
TreeManager(temp_path: Path)
```

#### Key Methods

##### `clean_newick_tree(tree_file: Path) -> Path`
Cleans and validates Newick format trees.

**Parameters:**
- `tree_file`: Input tree file path

**Returns:** Path to cleaned tree file.

##### `generate_constraint_trees(tree: str, taxa_list: List[str]) -> List[Dict[str, Any]]`
Generates constraint definitions for all internal branches.

**Parameters:**
- `tree`: Newick format tree string
- `taxa_list`: List of all taxon names

**Returns:** List of constraint dictionaries with clade definitions.

##### `annotate_tree_with_support(tree: str, support_data: Dict[str, float]) -> str`
Annotates tree with branch support values.

**Parameters:**
- `tree`: Input Newick tree
- `support_data`: Dictionary mapping clades to support values

**Returns:** Annotated Newick tree string.

### DatasetNormalizer

Performs cross-study normalization calculations.

#### Constructor
```python
DatasetNormalizer()
```

#### Key Methods

##### `apply_dataset_relative_normalizations(decay_indices: Dict[str, Any], calculate_dataset_relative: bool)`
Applies dataset-relative normalization to decay indices.

**Parameters:**
- `decay_indices`: Decay index data dictionary
- `calculate_dataset_relative`: Whether to calculate dataset-relative metrics

**Modifies:** Input dictionary with normalized values including percentile ranks, z-scores, and dataset-relative rankings.

##### `calculate_effect_sizes(site_data: Dict[str, float], ml_decay: float) -> Dict[str, float]`
Calculates effect size metrics for branch support.

**Parameters:**
- `site_data`: Per-site likelihood analysis data
- `ml_decay`: ML decay index value

**Returns:** Dictionary with effect size calculations.

## Configuration Models

### PanDecayConfig

Root configuration model with comprehensive validation.

```python
class PanDecayConfig(BaseModel):
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    computational: ComputationalConfig = Field(default_factory=ComputationalConfig)
    visualization: VisualizationConfig = Field(default_factory=VisualizationConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)
    external_tools: ExternalToolsConfig = Field(default_factory=ExternalToolsConfig)
```

### AnalysisConfig

Analysis-specific configuration parameters.

```python
class AnalysisConfig(BaseModel):
    analysis_type: Literal["ml", "bayesian", "parsimony", "ml+bayesian", "all"] = "ml"
    model: str = "GTR"
    gamma: bool = True
    invariant: bool = False
    constraint_mode: Literal["all", "specific"] = "all"
    normalization: Literal["none", "basic", "full"] = "basic"
    site_analysis: bool = False
    bootstrap_replicates: int = 0
```

### ComputationalConfig

Performance and resource configuration.

```python
class ComputationalConfig(BaseModel):
    threads: Union[int, str] = "auto"
    async_constraints: bool = False
    max_async_workers: int = 4
    constraint_timeout: int = 1800
    ml_timeout: int = 7200
    memory_limit: Optional[str] = None
```

### VisualizationConfig

Visualization and plotting configuration.

```python
class VisualizationConfig(BaseModel):
    enabled: bool = True
    format: Literal["static", "interactive", "both"] = "both"
    static_formats: List[str] = Field(default=["png", "pdf"])
    interactive_format: str = "html"
    dpi: int = 300
    theme: str = "publication"
```

## Async Processing

### AsyncConstraintProcessor

Manages parallel constraint analysis execution.

#### Constructor
```python
AsyncConstraintProcessor(
    max_workers: int = 4,
    timeout: int = 1800
)
```

#### Key Methods

##### `process_constraints_async(tasks: List[ConstraintTask]) -> List[ConstraintResult]`
Executes constraint analyses in parallel using ThreadPoolExecutor.

**Parameters:**
- `tasks`: List of constraint analysis tasks

**Returns:** List of analysis results with timing and error information.

### ConstraintTask

Represents a single constraint analysis task.

```python
@dataclass
class ConstraintTask:
    clade_id: str
    taxa: List[str]
    command_file: str
    log_file: str
    timeout: int = 1800
```

### ConstraintResult

Contains results from constraint analysis execution.

```python
@dataclass  
class ConstraintResult:
    clade_id: str
    success: bool
    execution_time: float
    result_data: Optional[Dict[str, Any]] = None
    error_message: Optional[str] = None
```

## Visualization Systems

### PlotManager

Static matplotlib-based visualization system for phylogenetic analysis results.

#### Constructor
```python
PlotManager(
    output_dir: Path,
    file_tracker: Optional[FileTracker] = None
)
```

#### Key Methods

##### `create_distribution_plot(data: Dict[str, Any], analysis_type: str) -> Path`
Creates publication-ready static distribution plots of support values.

**Parameters:**
- `data`: Analysis results data
- `analysis_type`: Type of analysis performed

**Returns:** Dictionary mapping visualization types to file paths.

##### `create_decay_scatter_plot(data: Dict[str, Any]) -> Tuple[Optional[Path], Optional[Path]]`
Creates scatter plot of decay indices vs. significance.

**Returns:** Tuple of (static_file_path, interactive_file_path).

##### `create_support_histogram(data: Dict[str, Any]) -> Tuple[Optional[Path], Optional[Path]]`
Creates histogram of support value distributions.

**Returns:** Tuple of (static_file_path, interactive_file_path).

## Utility Functions

### Configuration Loading

##### `load_configuration(config_path: Union[str, Path]) -> PanDecayConfig`
Loads and validates configuration from YAML, TOML, or INI files.

**Parameters:**
- `config_path`: Path to configuration file

**Returns:** Validated configuration object.

**Raises:** `ConfigurationError` for invalid configurations.

##### `detect_config_format(config_path: Path) -> str`
Automatically detects configuration file format.

**Returns:** Format string ("yaml", "toml", or "ini").

### Logging Setup

##### `setup_logging(debug_mode: bool = False, log_file: str = None) -> logging.Logger`
Configures application logging with appropriate handlers and formatting.

**Parameters:**
- `debug_mode`: Enable debug-level logging
- `log_file`: Optional log file path

**Returns:** Configured logger instance.

### Path Utilities

##### `get_display_path(path: Union[str, Path], max_length: int = 50) -> str`
Creates shortened display version of file paths.

**Parameters:**
- `path`: File path to shorten
- `max_length`: Maximum display length

**Returns:** Shortened path string with ellipsis if needed.

## Configuration Schema

### Complete YAML Configuration Example

```yaml
# panDecay Configuration File (YAML)
analysis:
  analysis_type: "ml"           # ml, bayesian, parsimony, ml+bayesian, all
  model: "GTR"                  # GTR, HKY, JC, AUTO (for proteins)
  gamma: true                   # Gamma rate variation
  invariant: false              # Invariant sites
  constraint_mode: "all"        # all, specific
  normalization: "basic"        # none, basic, full
  site_analysis: false          # Per-site likelihood analysis
  bootstrap_replicates: 0       # Bootstrap support (0 = disabled)

computational:
  threads: "auto"               # auto or integer
  async_constraints: true       # Enable parallel constraint processing
  max_async_workers: 4          # Maximum parallel workers
  constraint_timeout: 1800      # Timeout per constraint (seconds)
  ml_timeout: 7200             # ML search timeout (seconds)
  memory_limit: "8G"           # Optional memory limit

visualization:
  enabled: true                 # Enable visualization
  format: "both"               # static, interactive, both
  static_formats: ["png", "pdf"] # Static output formats
  interactive_format: "html"    # Interactive format
  dpi: 300                     # Resolution for static plots
  theme: "publication"         # Styling theme

output:
  output_dir: "."              # Output directory
  prefix: ""                   # File prefix
  keep_files: false            # Retain temporary files
  debug: false                 # Debug mode
  detailed_report: true        # Generate detailed markdown report

external_tools:
  paup_path: "paup"           # PAUP* executable path
  mrbayes_path: "mb"          # MrBayes executable path
  paup_block_file: null       # Optional PAUP* command block

mrbayes:
  ngen: 1000000               # MCMC generations
  chains: 4                   # Number of chains
  burnin: 0.25               # Burnin fraction
  sample_freq: 1000          # Sampling frequency
  use_mpi: false             # MPI support
  use_beagle: false          # BEAGLE acceleration
```

### Schema Validation

All configuration files are validated against Pydantic models that provide:

- **Type checking**: Automatic type conversion and validation
- **Range validation**: Numeric ranges and enum constraints  
- **Required fields**: Identification of missing required parameters
- **Default values**: Sensible defaults for optional parameters
- **Error messages**: Clear, actionable error descriptions

## Usage Examples

### Programmatic Usage

```python
from pathlib import Path
from src.panDecay import panDecayIndices
from src.config_loader import load_configuration

# Basic usage
calc = panDecayIndices(
    alignment_file="examples/data/alignment.fas",
    analysis_type="ml",
    model="GTR",
    gamma=True,
    threads=4
)

results = calc.run_analysis()
print(f"Analysis completed with {len(results['decay_indices'])} clades")

# With configuration file
config = load_configuration("examples/configs/example_config.yaml")
calc = panDecayIndices(
    alignment_file="examples/data/alignment.fas",
    **config.analysis.dict()
)

results = calc.run_analysis()
```

### Async Processing Example

```python
from src.async_constraint_processor import AsyncConstraintProcessor, ConstraintTask

# Create constraint tasks
tasks = [
    ConstraintTask(
        clade_id=f"Clade_{i}",
        taxa=clade_taxa,
        command_file=f"constraint_{i}.nex",
        log_file=f"constraint_{i}.log"
    )
    for i, clade_taxa in enumerate(constraint_list)
]

# Process in parallel
processor = AsyncConstraintProcessor(max_workers=4, timeout=1800)
results = processor.process_constraints_async(tasks)

# Check results
successful = [r for r in results if r.success]
failed = [r for r in results if not r.success]
print(f"Successful: {len(successful)}, Failed: {len(failed)}")
```

### Custom Visualization Example

```python
from src.visualization.plot_manager import PlotManager
from pathlib import Path

# Create plot manager with organized output directory
output_dir = Path("20241015_143022_panDecay_alignment/visualizations/")
plot_manager = PlotManager(output_dir)

# Create static distribution plot
distribution_plot = plot_manager.create_distribution_plot(results, "ml")

print(f"Created distribution plot: {distribution_plot}")

# Create site-specific alignment visualizations (if site analysis performed)
if site_analysis_data:
    alignment_plots = plot_manager.create_alignment_visualizations(
        site_analysis_data, 
        alignment_data
    )
    print(f"Created {len(alignment_plots)} site-specific plots")
```

---

## Error Handling

All API functions use consistent error handling patterns:

- **ConfigurationError**: Invalid configuration parameters
- **RuntimeError**: External tool execution failures  
- **FileNotFoundError**: Missing input files
- **TimeoutError**: Operation timeouts
- **ValidationError**: Pydantic schema validation failures

Always wrap API calls in appropriate try-catch blocks and check return values for optional operations.