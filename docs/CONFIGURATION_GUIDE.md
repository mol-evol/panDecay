# panDecay Configuration Guide

Comprehensive guide to configuring panDecay using YAML, TOML, and legacy INI formats.

## Table of Contents

- [Configuration Formats](#configuration-formats)
- [Complete Configuration Reference](#complete-configuration-reference)
- [Validation and Testing](#validation-and-testing)
- [Common Configurations](#common-configurations)
- [Advanced Features](#advanced-features)

## Configuration Formats

panDecay supports three configuration formats:

### YAML (Recommended)
```yaml
# Modern format with validation and better structure
analysis:
  type: ml+bayesian
  model: GTR
  gamma: true
```

### TOML
```toml
# Alternative modern format
[analysis]
type = "ml+bayesian"
model = "GTR"
gamma = true
```

### INI (Legacy)
```ini
# Backward compatibility
[analysis]
analysis_type = ml+bayesian
model = GTR
gamma = true
```

## Complete Configuration Reference

### Input/Output Section

```yaml
input_output:
  # Input alignment file (overridden by command line)
  alignment_file: "alignment.fas"
  
  # Output directory for results
  output_dir: "output"
  
  # Output file prefix
  output_prefix: "pandecay_results"
  
  # Keep temporary files for debugging
  keep_files: false
  
  # Debug mode - retain all intermediate files
  debug: false
  
  # Temporary directory (auto-detected if not specified)
  temp_dir: "/tmp"
```

### Analysis Configuration

```yaml
analysis:
  # Analysis type: ml, bayesian, parsimony, all, ml+bayesian, ml+parsimony, bayesian+parsimony
  type: "ml+bayesian"
  
  # Evolutionary model for ML and Bayesian analysis
  model: "GTR"
  
  # Enable gamma rate variation
  gamma: true
  
  # Number of gamma categories (if gamma enabled)
  gamma_categories: 4
  
  # Enable proportion of invariant sites
  invariant: true
  
  # Site-specific likelihood analysis
  site_specific: false
  
  # Generate effect size analysis for cross-study comparison
  effect_size: true
  
  # Custom PAUP block for advanced ML settings
  paup_block: null
```

### Computational Settings

```yaml
computational:
  # Number of threads (auto = detect optimal)
  threads: "auto"
  
  # Path to PAUP* executable
  paup_path: "paup"
  
  # Path to MrBayes executable  
  mrbayes_path: "mb"
  
  # Enable async constraint processing (recommended)
  async_constraints: true
  
  # Maximum parallel constraint workers
  max_async_workers: 4
  
  # Timeout per constraint analysis (seconds)
  constraint_timeout: 1800
  
  # Memory limit for analysis processes
  memory_limit: null
```

### Bayesian Analysis Settings

```yaml
bayesian:
  # Number of MCMC generations
  generations: 1000000
  
  # Sampling frequency
  sample_frequency: 1000
  
  # Number of runs
  nruns: 2
  
  # Number of chains per run
  nchains: 4
  
  # Burn-in fraction
  burnin: 0.25
  
  # Temperature for heated chains
  temperature: 0.2
  
  # Enable MPI for parallel processing
  use_mpi: false
  
  # Number of MPI processes (if use_mpi enabled)
  mpi_processes: 4
  
  # Enable BEAGLE library for GPU acceleration
  use_beagle: false
  
  # Marginal likelihood estimation method: stepping_stone, harmonic_mean
  marginal_likelihood_method: "stepping_stone"
  
  # Stepping stone sampling parameters
  stepping_stone:
    nsteps: 50
    step_generations: 10000
    alpha: 0.4
  
  # Convergence diagnostics
  diagnostics:
    # Minimum effective sample size
    min_ess: 200
    
    # Maximum potential scale reduction factor
    max_psrf: 1.1
    
    # Maximum average standard deviation of split frequencies
    max_asdsf: 0.01
```

### Visualization Settings

```yaml
visualization:
  # Output format: static (only supported format)
  format: "static"
  
  # Static plot formats
  static_formats: ["png"]
  
  # Plot resolution (DPI)
  dpi: 300
  
  # Visualization theme: publication, presentation, poster
  theme: "publication"
  
  # Color palette: viridis, plasma, inferno, magma, Set1, Set2, paired
  color_palette: "viridis"
  
  # Figure size (width, height in inches)
  figure_size: [12, 8]
  
  # Font family for plots
  font_family: ["DejaVu Sans", "Helvetica", "Arial"]
  
  # Font sizes
  font_sizes:
    title: 16
    labels: 12
    legend: 10
    ticks: 10
  
  # Generate annotated tree files
  generate_trees: true
  
  # Tree annotation types
  tree_annotations: ["support", "decay", "combined"]
```

### Constraint Settings

```yaml
constraints:
  # Constraint generation strategy: comprehensive, fast, custom
  strategy: "comprehensive"
  
  # Minimum branch length for constraint generation
  min_branch_length: 0.001
  
  # Maximum number of constraints (0 = unlimited)
  max_constraints: 0
  
  # Custom constraint definitions (advanced)
  custom_constraints: []
  
  # Bootstrap support threshold for constraint pruning
  bootstrap_threshold: 0.5
  
  # Skip constraints for very short branches
  skip_short_branches: true
```

## Validation and Testing

### Configuration Validation

```bash
# Validate YAML configuration
python3 config_loader.py validate config.yaml

# Test configuration without running analysis
python3 panDecay.py --config config.yaml --dry-run

# Generate and validate example configuration
python3 panDecay.py --generate-yaml-config test_config.yaml
python3 config_loader.py validate test_config.yaml
```

### Common Validation Errors

#### 1. Invalid Analysis Type
```yaml
# ERROR: Invalid analysis type
analysis:
  type: "invalid_type"

# FIX: Use valid types
analysis:
  type: "ml"  # or bayesian, parsimony, all, ml+bayesian, etc.
```

#### 2. Invalid Model
```yaml
# ERROR: Invalid model
analysis:
  model: "INVALID"

# FIX: Use supported models  
analysis:
  model: "GTR"  # or JC, F81, HKY, K2P, SYM, etc.
```

#### 3. Invalid Thread Count
```yaml
# ERROR: Negative threads
computational:
  threads: -1

# FIX: Use positive integer or "auto"
computational:
  threads: "auto"  # or 1, 2, 4, 8, etc.
```

#### 4. Invalid File Paths
```yaml
# ERROR: Non-existent executable
computational:
  paup_path: "/nonexistent/paup"

# FIX: Provide valid path or use default
computational:
  paup_path: "paup"  # Use PATH lookup
```

## Common Configurations

### Quick ML Analysis

```yaml
# Minimal configuration for fast ML analysis
analysis:
  type: ml
  model: GTR
  gamma: true

computational:
  threads: auto
  async_constraints: true

visualization:
  format: static
```

### Comprehensive Bayesian Analysis

```yaml
# Full Bayesian analysis with all features
analysis:
  type: bayesian
  model: GTR
  gamma: true
  invariant: true
  effect_size: true

computational:
  threads: auto
  async_constraints: true
  max_async_workers: 4

bayesian:
  generations: 2000000
  nruns: 4
  nchains: 4
  burnin: 0.25
  use_mpi: false
  marginal_likelihood_method: stepping_stone

visualization:
  format: both
  theme: publication
  dpi: 300
```

### High-Performance Configuration

```yaml
# Optimized for speed and parallel processing
analysis:
  type: ml
  model: GTR
  site_specific: false

computational:
  threads: auto
  async_constraints: true
  max_async_workers: 8
  constraint_timeout: 900

constraints:
  strategy: fast
  skip_short_branches: true

visualization:
  format: static
  generate_trees: false
```

### Publication-Ready Analysis

```yaml
# Configuration for publication-quality results
analysis:
  type: all
  model: GTR
  gamma: true
  invariant: true
  site_specific: true
  effect_size: true

computational:
  threads: auto
  async_constraints: true

bayesian:
  generations: 5000000
  nruns: 4
  burnin: 0.25
  marginal_likelihood_method: stepping_stone

visualization:
  format: both
  theme: publication
  dpi: 300
  static_formats: [png, pdf, svg]
  generate_trees: true
```

### Large Dataset Configuration

```yaml
# Optimized for large alignments and many taxa
analysis:
  type: ml+bayesian
  model: GTR
  gamma: true
  site_specific: false  # Disable for speed

computational:
  threads: auto
  async_constraints: true
  max_async_workers: 16
  constraint_timeout: 3600  # Longer timeout
  memory_limit: "32G"

constraints:
  strategy: fast
  max_constraints: 50  # Limit constraints
  skip_short_branches: true

bayesian:
  generations: 1000000  # Fewer generations
  nruns: 2
  sample_frequency: 2000  # Less frequent sampling

visualization:
  format: static  # Skip interactive for speed
  dpi: 150  # Lower resolution
```

## Advanced Features

### Custom PAUP Commands

```yaml
analysis:
  type: ml
  paup_block: |
    set criterion=likelihood;
    lscores 1/ userbrlen sitelike;
    savetrees file=custom_trees.tre brlens;
```

### Custom Constraints

```yaml
constraints:
  strategy: custom
  custom_constraints:
    - name: "clade_A"
      taxa: ["species1", "species2", "species3"]
      monophyletic: true
    - name: "clade_B"  
      taxa: ["species4", "species5"]
      monophyletic: false
```

### MPI Configuration

```yaml
bayesian:
  use_mpi: true
  mpi_processes: 16
  
computational:
  # MPI overrides local thread settings
  threads: 1  # Per MPI process
```

### Environment Variable Override

Configuration values can be overridden with environment variables:

```bash
# Override thread count
export PANDECAY_THREADS=16

# Override debug mode
export PANDECAY_DEBUG=true

# Override output directory
export PANDECAY_OUTPUT_DIR=/scratch/results

# Run with environment overrides
python3 panDecay.py --config config.yaml alignment.fas
```

### Configuration Validation Schema

The configuration uses Pydantic models for validation. You can extend validation:

```python
# Custom validator example
from config_models import PanDecayConfig

# Load and validate
config = PanDecayConfig.from_yaml("config.yaml")

# Check specific constraints
if config.computational.threads == "auto":
    import multiprocessing
    actual_threads = multiprocessing.cpu_count()
    print(f"Auto-detected {actual_threads} threads")
```

### Template Generation

```bash
# Generate configuration templates
python3 panDecay.py --generate-yaml-config example.yaml
python3 panDecay.py --generate-toml-config example.toml

# Generate minimal configuration
python3 panDecay.py --generate-yaml-config minimal.yaml --minimal

# Generate configuration with comments
python3 panDecay.py --generate-yaml-config commented.yaml --comments
```

This guide covers all aspects of panDecay configuration. For specific use cases or troubleshooting, refer to the main README.md.