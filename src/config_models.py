#!/usr/bin/env python3
"""
Configuration models for panDecay using Pydantic for validation.

This module defines the structure and validation rules for panDecay
configuration files, supporting both YAML and TOML formats.
"""

from typing import Dict, List, Optional, Union, Any, Literal
from pathlib import Path
from pydantic import BaseModel, Field, validator, ConfigDict
import os


class InputOutputConfig(BaseModel):
    """Input/Output configuration settings."""
    
    alignment_file: Path = Field(..., description="Path to input alignment file")
    alignment_format: Literal["fasta", "phylip", "nexus", "clustal", "stockholm"] = Field(
        default="fasta", description="Format of input alignment"
    )
    data_type: Literal["dna", "protein", "discrete"] = Field(
        default="dna", description="Type of sequence data"
    )
    output_prefix: str = Field(
        default="pan_decay_indices", description="Prefix for output files"
    )
    tree_prefix: str = Field(
        default="annotated_tree", description="Prefix for annotated tree files"
    )
    temp_directory: Optional[Path] = Field(
        default=None, description="Custom temporary directory path"
    )
    keep_files: bool = Field(
        default=False, description="Keep temporary files after analysis"
    )
    debug: bool = Field(
        default=False, description="Enable debug mode with detailed logging"
    )
    
    @validator('alignment_file')
    def validate_alignment_file(cls, v):
        """Validate that alignment file exists."""
        if not Path(v).exists():
            raise ValueError(f"Alignment file not found: {v}")
        return v
    
    @validator('temp_directory')
    def validate_temp_directory(cls, v):
        """Validate temp directory if specified."""
        if v is not None:
            v = Path(v)
            if not v.exists():
                try:
                    v.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    raise ValueError(f"Cannot create temp directory {v}: {e}")
        return v


class AnalysisConfig(BaseModel):
    """Analysis mode configuration settings."""
    
    analysis_types: List[Literal["ml", "bayesian", "parsimony"]] = Field(
        default=["ml"], description="Types of analyses to perform"
    )
    bootstrap: bool = Field(
        default=False, description="Perform bootstrap analysis"
    )
    bootstrap_reps: int = Field(
        default=100, ge=10, le=10000, description="Number of bootstrap replicates"
    )
    site_analysis: bool = Field(
        default=False, description="Perform site-specific analysis"
    )
    normalization: Literal["none", "basic", "full"] = Field(
        default="basic", description="Normalization level: none (no normalization), basic (per-site and relative metrics), full (includes dataset-relative rankings)"
    )
    
    @validator('analysis_types')
    def validate_analysis_types(cls, v):
        """Ensure analysis types are valid and not duplicated."""
        if not v:
            raise ValueError("At least one analysis type must be specified")
        return list(set(v))  # Remove duplicates


class ModelConfig(BaseModel):
    """Evolutionary model configuration settings."""
    
    # DNA models
    dna_model: Literal["JC", "K80", "HKY", "GTR"] = Field(
        default="GTR", description="DNA substitution model"
    )
    
    # Protein models  
    protein_model: Literal["JTT", "WAG", "LG", "Dayhoff"] = Field(
        default="WAG", description="Protein substitution model"
    )
    
    # Discrete models
    discrete_model: Literal["Mk"] = Field(
        default="Mk", description="Discrete morphology model"
    )
    
    # Rate variation
    gamma: bool = Field(
        default=False, description="Add gamma rate heterogeneity"
    )
    gamma_shape: Optional[float] = Field(
        default=None, gt=0, description="Fixed gamma shape parameter"
    )
    invariant: bool = Field(
        default=False, description="Add proportion of invariable sites"
    )
    prop_invar: Optional[float] = Field(
        default=None, ge=0, le=1, description="Fixed proportion of invariable sites"
    )
    
    # Base frequencies
    base_freq: Literal["equal", "estimate", "empirical"] = Field(
        default="estimate", description="Base frequency estimation method"
    )
    
    # NST for DNA
    nst: Optional[Literal[1, 2, 6]] = Field(
        default=None, description="Number of substitution types for DNA"
    )


class ComputationalConfig(BaseModel):
    """Computational settings configuration."""
    
    threads: Union[int, Literal["auto", "all"]] = Field(
        default="auto", description="Number of threads for PAUP*"
    )
    paup_path: str = Field(
        default="paup", description="Path to PAUP* executable"
    )
    starting_tree: Optional[Path] = Field(
        default=None, description="Starting tree file in Newick format"
    )
    paup_block: Optional[Path] = Field(
        default=None, description="Custom PAUP* commands file"
    )
    
    # Timeouts (seconds)
    ml_timeout: int = Field(
        default=3600, ge=60, description="Timeout for ML tree search"
    )
    constraint_timeout: int = Field(
        default=600, ge=60, description="Timeout for constraint analysis"
    )
    site_analysis_timeout: int = Field(
        default=600, ge=60, description="Timeout for site-specific analysis"
    )
    
    # Async constraint processing
    use_async_constraints: bool = Field(
        default=False, description="Enable parallel constraint processing"
    )
    max_async_workers: Optional[int] = Field(
        default=None, ge=1, le=32, description="Maximum number of parallel constraint workers"
    )
    
    @validator('threads')
    def validate_threads(cls, v):
        """Validate thread specification."""
        if isinstance(v, int):
            if v < 1:
                raise ValueError("Thread count must be positive")
            if v > os.cpu_count():
                raise ValueError(f"Thread count exceeds available cores ({os.cpu_count()})")
        return v
    
    @validator('starting_tree', 'paup_block')
    def validate_file_paths(cls, v):
        """Validate optional file paths exist."""
        if v is not None and not Path(v).exists():
            raise ValueError(f"File not found: {v}")
        return v


class BayesianConfig(BaseModel):
    """Bayesian analysis configuration settings."""
    
    software: Literal["mrbayes"] = Field(
        default="mrbayes", description="Bayesian software to use"
    )
    mrbayes_path: str = Field(
        default="mb", description="Path to MrBayes executable"
    )
    
    # MCMC settings
    ngen: int = Field(
        default=1000000, ge=10000, description="Number of MCMC generations"
    )
    burnin: float = Field(
        default=0.25, ge=0, le=0.9, description="Burnin fraction"
    )
    chains: int = Field(
        default=4, ge=1, le=32, description="Number of MCMC chains"
    )
    sample_freq: int = Field(
        default=1000, ge=1, description="Sample frequency"
    )
    print_freq: int = Field(
        default=10000, ge=1, description="Print frequency"
    )
    diag_freq: int = Field(
        default=50000, ge=1, description="Diagnostics frequency"
    )
    
    # Marginal likelihood
    marginal_likelihood: Literal["ss", "ps", "hm"] = Field(
        default="ss", description="Marginal likelihood estimation method"
    )
    ss_alpha: float = Field(
        default=0.4, gt=0, lt=1, description="Stepping-stone alpha parameter"
    )
    ss_nsteps: int = Field(
        default=50, ge=10, description="Number of stepping-stone steps"
    )
    
    # Convergence checking
    check_convergence: bool = Field(
        default=True, description="Check MCMC convergence diagnostics"
    )
    min_ess: float = Field(
        default=200, ge=50, description="Minimum ESS threshold"
    )
    max_psrf: float = Field(
        default=1.01, gt=1, description="Maximum PSRF threshold"
    )
    max_asdsf: float = Field(
        default=0.01, gt=0, description="Maximum ASDSF threshold"
    )
    convergence_strict: bool = Field(
        default=False, description="Fail analysis if convergence criteria not met"
    )
    
    # Parallel processing
    use_mpi: bool = Field(
        default=False, description="Use MPI version of MrBayes"
    )
    mpi_processors: Optional[int] = Field(
        default=None, ge=1, description="Number of MPI processors"
    )
    mpirun_path: str = Field(
        default="mpirun", description="Path to mpirun executable"
    )
    
    # BEAGLE acceleration
    use_beagle: bool = Field(
        default=False, description="Use BEAGLE library for acceleration"
    )
    beagle_device: Literal["cpu", "gpu", "auto"] = Field(
        default="auto", description="BEAGLE device"
    )
    beagle_precision: Literal["single", "double"] = Field(
        default="double", description="BEAGLE precision"
    )
    beagle_scaling: Literal["none", "dynamic", "always"] = Field(
        default="dynamic", description="BEAGLE scaling"
    )


class VisualizationConfig(BaseModel):
    """Visualization configuration settings."""
    
    enable: bool = Field(
        default=False, description="Generate visualizations"
    )
    format: Literal["static", "interactive", "both", "png", "pdf", "svg", "html"] = Field(
        default="both", description="Visualization format"
    )
    
    # Static plot settings
    static: Dict[str, Any] = Field(
        default={
            "dpi": 300,
            "formats": ["png", "pdf"],
            "style": "publication",
            "figsize": [10, 8],
            "font_size": 12,
            "color_palette": "viridis"
        },
        description="Static plot configuration"
    )
    
    # Interactive plot settings
    interactive: Dict[str, Any] = Field(
        default={
            "theme": "plotly_white",
            "export_html": True,
            "include_controls": True,
            "width": 1000,
            "height": 800,
            "color_scale": "viridis"
        },
        description="Interactive plot configuration"
    )
    
    annotation: Literal["au", "lnl", "decay"] = Field(
        default="decay", description="Annotation type for visualization"
    )


class ConstraintConfig(BaseModel):
    """Constraint analysis configuration settings."""
    
    mode: Literal["all", "specific", "exclude"] = Field(
        default="all", description="Constraint testing mode"
    )
    test_branches: Optional[List[str]] = Field(
        default=None, description="Specific branches to test"
    )
    constraint_file: Optional[Path] = Field(
        default=None, description="File containing constraint definitions"
    )
    custom_constraints: Optional[Dict[str, List[str]]] = Field(
        default=None, description="Custom constraint definitions"
    )
    
    @validator('constraint_file')
    def validate_constraint_file(cls, v):
        """Validate constraint file exists."""
        if v is not None and not Path(v).exists():
            raise ValueError(f"Constraint file not found: {v}")
        return v


class PanDecayConfig(BaseModel):
    """Main panDecay configuration model."""
    
    model_config = ConfigDict(
        extra='forbid',  # Don't allow extra fields
        validate_assignment=True,  # Validate on assignment
        use_enum_values=True  # Use enum values in serialization
    )
    
    # Configuration sections
    input_output: InputOutputConfig
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    model: ModelConfig = Field(default_factory=ModelConfig)
    computational: ComputationalConfig = Field(default_factory=ComputationalConfig)
    bayesian: BayesianConfig = Field(default_factory=BayesianConfig)
    visualization: VisualizationConfig = Field(default_factory=VisualizationConfig)
    constraints: ConstraintConfig = Field(default_factory=ConstraintConfig)
    
    @validator('analysis')
    def validate_bayesian_requirements(cls, v, values):
        """Ensure Bayesian config is complete if Bayesian analysis requested."""
        if 'bayesian' in v.analysis_types:
            # Additional validation can be added here
            pass
        return v
    
    def get_model_for_data_type(self, data_type: str) -> str:
        """Get the appropriate model for the specified data type."""
        if data_type == "dna":
            return self.model.dna_model
        elif data_type == "protein":
            return self.model.protein_model
        elif data_type == "discrete":
            return self.model.discrete_model
        else:
            raise ValueError(f"Unknown data type: {data_type}")
    
    def get_thread_count(self) -> int:
        """Get the actual thread count to use."""
        if self.computational.threads == "auto":
            return max(1, os.cpu_count() - 2)
        elif self.computational.threads == "all":
            return os.cpu_count()
        else:
            return self.computational.threads
    
    def to_legacy_args(self) -> Dict[str, Any]:
        """Convert configuration to legacy argument format for backward compatibility."""
        args = {}
        
        # Map configuration to legacy argument names
        args['alignment'] = str(self.input_output.alignment_file)
        args['format'] = self.input_output.alignment_format
        args['data_type'] = self.input_output.data_type
        args['output'] = self.input_output.output_prefix
        args['tree'] = self.input_output.tree_prefix
        args['keep_files'] = self.input_output.keep_files
        args['debug'] = self.input_output.debug
        
        # Analysis settings
        args['analysis'] = '+'.join(self.analysis.analysis_types)
        args['bootstrap'] = self.analysis.bootstrap
        args['bootstrap_reps'] = self.analysis.bootstrap_reps
        args['site_analysis'] = self.analysis.site_analysis
        args['normalization'] = self.analysis.normalization
        
        # Model settings
        args['model'] = self.get_model_for_data_type(self.input_output.data_type)
        args['gamma'] = self.model.gamma
        args['invariant'] = self.model.invariant
        
        # Computational settings
        args['threads'] = self.get_thread_count()
        args['paup'] = self.computational.paup_path
        args['async_constraints'] = self.computational.use_async_constraints
        args['max_async_workers'] = self.computational.max_async_workers
        args['constraint_timeout'] = self.computational.constraint_timeout
        
        # Visualization
        args['visualize'] = self.visualization.enable
        args['viz_format'] = self.visualization.format
        
        return args