#!/usr/bin/env python3
"""
Base classes and interfaces for panDecay analysis engines.

This module defines the common interface that all analysis engines must implement,
providing a consistent API for ML, Bayesian, and Parsimony analyses.
"""

from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


@dataclass
class AnalysisResult:
    """
    Standardized result structure for all analysis types.
    
    Attributes:
        clade_id: Unique identifier for the clade
        support_value: Primary support metric (ΔlnL, Bayes factor, decay steps)
        p_value: Statistical significance (AU p-value, posterior probability)
        tree_likelihood: Likelihood of the optimal tree
        constrained_likelihood: Likelihood of the constrained tree
        additional_metrics: Method-specific additional data
        success: Whether the analysis completed successfully
        error_message: Error description if analysis failed
    """
    clade_id: str
    support_value: float
    p_value: Optional[float] = None
    tree_likelihood: Optional[float] = None
    constrained_likelihood: Optional[float] = None
    additional_metrics: Optional[Dict[str, Any]] = None
    success: bool = True
    error_message: Optional[str] = None


class AnalysisEngine(ABC):
    """
    Abstract base class for all phylogenetic analysis engines.
    
    This class defines the common interface that ML, Bayesian, and Parsimony
    analysis engines must implement, ensuring consistent behavior and API.
    """
    
    def __init__(self, temp_dir: Path, external_runner, debug: bool = False):
        """
        Initialize the analysis engine.
        
        Args:
            temp_dir: Temporary directory for analysis files
            external_runner: ExternalToolRunner instance for running external tools
            debug: Enable debug logging and file retention
        """
        self.temp_dir = temp_dir
        self.external_runner = external_runner
        self.debug = debug
        self.logger = logging.getLogger(self.__class__.__name__)
    
    @abstractmethod
    def build_optimal_tree(self, alignment_file: Path, model_setup: str) -> Tuple[Optional[str], Optional[float]]:
        """
        Build the optimal tree using this analysis method.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: Model configuration string
            
        Returns:
            Tuple of (tree_string, likelihood_score) or (None, None) if failed
        """
        pass
    
    @abstractmethod
    def analyze_constraint(self, alignment_file: Path, model_setup: str, 
                          clade_taxa: List[str], clade_id: str) -> AnalysisResult:
        """
        Analyze support for a specific clade by constraint testing.
        
        Args:
            alignment_file: Path to the sequence alignment  
            model_setup: Model configuration string
            clade_taxa: List of taxa that should form the clade
            clade_id: Unique identifier for this clade
            
        Returns:
            AnalysisResult with support metrics
        """
        pass
    
    @abstractmethod
    def get_analysis_type(self) -> str:
        """
        Return the analysis type identifier.
        
        Returns:
            String identifier (e.g., 'ml', 'bayesian', 'parsimony')
        """
        pass
    
    def validate_setup(self) -> bool:
        """
        Validate that the analysis engine is properly configured.
        
        Returns:
            True if setup is valid, False otherwise
        """
        if not self.temp_dir.exists():
            self.logger.error(f"Temporary directory does not exist: {self.temp_dir}")
            return False
            
        if not self.external_runner:
            self.logger.error("External tool runner not provided")
            return False
            
        return True
    
    def cleanup_temp_files(self, file_patterns: List[str]) -> None:
        """
        Clean up temporary files created during analysis.
        
        Args:
            file_patterns: List of file patterns to remove
        """
        if self.debug:
            self.logger.debug("Debug mode: keeping temporary files")
            return
            
        for pattern in file_patterns:
            for file_path in self.temp_dir.glob(pattern):
                try:
                    file_path.unlink()
                    self.logger.debug(f"Removed temporary file: {file_path}")
                except Exception as e:
                    self.logger.warning(f"Failed to remove {file_path}: {e}")


class AnalysisConfig:
    """
    Configuration container for analysis engines.
    
    This class holds all the configuration parameters needed by analysis engines,
    reducing the parameter explosion in constructors.
    """
    
    def __init__(self, 
                 alignment_file: Path,
                 alignment_format: str = "fasta",
                 data_type: str = "dna",
                 model: str = "GTR+G",
                 threads: Union[str, int] = "auto",
                 debug: bool = False,
                 keep_files: bool = False,
                 **kwargs):
        """
        Initialize analysis configuration.
        
        Args:
            alignment_file: Path to input alignment
            alignment_format: Format of alignment file
            data_type: Type of sequence data (dna, protein, discrete)
            model: Evolutionary model specification
            threads: Number of threads to use
            debug: Enable debug mode
            keep_files: Keep temporary files after analysis
            **kwargs: Additional analysis-specific parameters
        """
        self.alignment_file = Path(alignment_file)
        self.alignment_format = alignment_format
        self.data_type = data_type
        self.model = model
        self.threads = threads
        self.debug = debug
        self.keep_files = keep_files or debug
        
        # Store additional parameters
        for key, value in kwargs.items():
            setattr(self, key, value)
    
    def get_model_setup_commands(self) -> str:
        """
        Generate model setup commands for external tools.
        
        Returns:
            String containing model setup commands
        """
        # This will be implemented by specific analysis engines
        raise NotImplementedError("Subclasses must implement get_model_setup_commands")