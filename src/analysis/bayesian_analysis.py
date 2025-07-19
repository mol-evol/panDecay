#!/usr/bin/env python3
"""
Bayesian analysis engine for panDecay.

This module handles Bayesian phylogenetic analysis using MrBayes, including
marginal likelihood estimation and posterior probability calculations.
"""

import re
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass

from .analysis_base import AnalysisEngine, AnalysisResult, AnalysisConfig

logger = logging.getLogger(__name__)

# Constants from original code
NEXUS_ALIGNMENT_FN = "alignment.nex"
DEFAULT_BAYESIAN_TIMEOUT = 7200
BURNIN_FRACTION = 0.25
STEPPING_STONE_ALPHA = 0.4
STEPPING_STONE_STEPS = 50
DEFAULT_MIN_ESS = 100
DEFAULT_MAX_PSRF = 1.10
DEFAULT_MAX_ASDSF = 0.10


class BayesianAnalysisEngine(AnalysisEngine):
    """
    Bayesian analysis engine using MrBayes.
    
    This engine handles Bayesian phylogenetic analysis including MCMC sampling,
    marginal likelihood estimation, and posterior probability calculations.
    """
    
    def __init__(self, config: AnalysisConfig, temp_dir: Path, external_runner, debug: bool = False):
        """
        Initialize the Bayesian analysis engine.
        
        Args:
            config: Analysis configuration object
            temp_dir: Temporary directory for analysis files
            external_runner: ExternalToolRunner instance for MrBayes
            debug: Enable debug logging and file retention
        """
        super().__init__(temp_dir, external_runner, debug)
        self.config = config
        self.consensus_tree = None
        self.marginal_likelihood = None
        
        # Bayesian-specific configuration
        self.ngen = getattr(config, 'bayes_ngen', 1000000)
        self.burnin = getattr(config, 'bayes_burnin', BURNIN_FRACTION)
        self.nchains = getattr(config, 'bayes_chains', 4)
        self.sample_freq = getattr(config, 'bayes_sample_freq', 1000)
        self.marginal_method = getattr(config, 'marginal_likelihood', 'ss')
        self.use_mpi = getattr(config, 'use_mpi', False)
        self.mpi_processors = getattr(config, 'mpi_processors', None)
        
    def get_analysis_type(self) -> str:
        """Return the analysis type identifier."""
        return "bayesian"
    
    def build_optimal_tree(self, alignment_file: Path, model_setup: str) -> Tuple[Optional[str], Optional[float]]:
        """
        Build the optimal Bayesian consensus tree using MrBayes.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: MrBayes model configuration commands
            
        Returns:
            Tuple of (consensus_tree_string, marginal_likelihood) or (None, None) if failed
        """
        logger.info("Running unconstrained Bayesian analysis...")
        
        try:
            # Generate MrBayes NEXUS file
            nexus_file = self._generate_mrbayes_nexus(model_setup=model_setup)
            if not nexus_file:
                logger.error("Failed to generate MrBayes NEXUS file")
                return None, None
            
            # Run MrBayes analysis
            if self.marginal_method == "harmonic":
                marginal_lnl, con_tree_path = self._run_mrbayes_with_harmonic_mean(nexus_file, "unconstrained")
            else:
                marginal_lnl, con_tree_path = self._run_mrbayes(nexus_file, "unconstrained")
            
            if marginal_lnl is not None and con_tree_path and con_tree_path.exists():
                # Load consensus tree
                self.consensus_tree = con_tree_path.read_text().strip()
                self.marginal_likelihood = marginal_lnl
                
                logger.info(f"Successfully completed Bayesian analysis. Marginal lnL: {marginal_lnl}")
                return self.consensus_tree, self.marginal_likelihood
            else:
                logger.error("Bayesian analysis failed to produce valid results")
                return None, None
                
        except Exception as e:
            logger.error(f"Failed to build Bayesian tree: {e}")
            return None, None
    
    def analyze_constraint(self, alignment_file: Path, model_setup: str, 
                          clade_taxa: List[str], clade_id: str) -> AnalysisResult:
        """
        Analyze support for a specific clade using Bayesian constraint testing.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: MrBayes model configuration commands
            clade_taxa: List of taxa that should form the clade
            clade_id: Unique identifier for this clade
            
        Returns:
            AnalysisResult with Bayesian support metrics
        """
        try:
            # Generate constrained analysis
            constrained_ml = self._run_constrained_bayesian_analysis(clade_taxa, clade_id, model_setup)
            
            if constrained_ml is None:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="Failed to run constrained Bayesian analysis"
                )
            
            # Calculate Bayes factor (marginal likelihood difference)
            if self.marginal_likelihood is not None:
                bayes_decay = self.marginal_likelihood - constrained_ml
                
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=bayes_decay,
                    tree_likelihood=self.marginal_likelihood,
                    constrained_likelihood=constrained_ml,
                    additional_metrics={
                        'bayes_decay': bayes_decay,
                        'unconstrained_ml': self.marginal_likelihood,
                        'constrained_ml': constrained_ml,
                        'bayes_factor': bayes_decay  # Same as bayes_decay for compatibility
                    }
                )
            else:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="Unconstrained marginal likelihood not available"
                )
                
        except Exception as e:
            logger.error(f"Failed to analyze Bayesian constraint for {clade_id}: {e}")
            return AnalysisResult(
                clade_id=clade_id,
                support_value=0.0,
                success=False,
                error_message=str(e)
            )
    
    def _generate_mrbayes_nexus(self, constraint_tree_file: Optional[Path] = None, 
                               clade_taxa: Optional[List[str]] = None, 
                               constraint_id: Optional[str] = None,
                               model_setup: Optional[str] = None) -> Optional[Path]:
        """
        Generate a MrBayes NEXUS file for analysis.
        
        Args:
            constraint_tree_file: Optional constraint tree file
            clade_taxa: Optional list of taxa for constraint definition
            constraint_id: Optional constraint identifier
            model_setup: Optional model setup commands
            
        Returns:
            Path to generated NEXUS file or None if failed
        """
        try:
            # Read alignment file
            if not self.config.alignment_file.exists():
                logger.error(f"Alignment file not found: {self.config.alignment_file}")
                return None
            
            alignment_content = self.config.alignment_file.read_text()
            
            # Build MrBayes block
            mb_commands = []
            
            # Basic settings
            mb_commands.extend([
                f"set autoclose=yes nowarn=yes;",
                f"mcmcp ngen={self.ngen} printfreq={self.sample_freq} samplefreq={self.sample_freq} nchains={self.nchains};",
                f"mcmcp burninfrac={self.burnin};"
            ])
            
            # Model settings based on data type
            if self.config.data_type == "dna":
                if hasattr(self.config, 'bayes_model'):
                    model = self.config.bayes_model.lower()
                    if model == "jc":
                        mb_commands.append("lset nst=1 rates=equal;")
                    elif model == "hky":
                        mb_commands.append("lset nst=2 rates=gamma;")
                    elif model == "gtr":
                        mb_commands.append("lset nst=6 rates=gamma;")
                    else:
                        mb_commands.append("lset nst=6 rates=gamma;")  # Default to GTR
                else:
                    mb_commands.append("lset nst=6 rates=gamma;")  # Default
            elif self.config.data_type == "protein":
                mb_commands.append("prset aamodelpr=mixed;")
            
            # Add constraint if specified
            if constraint_tree_file and constraint_tree_file.exists():
                constraint_tree_content = constraint_tree_file.read_text().strip()
                mb_commands.extend([
                    f"constraint backbone = {constraint_tree_content};",
                    "prset topologypr=constraints(backbone);"
                ])
            elif clade_taxa:
                # Define negative constraint (force non-monophyly)
                formatted_taxa = [self._format_taxon_for_mrbayes(taxon) for taxon in clade_taxa]
                constraint_def = " ".join(formatted_taxa)
                mb_commands.extend([
                    f"constraint negative_monophyly = -{constraint_def};",
                    "prset topologypr=constraints(negative_monophyly);"
                ])
            
            # Marginal likelihood estimation
            if self.marginal_method == "ss":
                mb_commands.extend([
                    f"ss alpha={STEPPING_STONE_ALPHA} nsteps={STEPPING_STONE_STEPS};",
                ])
            elif self.marginal_method == "harmonic":
                mb_commands.append("mcmc;")
            
            # Output filename
            output_prefix = constraint_id if constraint_id else "unconstrained"
            mb_commands.extend([
                f"mcmc filename={output_prefix};",
                "quit;"
            ])
            
            # Create complete NEXUS content
            nexus_content = alignment_content
            if not nexus_content.strip().endswith("end;"):
                nexus_content += "\\n"
            
            nexus_content += "\\nbegin mrbayes;\\n"
            nexus_content += "\\n".join(mb_commands)
            nexus_content += "\\nend;\\n"
            
            # Write NEXUS file
            nexus_filename = f"{output_prefix}.nex" if constraint_id else "bayesian_analysis.nex"
            nexus_path = self.temp_dir / nexus_filename
            nexus_path.write_text(nexus_content)
            
            if self.debug:
                logger.debug(f"Generated MrBayes NEXUS file: {nexus_path}")
                logger.debug(f"MrBayes commands:\\n{nexus_content}")
            
            return nexus_path
            
        except Exception as e:
            logger.error(f"Failed to generate MrBayes NEXUS file: {e}")
            return None
    
    def _run_mrbayes(self, nexus_file: Path, output_prefix: str) -> Tuple[Optional[float], Optional[Path]]:
        """
        Run MrBayes analysis with stepping-stone marginal likelihood estimation.
        
        Args:
            nexus_file: Path to MrBayes NEXUS file
            output_prefix: Prefix for output files
            
        Returns:
            Tuple of (marginal_likelihood, consensus_tree_path) or (None, None) if failed
        """
        try:
            logger.info(f"Running MrBayes analysis: {output_prefix}")
            
            # Prepare MrBayes command
            if self.use_mpi and self.mpi_processors:
                cmd_result = self.external_runner.run_mrbayes_mpi(
                    str(nexus_file.name), 
                    f"mrbayes_{output_prefix}.log",
                    processors=self.mpi_processors,
                    timeout_sec=DEFAULT_BAYESIAN_TIMEOUT
                )
            else:
                cmd_result = self.external_runner.run_mrbayes(
                    str(nexus_file.name),
                    f"mrbayes_{output_prefix}.log", 
                    timeout_sec=DEFAULT_BAYESIAN_TIMEOUT
                )
            
            # Parse marginal likelihood from output
            marginal_lnl = self._parse_marginal_likelihood(f"mrbayes_{output_prefix}.log")
            
            # Find consensus tree file
            con_tree_path = self.temp_dir / f"{output_prefix}.con.tre"
            if not con_tree_path.exists():
                # Try alternative naming patterns
                for pattern in [f"{output_prefix}*.con.tre", "*.con.tre"]:
                    matches = list(self.temp_dir.glob(pattern))
                    if matches:
                        con_tree_path = matches[0]
                        break
            
            if marginal_lnl is not None and con_tree_path.exists():
                logger.info(f"MrBayes analysis completed. Marginal lnL: {marginal_lnl}")
                return marginal_lnl, con_tree_path
            else:
                logger.warning(f"MrBayes analysis incomplete for {output_prefix}")
                return None, None
                
        except Exception as e:
            logger.error(f"MrBayes analysis failed for {output_prefix}: {e}")
            return None, None
    
    def _run_mrbayes_with_harmonic_mean(self, nexus_file: Path, output_prefix: str) -> Tuple[Optional[float], Optional[Path]]:
        """
        Run MrBayes analysis with harmonic mean marginal likelihood estimation.
        
        Args:
            nexus_file: Path to MrBayes NEXUS file
            output_prefix: Prefix for output files
            
        Returns:
            Tuple of (marginal_likelihood, consensus_tree_path) or (None, None) if failed
        """
        # Similar to _run_mrbayes but with harmonic mean specific handling
        return self._run_mrbayes(nexus_file, output_prefix)
    
    def _run_constrained_bayesian_analysis(self, clade_taxa: List[str], clade_id: str, model_setup: str) -> Optional[float]:
        """
        Run constrained Bayesian analysis for a specific clade.
        
        Args:
            clade_taxa: List of taxa that should NOT form a monophyletic group
            clade_id: Identifier for this constraint
            model_setup: Model setup commands
            
        Returns:
            Marginal likelihood of constrained analysis or None if failed
        """
        try:
            # Generate constrained NEXUS file
            nexus_file = self._generate_mrbayes_nexus(
                clade_taxa=clade_taxa, 
                constraint_id=f"constraint_{clade_id}",
                model_setup=model_setup
            )
            
            if not nexus_file:
                logger.error(f"Failed to generate constrained NEXUS for {clade_id}")
                return None
            
            # Run constrained analysis
            if self.marginal_method == "harmonic":
                marginal_lnl, _ = self._run_mrbayes_with_harmonic_mean(nexus_file, f"constraint_{clade_id}")
            else:
                marginal_lnl, _ = self._run_mrbayes(nexus_file, f"constraint_{clade_id}")
            
            if marginal_lnl is not None:
                logger.debug(f"Constrained analysis for {clade_id} completed. Marginal lnL: {marginal_lnl}")
                return marginal_lnl
            else:
                logger.warning(f"Constrained analysis failed for {clade_id}")
                return None
                
        except Exception as e:
            logger.error(f"Failed constrained Bayesian analysis for {clade_id}: {e}")
            return None
    
    def _parse_marginal_likelihood(self, log_filename: str) -> Optional[float]:
        """
        Parse marginal likelihood from MrBayes log file.
        
        Args:
            log_filename: Name of MrBayes log file
            
        Returns:
            Marginal likelihood value or None if not found
        """
        try:
            log_path = self.temp_dir / log_filename
            if not log_path.exists():
                logger.error(f"MrBayes log file not found: {log_path}")
                return None
            
            log_content = log_path.read_text()
            
            # Parse stepping-stone marginal likelihood
            if self.marginal_method == "ss":
                pattern = r'Stepping-stone\\s+marginal\\s+likelihood\\s*=\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)'
                matches = re.findall(pattern, log_content, re.IGNORECASE)
                if matches:
                    return float(matches[-1])  # Take the last (final) estimate
            
            # Parse harmonic mean marginal likelihood
            elif self.marginal_method == "harmonic":
                pattern = r'Harmonic\\s+mean\\s+marginal\\s+likelihood\\s*=\\s*([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)'
                matches = re.findall(pattern, log_content, re.IGNORECASE)
                if matches:
                    return float(matches[-1])
            
            logger.warning(f"No marginal likelihood found in {log_filename}")
            return None
            
        except Exception as e:
            logger.error(f"Failed to parse marginal likelihood from {log_filename}: {e}")
            return None
    
    def _format_taxon_for_mrbayes(self, taxon_name: str) -> str:
        """
        Format taxon name for MrBayes compatibility.
        
        Args:
            taxon_name: Original taxon name
            
        Returns:
            MrBayes-compatible taxon name
        """
        # MrBayes has similar requirements to PAUP* for taxon names
        formatted = re.sub(r'[^a-zA-Z0-9_]', '_', str(taxon_name))
        
        # Ensure name doesn't start with a digit
        if formatted and formatted[0].isdigit():
            formatted = 'T_' + formatted
            
        return formatted
    
    def validate_setup(self) -> bool:
        """
        Validate that the Bayesian analysis engine is properly configured.
        
        Returns:
            True if setup is valid, False otherwise
        """
        if not super().validate_setup():
            return False
        
        # Check MrBayes availability
        if not hasattr(self.external_runner, 'mrbayes_path'):
            logger.error("MrBayes path not configured in external runner")
            return False
        
        # Validate MCMC parameters
        if self.ngen <= 0:
            logger.error(f"Invalid number of generations: {self.ngen}")
            return False
        
        if not 0 < self.burnin < 1:
            logger.error(f"Invalid burnin fraction: {self.burnin}")
            return False
        
        if self.nchains < 1:
            logger.error(f"Invalid number of chains: {self.nchains}")
            return False
        
        return True