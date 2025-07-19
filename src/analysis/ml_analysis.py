#!/usr/bin/env python3
"""
Maximum Likelihood analysis engine for panDecay.

This module handles ML tree building, constraint testing, and AU (Approximately Unbiased)
statistical tests using PAUP*.
"""

import re
import shutil
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass

from .analysis_base import AnalysisEngine, AnalysisResult, AnalysisConfig

logger = logging.getLogger(__name__)

# Constants from original code
NEXUS_ALIGNMENT_FN = "alignment.nex"
ML_TREE_FN = "ml_tree.tre"
ML_SCORE_FN = "ml_score.txt"
ML_SEARCH_NEX_FN = "ml_search.nex"
ML_LOG_FN = "paup_ml.log"
AU_TEST_NEX_FN = "au_test.nex"
AU_LOG_FN = "paup_au.log"
DEFAULT_ML_TIMEOUT = 3600
DEFAULT_CONSTRAINT_TIMEOUT = 1800
PARSIMONY_SEARCH_REPS = 10


class MLAnalysisEngine(AnalysisEngine):
    """
    Maximum Likelihood analysis engine using PAUP*.
    
    This engine handles ML tree building, constraint testing, and AU tests
    for calculating phylogenetic decay indices.
    """
    
    def __init__(self, config: AnalysisConfig, temp_dir: Path, external_runner, debug: bool = False):
        """
        Initialize the ML analysis engine.
        
        Args:
            config: Analysis configuration object
            temp_dir: Temporary directory for analysis files
            external_runner: ExternalToolRunner instance for PAUP*
            debug: Enable debug logging and file retention
        """
        super().__init__(temp_dir, external_runner, debug)
        self.config = config
        self.ml_tree = None
        self.ml_likelihood = None
        
    def get_analysis_type(self) -> str:
        """Return the analysis type identifier."""
        return "ml"
    
    def build_optimal_tree(self, alignment_file: Path, model_setup: str) -> Tuple[Optional[str], Optional[float]]:
        """
        Build the optimal ML tree using PAUP*.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: PAUP* model configuration commands
            
        Returns:
            Tuple of (tree_string, likelihood_score) or (None, None) if failed
        """
        logger.info("Building maximum likelihood tree...")
        
        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", model_setup]
        
        # Add search commands based on configuration
        if hasattr(self.config, 'user_paup_block') and self.config.user_paup_block:
            # User-provided PAUP block
            script_cmds.append(self.config.user_paup_block)
            
            # Add defensive commands if not present in user block
            block_lower = self.config.user_paup_block.lower()
            if "savetrees" not in block_lower:
                script_cmds.append(f"savetrees file={ML_TREE_FN} format=newick brlens=yes replace=yes;")
            if "lscores" not in block_lower and "lscore" not in block_lower:
                script_cmds.append(f"lscores 1 / scorefile={ML_SCORE_FN} replace=yes;")
        else:
            # Standard model processing
            if hasattr(self.config, 'starting_tree') and self.config.starting_tree and Path(self.config.starting_tree).exists():
                start_tree_fn_temp = "start_tree.tre"
                shutil.copy(str(self.config.starting_tree), str(self.temp_dir / start_tree_fn_temp))
                script_cmds.extend([
                    f"gettrees file={start_tree_fn_temp};",
                    "lscores 1 / userbrlen=yes;",
                    "hsearch start=current;"
                ])
            else:
                if hasattr(self.config, 'starting_tree') and self.config.starting_tree:
                    logger.warning(f"Starting tree file not found: {self.config.starting_tree}. Performing standard search.")
                script_cmds.append(f"hsearch start=stepwise addseq=random nreps={PARSIMONY_SEARCH_REPS};")
            
            script_cmds.extend([
                f"savetrees file={ML_TREE_FN} format=newick brlens=yes replace=yes;",
                f"lscores 1 / scorefile={ML_SCORE_FN} replace=yes;"
            ])
        
        # Create PAUP* script
        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        ml_search_cmd_path = self.temp_dir / ML_SEARCH_NEX_FN
        ml_search_cmd_path.write_text(paup_script_content)
        
        if self.debug:
            logger.debug(f"ML search PAUP* script ({ml_search_cmd_path}):ֿֿ\\n{paup_script_content}")
        
        try:
            # Run PAUP* analysis
            paup_result = self.external_runner.run_paup_command_file(
                ML_SEARCH_NEX_FN, ML_LOG_FN, timeout_sec=DEFAULT_ML_TIMEOUT
            )
            
            # Parse likelihood from score file
            self.ml_likelihood = self.external_runner.parse_likelihood_from_score_file(
                self.temp_dir / ML_SCORE_FN
            )
            logger.debug(f"ML likelihood from score file: {self.ml_likelihood}")
            
            # Fallback to log parsing if score file parsing failed
            if self.ml_likelihood is None and paup_result.stdout:
                logger.info(f"Fallback: Parsing ML likelihood from PAUP* log {ML_LOG_FN}")
                patterns = [r'-ln\\s*L\\s*=\\s*([0-9.]+)', r'likelihood\\s*=\\s*([0-9.]+)', r'score\\s*=\\s*([0-9.]+)']
                for pattern in patterns:
                    matches = re.findall(pattern, paup_result.stdout, re.IGNORECASE)
                    if matches:
                        self.ml_likelihood = float(matches[-1])
                        logger.debug(f"ML likelihood from log: {self.ml_likelihood}")
                        break
            
            # Load tree string
            ml_tree_path = self.temp_dir / ML_TREE_FN
            if ml_tree_path.exists():
                self.ml_tree = ml_tree_path.read_text().strip()
                logger.info(f"Successfully built ML tree. Log-likelihood: {self.ml_likelihood}")
                return self.ml_tree, self.ml_likelihood
            else:
                logger.error(f"ML tree file not found: {ml_tree_path}")
                return None, None
                
        except Exception as e:
            logger.error(f"Failed to build ML tree: {e}")
            return None, None
    
    def analyze_constraint(self, alignment_file: Path, model_setup: str, 
                          clade_taxa: List[str], clade_id: str) -> AnalysisResult:
        """
        Analyze support for a specific clade using constraint testing.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: PAUP* model configuration commands
            clade_taxa: List of taxa that should form the clade
            clade_id: Unique identifier for this clade
            
        Returns:
            AnalysisResult with ML support metrics
        """
        try:
            # Generate constraint tree
            constrained_likelihood = self._generate_and_score_constraint_tree(
                model_setup, clade_taxa, clade_id
            )
            
            if constrained_likelihood is None:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="Failed to generate constraint tree"
                )
            
            # Calculate likelihood difference
            if self.ml_likelihood is not None:
                delta_lnl = constrained_likelihood - self.ml_likelihood
                
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=delta_lnl,
                    tree_likelihood=self.ml_likelihood,
                    constrained_likelihood=constrained_likelihood,
                    additional_metrics={
                        'delta_lnl': delta_lnl,
                        'ml_likelihood': self.ml_likelihood,
                        'constrained_likelihood': constrained_likelihood
                    }
                )
            else:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="ML likelihood not available"
                )
                
        except Exception as e:
            logger.error(f"Failed to analyze constraint for {clade_id}: {e}")
            return AnalysisResult(
                clade_id=clade_id,
                support_value=0.0,
                success=False,
                error_message=str(e)
            )
    
    def _generate_and_score_constraint_tree(self, model_setup: str, clade_taxa: List[str], clade_id: str) -> Optional[float]:
        """
        Generate a constraint tree where the specified clade is forced to be non-monophyletic.
        
        Args:
            model_setup: PAUP* model configuration commands
            clade_taxa: List of taxa that should NOT form a monophyletic group
            clade_id: Identifier for this constraint
            
        Returns:
            Likelihood of the best constrained tree, or None if failed
        """
        try:
            # Format taxa names for PAUP*
            formatted_taxa = [self._format_taxon_for_paup(taxon) for taxon in clade_taxa]
            constraint_def = " ".join(formatted_taxa)
            
            # Generate constraint commands
            constraint_cmds = [
                f"execute {NEXUS_ALIGNMENT_FN};",
                model_setup,
                f"constraints clade_constraint (MONOPHYLY) = (({constraint_def}));",
                f"hsearch start=stepwise addseq=random nreps={PARSIMONY_SEARCH_REPS} enforce=yes converse=yes constraints=clade_constraint;",
                f"lscores 1 / scorefile=constraint_{clade_id}_score.txt replace=yes;"
            ]
            
            # Create and run constraint script
            constraint_script = f"#NEXUS\nbegin paup;\n" + "\n".join(constraint_cmds) + "\nquit;\nend;\n"
            constraint_cmd_path = self.temp_dir / f"constraint_{clade_id}.nex"
            constraint_cmd_path.write_text(constraint_script)
            
            if self.debug:
                logger.debug(f"Constraint script for {clade_id}:\\n{constraint_script}")
            
            # Run PAUP* constraint analysis
            paup_result = self.external_runner.run_paup_command_file(
                f"constraint_{clade_id}.nex", 
                f"constraint_{clade_id}.log", 
                timeout_sec=DEFAULT_CONSTRAINT_TIMEOUT
            )
            
            # Parse likelihood from score file
            score_file_path = self.temp_dir / f"constraint_{clade_id}_score.txt"
            constrained_likelihood = self.external_runner.parse_likelihood_from_score_file(score_file_path)
            
            if constrained_likelihood is not None:
                logger.debug(f"Constraint {clade_id} likelihood: {constrained_likelihood}")
                return constrained_likelihood
            else:
                logger.warning(f"Failed to parse constraint likelihood for {clade_id}")
                return None
                
        except Exception as e:
            logger.error(f"Failed to generate constraint tree for {clade_id}: {e}")
            return None
    
    def run_au_test(self, tree_files: List[str], results_storage: Dict[str, Any]) -> Dict[int, Dict[str, Any]]:
        """
        Run AU (Approximately Unbiased) test on a set of trees.
        
        Args:
            tree_files: List of tree filenames (relative to temp_dir)
            results_storage: Dictionary to store additional results
            
        Returns:
            Dictionary mapping tree indices to AU test results
        """
        logger.info(f"Running AU test on {len(tree_files)} trees...")
        
        try:
            # Prepare AU test commands
            tree_list = " ".join(tree_files)
            au_cmds = [
                f"execute {NEXUS_ALIGNMENT_FN};",
                self._get_model_setup_commands(),
                f"gettrees file={tree_list} mode=3;",  # mode=3 allows multiple tree files
                f"lscores all / au=yes userbrlen=no scorefile=au_scores.txt replace=yes;"
            ]
            
            # Create AU test script
            au_script = f"#NEXUS\nbegin paup;\n" + "\n".join(au_cmds) + "\nquit;\nend;\n"
            au_test_path = self.temp_dir / AU_TEST_NEX_FN
            au_test_path.write_text(au_script)
            
            if self.debug:
                logger.debug(f"AU test script:\\n{au_script}")
            
            # Run AU test
            paup_result = self.external_runner.run_paup_command_file(
                AU_TEST_NEX_FN, AU_LOG_FN, timeout_sec=DEFAULT_ML_TIMEOUT
            )
            
            # Parse AU test results
            au_results = self._parse_au_results(self.temp_dir / AU_LOG_FN)
            
            if au_results:
                logger.info(f"AU test completed successfully for {len(au_results)} trees")
                return au_results
            else:
                logger.warning("AU test parsing returned no results")
                return {}
                
        except Exception as e:
            logger.error(f"AU test failed: {e}")
            return {}
    
    def _parse_au_results(self, au_log_path: Path) -> Dict[int, Dict[str, Any]]:
        """
        Parse AU test results from PAUP* log file.
        
        Args:
            au_log_path: Path to AU test log file
            
        Returns:
            Dictionary mapping tree indices to AU test results
        """
        results = {}
        
        try:
            if not au_log_path.exists():
                logger.error(f"AU log file not found: {au_log_path}")
                return results
            
            log_content = au_log_path.read_text()
            
            # Parse AU test results using regex patterns
            # This is a simplified version - the full implementation would need
            # the complete parsing logic from the original _parse_au_results method
            au_pattern = r'Tree\\s+(\\d+).*?AU\\s+([0-9.]+)'
            matches = re.findall(au_pattern, log_content, re.IGNORECASE | re.DOTALL)
            
            for match in matches:
                tree_idx = int(match[0])
                au_pvalue = float(match[1])
                
                results[tree_idx] = {
                    'AU_pvalue': au_pvalue,
                    'significant': au_pvalue < 0.05
                }
            
            logger.debug(f"Parsed AU results for {len(results)} trees")
            return results
            
        except Exception as e:
            logger.error(f"Failed to parse AU results from {au_log_path}: {e}")
            return results
    
    def _format_taxon_for_paup(self, taxon_name: str) -> str:
        """
        Format taxon name for PAUP* compatibility.
        
        Args:
            taxon_name: Original taxon name
            
        Returns:
            PAUP*-compatible taxon name
        """
        # Remove problematic characters and ensure valid format
        formatted = re.sub(r'[^a-zA-Z0-9_]', '_', str(taxon_name))
        
        # Ensure name doesn't start with a digit
        if formatted and formatted[0].isdigit():
            formatted = 'T_' + formatted
            
        return formatted
    
    def _get_model_setup_commands(self) -> str:
        """
        Generate PAUP* model setup commands from configuration.
        
        Returns:
            String containing PAUP* model setup commands
        """
        # This would implement the model setup logic from the original
        # _get_paup_model_setup_cmds method
        # For now, return a basic setup
        return f"set criterion=likelihood; lset nst=6 basefreq=estimate rates=gamma;"
    
    def validate_setup(self) -> bool:
        """
        Validate that the ML analysis engine is properly configured.
        
        Returns:
            True if setup is valid, False otherwise
        """
        if not super().validate_setup():
            return False
        
        # Check PAUP* availability through external runner
        if not hasattr(self.external_runner, 'paup_path'):
            logger.error("PAUP* path not configured in external runner")
            return False
        
        # Validate alignment file
        if not self.config.alignment_file.exists():
            logger.error(f"Alignment file not found: {self.config.alignment_file}")
            return False
        
        return True