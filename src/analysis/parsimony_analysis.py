#!/usr/bin/env python3
"""
Parsimony analysis engine for panDecay.

This module handles parsimony-based phylogenetic analysis including Bremer support
(traditional decay indices) using PAUP*.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass

from .analysis_base import AnalysisEngine, AnalysisResult, AnalysisConfig

logger = logging.getLogger(__name__)

# Constants from original code
NEXUS_ALIGNMENT_FN = "alignment.nex"
PARSIMONY_TREE_FN = "parsimony_tree.tre"
PARSIMONY_SEARCH_NEX_FN = "parsimony_search.nex"
PARSIMONY_LOG_FN = "paup_parsimony.log"
DEFAULT_PARSIMONY_TIMEOUT = 1800
PARSIMONY_SEARCH_REPS = 10


class ParsimonyAnalysisEngine(AnalysisEngine):
    """
    Parsimony analysis engine using PAUP*.
    
    This engine handles parsimony tree building and Bremer support calculation
    (traditional decay indices).
    """
    
    def __init__(self, config: AnalysisConfig, temp_dir: Path, external_runner, debug: bool = False):
        """
        Initialize the Parsimony analysis engine.
        
        Args:
            config: Analysis configuration object
            temp_dir: Temporary directory for analysis files
            external_runner: ExternalToolRunner instance for PAUP*
            debug: Enable debug logging and file retention
        """
        super().__init__(temp_dir, external_runner, debug)
        self.config = config
        self.parsimony_tree = None
        self.parsimony_score = None
        
    def get_analysis_type(self) -> str:
        """Return the analysis type identifier."""
        return "parsimony"
    
    def build_optimal_tree(self, alignment_file: Path, model_setup: str) -> Tuple[Optional[str], Optional[float]]:
        """
        Build the optimal parsimony tree using PAUP*.
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: PAUP* model configuration commands (not used for parsimony)
            
        Returns:
            Tuple of (tree_string, parsimony_score) or (None, None) if failed
        """
        logger.info("Building parsimony tree...")
        
        try:
            # Prepare parsimony search commands
            script_cmds = [
                f"execute {NEXUS_ALIGNMENT_FN};",
                "set criterion=parsimony;",
            ]
            
            # Add starting tree if available
            if hasattr(self.config, 'starting_tree') and self.config.starting_tree and Path(self.config.starting_tree).exists():
                start_tree_fn_temp = "start_tree.tre"
                shutil.copy(str(self.config.starting_tree), str(self.temp_dir / start_tree_fn_temp))
                script_cmds.extend([
                    f"gettrees file={start_tree_fn_temp};",
                    "pscore 1;",
                    "hsearch start=current;"
                ])
            else:
                if hasattr(self.config, 'starting_tree') and self.config.starting_tree:
                    logger.warning(f"Starting tree file not found: {self.config.starting_tree}. Performing standard search.")
                script_cmds.append(f"hsearch start=stepwise addseq=random nreps={PARSIMONY_SEARCH_REPS};")
            
            # Save tree and get score
            script_cmds.extend([
                f"savetrees file={PARSIMONY_TREE_FN} format=newick brlens=yes replace=yes;",
                "describetrees 1 / plot=phylogram;"
            ])
            
            # Create PAUP* script
            paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
            parsimony_search_cmd_path = self.temp_dir / PARSIMONY_SEARCH_NEX_FN
            parsimony_search_cmd_path.write_text(paup_script_content)
            
            if self.debug:
                logger.debug(f"Parsimony search PAUP* script ({parsimony_search_cmd_path}):ֿֿ\\n{paup_script_content}")
            
            # Run PAUP* analysis
            paup_result = self.external_runner.run_paup_command_file(
                PARSIMONY_SEARCH_NEX_FN, PARSIMONY_LOG_FN, timeout_sec=DEFAULT_PARSIMONY_TIMEOUT
            )
            
            # Parse parsimony score from log
            self.parsimony_score = self._parse_parsimony_score(self.temp_dir / PARSIMONY_LOG_FN)
            
            # Load tree string
            parsimony_tree_path = self.temp_dir / PARSIMONY_TREE_FN
            if parsimony_tree_path.exists():
                self.parsimony_tree = parsimony_tree_path.read_text().strip()
                logger.info(f"Successfully built parsimony tree. Score: {self.parsimony_score}")
                return self.parsimony_tree, self.parsimony_score
            else:
                logger.error(f"Parsimony tree file not found: {parsimony_tree_path}")
                return None, None
                
        except Exception as e:
            logger.error(f"Failed to build parsimony tree: {e}")
            return None, None
    
    def analyze_constraint(self, alignment_file: Path, model_setup: str, 
                          clade_taxa: List[str], clade_id: str) -> AnalysisResult:
        """
        Analyze support for a specific clade using parsimony constraint testing (Bremer support).
        
        Args:
            alignment_file: Path to the sequence alignment
            model_setup: PAUP* model configuration commands (not used for parsimony)
            clade_taxa: List of taxa that should form the clade
            clade_id: Unique identifier for this clade
            
        Returns:
            AnalysisResult with Bremer support metrics
        """
        try:
            # Generate constraint tree and calculate parsimony score
            constrained_score = self._generate_and_score_parsimony_constraint(clade_taxa, clade_id)
            
            if constrained_score is None:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="Failed to generate constraint tree"
                )
            
            # Calculate Bremer support (additional steps required)
            if self.parsimony_score is not None:
                bremer_support = constrained_score - self.parsimony_score
                
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=bremer_support,
                    tree_likelihood=self.parsimony_score,
                    constrained_likelihood=constrained_score,
                    additional_metrics={
                        'bremer_support': bremer_support,
                        'optimal_score': self.parsimony_score,
                        'constrained_score': constrained_score,
                        'additional_steps': bremer_support
                    }
                )
            else:
                return AnalysisResult(
                    clade_id=clade_id,
                    support_value=0.0,
                    success=False,
                    error_message="Optimal parsimony score not available"
                )
                
        except Exception as e:
            logger.error(f"Failed to analyze parsimony constraint for {clade_id}: {e}")
            return AnalysisResult(
                clade_id=clade_id,
                support_value=0.0,
                success=False,
                error_message=str(e)
            )
    
    def _generate_and_score_parsimony_constraint(self, clade_taxa: List[str], clade_id: str) -> Optional[float]:
        """
        Generate a parsimony constraint tree where the specified clade is forced to be non-monophyletic.
        
        Args:
            clade_taxa: List of taxa that should NOT form a monophyletic group
            clade_id: Identifier for this constraint
            
        Returns:
            Parsimony score of the best constrained tree, or None if failed
        """
        try:
            # Format taxa names for PAUP*
            formatted_taxa = [self._format_taxon_for_paup(taxon) for taxon in clade_taxa]
            constraint_def = " ".join(formatted_taxa)
            
            # Generate constraint commands
            constraint_cmds = [
                f"execute {NEXUS_ALIGNMENT_FN};",
                "set criterion=parsimony;",
                f"constraint negative_monophyly (NOT ({constraint_def}));",
                f"hsearch start=stepwise addseq=random nreps={PARSIMONY_SEARCH_REPS} enforce=yes constraints=negative_monophyly;",
                "describetrees 1 / plot=phylogram;"
            ]
            
            # Create and run constraint script
            constraint_script = f"#NEXUS\nbegin paup;\n" + "\n".join(constraint_cmds) + "\nquit;\nend;\n"
            constraint_cmd_path = self.temp_dir / f"parsimony_constraint_{clade_id}.nex"
            constraint_cmd_path.write_text(constraint_script)
            
            if self.debug:
                logger.debug(f"Parsimony constraint script for {clade_id}:\\n{constraint_script}")
            
            # Run PAUP* constraint analysis
            paup_result = self.external_runner.run_paup_command_file(
                f"parsimony_constraint_{clade_id}.nex", 
                f"parsimony_constraint_{clade_id}.log", 
                timeout_sec=DEFAULT_PARSIMONY_TIMEOUT
            )
            
            # Parse parsimony score from log
            constrained_score = self._parse_parsimony_score(
                self.temp_dir / f"parsimony_constraint_{clade_id}.log"
            )
            
            if constrained_score is not None:
                logger.debug(f"Parsimony constraint {clade_id} score: {constrained_score}")
                return constrained_score
            else:
                logger.warning(f"Failed to parse constraint parsimony score for {clade_id}")
                return None
                
        except Exception as e:
            logger.error(f"Failed to generate parsimony constraint tree for {clade_id}: {e}")
            return None
    
    def _parse_parsimony_score(self, log_path: Path) -> Optional[float]:
        """
        Parse parsimony score from PAUP* log file.
        
        Args:
            log_path: Path to PAUP* log file
            
        Returns:
            Parsimony score or None if not found
        """
        try:
            if not log_path.exists():
                logger.error(f"PAUP* log file not found: {log_path}")
                return None
            
            log_content = log_path.read_text()
            
            # Parse parsimony score using regex patterns
            patterns = [
                r'Tree\\s+length\\s*=\\s*([0-9]+(?:\\.[0-9]+)?)',
                r'Length\\s*=\\s*([0-9]+(?:\\.[0-9]+)?)',
                r'Score\\s*=\\s*([0-9]+(?:\\.[0-9]+)?)',
                r'Steps\\s*=\\s*([0-9]+(?:\\.[0-9]+)?)'
            ]
            
            for pattern in patterns:
                matches = re.findall(pattern, log_content, re.IGNORECASE)
                if matches:
                    score = float(matches[-1])  # Take the last (final) score
                    logger.debug(f"Parsed parsimony score: {score}")
                    return score
            
            logger.warning(f"No parsimony score found in {log_path}")
            return None
            
        except Exception as e:
            logger.error(f"Failed to parse parsimony score from {log_path}: {e}")
            return None
    
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
    
    def calculate_bremer_support_matrix(self, clade_list: List[Tuple[str, List[str]]]) -> Dict[str, float]:
        """
        Calculate Bremer support for multiple clades efficiently.
        
        Args:
            clade_list: List of tuples (clade_id, clade_taxa)
            
        Returns:
            Dictionary mapping clade IDs to Bremer support values
        """
        results = {}
        
        for clade_id, clade_taxa in clade_list:
            result = self.analyze_constraint(self.config.alignment_file, "", clade_taxa, clade_id)
            if result.success:
                results[clade_id] = result.support_value
            else:
                results[clade_id] = 0.0
                logger.warning(f"Failed to calculate Bremer support for {clade_id}: {result.error_message}")
        
        return results
    
    def validate_setup(self) -> bool:
        """
        Validate that the Parsimony analysis engine is properly configured.
        
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
        
        # Parsimony analysis works with all data types, so no specific validation needed
        return True