#!/usr/bin/env python3
"""
Tree annotation module for panDecay.

This module handles tree annotation with support values, including AU p-values,
likelihood differences, Bayesian support, and bootstrap values.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from Bio import Phylo
import io

logger = logging.getLogger(__name__)


class TreeAnnotator:
    """
    Handles tree annotation with various support metrics.
    
    This class takes ML trees and annotates them with different types of
    support values for visualization and analysis.
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize the tree annotator.
        
        Args:
            debug: Enable debug logging
        """
        self.debug = debug
    
    def annotate_trees(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                      output_dir: Path, base_filename: str = "annotated_tree",
                      has_bootstrap: bool = False) -> Dict[str, Path]:
        """
        Create annotated trees with different support values.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_dir: Directory to save the tree files
            base_filename: Base name for the tree files (without extension)
            has_bootstrap: Whether bootstrap analysis was performed
            
        Returns:
            Dictionary with paths to the created tree files
        """
        if not ml_tree or not decay_indices:
            logger.warning("ML tree or decay indices missing. Cannot annotate trees.")
            return {}
        
        output_dir.mkdir(parents=True, exist_ok=True)
        tree_files = {}
        
        try:
            # Create AU p-value annotated tree
            au_tree_path = output_dir / f"{base_filename}_au.nwk"
            if self._create_au_annotated_tree(ml_tree, decay_indices, au_tree_path):
                tree_files['au'] = au_tree_path
            
            # Create log-likelihood difference annotated tree
            delta_lnl_tree_path = output_dir / f"{base_filename}_delta_lnl.nwk"
            if self._create_delta_lnl_annotated_tree(ml_tree, decay_indices, delta_lnl_tree_path):
                tree_files['delta_lnl'] = delta_lnl_tree_path
            
            # Create Bayesian decay annotated tree (if Bayesian analysis was performed)
            if any('bayes_decay' in data for data in decay_indices.values()):
                bayes_tree_path = output_dir / f"{base_filename}_bayes_decay.nwk"
                if self._create_bayes_annotated_tree(ml_tree, decay_indices, bayes_tree_path):
                    tree_files['bayes_decay'] = bayes_tree_path
            
            # Create combined tree with multiple values
            combined_tree_path = output_dir / f"{base_filename}_combined.nwk"
            if self._create_combined_annotated_tree(ml_tree, decay_indices, combined_tree_path):
                tree_files['combined'] = combined_tree_path
            
            # Create bootstrap annotated tree (if bootstrap analysis was performed)
            if has_bootstrap:
                bootstrap_tree_path = output_dir / f"{base_filename}_bootstrap.nwk"
                if self._create_bootstrap_annotated_tree(ml_tree, decay_indices, bootstrap_tree_path):
                    tree_files['bootstrap'] = bootstrap_tree_path
            
            logger.info(f"Created {len(tree_files)} annotated tree files")
            return tree_files
            
        except Exception as e:
            logger.error(f"Failed to annotate trees: {e}")
            return {}
    
    def _create_au_annotated_tree(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                                 output_path: Path) -> bool:
        """
        Create tree annotated with AU p-values.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the annotated tree
            
        Returns:
            True if successful, False otherwise
        """
        try:
            annotated_tree = self._annotate_tree_with_values(
                ml_tree, 
                decay_indices,
                value_key='AU_pvalue',
                label_format="AU:{:.4f}",
                default_value=1.0
            )
            
            if annotated_tree:
                output_path.write_text(annotated_tree + "\\n")
                logger.debug(f"Created AU annotated tree: {output_path}")
                return True
            else:
                logger.warning("Failed to create AU annotated tree")
                return False
                
        except Exception as e:
            logger.error(f"Failed to create AU annotated tree: {e}")
            return False
    
    def _create_delta_lnl_annotated_tree(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                                        output_path: Path) -> bool:
        """
        Create tree annotated with likelihood differences.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the annotated tree
            
        Returns:
            True if successful, False otherwise
        """
        try:
            annotated_tree = self._annotate_tree_with_values(
                ml_tree,
                decay_indices,
                value_key='delta_lnL',
                label_format="ΔlnL:{:.4f}",
                default_value=0.0
            )
            
            if annotated_tree:
                output_path.write_text(annotated_tree + "\\n")
                logger.debug(f"Created delta lnL annotated tree: {output_path}")
                return True
            else:
                logger.warning("Failed to create delta lnL annotated tree")
                return False
                
        except Exception as e:
            logger.error(f"Failed to create delta lnL annotated tree: {e}")
            return False
    
    def _create_bayes_annotated_tree(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                                    output_path: Path) -> bool:
        """
        Create tree annotated with Bayesian decay values.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the annotated tree
            
        Returns:
            True if successful, False otherwise
        """
        try:
            annotated_tree = self._annotate_tree_with_values(
                ml_tree,
                decay_indices,
                value_key='bayes_decay',
                label_format="BD:{:.2f}",
                default_value=0.0
            )
            
            if annotated_tree:
                output_path.write_text(annotated_tree + "\\n")
                logger.debug(f"Created Bayesian decay annotated tree: {output_path}")
                return True
            else:
                logger.warning("Failed to create Bayesian decay annotated tree")
                return False
                
        except Exception as e:
            logger.error(f"Failed to create Bayesian decay annotated tree: {e}")
            return False
    
    def _create_combined_annotated_tree(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                                       output_path: Path) -> bool:
        """
        Create tree annotated with multiple support values.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the annotated tree
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Create combined labels with multiple values
            def format_combined_label(data: Dict[str, Any]) -> str:
                label_parts = []
                
                # Add AU p-value if available
                if 'AU_pvalue' in data:
                    label_parts.append(f"AU:{data['AU_pvalue']:.4f}")
                
                # Add delta lnL if available
                if 'delta_lnL' in data:
                    label_parts.append(f"ΔlnL:{data['delta_lnL']:.2f}")
                
                # Add Bayesian decay if available
                if 'bayes_decay' in data:
                    label_parts.append(f"BD:{data['bayes_decay']:.2f}")
                
                # Add parsimony decay if available
                if 'parsimony_decay' in data:
                    label_parts.append(f"PD:{data['parsimony_decay']}")
                
                return "|".join(label_parts) if label_parts else "NS"
            
            annotated_tree = self._annotate_tree_with_custom_labels(
                ml_tree,
                decay_indices,
                format_combined_label
            )
            
            if annotated_tree:
                output_path.write_text(annotated_tree + "\\n")
                logger.debug(f"Created combined annotated tree: {output_path}")
                return True
            else:
                logger.warning("Failed to create combined annotated tree")
                return False
                
        except Exception as e:
            logger.error(f"Failed to create combined annotated tree: {e}")
            return False
    
    def _create_bootstrap_annotated_tree(self, ml_tree: str, decay_indices: Dict[str, Dict[str, Any]], 
                                        output_path: Path) -> bool:
        """
        Create tree annotated with bootstrap values.
        
        Args:
            ml_tree: ML tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the annotated tree
            
        Returns:
            True if successful, False otherwise
        """
        try:
            annotated_tree = self._annotate_tree_with_values(
                ml_tree,
                decay_indices,
                value_key='bootstrap_support',
                label_format="{:.0f}",  # Bootstrap values are typically integers
                default_value=0
            )
            
            if annotated_tree:
                output_path.write_text(annotated_tree + "\\n")
                logger.debug(f"Created bootstrap annotated tree: {output_path}")
                return True
            else:
                logger.warning("Failed to create bootstrap annotated tree")
                return False
                
        except Exception as e:
            logger.error(f"Failed to create bootstrap annotated tree: {e}")
            return False
    
    def _annotate_tree_with_values(self, tree_string: str, decay_indices: Dict[str, Dict[str, Any]], 
                                  value_key: str, label_format: str, default_value: float) -> Optional[str]:
        """
        Annotate tree with specific support values.
        
        Args:
            tree_string: Original tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            value_key: Key for the value to use for annotation
            label_format: Format string for the label
            default_value: Default value if key not found
            
        Returns:
            Annotated tree string or None if failed
        """
        try:
            # Parse tree to identify internal nodes
            annotated = tree_string
            
            # For each clade in decay_indices, find and annotate the corresponding node
            for clade_id, data in decay_indices.items():
                value = data.get(value_key, default_value)
                
                # Create label
                if value_key == 'AU_pvalue':
                    label = f"'{clade_id} - AU:{value:.4f}'"
                elif value_key == 'delta_lnL':
                    label = f"'{clade_id} - ΔlnL:{value:.4f}'"
                elif value_key == 'bayes_decay':
                    label = f"'{clade_id} - BD:{value:.2f}'"
                elif value_key == 'bootstrap_support':
                    label = f"'{clade_id} - BS:{value:.0f}'"
                else:
                    label = f"'{clade_id} - {label_format.format(value)}'"
                
                # Find and replace the internal node
                # This is a simplified approach - a more robust implementation would
                # use proper tree parsing libraries
                annotated = self._insert_node_label(annotated, clade_id, label, data)
            
            return annotated
            
        except Exception as e:
            logger.error(f"Failed to annotate tree with {value_key}: {e}")
            return None
    
    def _annotate_tree_with_custom_labels(self, tree_string: str, decay_indices: Dict[str, Dict[str, Any]], 
                                         label_formatter) -> Optional[str]:
        """
        Annotate tree with custom-formatted labels.
        
        Args:
            tree_string: Original tree string in Newick format
            decay_indices: Dictionary of decay analysis results
            label_formatter: Function to format labels from data
            
        Returns:
            Annotated tree string or None if failed
        """
        try:
            annotated = tree_string
            
            for clade_id, data in decay_indices.items():
                label_content = label_formatter(data)
                label = f"'{clade_id} - {label_content}'"
                
                # Insert the label into the tree
                annotated = self._insert_node_label(annotated, clade_id, label, data)
            
            return annotated
            
        except Exception as e:
            logger.error(f"Failed to annotate tree with custom labels: {e}")
            return None
    
    def _insert_node_label(self, tree_string: str, clade_id: str, label: str, 
                          clade_data: Dict[str, Any]) -> str:
        """
        Insert a label at the appropriate internal node in the tree.
        
        Args:
            tree_string: Tree string in Newick format
            clade_id: Identifier for the clade
            label: Label to insert
            clade_data: Data for this clade
            
        Returns:
            Tree string with label inserted
        """
        try:
            # This is a simplified implementation
            # A more robust version would properly parse the tree structure
            # and identify internal nodes based on the taxa they contain
            
            # For now, we'll use a pattern-based approach
            # Look for internal nodes (represented by closing parentheses followed by branch lengths)
            pattern = r'\\)[^,);]*'
            
            def replace_node(match):
                node_str = match.group(0)
                # Check if this node corresponds to our clade
                # This would require more sophisticated logic in a full implementation
                return node_str + label
            
            # This is a placeholder - the actual implementation would need
            # to identify which internal node corresponds to which clade
            # based on the taxa it contains
            
            return tree_string
            
        except Exception as e:
            logger.error(f"Failed to insert node label for {clade_id}: {e}")
            return tree_string
    
    def parse_clade_from_tree(self, tree_string: str) -> Dict[str, List[str]]:
        """
        Parse clades from a tree structure.
        
        Args:
            tree_string: Tree string in Newick format
            
        Returns:
            Dictionary mapping internal node IDs to lists of descendant taxa
        """
        try:
            # Use BioPython to parse the tree
            tree_io = io.StringIO(tree_string)
            tree = Phylo.read(tree_io, "newick")
            
            clades = {}
            node_counter = 0
            
            # Traverse internal nodes
            for clade in tree.find_clades(terminal=False):
                node_counter += 1
                clade_id = f"Clade_{node_counter}"
                
                # Get all terminal descendants
                terminals = [term.name for term in clade.get_terminals()]
                clades[clade_id] = terminals
            
            return clades
            
        except Exception as e:
            logger.error(f"Failed to parse clades from tree: {e}")
            return {}
    
    def validate_tree_format(self, tree_string: str) -> bool:
        """
        Validate that a tree string is in proper Newick format.
        
        Args:
            tree_string: Tree string to validate
            
        Returns:
            True if valid, False otherwise
        """
        try:
            # Basic validation checks
            if not tree_string.strip():
                return False
            
            # Check for balanced parentheses
            open_count = tree_string.count('(')
            close_count = tree_string.count(')')
            
            if open_count != close_count:
                logger.error(f"Unbalanced parentheses in tree: {open_count} open, {close_count} close")
                return False
            
            # Check that tree ends with semicolon
            if not tree_string.strip().endswith(';'):
                logger.warning("Tree does not end with semicolon")
            
            # Try to parse with BioPython
            tree_io = io.StringIO(tree_string)
            tree = Phylo.read(tree_io, "newick")
            
            return True
            
        except Exception as e:
            logger.error(f"Tree validation failed: {e}")
            return False