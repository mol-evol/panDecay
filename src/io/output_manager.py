#!/usr/bin/env python3
"""
Output management for panDecay.

This module handles all output formatting and writing operations including
results tables, support summaries, and file management.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union, IO
import datetime

logger = logging.getLogger(__name__)


class OutputManager:
    """
    Manages output generation including reports, tables, and visualizations.
    
    This class handles all the various output formats that panDecay produces,
    including the main results table, detailed reports, and summary statistics.
    """
    
    def __init__(self, temp_path: Optional[Path] = None, debug: bool = False, output_style: str = "unicode"):
        """
        Initialize the output manager.
        
        Args:
            temp_path: Temporary directory path
            debug: Enable debug mode
            output_style: Output style (unicode, ascii, minimal)
        """
        self.temp_path = temp_path or Path.cwd()
        self.debug = debug
        self.output_style = output_style
        
    def get_box_chars(self) -> Dict[str, str]:
        """
        Return appropriate box drawing characters based on output style.
        
        Returns:
            Dictionary mapping box element names to characters
        """
        if self.output_style == "unicode":
            return {
                'horizontal': '─', 'vertical': '│', 'top_left': '┌', 'top_right': '┐',
                'bottom_left': '└', 'bottom_right': '┘', 'cross': '┼', 'top_tee': '┬',
                'bottom_tee': '┴', 'left_tee': '├', 'right_tee': '┤'
            }
        else:
            return {
                'horizontal': '-', 'vertical': '|', 'top_left': '+', 'top_right': '+',
                'bottom_left': '+', 'bottom_right': '+', 'cross': '+', 'top_tee': '+',
                'bottom_tee': '+', 'left_tee': '+', 'right_tee': '+'
            }
    
    def write_support_table(self, f: IO, decay_indices: Dict[str, Dict[str, Any]], 
                           has_ml: bool, has_bayesian: bool, has_parsimony: bool, 
                           has_posterior: bool, has_bootstrap: bool) -> None:
        """
        Write formatted support values table.
        
        Args:
            f: File handle to write to
            decay_indices: Dictionary of decay analysis results
            has_ml: Whether ML analysis was performed
            has_bayesian: Whether Bayesian analysis was performed
            has_parsimony: Whether parsimony analysis was performed
            has_posterior: Whether posterior probabilities are available
            has_bootstrap: Whether bootstrap analysis was performed
        """
        box = self.get_box_chars()
        
        # Calculate column widths
        clade_width = max(10, max(len(str(clade_id)) for clade_id in decay_indices.keys()) if decay_indices else 10)
        taxa_width = 6
        
        # Determine which columns to show
        ml_cols = has_ml
        bayes_cols = has_bayesian
        pars_cols = has_parsimony
        
        # Build header
        header_parts = [f"{'Clade ID':<{clade_width}}", f"{'Taxa':<{taxa_width}}"]
        
        if ml_cols:
            header_parts.extend([f"{'ΔlnL':<10}", f"{'AU p-val':<10}", f"{'Support':<10}"])
        
        if bayes_cols:
            header_parts.extend([f"{'BD':<8}", f"{'LD/site':<8}", f"{'LD%':<8}"])
        
        if pars_cols:
            header_parts.append(f"{'Decay':<23}")
        
        # Write table header
        header_line = box['vertical'] + (' ' + box['vertical'] + ' ').join(header_parts) + box['vertical']
        separator_line = self._create_separator_line(header_line, box)
        
        f.write(separator_line + '\\n')
        f.write(header_line + '\\n')
        
        # Write sub-header for complex columns
        if ml_cols and bayes_cols:
            sub_header_parts = [' ' * clade_width, ' ' * taxa_width]
            
            if ml_cols:
                sub_header_parts.extend([
                    f"{'':─<10}",
                    f"{'':─<10}",  
                    f"{'':─<10}"
                ])
            
            if bayes_cols:
                sub_header_parts.extend([
                    f"{'(Normalized)':─<24}"
                ])
            
            sub_header_line = box['vertical'] + box['horizontal'].join(sub_header_parts) + box['vertical']
            f.write(sub_header_line + '\\n')
        
        f.write(separator_line + '\\n')
        
        # Write data rows
        for clade_id, data in decay_indices.items():
            row_parts = [f"{clade_id:<{clade_width}}", f"{data.get('taxa_count', 0):<{taxa_width}}"]
            
            if ml_cols:
                delta_lnl = data.get('delta_lnL', 0.0)
                au_pval = data.get('AU_pvalue', 1.0)
                support = self._format_significance(au_pval)
                
                row_parts.extend([
                    f"{delta_lnl:>10.3f}",
                    f"{au_pval:>10.4f}",
                    f"{support:>10}"
                ])
            
            if bayes_cols:
                bd = data.get('bayes_decay', 0.0)
                ld_site = data.get('LD_per_site', 0.0)
                ld_pct = data.get('LD_percent', 0.0)
                
                row_parts.extend([
                    f"{bd:>8.2f}",
                    f"{ld_site:>8.6f}",
                    f"{ld_pct:>8.3f}"
                ])
            
            if pars_cols:
                decay = data.get('parsimony_decay', 0)
                row_parts.append(f"{decay:<23}")
            
            row_line = box['vertical'] + (' ' + box['vertical'] + ' ').join(row_parts) + box['vertical']
            f.write(row_line + '\\n')
        
        f.write(separator_line + '\\n')
    
    def write_dataset_relative_rankings(self, f: IO, decay_indices: Dict[str, Dict[str, Any]], 
                                       has_bayesian: bool, has_ml: bool, has_parsimony: bool) -> None:
        """
        Write dataset-relative support rankings.
        
        Args:
            f: File handle to write to
            decay_indices: Dictionary of decay analysis results
            has_bayesian: Whether Bayesian analysis was performed
            has_ml: Whether ML analysis was performed
            has_parsimony: Whether parsimony analysis was performed
        """
        f.write("\\n### Dataset-Relative Support Rankings\\n\\n")
        f.write("Rankings based on ML likelihood decay values within this dataset. \\n\\n")
        f.write("**Important**: \\n")
        f.write("These rankings show relative support within this dataset only and do not indicate absolute phylogenetic reliability.\\n\\n")
        
        # Sort by support values (using ML delta_lnL as primary metric)
        if has_ml:
            sorted_clades = sorted(
                decay_indices.items(), 
                key=lambda x: x[1].get('delta_lnL', 0.0), 
                reverse=True
            )
        elif has_bayesian:
            sorted_clades = sorted(
                decay_indices.items(),
                key=lambda x: x[1].get('bayes_decay', 0.0),
                reverse=True
            )
        else:
            sorted_clades = sorted(
                decay_indices.items(),
                key=lambda x: x[1].get('parsimony_decay', 0),
                reverse=True
            )
        
        # Write highest support
        f.write("**Highest relative support within dataset:**\\n\\n")
        for i, (clade_id, data) in enumerate(sorted_clades[:5]):
            rank = i + 1
            if has_ml:
                raw_value = data.get('delta_lnL', 0.0)
                percentile = data.get('percentile_rank', 0.0)
                z_score = data.get('z_score', 0.0)
                f.write(f"{rank}. {clade_id}: {percentile:.1f}th percentile, Z-score={z_score:+.2f}, Raw ML={raw_value:.2f}\\n")
            elif has_bayesian:
                raw_value = data.get('bayes_decay', 0.0)
                f.write(f"{rank}. {clade_id}: Raw BD={raw_value:.2f}\\n")
            else:
                raw_value = data.get('parsimony_decay', 0)
                f.write(f"{rank}. {clade_id}: Raw Parsimony={raw_value}\\n")
        
        # Write lowest support
        f.write("\\n**Lowest relative support within dataset:**\\n\\n")
        for i, (clade_id, data) in enumerate(sorted_clades[-3:]):
            rank = len(sorted_clades) - 2 + i
            if has_ml:
                raw_value = data.get('delta_lnL', 0.0)
                percentile = data.get('percentile_rank', 0.0)
                z_score = data.get('z_score', 0.0)
                f.write(f"{rank}. {clade_id}: {percentile:.1f}th percentile, Z-score={z_score:+.2f}, Raw ML={raw_value:.2f}\\n")
            elif has_bayesian:
                raw_value = data.get('bayes_decay', 0.0)
                f.write(f"{rank}. {clade_id}: Raw BD={raw_value:.2f}\\n")
            else:
                raw_value = data.get('parsimony_decay', 0)
                f.write(f"{rank}. {clade_id}: Raw Parsimony={raw_value}\\n")
    
    def write_analysis_summary(self, f: IO, analysis_config: Dict[str, Any], 
                              decay_indices: Dict[str, Dict[str, Any]]) -> None:
        """
        Write analysis summary section.
        
        Args:
            f: File handle to write to
            analysis_config: Analysis configuration information
            decay_indices: Dictionary of decay analysis results
        """
        # Title and separator
        title = "panDecay Branch Support Analysis Results"
        separator = "=" * len(title)
        
        f.write(separator + "\\n")
        f.write(f"{title:^{len(separator)}}\\n")
        f.write(separator + "\\n\\n")
        
        # Analysis summary
        f.write("Analysis Summary\\n")
        f.write("─" * 16 + "\\n")
        
        # ML tree likelihood if available
        if 'ml_likelihood' in analysis_config:
            f.write(f"• ML tree log-likelihood: {analysis_config['ml_likelihood']:.3f}\\n")
        
        # Analysis types
        analysis_types = []
        if analysis_config.get('has_ml'):
            analysis_types.append("Maximum Likelihood")
        if analysis_config.get('has_bayesian'):
            analysis_types.append("Bayesian")
        if analysis_config.get('has_parsimony'):
            analysis_types.append("Parsimony")
        
        f.write(f"• Analysis types: {' + '.join(analysis_types)}\\n")
        f.write(f"• Total clades analyzed: {len(decay_indices)}\\n\\n")
    
    def write_clade_details(self, f: IO, decay_indices: Dict[str, Dict[str, Any]], 
                           has_ml: bool, has_bayesian: bool, has_parsimony: bool) -> None:
        """
        Write detailed clade information.
        
        Args:
            f: File handle to write to
            decay_indices: Dictionary of decay analysis results  
            has_ml: Whether ML analysis was performed
            has_bayesian: Whether Bayesian analysis was performed
            has_parsimony: Whether parsimony analysis was performed
        """
        f.write("Clade Details\\n")
        f.write("─" * 13 + "\\n")
        
        for clade_id, data in decay_indices.items():
            # Determine support level
            support_methods = []
            
            if has_ml and data.get('AU_pvalue', 1.0) < 0.05:
                support_methods.append("ML")
            
            if has_bayesian and data.get('bayes_decay', 0.0) > 0:
                support_methods.append("Bayes")
                
            if has_parsimony and data.get('parsimony_decay', 0) > 0:
                support_methods.append("Parsimony")
            
            if support_methods:
                support_desc = f"Strong support: {'/'.join(support_methods)}"
            else:
                support_desc = "Weak support"
            
            f.write(f"→ {clade_id} ({support_desc})\\n")
            
            # Write taxa list
            taxa_list = data.get('taxa_names', [])
            if len(taxa_list) <= 4:
                f.write(", ".join(taxa_list) + "\\n")
            else:
                # Write first few taxa, then wrap
                line = ""
                for i, taxon in enumerate(taxa_list):
                    if i == 0:
                        line = taxon
                    elif len(line) + len(taxon) + 2 <= 72:  # 72 char wrap
                        line += f", {taxon}"
                    else:
                        f.write(line + ",\\n")
                        line = f"        {taxon}"  # Indent continuation
                
                if line:
                    f.write(line + "\\n")
            
            f.write("\\n")
    
    def _create_separator_line(self, header_line: str, box: Dict[str, str]) -> str:
        """
        Create a separator line matching the header structure.
        
        Args:
            header_line: The header line to match
            box: Box drawing characters
            
        Returns:
            Separator line string
        """
        # Replace characters to create separator
        separator = ""
        for char in header_line:
            if char == box['vertical']:
                separator += box['cross']
            elif char == ' ':
                separator += box['horizontal']
            else:
                separator += box['horizontal']
        
        # Fix the corners
        separator = box['top_left'] + separator[1:-1] + box['top_right']
        return separator
    
    def _format_significance(self, p_value: float) -> str:
        """
        Format statistical significance indicators.
        
        Args:
            p_value: P-value to format
            
        Returns:
            Significance indicator string
        """
        if p_value < 0.001:
            return "***"
        elif p_value < 0.01:
            return "**"
        elif p_value < 0.05:
            return "*"
        else:
            return "ns"
    
    def create_results_summary(self, decay_indices: Dict[str, Dict[str, Any]], 
                              analysis_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Create a summary of analysis results.
        
        Args:
            decay_indices: Dictionary of decay analysis results
            analysis_config: Analysis configuration
            
        Returns:
            Dictionary containing summary statistics
        """
        summary = {
            'total_clades': len(decay_indices),
            'timestamp': datetime.datetime.now().isoformat()
        }
        
        # ML summary statistics
        if analysis_config.get('has_ml'):
            ml_values = [data.get('delta_lnL', 0.0) for data in decay_indices.values()]
            au_significant = sum(1 for data in decay_indices.values() 
                               if data.get('AU_pvalue', 1.0) < 0.05)
            
            summary.update({
                'ml_tree_likelihood': analysis_config.get('ml_likelihood'),
                'avg_delta_lnl': sum(ml_values) / len(ml_values) if ml_values else 0,
                'min_delta_lnl': min(ml_values) if ml_values else 0,
                'max_delta_lnl': max(ml_values) if ml_values else 0,
                'au_significant_count': au_significant,
                'au_significant_fraction': au_significant / len(decay_indices) if decay_indices else 0
            })
        
        # Bayesian summary statistics
        if analysis_config.get('has_bayesian'):
            bd_values = [data.get('bayes_decay', 0.0) for data in decay_indices.values()]
            
            summary.update({
                'avg_bayes_decay': sum(bd_values) / len(bd_values) if bd_values else 0,
                'min_bayes_decay': min(bd_values) if bd_values else 0,
                'max_bayes_decay': max(bd_values) if bd_values else 0
            })
        
        # Parsimony summary statistics  
        if analysis_config.get('has_parsimony'):
            pars_values = [data.get('parsimony_decay', 0) for data in decay_indices.values()]
            
            summary.update({
                'avg_parsimony_decay': sum(pars_values) / len(pars_values) if pars_values else 0,
                'min_parsimony_decay': min(pars_values) if pars_values else 0,
                'max_parsimony_decay': max(pars_values) if pars_values else 0
            })
        
        return summary