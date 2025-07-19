#!/usr/bin/env python3
"""
Report generation module for panDecay.

This module handles the creation of detailed markdown reports with
analysis results, interpretation guides, and statistical summaries.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
import datetime

logger = logging.getLogger(__name__)


class ReportGenerator:
    """
    Generates comprehensive markdown reports for panDecay analyses.
    
    This class creates detailed reports that include analysis parameters,
    results interpretation, and statistical summaries.
    """
    
    def __init__(self, version: str = "1.1.0"):
        """
        Initialize the report generator.
        
        Args:
            version: panDecay version string
        """
        self.version = version
    
    def generate_detailed_report(self, output_path: Path, analysis_config: Dict[str, Any], 
                               decay_indices: Dict[str, Dict[str, Any]], 
                               summary_stats: Dict[str, Any]) -> bool:
        """
        Generate a comprehensive markdown report.
        
        Args:
            output_path: Path to save the report
            analysis_config: Analysis configuration and parameters
            decay_indices: Dictionary of decay analysis results
            summary_stats: Summary statistics
            
        Returns:
            True if report generated successfully, False otherwise
        """
        try:
            with open(output_path, 'w') as f:
                self._write_header(f, analysis_config)
                self._write_parameters_section(f, analysis_config)
                self._write_summary_statistics(f, summary_stats, analysis_config)
                
                # Add Bayesian convergence diagnostics if available
                if analysis_config.get('has_bayesian'):
                    self._write_bayesian_diagnostics(f, analysis_config)
                
                self._write_detailed_results(f, decay_indices, analysis_config)
                self._write_interpretation_guide(f, analysis_config)
                
                # Add dataset-relative rankings if normalization was performed
                if analysis_config.get('normalization') in ['basic', 'full']:
                    self._write_dataset_relative_section(f, decay_indices, analysis_config)
            
            logger.info(f"Generated detailed report: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to generate report: {e}")
            return False
    
    def _write_header(self, f, analysis_config: Dict[str, Any]) -> None:
        """Write the report header."""
        # Determine analysis types for title
        analysis_types = []
        if analysis_config.get('has_ml'):
            analysis_types.append("ML")
        if analysis_config.get('has_bayesian'):
            analysis_types.append("Bayesian")
        if analysis_config.get('has_parsimony'):
            analysis_types.append("Parsimony")
        
        title_suffix = " + ".join(analysis_types) if analysis_types else ""
        
        f.write(f"# panDecay {title_suffix} Branch Support Analysis Report (v{self.version})\\n\\n")
        f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")
    
    def _write_parameters_section(self, f, analysis_config: Dict[str, Any]) -> None:
        """Write the analysis parameters section."""
        f.write("## Analysis Parameters\\n\\n")
        
        # Basic parameters
        f.write(f"- Alignment file: `{analysis_config.get('alignment_file', 'Unknown')}`\\n")
        f.write(f"- Data type: `{analysis_config.get('data_type', 'Unknown')}`\\n")
        f.write(f"- Analysis mode: `{analysis_config.get('analysis_mode', 'Unknown')}`\\n")
        
        # ML parameters
        if analysis_config.get('has_ml'):
            f.write(f"- ML Model string: `{analysis_config.get('model', 'Unknown')}`\\n")
            f.write(f"- PAUP* `lset` command: `{analysis_config.get('paup_lset', 'Unknown')}`\\n")
        
        # Bayesian parameters
        if analysis_config.get('has_bayesian'):
            f.write(f"- Bayesian software: `{analysis_config.get('bayesian_software', 'mrbayes')}`\\n")
            f.write(f"- Bayesian model: `{analysis_config.get('bayes_model', 'Unknown')}`\\n")
            f.write(f"- MCMC generations: `{analysis_config.get('bayes_ngen', 'Unknown')}`\\n")
            f.write(f"- Burnin fraction: `{analysis_config.get('bayes_burnin', 'Unknown')}`\\n")
            f.write(f"- Marginal likelihood method: `{analysis_config.get('marginal_likelihood', 'Unknown')}`\\n")
        
        # Parsimony parameters
        if analysis_config.get('has_parsimony'):
            f.write(f"- Parsimony software: `PAUP*`\\n")
            f.write(f"- Parsimony method: `Traditional Bremer support`\\n")
            f.write(f"- Tree search: `Heuristic search with TBR swapping`\\n")
        
        f.write("\\n")
    
    def _write_summary_statistics(self, f, summary_stats: Dict[str, Any], 
                                 analysis_config: Dict[str, Any]) -> None:
        """Write the summary statistics section."""
        f.write("## Summary Statistics\\n\\n")
        
        # ML statistics
        if analysis_config.get('has_ml'):
            f.write(f"- ML tree log-likelihood: **{summary_stats.get('ml_tree_likelihood', 'Unknown')}**\\n")
            f.write(f"- Number of internal branches tested: {summary_stats.get('total_clades', 0)}\\n")
            f.write(f"- Avg ML log-likelihood difference (constrained vs ML): {summary_stats.get('avg_delta_lnl', 0.0):.4f}\\n")
            f.write(f"- Min ML log-likelihood difference: {summary_stats.get('min_delta_lnl', 0.0):.4f}\\n")
            f.write(f"- Max ML log-likelihood difference: {summary_stats.get('max_delta_lnl', 0.0):.4f}\\n")
            f.write(f"- Branches with significant AU support (p < 0.05): {summary_stats.get('au_significant_count', 0)} / {summary_stats.get('total_clades', 0)} evaluated\\n")
        
        # Bayesian statistics
        if analysis_config.get('has_bayesian'):
            f.write(f"- Avg Bayesian decay (marginal lnL difference): {summary_stats.get('avg_bayes_decay', 0.0):.4f}\\n")
            f.write(f"- Min Bayesian decay: {summary_stats.get('min_bayes_decay', 0.0):.4f}\\n")
            f.write(f"- Max Bayesian decay: {summary_stats.get('max_bayes_decay', 0.0):.4f}\\n")
        
        # Parsimony statistics
        if analysis_config.get('has_parsimony'):
            f.write(f"- Avg Parsimony decay (Bremer support): {summary_stats.get('avg_parsimony_decay', 0.0):.4f}\\n")
            f.write(f"- Min Parsimony decay: {summary_stats.get('min_parsimony_decay', 0.0):.4f}\\n")
            f.write(f"- Max Parsimony decay: {summary_stats.get('max_parsimony_decay', 0.0):.4f}\\n")
        
        f.write("\\n")
    
    def _write_bayesian_diagnostics(self, f, analysis_config: Dict[str, Any]) -> None:
        """Write Bayesian convergence diagnostics section."""
        f.write("## Bayesian Convergence Diagnostics\\n\\n")
        f.write("*Note: Detailed convergence diagnostics would be extracted from MrBayes output*\\n\\n")
        # This would be expanded with actual convergence diagnostic parsing
    
    def _write_detailed_results(self, f, decay_indices: Dict[str, Dict[str, Any]], 
                               analysis_config: Dict[str, Any]) -> None:
        """Write the detailed results table."""
        f.write("## Detailed Branch Support Results\\n\\n")
        
        # Create table header
        headers = ["Clade ID", "Taxa Count"]
        
        if analysis_config.get('has_ml'):
            headers.extend(["Constrained lnL", "ΔlnL (from ML)", "AU p-value", "Significant (AU)"])
        
        if analysis_config.get('has_bayesian'):
            headers.extend(["Bayes Decay", "BD/site", "BD%"])
        
        if analysis_config.get('has_parsimony'):
            headers.extend(["Parsimony Decay"])
        
        headers.append("Included Taxa (sample)")
        
        # Write table header
        f.write("| " + " | ".join(headers) + " |\\n")
        f.write("|" + "".join("------------ |" for _ in headers) + "\\n")
        
        # Write data rows
        for clade_id, data in decay_indices.items():
            row = [clade_id, str(data.get('taxa_count', 0))]
            
            if analysis_config.get('has_ml'):
                constrained_lnl = data.get('constrained_likelihood', 0.0)
                delta_lnl = data.get('delta_lnL', 0.0)
                au_pval = data.get('AU_pvalue', 1.0)
                significant = "**Yes**" if au_pval < 0.05 else "No"
                
                row.extend([
                    f"{constrained_lnl:.4f}",
                    f"{delta_lnl:.4f}",
                    f"{au_pval:.4f}",
                    significant
                ])
            
            if analysis_config.get('has_bayesian'):
                bd = data.get('bayes_decay', 0.0)
                bd_site = data.get('LD_per_site', 0.0)
                bd_pct = data.get('LD_percent', 0.0)
                
                row.extend([
                    f"{bd:.4f}",
                    f"{bd_site:.6f}",
                    f"{bd_pct:.3f}%"
                ])
            
            if analysis_config.get('has_parsimony'):
                pars_decay = data.get('parsimony_decay', 0)
                row.append(str(pars_decay))
            
            # Add taxa sample
            taxa_names = data.get('taxa_names', [])
            if len(taxa_names) <= 3:
                taxa_sample = ", ".join(taxa_names)
            else:
                taxa_sample = ", ".join(taxa_names[:3]) + "..."
            row.append(taxa_sample)
            
            f.write("| " + " | ".join(row) + " |\\n")
        
        f.write("\\n")
    
    def _write_interpretation_guide(self, f, analysis_config: Dict[str, Any]) -> None:
        """Write the interpretation guide section."""
        f.write("## Interpretation Guide\\n\\n")
        
        # ML interpretation
        if analysis_config.get('has_ml'):
            f.write("### ML Analysis\\n")
            f.write("- **ΔlnL (from ML)**: Log-likelihood difference between the constrained tree (without the clade) and the ML tree. Calculated as: constrained_lnL - ML_lnL. Larger positive values indicate stronger support for the clade's presence in the ML tree.\\n")
            f.write("- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.\\n\\n")
        
        # Bayesian interpretation
        if analysis_config.get('has_bayesian'):
            f.write("### Bayesian Analysis\\n")
            f.write("- **Bayes Decay (BD)**: Marginal log-likelihood difference (unconstrained - constrained). This is the primary metric for Bayesian support.\\n")
            f.write("  - **Key insight**: In phylogenetic topology testing, BD values typically closely approximate ML log-likelihood differences\\n")
            f.write("  - **Important**: BD values scale with dataset characteristics (alignment length, sequence diversity, substitution rates)\\n")
            f.write("  - **Interpretation**: Compare BD values only within this analysis. Higher BD values indicate stronger relative support\\n")
            f.write("  - **Cross-study comparisons**: Use normalized metrics (effect sizes) rather than raw BD values\\n")
            f.write("  - **Why BD ≈ ΔlnL**: When comparing models that differ only by a topological constraint, the marginal likelihood is dominated by the likelihood component\\n")
            f.write("  - **Normalized LD Metrics** (for cross-study comparisons):\\n")
            f.write("    - **LD/site**: Per-site Likelihood Decay (LD ÷ alignment length) - enables comparison across different alignment lengths\\n")
            f.write("    - **LD%**: Relative LD as percentage of unconstrained marginal likelihood\\n")
            f.write("      - Provides context relative to total likelihood magnitude\\n")
            f.write("      - Useful for comparing datasets with different likelihood scales\\n")
            f.write("    - **Why normalize?**: Raw LD values scale with dataset size, making cross-study comparisons misleading\\n")
            f.write("    - **Recommendation**: Use effect sizes for cross-study comparisons and meta-analyses\\n")
            f.write("  - **Negative values** suggest the constrained analysis had higher marginal likelihood, which may indicate:\\n")
            f.write("    - Poor chain convergence or insufficient MCMC sampling\\n")
            f.write("    - Complex posterior distribution requiring more steps\\n")
            f.write("    - Genuine lack of support for the clade\\n")
        
        # Parsimony interpretation
        if analysis_config.get('has_parsimony'):
            f.write("### Parsimony Analysis\\n")
            f.write("- **Parsimony Decay**: Traditional Bremer support value - the increase in parsimony tree length required to break the clade (make it non-monophyletic)\\n")
            f.write("  - **Calculation**: Constrained tree length - Optimal tree length\\n")
            f.write("  - **Interpretation**: Higher values indicate stronger support for the clade\\n")
            f.write("  - **Units**: Additional parsimony steps required to reject the clade\\n")
            f.write("  - **Comparison**: Values are comparable within the same dataset but not across different datasets\\n")
            f.write("  - **Zero values**: Indicate very weak support (clade can be broken with no additional steps)\\n")
            f.write("  - **Negative values**: Should not occur in proper Bremer support analysis (indicates analysis error)\\n\\n")
        
        # Understanding BD ≈ ΔlnL section (if both ML and Bayesian)
        if analysis_config.get('has_ml') and analysis_config.get('has_bayesian'):
            f.write("## Understanding BD ≈ ΔlnL in Phylogenetics\\n\\n")
            f.write("You may notice that Bayesian Decay (BD) values closely approximate the ML log-likelihood differences (ΔlnL). This is **expected behavior** in phylogenetic topology testing, not an anomaly. Here's why:\\n\\n")
            f.write("1. **Identical Models**: The constrained and unconstrained analyses use the same substitution model, differing only in whether a specific clade is allowed to exist.\\n\\n")
            f.write("2. **Likelihood Dominance**: When data strongly support a topology, the marginal likelihood (which integrates over all parameters) becomes dominated by the likelihood at the optimal parameter values.\\n\\n")
            f.write("3. **Minimal Prior Effects**: Since both analyses explore nearly identical parameter spaces (same model parameters, only different tree topologies), the prior's influence is minimal.\\n\\n")
            f.write("**Practical Implications**:\\n")
            f.write("- Similar BD and ΔlnL values confirm your analyses are working correctly\\n")
            f.write("- LD values scale with dataset characteristics - do not compare absolute values across studies\\n")
            f.write("- For cross-study comparisons, use normalized metrics (effect sizes) instead of raw LD values\\n")
            f.write("- Compare relative LD values across branches within this analysis to identify well-supported clades\\n\\n")
    
    def _write_dataset_relative_section(self, f, decay_indices: Dict[str, Dict[str, Any]], 
                                       analysis_config: Dict[str, Any]) -> None:
        """Write dataset-relative rankings section."""
        f.write("### Dataset-Relative Support Rankings\\n\\n")
        f.write("Rankings based on ML likelihood decay values within this dataset. \\n\\n")
        f.write("**Important**: \\n")
        f.write("These rankings show relative support within this dataset only and do not indicate absolute phylogenetic reliability.\\n\\n")
        
        # Sort clades by support values
        if analysis_config.get('has_ml'):
            sorted_clades = sorted(
                decay_indices.items(),
                key=lambda x: x[1].get('delta_lnL', 0.0),
                reverse=True
            )
            value_key = 'delta_lnL'
            value_label = 'Raw ML'
        elif analysis_config.get('has_bayesian'):
            sorted_clades = sorted(
                decay_indices.items(),
                key=lambda x: x[1].get('bayes_decay', 0.0),
                reverse=True
            )
            value_key = 'bayes_decay'
            value_label = 'Raw BD'
        else:
            sorted_clades = sorted(
                decay_indices.items(),
                key=lambda x: x[1].get('parsimony_decay', 0),
                reverse=True
            )
            value_key = 'parsimony_decay'
            value_label = 'Raw Parsimony'
        
        # Write highest support
        f.write("**Highest relative support within dataset:**\\n\\n")
        for i, (clade_id, data) in enumerate(sorted_clades[:5]):
            rank = i + 1
            raw_value = data.get(value_key, 0.0)
            
            if 'percentile_rank' in data and 'z_score' in data:
                percentile = data['percentile_rank']
                z_score = data['z_score']
                f.write(f"{rank}. {clade_id}: {percentile:.1f}th percentile, Z-score={z_score:+.2f}, {value_label}={raw_value:.2f}\\n")
            else:
                f.write(f"{rank}. {clade_id}: {value_label}={raw_value:.2f}\\n")
        
        # Write lowest support
        f.write("\\n**Lowest relative support within dataset:**\\n\\n")
        for i, (clade_id, data) in enumerate(sorted_clades[-3:]):
            rank = len(sorted_clades) - 2 + i
            raw_value = data.get(value_key, 0.0)
            
            if 'percentile_rank' in data and 'z_score' in data:
                percentile = data['percentile_rank']
                z_score = data['z_score']
                f.write(f"{rank}. {clade_id}: {percentile:.1f}th percentile, Z-score={z_score:+.2f}, {value_label}={raw_value:.2f}\\n")
            else:
                f.write(f"{rank}. {clade_id}: {value_label}={raw_value:.2f}\\n")
    
    def generate_quick_summary(self, decay_indices: Dict[str, Dict[str, Any]], 
                              analysis_config: Dict[str, Any]) -> str:
        """
        Generate a quick text summary of results.
        
        Args:
            decay_indices: Dictionary of decay analysis results
            analysis_config: Analysis configuration
            
        Returns:
            Quick summary string
        """
        summary_lines = []
        
        # Basic info
        total_clades = len(decay_indices)
        summary_lines.append(f"panDecay Analysis Summary ({total_clades} clades analyzed)")
        summary_lines.append("=" * 50)
        
        # ML summary
        if analysis_config.get('has_ml'):
            significant_au = sum(1 for data in decay_indices.values() 
                               if data.get('AU_pvalue', 1.0) < 0.05)
            summary_lines.append(f"ML Analysis: {significant_au}/{total_clades} clades with significant AU support (p<0.05)")
        
        # Bayesian summary
        if analysis_config.get('has_bayesian'):
            avg_bd = sum(data.get('bayes_decay', 0.0) for data in decay_indices.values()) / total_clades
            summary_lines.append(f"Bayesian Analysis: Average Bayes Decay = {avg_bd:.3f}")
        
        # Parsimony summary
        if analysis_config.get('has_parsimony'):
            avg_bremer = sum(data.get('parsimony_decay', 0) for data in decay_indices.values()) / total_clades
            summary_lines.append(f"Parsimony Analysis: Average Bremer Support = {avg_bremer:.1f}")
        
        return "\\n".join(summary_lines)