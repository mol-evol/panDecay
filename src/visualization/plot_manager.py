#!/usr/bin/env python3
"""
Plot management module for panDecay visualizations.

This module centralizes all matplotlib imports and provides a unified
interface for creating phylogenetic decay analysis plots.
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
import numpy as np

# Centralized matplotlib imports
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None
    sns = None

logger = logging.getLogger(__name__)


class PlotManager:
    """
    Manages all plotting operations for panDecay.
    
    This class centralizes matplotlib usage and provides a clean interface
    for creating various types of phylogenetic analysis plots.
    """
    
    def __init__(self, style: str = "seaborn", dpi: int = 150, figsize: Tuple[int, int] = (10, 6)):
        """
        Initialize the plot manager.
        
        Args:
            style: Matplotlib style to use
            dpi: Resolution for saved plots
            figsize: Default figure size (width, height)
        """
        self.style = style
        self.dpi = dpi
        self.figsize = figsize
        
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available. Install matplotlib and seaborn for visualization support.")
            return
        
        # Configure matplotlib style
        try:
            plt.style.use(style)
        except OSError:
            logger.warning(f"Style '{style}' not available, using default")
            plt.style.use('default')
        
        # Set seaborn style if available
        if sns:
            sns.set_context("notebook")
            sns.set_palette("viridis")
    
    def create_support_distribution_plot(self, decay_indices: Dict[str, Dict[str, Any]], 
                                       output_path: Path, value_type: str = "lnl",
                                       format: str = "png") -> bool:
        """
        Create a distribution plot of support values.
        
        Args:
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the plot
            value_type: Type of values to plot ('lnl', 'au', 'bayes', 'parsimony')
            format: Output format (png, pdf, svg)
            
        Returns:
            True if plot created successfully, False otherwise
        """
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available for plotting")
            return False
        
        try:
            # Extract values based on type
            if value_type == "lnl":
                values = [data.get('delta_lnL', 0.0) for data in decay_indices.values()]
                title = "Distribution of ML Log-Likelihood Differences (ΔlnL)"
                xlabel = "ΔlnL"
            elif value_type == "au":
                values = [data.get('AU_pvalue', 1.0) for data in decay_indices.values()]
                title = "Distribution of AU p-values"
                xlabel = "AU p-value"
            elif value_type == "bayes":
                values = [data.get('bayes_decay', 0.0) for data in decay_indices.values()]
                title = "Distribution of Bayesian Decay Values"
                xlabel = "Bayes Decay"
            elif value_type == "parsimony":
                values = [data.get('parsimony_decay', 0) for data in decay_indices.values()]
                title = "Distribution of Parsimony Decay (Bremer Support)"
                xlabel = "Parsimony Decay"
            else:
                logger.error(f"Unknown value type: {value_type}")
                return False
            
            # Create figure
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=self.figsize)
            
            # Histogram
            ax1.hist(values, bins=min(len(values), 20), alpha=0.7, edgecolor='black')
            ax1.set_xlabel(xlabel)
            ax1.set_ylabel("Frequency")
            ax1.set_title(f"Histogram of {xlabel}")
            ax1.grid(True, alpha=0.3)
            
            # Box plot
            ax2.boxplot(values)
            ax2.set_ylabel(xlabel)
            ax2.set_title(f"Box Plot of {xlabel}")
            ax2.grid(True, alpha=0.3)
            
            # Overall title
            fig.suptitle(title, fontsize=14, fontweight='bold')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save plot
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Created support distribution plot: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to create support distribution plot: {e}")
            return False
    
    def create_scatter_plot(self, decay_indices: Dict[str, Dict[str, Any]], 
                           output_path: Path, x_key: str, y_key: str,
                           format: str = "png") -> bool:
        """
        Create a scatter plot comparing two support metrics.
        
        Args:
            decay_indices: Dictionary of decay analysis results
            output_path: Path to save the plot
            x_key: Key for x-axis values
            y_key: Key for y-axis values
            format: Output format (png, pdf, svg)
            
        Returns:
            True if plot created successfully, False otherwise
        """
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available for plotting")
            return False
        
        try:
            # Extract values
            x_values = []
            y_values = []
            labels = []
            
            for clade_id, data in decay_indices.items():
                x_val = data.get(x_key)
                y_val = data.get(y_key)
                
                if x_val is not None and y_val is not None:
                    x_values.append(x_val)
                    y_values.append(y_val)
                    labels.append(clade_id)
            
            if not x_values or not y_values:
                logger.warning(f"No data available for scatter plot ({x_key} vs {y_key})")
                return False
            
            # Create figure
            fig, ax = plt.subplots(figsize=self.figsize)
            
            # Create scatter plot
            scatter = ax.scatter(x_values, y_values, alpha=0.7, s=50)
            
            # Add labels for points
            for i, label in enumerate(labels):
                ax.annotate(label, (x_values[i], y_values[i]), 
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=8, alpha=0.7)
            
            # Set labels and title
            ax.set_xlabel(self._get_axis_label(x_key))
            ax.set_ylabel(self._get_axis_label(y_key))
            ax.set_title(f"{self._get_axis_label(y_key)} vs {self._get_axis_label(x_key)}")
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
            # Calculate and add correlation if both are numeric
            if len(x_values) > 2:
                correlation = np.corrcoef(x_values, y_values)[0, 1]
                ax.text(0.05, 0.95, f"r = {correlation:.3f}", 
                       transform=ax.transAxes, fontsize=10,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            # Save plot
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Created scatter plot: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to create scatter plot: {e}")
            return False
    
    def create_site_analysis_plot(self, clade_id: str, site_data: Dict[int, Dict[str, Any]], 
                                 output_path: Path, format: str = "png") -> bool:
        """
        Create a site-specific analysis plot.
        
        Args:
            clade_id: Identifier for the clade
            site_data: Dictionary mapping site numbers to site-specific data
            output_path: Path to save the plot
            format: Output format (png, pdf, svg)
            
        Returns:
            True if plot created successfully, False otherwise
        """
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available for plotting")
            return False
        
        try:
            # Extract site numbers and delta values
            sites = sorted(site_data.keys())
            deltas = [site_data[site].get('delta_lnL', 0.0) for site in sites]
            
            if not sites or not deltas:
                logger.warning(f"No site data available for {clade_id}")
                return False
            
            # Create figure
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Create line plot
            ax.plot(sites, deltas, 'o-', alpha=0.7, markersize=3)
            
            # Add horizontal line at y=0
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
            
            # Set labels and title
            ax.set_xlabel("Site Position")
            ax.set_ylabel("Site-specific ΔlnL")
            ax.set_title(f"Site-specific Support for {clade_id}")
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
            # Add summary statistics
            mean_delta = np.mean(deltas)
            std_delta = np.std(deltas)
            ax.text(0.02, 0.98, f"Mean: {mean_delta:.3f}\\nStd: {std_delta:.3f}", 
                   transform=ax.transAxes, fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            # Save plot
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Created site analysis plot for {clade_id}: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to create site analysis plot for {clade_id}: {e}")
            return False
    
    def create_site_histogram(self, clade_id: str, delta_values: List[float], 
                             output_path: Path, format: str = "png") -> bool:
        """
        Create a histogram of site-specific delta values.
        
        Args:
            clade_id: Identifier for the clade
            delta_values: List of delta lnL values
            output_path: Path to save the plot
            format: Output format (png, pdf, svg)
            
        Returns:
            True if plot created successfully, False otherwise
        """
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available for plotting")
            return False
        
        try:
            if not delta_values:
                logger.warning(f"No delta values for histogram of {clade_id}")
                return False
            
            # Create figure
            fig, ax = plt.subplots(figsize=self.figsize)
            
            # Create histogram
            n_bins = min(len(delta_values) // 5, 50)  # Adaptive bin count
            n, bins, patches = ax.hist(delta_values, bins=n_bins, alpha=0.7, edgecolor='black')
            
            # Color bars based on sign (positive = supporting, negative = conflicting)
            for i, (patch, bin_left, bin_right) in enumerate(zip(patches, bins[:-1], bins[1:])):
                if (bin_left + bin_right) / 2 > 0:
                    patch.set_facecolor('green')
                    patch.set_alpha(0.6)
                else:
                    patch.set_facecolor('red')
                    patch.set_alpha(0.6)
            
            # Add vertical line at x=0
            ax.axvline(x=0, color='black', linestyle='--', alpha=0.7)
            
            # Set labels and title
            ax.set_xlabel("Site-specific ΔlnL")
            ax.set_ylabel("Number of Sites")
            ax.set_title(f"Distribution of Site Support for {clade_id}")
            
            # Add grid
            ax.grid(True, alpha=0.3)
            
            # Add summary statistics
            supporting_sites = sum(1 for d in delta_values if d > 0)
            conflicting_sites = sum(1 for d in delta_values if d < 0)
            neutral_sites = len(delta_values) - supporting_sites - conflicting_sites
            
            stats_text = f"Supporting: {supporting_sites}\\nConflicting: {conflicting_sites}\\nNeutral: {neutral_sites}"
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10, 
                   verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            # Save plot
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            logger.info(f"Created site histogram for {clade_id}: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Failed to create site histogram for {clade_id}: {e}")
            return False
    
    def _get_axis_label(self, key: str) -> str:
        """
        Get a human-readable axis label for a data key.
        
        Args:
            key: Data key
            
        Returns:
            Human-readable label
        """
        label_map = {
            'delta_lnL': 'ΔlnL',
            'AU_pvalue': 'AU p-value',
            'bayes_decay': 'Bayes Decay',
            'parsimony_decay': 'Parsimony Decay',
            'bootstrap_support': 'Bootstrap Support',
            'LD_per_site': 'LD per site',
            'LD_percent': 'LD %'
        }
        
        return label_map.get(key, key.replace('_', ' ').title())
    
    def check_matplotlib_availability(self) -> bool:
        """
        Check if matplotlib is available for plotting.
        
        Returns:
            True if matplotlib is available, False otherwise
        """
        return HAS_MATPLOTLIB