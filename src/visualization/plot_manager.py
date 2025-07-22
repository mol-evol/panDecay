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
            
            # Color bars based on sign (negative = supporting, positive = conflicting)
            for i, (patch, bin_left, bin_right) in enumerate(zip(patches, bins[:-1], bins[1:])):
                if (bin_left + bin_right) / 2 < 0:  # Negative delta = supporting = green
                    patch.set_facecolor('green')
                    patch.set_alpha(0.6)
                else:  # Positive delta = conflicting = red
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
    
    def create_alignment_visualization(self, alignment, site_data_by_clade: Dict[str, Dict[int, Dict[str, Any]]], 
                                     output_dir: Path, chunk_size: int = 2000, 
                                     layout: str = "auto", format: str = "png", 
                                     file_tracker=None) -> bool:
        """
        Create combined alignment and site-specific support visualizations.
        Each clade gets one image with alignment on top and support bars below.
        
        Args:
            alignment: BioPython alignment object
            site_data_by_clade: Dictionary mapping clade IDs to site-specific data
            output_dir: Directory to save visualization files
            chunk_size: Maximum sites per visualization chunk
            layout: Layout mode (maintained for compatibility, always creates combined views)
            format: Output format (png, pdf, svg)
            
        Returns:
            True if visualizations created successfully, False otherwise
        """
        if not HAS_MATPLOTLIB:
            logger.error("Matplotlib not available for alignment visualization")
            return False
        
        try:
            # Create output directory
            viz_dir = output_dir / "alignment_visualization"
            viz_dir.mkdir(exist_ok=True)
            
            # Get alignment dimensions
            num_sites = alignment.get_alignment_length()
            num_clades = len(site_data_by_clade)
            
            # Create alignment chunks if needed
            chunks = self._chunk_alignment(alignment, chunk_size)
            
            success = True
            for i, (chunk_start, chunk_end, chunk_alignment) in enumerate(chunks):
                chunk_suffix = f"_chunk_{i+1}_sites_{chunk_start+1}-{chunk_end}" if len(chunks) > 1 else ""
                
                # Create one combined visualization per clade
                for clade_id, site_data in site_data_by_clade.items():
                    file_name = f"alignment_with_support_{clade_id}{chunk_suffix}.{format}"
                    output_path = viz_dir / file_name
                    
                    clade_success = self._create_combined_alignment_and_support(
                        chunk_alignment, chunk_start, clade_id, site_data, 
                        chunk_start, chunk_end, output_path, format
                    )
                    
                    # Track file with FileTracker if available
                    if clade_success and file_tracker:
                        try:
                            file_tracker.track_file('visualizations', output_path, f'Alignment visualization with support for {clade_id}')
                        except Exception as e:
                            logger.debug(f"Could not track alignment visualization file: {e}")
                    
                    if not clade_success:
                        success = False
            
            # Create summary file
            summary_path = self._create_visualization_summary(viz_dir, num_sites, num_clades, len(chunks), "combined")
            
            # Track summary file with FileTracker if available
            if summary_path and file_tracker:
                try:
                    file_tracker.track_file('visualizations', summary_path, 'Alignment visualization summary')
                except Exception as e:
                    logger.debug(f"Could not track alignment visualization summary: {e}")
            
            if success:
                logger.info(f"Created alignment visualizations in {viz_dir}")
            
            return success
            
        except Exception as e:
            logger.error(f"Failed to create alignment visualization: {e}")
            return False
    
    def _determine_layout_strategy(self, num_sites: int, num_clades: int) -> str:
        """
        Determine optimal layout strategy based on dataset size.
        
        Args:
            num_sites: Number of alignment sites
            num_clades: Number of clades to visualize
            
        Returns:
            Layout strategy ("single" or "separate")
        """
        if num_sites <= 500 and num_clades <= 5:
            return "single"
        else:
            return "separate"
    
    def _chunk_alignment(self, alignment, chunk_size: int) -> List[Tuple[int, int, Any]]:
        """
        Split alignment into chunks for large datasets.
        
        Args:
            alignment: BioPython alignment object
            chunk_size: Maximum sites per chunk
            
        Returns:
            List of tuples (start_pos, end_pos, chunk_alignment)
        """
        alignment_length = alignment.get_alignment_length()
        
        if alignment_length <= chunk_size:
            return [(0, alignment_length, alignment)]
        
        chunks = []
        overlap = 100  # Sites to overlap between chunks
        
        start = 0
        while start < alignment_length:
            end = min(start + chunk_size, alignment_length)
            
            # Create chunk alignment
            chunk_records = []
            for record in alignment:
                chunk_seq = record.seq[start:end]
                chunk_record = record[:]  # Copy record
                chunk_record.seq = chunk_seq
                chunk_records.append(chunk_record)
            
            # Create alignment object for chunk
            from Bio.Align import MultipleSeqAlignment
            chunk_alignment = MultipleSeqAlignment(chunk_records)
            
            chunks.append((start, end, chunk_alignment))
            
            # Move start position with overlap
            start = end - overlap
            if start >= alignment_length - overlap:
                break
        
        return chunks
    
    def _create_combined_alignment_and_support(self, alignment, start_pos: int, clade_id: str,
                                            site_data: Dict[int, Dict[str, Any]], 
                                            chunk_start: int, chunk_end: int, 
                                            output_path: Path, format: str) -> bool:
        """
        Create combined visualization with alignment on top and site support below.
        
        Args:
            alignment: BioPython alignment object (chunk)
            start_pos: Starting position in original alignment
            clade_id: Clade identifier
            site_data: Site-specific data for the clade
            chunk_start: Start position of chunk
            chunk_end: End position of chunk
            output_path: Path to save combined image
            format: Output format
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Filter sites in chunk range
            chunk_sites = {pos: data for pos, data in site_data.items() 
                          if chunk_start <= pos < chunk_end}
            
            if not chunk_sites:
                logger.warning(f"No site data for {clade_id} in range {chunk_start}-{chunk_end}")
                return False
            
            # Determine data type and color scheme for alignment
            color_scheme = self._get_color_scheme(alignment)
            
            # Convert alignment to numeric matrix
            matrix = self._alignment_to_matrix(alignment, color_scheme)
            
            # Calculate figure dimensions
            height = len(alignment)
            width = alignment.get_alignment_length()
            
            # Scale figure size based on alignment dimensions
            fig_width = min(max(width / 100, 10), 24)  # 10-24 inches wide
            
            # Create figure with two subplots: alignment on top, support below
            fig, (ax_align, ax_support) = plt.subplots(2, 1, figsize=(fig_width, 8), 
                                                      height_ratios=[height, 3],
                                                      sharex=True)
            
            # === TOP PANEL: ALIGNMENT ===
            im = ax_align.imshow(matrix, aspect='auto', interpolation='nearest')
            
            # Set colormap for alignment based on data type
            if 'dna' in color_scheme:
                colors = ['white', 'red', 'blue', 'green', 'yellow']  # Gap, A, T, G, C
                cmap = plt.matplotlib.colors.ListedColormap(colors)
                im.set_cmap(cmap)
            
            # Set alignment labels
            ax_align.set_ylabel("Sequences")
            ax_align.set_title(f"Alignment with Site Support for {clade_id}")
            
            # Set sequence labels
            seq_names = [record.id for record in alignment]
            ax_align.set_yticks(range(len(seq_names)))
            ax_align.set_yticklabels(seq_names, fontsize=8)
            
            # === BOTTOM PANEL: SITE SUPPORT ===
            sites = sorted(chunk_sites.keys())
            delta_values = [chunk_sites[site].get('delta_lnL', 0.0) for site in sites]
            
            # Adjust site positions relative to chunk start
            plot_positions = [site - chunk_start for site in sites]
            
            # Create bar colors (negative delta = supporting = green)
            colors = ['green' if delta < 0 else 'red' for delta in delta_values]
            
            # Plot support bars
            ax_support.bar(plot_positions, delta_values, color=colors, alpha=0.7, width=1.0)
            
            # Add horizontal line at y=0
            ax_support.axhline(y=0, color='black', linestyle='-', alpha=0.3)
            
            # Set support panel labels
            ax_support.set_ylabel("ΔlnL")
            ax_support.set_xlabel(f"Alignment Position (starting from {chunk_start + 1})")
            ax_support.grid(True, alpha=0.3)
            
            # Add summary statistics
            supporting = sum(1 for d in delta_values if d < 0)
            conflicting = sum(1 for d in delta_values if d >= 0)
            ax_support.text(0.02, 0.98, f"Supporting: {supporting} | Conflicting: {conflicting}", 
                           transform=ax_support.transAxes, fontsize=10, verticalalignment='top',
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            # === SHARED X-AXIS FORMATTING ===
            # Adjust tick frequency for readability
            if width > 100:
                tick_interval = max(width // 20, 10)
                tick_positions = range(0, width, tick_interval)
                ax_support.set_xticks(tick_positions)
                ax_support.set_xticklabels([str(chunk_start + pos + 1) for pos in tick_positions])
            
            # Save plot
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create combined alignment visualization for {clade_id}: {e}")
            return False
    
    def _create_alignment_image(self, alignment, start_pos: int, output_path: Path, format: str) -> bool:
        """
        Create colored alignment image.
        
        Args:
            alignment: BioPython alignment object (chunk)
            start_pos: Starting position in original alignment
            output_path: Path to save image
            format: Output format
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Determine data type and color scheme
            color_scheme = self._get_color_scheme(alignment)
            
            # Convert alignment to numeric matrix
            matrix = self._alignment_to_matrix(alignment, color_scheme)
            
            # Create figure
            height = len(alignment)
            width = alignment.get_alignment_length()
            
            # Scale figure size based on alignment dimensions
            fig_width = min(max(width / 100, 8), 20)  # 8-20 inches wide
            fig_height = min(max(height / 10, 4), 12)  # 4-12 inches tall
            
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            
            # Create image
            im = ax.imshow(matrix, aspect='auto', interpolation='nearest')
            
            # Set colormap based on data type
            if 'dna' in color_scheme:
                colors = ['white', 'red', 'blue', 'green', 'yellow']  # Gap, A, T, G, C
                cmap = plt.matplotlib.colors.ListedColormap(colors)
                im.set_cmap(cmap)
            
            # Set labels
            ax.set_xlabel(f"Alignment Position (starting from {start_pos + 1})")
            ax.set_ylabel("Sequences")
            ax.set_title(f"Multiple Sequence Alignment Visualization")
            
            # Set sequence labels
            seq_names = [record.id for record in alignment]
            ax.set_yticks(range(len(seq_names)))
            ax.set_yticklabels(seq_names, fontsize=8)
            
            # Adjust tick frequency for readability
            if width > 100:
                tick_interval = max(width // 20, 10)
                tick_positions = range(0, width, tick_interval)
                ax.set_xticks(tick_positions)
                ax.set_xticklabels([str(start_pos + pos + 1) for pos in tick_positions])
            
            # Save plot
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create alignment image: {e}")
            return False
    
    def _create_combined_site_overlay(self, site_data_by_clade: Dict[str, Dict[int, Dict[str, Any]]], 
                                    start_pos: int, end_pos: int, output_path: Path, format: str) -> bool:
        """
        Create combined site support overlay for all clades.
        
        Args:
            site_data_by_clade: Site data for all clades
            start_pos: Start position in alignment
            end_pos: End position in alignment
            output_path: Path to save overlay
            format: Output format
            
        Returns:
            True if successful, False otherwise
        """
        try:
            num_clades = len(site_data_by_clade)
            fig, axes = plt.subplots(num_clades, 1, figsize=(12, 2 * num_clades), sharex=True)
            
            if num_clades == 1:
                axes = [axes]
            
            for i, (clade_id, site_data) in enumerate(site_data_by_clade.items()):
                # Filter sites in chunk range
                chunk_sites = {pos: data for pos, data in site_data.items() 
                              if start_pos <= pos < end_pos}
                
                if chunk_sites:
                    self._plot_site_support_bars(axes[i], chunk_sites, start_pos, clade_id)
                else:
                    axes[i].text(0.5, 0.5, f"No site data for {clade_id}", 
                               transform=axes[i].transAxes, ha='center', va='center')
            
            plt.xlabel(f"Alignment Position (starting from {start_pos + 1})")
            plt.suptitle("Site-specific Branch Support")
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create combined site overlay: {e}")
            return False
    
    def _create_single_clade_overlay(self, clade_id: str, site_data: Dict[int, Dict[str, Any]], 
                                   start_pos: int, end_pos: int, output_path: Path, format: str) -> bool:
        """
        Create site support overlay for single clade.
        
        Args:
            clade_id: Clade identifier
            site_data: Site-specific data for the clade
            start_pos: Start position in alignment
            end_pos: End position in alignment
            output_path: Path to save overlay
            format: Output format
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Filter sites in chunk range
            chunk_sites = {pos: data for pos, data in site_data.items() 
                          if start_pos <= pos < end_pos}
            
            if not chunk_sites:
                logger.warning(f"No site data for {clade_id} in range {start_pos}-{end_pos}")
                return False
            
            fig, ax = plt.subplots(figsize=(12, 3))
            self._plot_site_support_bars(ax, chunk_sites, start_pos, clade_id)
            
            plt.xlabel(f"Alignment Position (starting from {start_pos + 1})")
            plt.title(f"Site-specific Support for {clade_id}")
            plt.tight_layout()
            plt.savefig(str(output_path), dpi=self.dpi, format=format, bbox_inches='tight')
            plt.close()
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to create single clade overlay for {clade_id}: {e}")
            return False
    
    def _plot_site_support_bars(self, ax, site_data: Dict[int, Dict[str, Any]], 
                               start_pos: int, clade_id: str):
        """
        Plot site support bars on given axis.
        
        Args:
            ax: Matplotlib axis
            site_data: Site-specific data
            start_pos: Starting position offset
            clade_id: Clade identifier
        """
        sites = sorted(site_data.keys())
        delta_values = [site_data[site].get('delta_lnL', 0.0) for site in sites]
        
        # Adjust site positions relative to chunk start
        plot_positions = [site - start_pos for site in sites]
        
        # Create bar colors based on support
        # Negative delta = ML tree better than constrained = supports branch = green
        # Positive delta = constrained tree better than ML = conflicts with branch = red
        colors = ['green' if delta < 0 else 'red' for delta in delta_values]
        
        # Plot bars
        bars = ax.bar(plot_positions, delta_values, color=colors, alpha=0.7, width=1.0)
        
        # Add horizontal line at y=0
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
        
        # Set labels
        ax.set_ylabel("ΔlnL")
        ax.set_title(f"{clade_id}")
        ax.grid(True, alpha=0.3)
        
        # Add summary statistics
        supporting = sum(1 for d in delta_values if d > 0)
        conflicting = sum(1 for d in delta_values if d < 0)
        ax.text(0.02, 0.98, f"Support: {supporting} | Conflict: {conflicting}", 
               transform=ax.transAxes, fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    def _get_color_scheme(self, alignment) -> Dict[str, int]:
        """
        Determine appropriate color scheme based on alignment content.
        
        Args:
            alignment: BioPython alignment object
            
        Returns:
            Dictionary mapping characters to numeric codes
        """
        # Sample first sequence to determine data type
        sample_seq = str(alignment[0].seq).upper()
        
        # Check for DNA
        dna_chars = set('ATGC')
        if set(sample_seq.replace('-', '').replace('N', '')) <= dna_chars:
            return {
                'dna': True,
                '-': 0, 'N': 0,  # Gaps and ambiguous
                'A': 1, 'T': 2, 'G': 3, 'C': 4
            }
        
        # Default to discrete character mapping
        all_chars = set()
        for record in alignment:
            all_chars.update(str(record.seq).upper())
        
        color_map = {'discrete': True}
        for i, char in enumerate(sorted(all_chars)):
            color_map[char] = i
        
        return color_map
    
    def _alignment_to_matrix(self, alignment, color_scheme: Dict[str, int]) -> np.ndarray:
        """
        Convert alignment to numeric matrix for visualization.
        
        Args:
            alignment: BioPython alignment object
            color_scheme: Character to number mapping
            
        Returns:
            Numeric matrix representing alignment
        """
        height = len(alignment)
        width = alignment.get_alignment_length()
        matrix = np.zeros((height, width), dtype=int)
        
        for i, record in enumerate(alignment):
            seq_str = str(record.seq).upper()
            for j, char in enumerate(seq_str):
                matrix[i, j] = color_scheme.get(char, 0)  # Default to 0 for unknown chars
        
        return matrix
    
    def _create_visualization_summary(self, viz_dir: Path, num_sites: int, 
                                    num_clades: int, num_chunks: int, layout: str) -> Path:
        """
        Create summary file describing the visualization outputs.
        
        Args:
            viz_dir: Visualization directory
            num_sites: Total number of sites
            num_clades: Number of clades
            num_chunks: Number of chunks created
            layout: Layout strategy used
            
        Returns:
            Path to created summary file
        """
        summary_path = viz_dir / "visualization_summary.txt"
        
        with open(summary_path, 'w') as f:
            f.write("Alignment Visualization Summary\n")
            f.write("=" * 35 + "\n\n")
            f.write(f"Total alignment sites: {num_sites}\n")
            f.write(f"Number of clades analyzed: {num_clades}\n")
            f.write(f"Number of visualization chunks: {num_chunks}\n")
            f.write(f"Layout strategy: {layout}\n\n")
            
            f.write("Files created:\n")
            f.write("-" * 15 + "\n")
            
            if num_chunks == 1:
                f.write("• alignment_with_support_<clade>.png - Combined alignment + support visualization per clade\n")
                f.write("  (Top panel: colored sequence alignment, Bottom panel: site support bars)\n")
            else:
                f.write("• alignment_with_support_<clade>_chunk_N_sites_X-Y.png - Combined visualizations per clade and chunk\n")
                f.write("  (Top panel: colored sequence alignment, Bottom panel: site support bars)\n")
            
            f.write("\nVisualization layout:\n")
            f.write("-" * 20 + "\n")
            f.write("• Top panel: Colored sequence alignment\n")
            f.write("  - DNA: A=red, T=blue, G=green, C=yellow, gaps=white\n")
            f.write("• Bottom panel: Site-specific branch support bars\n")
            f.write("  - Perfectly aligned columns between panels\n")
            
            f.write("\nSupport color scheme:\n")
            f.write("-" * 20 + "\n")
            f.write("• Green bars: Sites supporting the branch (negative ΔlnL)\n")
            f.write("  (ML tree has better likelihood than constrained tree)\n")
            f.write("• Red bars: Sites conflicting with the branch (positive ΔlnL)\n")
            f.write("  (Constrained tree has better likelihood than ML tree)\n")
            f.write("• Bar height: Proportional to |ΔlnL| magnitude\n")
        
        return summary_path

    def check_matplotlib_availability(self) -> bool:
        """
        Check if matplotlib is available for plotting.
        
        Returns:
            True if matplotlib is available, False otherwise
        """
        return HAS_MATPLOTLIB