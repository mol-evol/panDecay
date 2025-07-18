#!/usr/bin/env python3
"""
Dual visualization system for panDecay.

This module provides both static (matplotlib/seaborn) and interactive (Plotly)
visualizations for phylogenetic decay analysis results.
"""

import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
import warnings

# Static plotting imports
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    plt = None
    sns = None

# Interactive plotting imports
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.offline as pyo
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None

logger = logging.getLogger(__name__)


class VisualizationError(Exception):
    """Exception raised for visualization-related errors."""
    pass


class DualVisualizationSystem:
    """
    Comprehensive visualization system supporting both static and interactive plots.
    """
    
    def __init__(self, 
                 output_format: str = "both",
                 static_config: Optional[Dict[str, Any]] = None,
                 interactive_config: Optional[Dict[str, Any]] = None):
        """
        Initialize dual visualization system.
        
        Args:
            output_format: Format to generate ("static", "interactive", "both")
            static_config: Configuration for static plots
            interactive_config: Configuration for interactive plots
        """
        self.output_format = output_format
        
        # Default static configuration
        self.static_config = {
            'dpi': 300,
            'formats': ['png', 'pdf'],
            'style': 'publication',
            'figsize': (10, 8),
            'color_palette': 'viridis',
            'font_size': 12
        }
        if static_config:
            self.static_config.update(static_config)
        
        # Default interactive configuration
        self.interactive_config = {
            'theme': 'plotly_white',
            'export_html': True,
            'include_controls': True,
            'width': 1000,
            'height': 800,
            'color_scale': 'viridis'
        }
        if interactive_config:
            self.interactive_config.update(interactive_config)
        
        # Check availability
        self.has_matplotlib = HAS_MATPLOTLIB
        self.has_plotly = HAS_PLOTLY
        
        # Validate format choice
        if output_format == "static" and not self.has_matplotlib:
            raise VisualizationError("Matplotlib not available for static plots")
        elif output_format == "interactive" and not self.has_plotly:
            raise VisualizationError("Plotly not available for interactive plots")
        elif output_format == "both" and not (self.has_matplotlib or self.has_plotly):
            raise VisualizationError("Neither Matplotlib nor Plotly available")
        
        logger.info(f"DualVisualizationSystem initialized: format={output_format}")
        logger.info(f"  Matplotlib available: {self.has_matplotlib}")
        logger.info(f"  Plotly available: {self.has_plotly}")
    
    def create_visualizations(self, 
                            decay_data: Dict[str, Dict[str, Any]], 
                            output_dir: Path,
                            viz_type: str = "decay_summary") -> Dict[str, List[str]]:
        """
        Create visualizations based on the specified format.
        
        Args:
            decay_data: Decay indices data
            output_dir: Output directory for plots
            viz_type: Type of visualization to create
            
        Returns:
            Dictionary mapping format type to list of created files
        """
        created_files = {'static': [], 'interactive': []}
        
        # Ensure output directories exist
        static_dir = output_dir / "static_plots"
        interactive_dir = output_dir / "interactive_plots"
        
        if self.output_format in ["static", "both"] and self.has_matplotlib:
            static_dir.mkdir(exist_ok=True)
            static_files = self._create_static_visualizations(decay_data, static_dir, viz_type)
            created_files['static'] = static_files
        
        if self.output_format in ["interactive", "both"] and self.has_plotly:
            interactive_dir.mkdir(exist_ok=True)
            interactive_files = self._create_interactive_visualizations(decay_data, interactive_dir, viz_type)
            created_files['interactive'] = interactive_files
        
        # Create combined report if both formats are available
        if self.output_format == "both" and created_files['static'] and created_files['interactive']:
            report_file = self._create_combined_report(created_files, output_dir)
            created_files['report'] = [report_file]
        
        return created_files
    
    def _create_static_visualizations(self, 
                                    decay_data: Dict[str, Dict[str, Any]], 
                                    output_dir: Path,
                                    viz_type: str) -> List[str]:
        """Create static matplotlib/seaborn visualizations."""
        if not self.has_matplotlib:
            logger.warning("Matplotlib not available, skipping static plots")
            return []
        
        created_files = []
        
        # Set up matplotlib style
        self._setup_matplotlib_style()
        
        if viz_type == "decay_summary":
            created_files.extend(self._create_static_decay_summary(decay_data, output_dir))
        elif viz_type == "site_analysis":
            created_files.extend(self._create_static_site_analysis(decay_data, output_dir))
        elif viz_type == "tree_support":
            created_files.extend(self._create_static_tree_support(decay_data, output_dir))
        
        return created_files
    
    def _create_interactive_visualizations(self, 
                                         decay_data: Dict[str, Dict[str, Any]], 
                                         output_dir: Path,
                                         viz_type: str) -> List[str]:
        """Create interactive Plotly visualizations."""
        if not self.has_plotly:
            logger.warning("Plotly not available, skipping interactive plots")
            return []
        
        created_files = []
        
        if viz_type == "decay_summary":
            created_files.extend(self._create_interactive_decay_summary(decay_data, output_dir))
        elif viz_type == "site_analysis":
            created_files.extend(self._create_interactive_site_analysis(decay_data, output_dir))
        elif viz_type == "tree_support":
            created_files.extend(self._create_interactive_tree_support(decay_data, output_dir))
        
        return created_files
    
    def _setup_matplotlib_style(self) -> None:
        """Set up matplotlib style for publication-quality plots."""
        style = self.static_config.get('style', 'publication')
        
        if style == 'publication':
            # Use available seaborn style or default
            available_styles = plt.style.available
            if 'seaborn-v0_8-whitegrid' in available_styles:
                plt.style.use('seaborn-v0_8-whitegrid')
            elif 'seaborn-whitegrid' in available_styles:
                plt.style.use('seaborn-whitegrid')
            else:
                plt.style.use('default')
                
            plt.rcParams.update({
                'font.size': self.static_config.get('font_size', 12),
                'font.family': 'sans-serif',
                'font.sans-serif': ['DejaVu Sans', 'Arial', 'Liberation Sans', 'sans-serif'],
                'axes.linewidth': 1.2,
                'axes.spines.top': False,
                'axes.spines.right': False,
                'xtick.direction': 'out',
                'ytick.direction': 'out',
                'figure.dpi': self.static_config.get('dpi', 300),
                'savefig.dpi': self.static_config.get('dpi', 300),
                'savefig.bbox': 'tight',
                'savefig.pad_inches': 0.1
            })
    
    def _create_static_decay_summary(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create static decay summary visualizations."""
        created_files = []
        
        if not decay_data:
            logger.warning("No decay data available for static visualization")
            return created_files
        
        # Extract data for plotting
        plot_data = self._extract_plot_data(decay_data)
        
        if not plot_data['clades']:
            logger.warning("No valid decay data for plotting")
            return created_files
        
        # Create multi-panel figure
        fig, axes = plt.subplots(2, 2, figsize=self.static_config['figsize'])
        fig.suptitle('Phylogenetic Decay Analysis Summary', fontsize=16, fontweight='bold')
        
        # 1. Decay values distribution
        if plot_data['ml_decay']:
            axes[0, 0].hist(plot_data['ml_decay'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
            axes[0, 0].set_title('ML Decay Distribution')
            axes[0, 0].set_xlabel('ML Decay Index')
            axes[0, 0].set_ylabel('Frequency')
        
        # 2. AU p-values distribution
        if plot_data['au_pvalues']:
            axes[0, 1].hist(plot_data['au_pvalues'], bins=20, alpha=0.7, color='lightcoral', edgecolor='black')
            axes[0, 1].axvline(x=0.05, color='red', linestyle='--', alpha=0.8, label='α = 0.05')
            axes[0, 1].set_title('AU P-values Distribution')
            axes[0, 1].set_xlabel('AU P-value')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].legend()
        
        # 3. ML vs Bayesian decay (if available)
        if plot_data['ml_decay'] and plot_data['bayesian_decay']:
            axes[1, 0].scatter(plot_data['ml_decay'], plot_data['bayesian_decay'], alpha=0.7, color='green')
            axes[1, 0].plot([0, max(plot_data['ml_decay'])], [0, max(plot_data['ml_decay'])], 'r--', alpha=0.5)
            axes[1, 0].set_title('ML vs Bayesian Decay')
            axes[1, 0].set_xlabel('ML Decay Index')
            axes[1, 0].set_ylabel('Bayesian Decay Index')
        
        # 4. Support summary by clade size
        if plot_data['clade_sizes'] and plot_data['ml_decay']:
            scatter = axes[1, 1].scatter(plot_data['clade_sizes'], plot_data['ml_decay'], 
                                       c=plot_data['au_pvalues'] if plot_data['au_pvalues'] else 'blue',
                                       cmap='RdYlBu', alpha=0.7)
            axes[1, 1].set_title('Support vs Clade Size')
            axes[1, 1].set_xlabel('Clade Size (# taxa)')
            axes[1, 1].set_ylabel('ML Decay Index')
            if plot_data['au_pvalues']:
                plt.colorbar(scatter, ax=axes[1, 1], label='AU P-value')
        
        plt.tight_layout()
        
        # Save in multiple formats
        for fmt in self.static_config['formats']:
            filename = output_dir / f"decay_summary.{fmt}"
            plt.savefig(filename, format=fmt, dpi=self.static_config['dpi'])
            created_files.append(str(filename))
        
        plt.close()
        logger.info(f"Created static decay summary plots: {created_files}")
        return created_files
    
    def _create_interactive_decay_summary(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create interactive decay summary visualizations."""
        created_files = []
        
        if not decay_data:
            logger.warning("No decay data available for interactive visualization")
            return created_files
        
        # Extract data for plotting
        plot_data = self._extract_plot_data(decay_data)
        
        if not plot_data['clades']:
            logger.warning("No valid decay data for interactive plotting")
            return created_files
        
        # Create interactive dashboard
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=['ML Decay Distribution', 'AU P-values Distribution',
                           'ML vs Bayesian Decay', 'Support vs Clade Size'],
            specs=[[{"type": "histogram"}, {"type": "histogram"}],
                   [{"type": "scatter"}, {"type": "scatter"}]]
        )
        
        # 1. ML Decay distribution
        if plot_data['ml_decay']:
            fig.add_trace(
                go.Histogram(
                    x=plot_data['ml_decay'],
                    name='ML Decay',
                    nbinsx=20,
                    marker_color='skyblue',
                    hovertemplate='ML Decay: %{x:.3f}<br>Count: %{y}<extra></extra>'
                ),
                row=1, col=1
            )
        
        # 2. AU P-values distribution
        if plot_data['au_pvalues']:
            fig.add_trace(
                go.Histogram(
                    x=plot_data['au_pvalues'],
                    name='AU P-values',
                    nbinsx=20,
                    marker_color='lightcoral',
                    hovertemplate='AU P-value: %{x:.3f}<br>Count: %{y}<extra></extra>'
                ),
                row=1, col=2
            )
            # Add significance line
            fig.add_vline(x=0.05, line_dash="dash", line_color="red", row=1, col=2)
        
        # 3. ML vs Bayesian decay
        if plot_data['ml_decay'] and plot_data['bayesian_decay']:
            fig.add_trace(
                go.Scatter(
                    x=plot_data['ml_decay'],
                    y=plot_data['bayesian_decay'],
                    mode='markers',
                    name='ML vs Bayesian',
                    text=plot_data['clades'],
                    marker=dict(color='green', size=8, opacity=0.7),
                    hovertemplate='<b>%{text}</b><br>ML: %{x:.3f}<br>Bayesian: %{y:.3f}<extra></extra>'
                ),
                row=2, col=1
            )
        
        # 4. Support vs Clade Size
        if plot_data['clade_sizes'] and plot_data['ml_decay']:
            fig.add_trace(
                go.Scatter(
                    x=plot_data['clade_sizes'],
                    y=plot_data['ml_decay'],
                    mode='markers',
                    name='Support vs Size',
                    text=plot_data['clades'],
                    marker=dict(
                        color=plot_data['au_pvalues'] if plot_data['au_pvalues'] else 'blue',
                        colorscale='RdYlBu',
                        size=10,
                        opacity=0.7,
                        colorbar=dict(title="AU P-value", x=1.05)
                    ),
                    hovertemplate='<b>%{text}</b><br>Size: %{x}<br>ML Decay: %{y:.3f}<extra></extra>'
                ),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            title={
                'text': 'Interactive Phylogenetic Decay Analysis Dashboard',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 18}
            },
            template=self.interactive_config['theme'],
            width=self.interactive_config['width'],
            height=self.interactive_config['height'],
            showlegend=True
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="ML Decay Index", row=1, col=1)
        fig.update_yaxes(title_text="Frequency", row=1, col=1)
        fig.update_xaxes(title_text="AU P-value", row=1, col=2)
        fig.update_yaxes(title_text="Frequency", row=1, col=2)
        fig.update_xaxes(title_text="ML Decay Index", row=2, col=1)
        fig.update_yaxes(title_text="Bayesian Decay Index", row=2, col=1)
        fig.update_xaxes(title_text="Clade Size (# taxa)", row=2, col=2)
        fig.update_yaxes(title_text="ML Decay Index", row=2, col=2)
        
        # Save as HTML
        if self.interactive_config['export_html']:
            html_file = output_dir / "decay_dashboard.html"
            fig.write_html(str(html_file))
            created_files.append(str(html_file))
        
        logger.info(f"Created interactive decay dashboard: {created_files}")
        return created_files
    
    def _create_static_site_analysis(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create static site analysis visualizations."""
        created_files = []
        
        for clade_id, data in decay_data.items():
            if 'site_data' not in data:
                continue
            
            site_data = data['site_data']
            if not site_data:
                continue
            
            # Extract site-specific data
            site_nums = sorted(site_data.keys())
            deltas = [site_data[site]['delta_lnL'] for site in site_nums if 'delta_lnL' in site_data[site]]
            
            if not deltas:
                continue
            
            # Create site analysis plot
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            fig.suptitle(f'Site Analysis for Clade {clade_id}', fontsize=14, fontweight='bold')
            
            # Site-by-site likelihood differences
            ax1.plot(site_nums, deltas, 'o-', markersize=3, linewidth=1, alpha=0.7)
            ax1.axhline(y=0, color='red', linestyle='--', alpha=0.5)
            ax1.set_xlabel('Site Position')
            ax1.set_ylabel('ΔlnL (ML - Constrained)')
            ax1.set_title('Site-specific Likelihood Differences')
            ax1.grid(True, alpha=0.3)
            
            # Histogram of likelihood differences
            ax2.hist(deltas, bins=30, alpha=0.7, color='lightblue', edgecolor='black')
            ax2.axvline(x=0, color='red', linestyle='--', alpha=0.5)
            ax2.set_xlabel('ΔlnL')
            ax2.set_ylabel('Frequency')
            ax2.set_title('Distribution of Likelihood Differences')
            
            plt.tight_layout()
            
            # Save plot
            for fmt in self.static_config['formats']:
                filename = output_dir / f"site_analysis_clade_{clade_id}.{fmt}"
                plt.savefig(filename, format=fmt, dpi=self.static_config['dpi'])
                created_files.append(str(filename))
            
            plt.close()
        
        logger.info(f"Created {len(created_files)} static site analysis plots")
        return created_files
    
    def _create_interactive_site_analysis(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create interactive site analysis visualizations."""
        created_files = []
        
        for clade_id, data in decay_data.items():
            if 'site_data' not in data:
                continue
            
            site_data = data['site_data']
            if not site_data:
                continue
            
            # Extract site-specific data
            site_nums = sorted(site_data.keys())
            deltas = [site_data[site]['delta_lnL'] for site in site_nums if 'delta_lnL' in site_data[site]]
            
            if not deltas:
                continue
            
            # Create interactive site analysis
            fig = make_subplots(
                rows=1, cols=2,
                subplot_titles=['Site-specific Likelihood Differences', 'Distribution of ΔlnL'],
                specs=[[{"type": "scatter"}, {"type": "histogram"}]]
            )
            
            # Site-by-site plot
            fig.add_trace(
                go.Scatter(
                    x=site_nums,
                    y=deltas,
                    mode='lines+markers',
                    name=f'Clade {clade_id}',
                    line=dict(width=2),
                    marker=dict(size=4),
                    hovertemplate='Site: %{x}<br>ΔlnL: %{y:.4f}<extra></extra>'
                ),
                row=1, col=1
            )
            
            # Add zero line
            fig.add_hline(y=0, line_dash="dash", line_color="red", row=1, col=1)
            
            # Histogram
            fig.add_trace(
                go.Histogram(
                    x=deltas,
                    name='ΔlnL Distribution',
                    nbinsx=30,
                    marker_color='lightblue',
                    hovertemplate='ΔlnL: %{x:.4f}<br>Count: %{y}<extra></extra>'
                ),
                row=1, col=2
            )
            
            # Add zero line to histogram
            fig.add_vline(x=0, line_dash="dash", line_color="red", row=1, col=2)
            
            # Update layout
            fig.update_layout(
                title=f'Interactive Site Analysis for Clade {clade_id}',
                template=self.interactive_config['theme'],
                width=self.interactive_config['width'],
                height=600,
                showlegend=True
            )
            
            fig.update_xaxes(title_text="Site Position", row=1, col=1)
            fig.update_yaxes(title_text="ΔlnL (ML - Constrained)", row=1, col=1)
            fig.update_xaxes(title_text="ΔlnL", row=1, col=2)
            fig.update_yaxes(title_text="Frequency", row=1, col=2)
            
            # Save as HTML
            html_file = output_dir / f"site_analysis_clade_{clade_id}.html"
            fig.write_html(str(html_file))
            created_files.append(str(html_file))
        
        logger.info(f"Created {len(created_files)} interactive site analysis plots")
        return created_files
    
    def _create_static_tree_support(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create static tree support visualizations."""
        created_files = []
        
        # Implementation for tree support visualization would go here
        # This would involve creating tree diagrams with support values
        logger.info("Static tree support visualization not yet implemented")
        return created_files
    
    def _create_interactive_tree_support(self, decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
        """Create interactive tree support visualizations."""
        created_files = []
        
        # Implementation for interactive tree support visualization would go here
        # This would involve creating interactive tree diagrams
        logger.info("Interactive tree support visualization not yet implemented")
        return created_files
    
    def _extract_plot_data(self, decay_data: Dict[str, Dict[str, Any]]) -> Dict[str, List]:
        """Extract and organize data for plotting."""
        plot_data = {
            'clades': [],
            'ml_decay': [],
            'bayesian_decay': [],
            'au_pvalues': [],
            'clade_sizes': []
        }
        
        for clade_id, data in decay_data.items():
            plot_data['clades'].append(clade_id)
            
            # ML decay values
            ml_decay = data.get('ml_decay_index', data.get('ml_decay', None))
            plot_data['ml_decay'].append(ml_decay if ml_decay is not None else 0)
            
            # Bayesian decay values
            bayesian_decay = data.get('bayesian_decay_index', data.get('bayesian_decay', None))
            plot_data['bayesian_decay'].append(bayesian_decay if bayesian_decay is not None else 0)
            
            # AU p-values
            au_pvalue = data.get('au_pvalue', data.get('au_p', None))
            plot_data['au_pvalues'].append(au_pvalue if au_pvalue is not None else 1.0)
            
            # Clade sizes
            clade_size = data.get('clade_size', len(data.get('taxa', [])))
            plot_data['clade_sizes'].append(clade_size if clade_size > 0 else 1)
        
        # Filter out None values
        for key in plot_data:
            if key != 'clades':
                plot_data[key] = [x for x in plot_data[key] if x is not None]
        
        return plot_data
    
    def _create_combined_report(self, created_files: Dict[str, List[str]], output_dir: Path) -> str:
        """Create a combined HTML report with both static and interactive visualizations."""
        report_file = output_dir / "visualization_report.html"
        
        html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>panDecay Visualization Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1 { color: #2c3e50; text-align: center; }
        h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
        .section { margin: 30px 0; }
        .file-list { margin: 15px 0; }
        .file-list ul { list-style-type: none; padding: 0; }
        .file-list li { margin: 8px 0; padding: 8px; background: #f8f9fa; border-radius: 4px; }
        .file-list a { color: #2980b9; text-decoration: none; }
        .file-list a:hover { text-decoration: underline; }
        .info { background: #e8f4f8; padding: 15px; border-radius: 5px; margin: 15px 0; }
    </style>
</head>
<body>
    <h1>panDecay Visualization Report</h1>
    
    <div class="info">
        <p><strong>Generated by:</strong> panDecay Dual Visualization System</p>
        <p><strong>Date:</strong> {date}</p>
        <p><strong>Output Format:</strong> Static + Interactive</p>
    </div>
    
    <div class="section">
        <h2>Static Visualizations (Publication Ready)</h2>
        <p>High-resolution static plots suitable for publications and reports.</p>
        <div class="file-list">
            <ul>
""".format(date=__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
        
        # Add static files
        for file_path in created_files.get('static', []):
            filename = Path(file_path).name
            html_content += f'                <li><a href="static_plots/{filename}">{filename}</a></li>\n'
        
        html_content += """
            </ul>
        </div>
    </div>
    
    <div class="section">
        <h2>Interactive Visualizations (Exploratory)</h2>
        <p>Interactive plots for data exploration with zoom, hover, and filtering capabilities.</p>
        <div class="file-list">
            <ul>
"""
        
        # Add interactive files
        for file_path in created_files.get('interactive', []):
            filename = Path(file_path).name
            html_content += f'                <li><a href="interactive_plots/{filename}">{filename}</a></li>\n'
        
        html_content += """
            </ul>
        </div>
    </div>
    
    <div class="section">
        <h2>Usage Guidelines</h2>
        <ul>
            <li><strong>Static plots:</strong> Use for manuscripts, presentations, and formal reports</li>
            <li><strong>Interactive plots:</strong> Use for data exploration, detailed analysis, and presentations</li>
            <li><strong>File formats:</strong> PNG/PDF for static, HTML for interactive</li>
        </ul>
    </div>
</body>
</html>
"""
        
        with open(report_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"Created combined visualization report: {report_file}")
        return str(report_file)


def create_visualizations(decay_data: Dict[str, Dict[str, Any]], 
                         output_dir: Path,
                         viz_format: str = "both",
                         static_config: Optional[Dict[str, Any]] = None,
                         interactive_config: Optional[Dict[str, Any]] = None) -> Dict[str, List[str]]:
    """
    Convenience function to create visualizations.
    
    Args:
        decay_data: Decay analysis results
        output_dir: Output directory
        viz_format: Format ("static", "interactive", "both")
        static_config: Static plot configuration
        interactive_config: Interactive plot configuration
        
    Returns:
        Dictionary of created files by format
    """
    try:
        viz_system = DualVisualizationSystem(
            output_format=viz_format,
            static_config=static_config,
            interactive_config=interactive_config
        )
        
        return viz_system.create_visualizations(decay_data, output_dir, "decay_summary")
        
    except VisualizationError as e:
        logger.error(f"Visualization error: {e}")
        return {'static': [], 'interactive': []}
    except Exception as e:
        logger.error(f"Unexpected visualization error: {e}")
        return {'static': [], 'interactive': []}


# Backward compatibility functions
def create_static_plots(decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
    """Create static plots only (backward compatibility)."""
    result = create_visualizations(decay_data, output_dir, "static")
    return result.get('static', [])


def create_interactive_plots(decay_data: Dict[str, Dict[str, Any]], output_dir: Path) -> List[str]:
    """Create interactive plots only."""
    result = create_visualizations(decay_data, output_dir, "interactive")
    return result.get('interactive', [])