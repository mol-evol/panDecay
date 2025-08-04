#!/usr/bin/env python3
"""
Utility functions for panDecay.

Common utility functions used across the panDecay modules.
"""

import argparse
import re
from pathlib import Path
from typing import Union

from pandecay.core.constants import VERSION


def get_display_path(path: Union[str, Path]) -> str:
    """Get a display-friendly path representation."""
    if isinstance(path, Path):
        path_obj = path
    else:
        path_obj = Path(path)
    
    if path_obj.is_absolute():
        try:
            # Try to get relative path from current directory
            rel_path = path_obj.relative_to(Path.cwd())
            # Use relative path if it's shorter and doesn't go up too many levels
            if len(str(rel_path)) < len(str(path_obj)) and not str(rel_path).startswith('../../../'):
                return str(rel_path)
        except ValueError:
            pass  # Path is not relative to current directory
    
    return str(path_obj)


def print_runtime_parameters(args_ns: argparse.Namespace, model_str_for_print: str) -> None:
    """Print runtime parameters in a clean banner format."""
    
    # Fixed width for consistency (80 characters)
    banner_width = 80
    separator = "═" * banner_width
    
    # Banner header with citation
    print(separator)
    print(f"{'panDecay v' + VERSION:^80}")
    print(f"{'McInerney, J.O. (2025) panDecay: Phylogenetic Analysis with Decay Indices':^80}")
    print(f"{'http://github.com/mol-evol/panDecay/':^80}")
    print(separator)
    
    # Analysis Configuration Section
    print("\nANALYSIS CONFIGURATION")
    alignment_info = f"{get_display_path(args_ns.alignment)} ({args_ns.format}, {args_ns.data_type})"
    print(f"  Alignment: {alignment_info}")
    print(f"  Model: {model_str_for_print}")
    print(f"  Mode: {args_ns.analysis}")
    
    # Advanced Model Parameters (if relevant)
    model_details = []
    if hasattr(args_ns, 'gamma') and args_ns.gamma:
        model_details.append("Gamma rate variation")
    if hasattr(args_ns, 'invariable') and args_ns.invariable:
        model_details.append("Invariable sites")
    if hasattr(args_ns, 'base_freq') and args_ns.base_freq and args_ns.base_freq != 'empirical':
        model_details.append(f"Base freq: {args_ns.base_freq}")
    
    if model_details:
        print(f"  Advanced: {', '.join(model_details)}")
    
    # Runtime Settings Section
    print("\nRUNTIME SETTINGS")
    print(f"  Threads: {args_ns.threads}")
    
    if hasattr(args_ns, 'project_name') and args_ns.project_name:
        print(f"  Project: {args_ns.project_name}")
    
    if hasattr(args_ns, 'output_dir') and args_ns.output_dir:
        print(f"  Output: {get_display_path(args_ns.output_dir)}")
    else:
        print(f"  Output: {get_display_path(args_ns.output)}")
    
    # Runtime flags and special settings
    special_settings = []
    if hasattr(args_ns, 'temp') and args_ns.temp:
        special_settings.append(f"Custom temp: {get_display_path(args_ns.temp)}")
    if hasattr(args_ns, 'debug') and args_ns.debug:
        special_settings.append("Debug mode")
    if hasattr(args_ns, 'keep_files') and args_ns.keep_files:
        special_settings.append("Keep temp files")
    
    if special_settings:
        for setting in special_settings:
            print(f"  {setting}")
    
    # Analysis Options Section
    analysis_options = []
    if args_ns.starting_tree:
        analysis_options.append(f"Starting tree: {get_display_path(args_ns.starting_tree)}")
    if args_ns.site_analysis:
        analysis_options.append("Site analysis: Enabled")
    if args_ns.bootstrap:
        analysis_options.append(f"Bootstrap: {args_ns.bootstrap_reps} replicates")
    if hasattr(args_ns, 'visualize') and args_ns.visualize:
        analysis_options.append(f"Visualization: {args_ns.viz_format.upper()} format")
    
    if analysis_options:
        print("\nANALYSIS OPTIONS")
        for option in analysis_options:
            print(f"  {option}")
    
    # Bayesian Analysis Section (if bayesian analysis is included)
    if 'bayesian' in args_ns.analysis or args_ns.analysis == 'all':
        print("\nBAYESIAN PARAMETERS")
        print(f"  MCMC generations: {args_ns.bayes_ngen:,}")
        
        # Combine burnin and chains on one line
        print(f"  Burnin: {args_ns.bayes_burnin:.2f}, Chains: {args_ns.bayes_chains}")
        
        ml_method_names = {
            'ss': 'Stepping-stone sampling',
            'ps': 'Path sampling', 
            'hm': 'Harmonic mean'
        }
        ml_method = ml_method_names.get(args_ns.marginal_likelihood, args_ns.marginal_likelihood)
        
        if args_ns.marginal_likelihood == 'ss':
            print(f"  Marginal likelihood: {ml_method} ({args_ns.ss_nsteps} steps, α={args_ns.ss_alpha:.2f})")
        else:
            print(f"  Marginal likelihood: {ml_method}")
    
    # Performance Settings Section (if any performance options are enabled)
    performance_options = []
    if hasattr(args_ns, 'use_mpi') and args_ns.use_mpi:
        mpi_procs = getattr(args_ns, 'mpi_processors', args_ns.bayes_chains)
        performance_options.append(f"MPI: {mpi_procs} processors")
    
    if hasattr(args_ns, 'use_beagle') and args_ns.use_beagle:
        beagle_device = getattr(args_ns, 'beagle_device', 'auto')
        performance_options.append(f"BEAGLE: {beagle_device}")
    
    if performance_options:
        print("\nPERFORMANCE")
        for option in performance_options:
            print(f"  {option}")
    
    # Banner footer
    print(separator)


def format_tree_annotation(clade_id: str, annotation_dict: dict, style: str = "compact") -> str:
    """Format tree annotation with support values and other metrics."""
    if style == "compact":
        # Compact format: AU=0.95;DeltaLnL=-1.23
        parts = []
        if "au_pvalue" in annotation_dict:
            au_val = annotation_dict["au_pvalue"]
            if isinstance(au_val, (int, float)):
                parts.append(f"AU={au_val:.3f}")
            else:
                parts.append(f"AU={au_val}")
        
        if "delta_lnl" in annotation_dict:
            delta_val = annotation_dict["delta_lnl"]
            if isinstance(delta_val, (int, float)):
                parts.append(f"DeltaLnL={delta_val:.3f}")
            else:
                parts.append(f"DeltaLnL={delta_val}")
        
        if "bootstrap_support" in annotation_dict:
            bs_val = annotation_dict["bootstrap_support"]
            if isinstance(bs_val, (int, float)):
                parts.append(f"BS={bs_val:.1f}")
            else:
                parts.append(f"BS={bs_val}")
        
        return ";".join(parts)
    
    elif style == "verbose":
        # Verbose format with labels
        parts = []
        if "au_pvalue" in annotation_dict:
            au_val = annotation_dict["au_pvalue"]
            parts.append(f"AU_pvalue={au_val}")
        
        if "delta_lnl" in annotation_dict:
            delta_val = annotation_dict["delta_lnl"]
            parts.append(f"Delta_LnL={delta_val}")
        
        if "bootstrap_support" in annotation_dict:
            bs_val = annotation_dict["bootstrap_support"]
            parts.append(f"Bootstrap_support={bs_val}")
        
        return "; ".join(parts)
    
    else:
        # Default to compact if unknown style
        return format_tree_annotation(clade_id, annotation_dict, "compact")


def validate_file_path(path: Union[str, Path], must_exist: bool = True) -> Path:
    """Validate and convert a file path, optionally checking existence."""
    path_obj = Path(path)
    
    if must_exist and not path_obj.exists():
        raise FileNotFoundError(f"File does not exist: {path}")
    
    return path_obj


def truncate_string(text: str, max_length: int = 50, suffix: str = "...") -> str:
    """Truncate a string to a maximum length with optional suffix."""
    if len(text) <= max_length:
        return text
    
    truncate_length = max_length - len(suffix)
    if truncate_length <= 0:
        return suffix[:max_length]
    
    return text[:truncate_length] + suffix


def format_taxon_for_paup(taxon_name: str) -> str:
    """Format a taxon name for PAUP* (handles spaces, special chars by quoting)."""
    if not isinstance(taxon_name, str): 
        taxon_name = str(taxon_name)
    
    # PAUP* needs quotes if name contains whitespace or NEXUS special chars
    special_chars = r'[\s\(\)\[\]\{\}/\\,;=\*`"\'<>]'
    if re.search(special_chars, taxon_name) or ':' in taxon_name:
        # Replace single quotes with underscores and wrap in quotes
        clean_name = taxon_name.replace("'", "_")
        return f"'{clean_name}'"
    
    return taxon_name


def format_support_symbol(pvalue: float) -> str:
    """Format support value with appropriate symbol or text."""
    try:
        if pvalue is None:
            return 'ns'
        elif pvalue >= 0.95:
            return f'{pvalue:.3f}'
        elif pvalue >= 0.05:
            return f'{pvalue:.3f}'
        else:  # p < 0.05, significant
            p_val = float(pvalue)
            if p_val < 0.001:
                return '<0.001*'
            else:
                return f'{p_val:.3f}*'
    except:
        return 'N/A'


def build_effective_model_display(args: argparse.Namespace) -> str:
    """Build an accurate model display string that reflects parameter overrides."""
    # Start with base model
    base_model = args.model.split("+")[0].upper()
    has_gamma = args.gamma or "+G" in args.model.upper()
    has_invar = args.invariable or "+I" in args.model.upper()
    
    # Collect override indicators
    overrides = []
    effective_model = base_model
    
    # Handle data type specific overrides
    if args.data_type == "dna":
        # NST override - this is the key fix for the reported issue
        if args.nst is not None:
            if args.nst == 1:
                effective_model = "JC"  # JC-like model
            elif args.nst == 2:
                effective_model = "HKY"  # HKY-like model  
            elif args.nst == 6:
                effective_model = "GTR"  # GTR-like model
            overrides.append(f"nst={args.nst}")
        
        # Base frequency override
        if args.base_freq:
            overrides.append(f"basefreq={args.base_freq}")
        
        # Rates override (overrides gamma flag)
        if args.rates:
            overrides.append(f"rates={args.rates}")
            if args.rates == "gamma":
                has_gamma = True
            elif args.rates == "equal":
                has_gamma = False
    
    elif args.data_type == "protein":
        # Protein model override
        if args.protein_model:
            effective_model = args.protein_model.upper()
            overrides.append(f"protein={args.protein_model}")
        elif base_model not in ["JTT", "WAG", "LG", "DAYHOFF", "MTREV", "CPREV", "BLOSUM62", "HIVB", "HIVW"]:
            # Generic protein model defaults to JTT in analysis
            effective_model = "JTT"
            overrides.append("protein=jtt")
    
    elif args.data_type == "discrete":
        effective_model = "Mk"
        overrides.append("nst=1")
        if args.base_freq:
            overrides.append(f"basefreq={args.base_freq}")
        else:
            overrides.append("basefreq=equal")
        
        # Parsimony model override
        if args.parsmodel is not None:
            overrides.append(f"parsmodel={str(args.parsmodel).lower()}")
    
    # Add gamma and invariable sites
    if has_gamma:
        effective_model += "+G"
        if args.gamma_shape is not None:
            overrides.append(f"shape={args.gamma_shape}")
    
    if has_invar:
        effective_model += "+I"
        if args.prop_invar is not None:
            overrides.append(f"pinvar={args.prop_invar}")
    
    # Build display string with override indicators
    if overrides:
        override_str = ", ".join(overrides)
        return f"{effective_model} ({override_str})"
    else:
        return effective_model


class ProgressIndicator:
    """
    Progress indicator with spinner animations for cleaner output.
    
    Provides visual feedback for long-running operations without verbose logging.
    """
    
    def __init__(self, style: str = "unicode"):
        """
        Initialize progress indicator.
        
        Args:
            style: Style for spinner ("unicode" or "ascii")
        """
        import sys
        import time
        import threading
        
        self.sys = sys
        self.time = time
        self.threading = threading
        
        if style == "unicode":
            self.spinner_chars = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        else:
            self.spinner_chars = ["|", "/", "-", "\\"]
        
        self.current_message = ""
        self.is_spinning = False
        self.spinner_thread = None
        self.current_char_index = 0
    
    def start(self, message: str):
        """
        Start spinner with message.
        
        Args:
            message: Status message to display
        """
        if self.is_spinning:
            self.stop()
        
        self.current_message = message
        self.is_spinning = True
        self.current_char_index = 0
        
        self.spinner_thread = self.threading.Thread(target=self._spin, daemon=True)
        self.spinner_thread.start()
    
    def update(self, message: str):
        """
        Update the spinner message.
        
        Args:
            message: New status message
        """
        self.current_message = message
    
    def stop(self, final_message: str = None):
        """
        Stop spinner and optionally show final message.
        
        Args:
            final_message: Optional final message to display
        """
        if not self.is_spinning:
            return
        
        self.is_spinning = False
        
        if self.spinner_thread and self.spinner_thread.is_alive():
            self.spinner_thread.join(timeout=0.1)
        
        # Clear current line completely
        try:
            # Try to get terminal width for complete line clearing
            import os
            terminal_width = os.get_terminal_size().columns
            self.sys.stdout.write('\r' + ' ' * terminal_width + '\r')
        except (OSError, AttributeError):
            # Fallback: clear based on current message length + safety margin
            full_line_length = len(self.current_message) + 10
            self.sys.stdout.write('\r' + ' ' * full_line_length + '\r')
        
        self.sys.stdout.flush()
        
        if final_message:
            print(final_message)
    
    def _spin(self):
        """Internal method to handle spinner animation."""
        while self.is_spinning:
            char = self.spinner_chars[self.current_char_index]
            self.sys.stdout.write(f'\r{char} {self.current_message}')
            self.sys.stdout.flush()
            
            self.current_char_index = (self.current_char_index + 1) % len(self.spinner_chars)
            self.time.sleep(0.1)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()


class OverwritingProgress:
    """
    Progress indicator that overwrites the same line for repetitive operations.
    
    Shows current item being processed and overwrites on the same line,
    then displays a final summary when complete.
    """
    
    def __init__(self):
        """Initialize overwriting progress indicator."""
        import sys
        self.sys = sys
        self.is_active = False
        self.current_line_length = 0
    
    def update(self, message: str):
        """
        Update the current progress message (overwrites previous line).
        
        Args:
            message: Current status message
        """
        # Clear previous line if it was longer
        if self.current_line_length > 0:
            self.sys.stdout.write('\r' + ' ' * self.current_line_length + '\r')
        
        # Write new message
        self.sys.stdout.write(f'\r{message}')
        self.sys.stdout.flush()
        
        self.current_line_length = len(message)
        self.is_active = True
    
    def finish(self, final_message: str):
        """
        Complete the progress and show final summary.
        
        Args:
            final_message: Final summary message to display
        """
        if self.is_active:
            # Clear current line
            self.sys.stdout.write('\r' + ' ' * self.current_line_length + '\r')
            self.sys.stdout.flush()
        
        # Print final message
        print(final_message)
        
        self.is_active = False
        self.current_line_length = 0
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.is_active:
            self.finish("")


def show_progress(message: str, style: str = "unicode"):
    """
    Convenience function to show progress with spinner.
    
    Args:
        message: Status message to display
        style: Spinner style ("unicode" or "ascii")
        
    Returns:
        ProgressIndicator instance
    """
    indicator = ProgressIndicator(style)
    indicator.start(message)
    return indicator