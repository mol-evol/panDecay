#!/usr/bin/env python3
"""
Utility functions for panDecay.

Common utility functions used across the panDecay modules.
"""

import argparse
import re
from pathlib import Path
from typing import Union

from .constants import VERSION


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
    """Print runtime parameters in a formatted box."""
    print("┌─" + "─" * 76 + "─┐")
    print("│" + f" panDecay v{VERSION} - Runtime Parameters".center(76) + " │")
    print("├─" + "─" * 76 + "─┤")
    print(f"│ {'Alignment:':<20} {get_display_path(args_ns.alignment):<54} │")
    print(f"│ {'Format:':<20} {args_ns.format:<54} │")
    print(f"│ {'Data type:':<20} {args_ns.data_type:<54} │")
    print(f"│ {'Model:':<20} {model_str_for_print:<54} │")
    print(f"│ {'Analysis mode:':<20} {args_ns.analysis:<54} │")
    print(f"│ {'Threads:':<20} {str(args_ns.threads):<54} │")
    print(f"│ {'Output file:':<20} {get_display_path(args_ns.output):<54} │")
    
    if args_ns.starting_tree:
        print(f"│ {'Starting tree:':<20} {get_display_path(args_ns.starting_tree):<54} │")
    
    if args_ns.site_analysis:
        print(f"│ {'Site analysis:':<20} {'Enabled':<54} │")
    
    if args_ns.bootstrap:
        print(f"│ {'Bootstrap:':<20} {f'{args_ns.bootstrap_reps} replicates':<54} │")
    
    if args_ns.visualize:
        print(f"│ {'Visualization:':<20} {f'{args_ns.viz_format.upper()} format':<54} │")
    
    print("└─" + "─" * 76 + "─┘")


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