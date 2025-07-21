#!/usr/bin/env python3
"""
Utility functions for panDecay core operations.
"""

import re
from typing import List


def format_taxon_for_paup(taxon_name: str) -> str:
    """
    Format taxon name for PAUP* compatibility.
    
    Args:
        taxon_name: Original taxon name
        
    Returns:
        PAUP*-compatible taxon name
    """
    # Replace spaces and special characters
    formatted = re.sub(r'[^\w]', '_', taxon_name)
    return formatted


def get_representative_taxa_sample(taxa_list: List[str], max_display: int = 3) -> str:
    """
    Get a representative sample of taxa for display purposes.
    
    Args:
        taxa_list: List of taxa names
        max_display: Maximum number of taxa to display
        
    Returns:
        Formatted string with representative taxa
    """
    if len(taxa_list) <= max_display:
        return ", ".join(taxa_list)
    else:
        sample = taxa_list[:max_display]
        return ", ".join(sample) + f"... (and {len(taxa_list) - max_display} more)"