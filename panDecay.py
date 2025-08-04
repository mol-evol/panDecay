#!/usr/bin/env python3
"""
panDecay - Phylogenetic Decay Indices Calculator

DEPRECATION WARNING: This script is deprecated. Please use the new package installation:

    pip install pandecay
    pandecay [arguments...]

This wrapper is provided for backward compatibility and will be removed in a future version.
"""

import sys
import warnings
from pathlib import Path

# Show deprecation warning
warnings.warn(
    "panDecay.py is deprecated. Please install with 'pip install pandecay' and use 'pandecay' command instead.",
    DeprecationWarning,
    stacklevel=2
)

# Add the current directory to Python path for package import
sys.path.insert(0, str(Path(__file__).parent))

try:
    from pandecay.cli import main
    if __name__ == "__main__":
        main()
except ImportError as e:
    print("Error: Could not import pandecay package.", file=sys.stderr)
    print("Please install pandecay with: pip install pandecay", file=sys.stderr)
    print(f"Import error: {e}", file=sys.stderr)
    sys.exit(1)