#!/usr/bin/env python3
"""
panDecay entry point script.
This script imports and runs the main panDecay module from the src directory.
"""

import sys
import os

# Add src directory to path to import panDecay modules
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

# Import and run the main panDecay module
from panDecay import main

if __name__ == "__main__":
    main()