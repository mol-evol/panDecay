#!/usr/bin/env python3
"""
Analysis coordinator for panDecay.

This module coordinates the different analysis engines and manages
the overall analysis workflow.
"""

import logging
from typing import Dict, List, Optional, Any
from pathlib import Path

logger = logging.getLogger(__name__)


class AnalysisCoordinator:
    """
    Coordinates multiple analysis engines and manages workflow.
    """
    
    def __init__(self, debug: bool = False):
        """
        Initialize analysis coordinator.
        
        Args:
            debug: Enable debug logging
        """
        self.debug = debug
        self.engines = {}
        
    def register_engine(self, name: str, engine):
        """
        Register an analysis engine.
        
        Args:
            name: Engine name
            engine: Analysis engine instance
        """
        self.engines[name] = engine
        
    def run_analysis(self, config) -> Dict[str, Any]:
        """
        Run coordinated analysis across all registered engines.
        
        Args:
            config: Analysis configuration
            
        Returns:
            Dictionary of results from all engines
        """
        results = {}
        
        for name, engine in self.engines.items():
            try:
                logger.debug(f"Running {name} analysis")
                results[name] = engine.run(config)
            except Exception as e:
                logger.error(f"Failed to run {name} analysis: {e}")
                results[name] = None
                
        return results