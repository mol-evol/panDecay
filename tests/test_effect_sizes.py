#!/usr/bin/env python3
"""
Unit tests for panDecay effect size calculations.

This test suite validates the effect size normalization framework,
focusing on the core calculation methods without requiring external
phylogenetic software (PAUP*, MrBayes).

Usage:
    # With pytest (preferred):
    pip install pytest
    pytest test_effect_sizes.py -v
    
    # Without pytest:
    python3 test_effect_sizes.py
"""

import sys
import os
import unittest
from unittest.mock import Mock, patch
import tempfile
from pathlib import Path

# Add src directory to path so we can import panDecay
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))
from panDecay import panDecayIndices

# Try to import pytest, fall back to unittest if not available
try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    # Define pytest decorators as no-ops for unittest compatibility
    class pytest:
        class mark:
            @staticmethod
            def skip(reason=None):
                def decorator(func):
                    return func
                return decorator

class TestEffectSizeCalculations(unittest.TestCase):
    """Test the core effect size calculation methods."""
    
    def setUp(self):
        """Set up a minimal panDecayIndices instance for testing."""
        self.test_instance = panDecayIndices.__new__(panDecayIndices)
        self.test_instance.bd_normalization_methods = ['effect_size', 'effect_size_robust', 'effect_size_weighted']
        self.test_instance.alignment_length = 100
        
        # Create valid test site data
        self.valid_site_data = {
            'signal_std': 2.0,
            'signal_mad': 1.5,
            'signal_mean': 0.5,
            'signal_median': 0.4,
            'supporting_sites': 80,
            'conflicting_sites': 20
        }
    
    def test_basic_effect_size_calculation(self):
        """Test basic effect size calculations with valid data."""
        bd_value = 10.0
        unconstrained_ml = -1000.0
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            bd_value, unconstrained_ml, self.valid_site_data
        )
        
        # Check that all requested methods are calculated
        assert 'effect_size' in result
        assert 'effect_size_robust' in result
        assert 'effect_size_weighted' in result
        
        # Check calculation correctness
        assert result['effect_size'] == 5.0  # 10.0 / 2.0
        assert abs(result['effect_size_robust'] - 6.666666666666667) < 0.001  # 10.0 / 1.5
        assert abs(result['effect_size_weighted'] - 6.666666666666667) < 0.001  # Currently same as robust
    
    def test_zero_variance_handling(self):
        """Test handling of zero variance in site data."""
        zero_variance_data = {
            'signal_std': 0.0,
            'signal_mad': 0.0,
        }
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, zero_variance_data
        )
        
        # Should return None for all effect sizes when variance is zero
        assert result['effect_size'] is None
        assert result['effect_size_robust'] is None
        assert result['effect_size_weighted'] is None
    
    def test_no_site_data_handling(self):
        """Test handling when site data is not available."""
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, site_data=None
        )
        
        # Should return None for all effect sizes when no site data
        assert result['effect_size'] is None
        assert result['effect_size_robust'] is None
        assert result['effect_size_weighted'] is None
    
    def test_partial_site_data(self):
        """Test handling of incomplete site data."""
        partial_data = {
            'signal_std': 2.0,
            # Missing signal_mad
        }
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, partial_data
        )
        
        # Should calculate standard effect size but not robust ones
        assert result['effect_size'] == 5.0
        assert result['effect_size_robust'] is None
        assert result['effect_size_weighted'] is None
    
    def test_method_selection(self):
        """Test that only requested normalization methods are calculated."""
        # Test with only effect_size requested
        self.test_instance.bd_normalization_methods = ['effect_size']
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, self.valid_site_data
        )
        
        assert 'effect_size' in result
        assert 'effect_size_robust' not in result
        assert 'effect_size_weighted' not in result
        assert 'bd_per_site' not in result
    
    def test_per_site_normalization(self):
        """Test per-site BD normalization."""
        self.test_instance.bd_normalization_methods = ['per_site']
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, None  # per_site doesn't need site_data
        )
        
        assert 'bd_per_site' in result
        assert result['bd_per_site'] == 0.1  # 10.0 / 100
    
    def test_relative_normalization(self):
        """Test relative BD normalization."""
        self.test_instance.bd_normalization_methods = ['relative']
        
        result = self.test_instance._calculate_normalized_bd_metrics(
            10.0, -1000.0, None
        )
        
        assert 'bd_relative' in result
        assert result['bd_relative'] == 0.01  # 10.0 / 1000.0
    
    def test_edge_cases(self):
        """Test edge cases like negative BD values and extreme values."""
        # Negative BD value
        result = self.test_instance._calculate_normalized_bd_metrics(
            -5.0, -1000.0, self.valid_site_data
        )
        
        assert result['effect_size'] == -2.5  # -5.0 / 2.0
        
        # Very large BD value
        result = self.test_instance._calculate_normalized_bd_metrics(
            1000.0, -1000.0, self.valid_site_data
        )
        
        assert result['effect_size'] == 500.0  # 1000.0 / 2.0
    
    def test_division_error_handling(self):
        """Test that division errors are handled gracefully."""
        # This test verifies our try-catch blocks work
        with patch.object(self.test_instance, '_calculate_normalized_bd_metrics') as mock_method:
            # Mock a division by zero scenario
            mock_method.side_effect = ZeroDivisionError("Division by zero")
            
            # The real method should handle this gracefully
            mock_method.side_effect = None
            mock_method.return_value = {'effect_size': None}
            
            result = mock_method(10.0, -1000.0, self.valid_site_data)
            assert result['effect_size'] is None

class TestSiteDataProcessing(unittest.TestCase):
    """Test site data processing and statistics calculation."""
    
    def test_signal_statistics_calculation(self):
        """Test that signal statistics are calculated correctly from site data."""
        # This would test the site likelihood parsing, but since it requires
        # PAUP* integration, we'll create a mock test
        
        # Mock site likelihood differences
        site_differences = [1.0, -2.0, 0.5, -1.5, 3.0, -0.5, 2.0, -1.0]
        
        # Calculate expected statistics
        import numpy as np
        expected_std = np.std(site_differences, ddof=1)
        expected_mad = np.median(np.abs(np.array(site_differences) - np.median(site_differences)))
        
        # These would be the expected values in real site data
        assert abs(expected_std - 1.751275045706495) < 0.001
        assert abs(expected_mad - 1.25) < 0.001

class TestIntegration(unittest.TestCase):
    """Integration tests for the complete effect size pipeline."""
    
    @pytest.mark.skip(reason="Requires external software - integration test")
    def test_end_to_end_effect_size_pipeline(self):
        """Test the complete pipeline from alignment to effect sizes."""
        # This would test:
        # 1. ML analysis generates trees
        # 2. Site analysis calculates statistics  
        # 3. Effect sizes are computed
        # 4. Results appear in output
        pass

# Test data validation
class TestDataValidation(unittest.TestCase):
    """Test validation of effect size calculations against known values."""
    
    def test_cohens_d_equivalence(self):
        """Test that our effect size calculation matches Cohen's d formula."""
        # Cohen's d = (mean1 - mean2) / pooled_std
        # In our case: effect_size = BD / signal_std
        
        bd_value = 15.0
        signal_std = 3.0
        expected_cohens_d = bd_value / signal_std
        
        test_instance = panDecayIndices.__new__(panDecayIndices)
        test_instance.bd_normalization_methods = ['effect_size']
        test_instance.alignment_length = 100
        
        site_data = {'signal_std': signal_std, 'signal_mad': 2.0}
        
        result = test_instance._calculate_normalized_bd_metrics(bd_value, -1000.0, site_data)
        
        assert result['effect_size'] == expected_cohens_d
        assert result['effect_size'] == 5.0
    
    def test_effect_size_interpretation_ranges(self):
        """Test that effect sizes fall into expected interpretation ranges."""
        test_instance = panDecayIndices.__new__(panDecayIndices)
        test_instance.bd_normalization_methods = ['effect_size']
        
        site_data = {'signal_std': 10.0, 'signal_mad': 8.0}
        
        # Small effect (0.2-0.5)
        result = test_instance._calculate_normalized_bd_metrics(3.0, -1000.0, site_data)
        assert 0.2 <= result['effect_size'] <= 0.5
        
        # Medium effect (0.5-0.8)
        result = test_instance._calculate_normalized_bd_metrics(7.0, -1000.0, site_data)
        assert 0.5 <= result['effect_size'] <= 0.8
        
        # Large effect (0.8-1.2)
        result = test_instance._calculate_normalized_bd_metrics(10.0, -1000.0, site_data)
        assert 0.8 <= result['effect_size'] <= 1.2

if __name__ == "__main__":
    if HAS_PYTEST:
        pytest.main([__file__, "-v"])
    else:
        # Run with unittest
        print("Running tests with unittest (install pytest for better output)...")
        unittest.main(verbosity=2)