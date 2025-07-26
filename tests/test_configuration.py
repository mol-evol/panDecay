"""
Comprehensive tests for configuration management.

Tests the centralized constants, configuration loading, and
validation functionality.
"""

import pytest
import tempfile
import shutil
from pathlib import Path

from src.config.constants import (
    FileNames, AnalysisThresholds, AnalysisTimeouts, ResourceLimits,
    ExternalTools, ModelDefaults, AnalysisModes, OutputFormats,
    ValidationLimits, StringPatterns, get_default_timeouts,
    get_default_resource_limits, get_file_naming_config
)


class TestFileNames:
    """Test the FileNames constants class."""
    
    def test_basic_file_names(self):
        """Test basic file name constants."""
        assert FileNames.NEXUS_ALIGNMENT == "alignment.nex"
        assert FileNames.ML_TREE == "ml_tree.tre"
        assert FileNames.ML_SCORE == "ml_score.txt"
        assert FileNames.ML_LOG == "paup_ml.log"
    
    def test_au_test_file_names(self):
        """Test AU test related file names."""
        assert FileNames.AU_TEST_NEX == "au_test.nex"
        assert FileNames.AU_TEST_LOG == "paup_au.log"
        assert FileNames.AU_TEST_RESULTS == "au_test_results.txt"
        assert FileNames.ALL_TREE_SCORES == "all_tree_scores.txt"
    
    def test_constraint_file_names(self):
        """Test constraint related file names."""
        assert FileNames.CONSTRAINT_PREFIX == "constraint"
        assert FileNames.CONSTRAINT_TREE_PREFIX == "constraint_tree"
        assert FileNames.CONSTRAINT_SCORE_PREFIX == "constraint_score"
    
    def test_bayesian_file_names(self):
        """Test Bayesian analysis file names."""
        assert FileNames.BAYESIAN_LOG == "mrbayes.log"
        assert FileNames.BAYESIAN_PSTAT == "bayesian.pstat"
        assert FileNames.BAYESIAN_TSTAT == "bayesian.tstat"
        assert FileNames.BAYESIAN_CON_TREE == "bayesian.con.tre"
    
    def test_parsimony_file_names(self):
        """Test parsimony analysis file names."""
        assert FileNames.PARSIMONY_LOG == "paup_parsimony.log"
        assert FileNames.PARSIMONY_TREES == "parsimony_trees.tre"
        assert FileNames.PARSIMONY_SCORE == "parsimony_score.txt"
        assert FileNames.CONSTRAINT_PARSIMONY_SCORE == "constraint_parsimony_score"


class TestAnalysisThresholds:
    """Test the AnalysisThresholds constants class."""
    
    def test_statistical_thresholds(self):
        """Test statistical threshold constants."""
        assert AnalysisThresholds.AU_SIGNIFICANCE == 0.05
        assert AnalysisThresholds.BURNIN_FRACTION == 0.25
        assert AnalysisThresholds.MIN_ESS == 200
        assert AnalysisThresholds.MAX_PSRF == 1.01
        assert AnalysisThresholds.MAX_ASDSF == 0.01
    
    def test_stepping_stone_parameters(self):
        """Test stepping stone sampling parameters."""
        assert AnalysisThresholds.STEPPING_STONE_ALPHA == 0.4
        assert AnalysisThresholds.STEPPING_STONE_STEPS == 50
    
    def test_parsimony_threshold(self):
        """Test parsimony-specific threshold."""
        assert AnalysisThresholds.PARSIMONY_DECAY_THRESHOLD == 1


class TestAnalysisTimeouts:
    """Test the AnalysisTimeouts dataclass."""
    
    def test_default_timeouts(self):
        """Test default timeout values."""
        timeouts = AnalysisTimeouts()
        
        assert timeouts.ml_timeout == 3600  # 1 hour
        assert timeouts.constraint_timeout == 600  # 10 minutes
        assert timeouts.site_analysis_timeout == 600  # 10 minutes
        assert timeouts.bayesian_timeout == 7200  # 2 hours
        assert timeouts.parsimony_timeout == 1800  # 30 minutes
    
    def test_custom_timeouts(self):
        """Test creating custom timeout values."""
        timeouts = AnalysisTimeouts(
            ml_timeout=7200,
            constraint_timeout=1200,
            bayesian_timeout=14400
        )
        
        assert timeouts.ml_timeout == 7200
        assert timeouts.constraint_timeout == 1200
        assert timeouts.bayesian_timeout == 14400
        # Defaults should still apply
        assert timeouts.site_analysis_timeout == 600
        assert timeouts.parsimony_timeout == 1800


class TestResourceLimits:
    """Test the ResourceLimits dataclass."""
    
    def test_default_resource_limits(self):
        """Test default resource limit values."""
        limits = ResourceLimits()
        
        assert limits.max_memory_mb == 1000
        assert limits.max_processes == 4
        assert limits.chunk_size == 8192
        assert limits.max_alignment_size_mb == 100
        assert limits.max_tree_count == 1000
    
    def test_custom_resource_limits(self):
        """Test creating custom resource limit values."""
        limits = ResourceLimits(
            max_memory_mb=2000,
            max_processes=8,
            chunk_size=16384
        )
        
        assert limits.max_memory_mb == 2000
        assert limits.max_processes == 8
        assert limits.chunk_size == 16384
        # Defaults should still apply
        assert limits.max_alignment_size_mb == 100
        assert limits.max_tree_count == 1000


class TestExternalTools:
    """Test the ExternalTools constants class."""
    
    def test_tool_defaults(self):
        """Test external tool default values."""
        assert ExternalTools.PAUP_DEFAULT == "paup"
        assert ExternalTools.MRBAYES_DEFAULT == "mb"
        assert ExternalTools.MPIRUN_DEFAULT == "mpirun"


class TestModelDefaults:
    """Test the ModelDefaults constants class."""
    
    def test_model_defaults(self):
        """Test model default values."""
        assert ModelDefaults.DNA_MODEL == "GTR"
        assert ModelDefaults.PROTEIN_MODEL == "WAG"
        assert ModelDefaults.DISCRETE_MODEL == "Mk"
        assert ModelDefaults.NST_DEFAULT == 2
        assert ModelDefaults.THREADS_AUTO == "auto"
        assert ModelDefaults.THREADS_ALL == "all"


class TestOutputFormats:
    """Test the OutputFormats enum classes."""
    
    def test_alignment_formats(self):
        """Test alignment format enums."""
        assert OutputFormats.Alignment.FASTA.value == "fasta"
        assert OutputFormats.Alignment.PHYLIP.value == "phylip"
        assert OutputFormats.Alignment.NEXUS.value == "nexus"
        assert OutputFormats.Alignment.CLUSTAL.value == "clustal"
        assert OutputFormats.Alignment.STOCKHOLM.value == "stockholm"
    
    def test_visualization_formats(self):
        """Test visualization format enums."""
        assert OutputFormats.Visualization.PNG.value == "png"
        assert OutputFormats.Visualization.PDF.value == "pdf"
        assert OutputFormats.Visualization.SVG.value == "svg"
    
    def test_tree_formats(self):
        """Test tree format enums."""
        assert OutputFormats.Tree.NEWICK.value == "newick"
        assert OutputFormats.Tree.NEXUS.value == "nexus"


class TestAnalysisModes:
    """Test the AnalysisModes constants class."""
    
    def test_single_analysis_modes(self):
        """Test single analysis mode constants."""
        assert AnalysisModes.ML == "ml"
        assert AnalysisModes.BAYESIAN == "bayesian"
        assert AnalysisModes.PARSIMONY == "parsimony"
    
    def test_combined_analysis_modes(self):
        """Test combined analysis mode constants."""
        assert AnalysisModes.ML_BAYESIAN == "ml+bayesian"
        assert AnalysisModes.ML_PARSIMONY == "ml+parsimony"
        assert AnalysisModes.BAYESIAN_PARSIMONY == "bayesian+parsimony"
        assert AnalysisModes.ALL == "all"


class TestValidationLimits:
    """Test the ValidationLimits constants class."""
    
    def test_sequence_limits(self):
        """Test sequence validation limits."""
        assert ValidationLimits.MIN_SEQUENCES == 3
        assert ValidationLimits.MAX_SEQUENCES == 10000
        assert ValidationLimits.MIN_SITES == 1
        assert ValidationLimits.MAX_SITES == 100000
    
    def test_analysis_limits(self):
        """Test analysis validation limits."""
        assert ValidationLimits.MIN_BOOTSTRAP_REPS == 10
        assert ValidationLimits.MAX_BOOTSTRAP_REPS == 10000
        assert ValidationLimits.MIN_THREADS == 1
        assert ValidationLimits.MAX_THREADS == 64


class TestStringPatterns:
    """Test the StringPatterns constants class."""
    
    def test_nexus_patterns(self):
        """Test NEXUS format patterns."""
        assert StringPatterns.NEXUS_HEADER == "#NEXUS"
        assert StringPatterns.PAUP_BEGIN == "begin paup;"
        assert StringPatterns.PAUP_END == "quit;\nend;"
        assert StringPatterns.MRBAYES_BEGIN == "begin mrbayes;"
        assert StringPatterns.MRBAYES_END == "end;"
    
    def test_nexus_template(self):
        """Test NEXUS template format."""
        template = StringPatterns.NEXUS_TEMPLATE
        
        assert "#NEXUS" in template
        assert "BEGIN DATA;" in template
        assert "DIMENSIONS NTAX={ntax} NCHAR={nchar};" in template
        assert "FORMAT DATATYPE={datatype}" in template
        assert "{matrix}" in template
        assert "END;" in template


class TestHelperFunctions:
    """Test helper functions for getting default configurations."""
    
    def test_get_default_timeouts(self):
        """Test getting default timeout configuration."""
        timeouts = get_default_timeouts()
        
        assert isinstance(timeouts, AnalysisTimeouts)
        assert timeouts.ml_timeout == 3600
        assert timeouts.bayesian_timeout == 7200
    
    def test_get_default_resource_limits(self):
        """Test getting default resource limits."""
        limits = get_default_resource_limits()
        
        assert isinstance(limits, ResourceLimits)
        assert limits.max_memory_mb == 1000
        assert limits.max_processes == 4
    
    def test_get_file_naming_config(self):
        """Test getting file naming configuration as dictionary."""
        config = get_file_naming_config()
        
        assert isinstance(config, dict)
        assert config['nexus_alignment'] == FileNames.NEXUS_ALIGNMENT
        assert config['ml_tree'] == FileNames.ML_TREE
        assert config['ml_score'] == FileNames.ML_SCORE
        assert config['au_test_nex'] == FileNames.AU_TEST_NEX
        
        # Test that all expected keys are present
        expected_keys = [
            'nexus_alignment', 'ml_tree', 'ml_score', 'ml_log',
            'au_test_nex', 'au_test_log', 'au_test_results', 'all_tree_scores'
        ]
        for key in expected_keys:
            assert key in config


class TestConstantsImmutability:
    """Test that constants cannot be accidentally modified."""
    
    def test_file_names_immutable(self):
        """Test that FileNames attributes cannot be modified."""
        original_value = FileNames.NEXUS_ALIGNMENT
        
        # Attempting to modify should not affect the original
        try:
            FileNames.NEXUS_ALIGNMENT = "modified.nex"
        except (AttributeError, TypeError):
            pass  # Expected behavior
        
        # Value should remain unchanged
        assert FileNames.NEXUS_ALIGNMENT == original_value
    
    def test_thresholds_immutable(self):
        """Test that AnalysisThresholds attributes cannot be modified."""
        original_value = AnalysisThresholds.AU_SIGNIFICANCE
        
        # Attempting to modify should not affect the original
        try:
            AnalysisThresholds.AU_SIGNIFICANCE = 0.01
        except (AttributeError, TypeError):
            pass  # Expected behavior
        
        # Value should remain unchanged
        assert AnalysisThresholds.AU_SIGNIFICANCE == original_value


class TestConstantsIntegration:
    """Test integration between different constant classes."""
    
    def test_timeout_consistency(self):
        """Test that timeout values are consistent with resource limits."""
        timeouts = AnalysisTimeouts()
        limits = ResourceLimits()
        
        # All timeouts should be positive
        assert timeouts.ml_timeout > 0
        assert timeouts.constraint_timeout > 0
        assert timeouts.bayesian_timeout > 0
        assert timeouts.parsimony_timeout > 0
        
        # Bayesian timeout should be longest (most computationally intensive)
        assert timeouts.bayesian_timeout >= timeouts.ml_timeout
        assert timeouts.bayesian_timeout >= timeouts.parsimony_timeout
    
    def test_validation_limits_consistency(self):
        """Test that validation limits are logically consistent."""
        # Minimum values should be less than maximum values
        assert ValidationLimits.MIN_SEQUENCES < ValidationLimits.MAX_SEQUENCES
        assert ValidationLimits.MIN_SITES < ValidationLimits.MAX_SITES
        assert ValidationLimits.MIN_BOOTSTRAP_REPS < ValidationLimits.MAX_BOOTSTRAP_REPS
        assert ValidationLimits.MIN_THREADS < ValidationLimits.MAX_THREADS
        
        # Minimum sequences should be at least 3 (required for phylogeny)
        assert ValidationLimits.MIN_SEQUENCES >= 3
        
        # Thread limits should be reasonable
        assert ValidationLimits.MIN_THREADS >= 1
        assert ValidationLimits.MAX_THREADS <= 64
    
    def test_threshold_value_ranges(self):
        """Test that threshold values are in reasonable ranges."""
        # Probabilities should be between 0 and 1
        assert 0 < AnalysisThresholds.AU_SIGNIFICANCE < 1
        assert 0 < AnalysisThresholds.BURNIN_FRACTION < 1
        
        # PSRF should be close to 1 (convergence indicator)
        assert 1 < AnalysisThresholds.MAX_PSRF < 2
        
        # ESS should be positive
        assert AnalysisThresholds.MIN_ESS > 0


if __name__ == "__main__":
    pytest.main([__file__])