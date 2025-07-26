"""
Comprehensive tests for utility modules.

Tests format detectors, config converters, thread calculators,
memory managers, and other utility functions.
"""

import pytest
import tempfile
import shutil
import os
from pathlib import Path
from unittest.mock import Mock, patch

from src.utils.format_detectors import (
    ConfigurationFormatDetector, AlignmentFormatMapper, AnalysisTypeConfiguration,
    YamlFormatDetector, TomlFormatDetector, IniFormatDetector, JsonFormatDetector
)
from src.utils.config_converters import (
    ConfigurationValueConverter, BooleanConverter, FloatConverter,
    IntegerConverter, NormalizationConverter, ChoiceConverter, RangeConverter
)
from src.utils.thread_calculator import (
    ThreadCalculator, AdaptiveThreadCalculator, HighCoreStrategy,
    MediumCoreStrategy, SingleCoreStrategy
)
from src.utils.memory_manager import (
    MemoryMonitor, ChunkedFileReader, MemoryEfficientCache,
    MemoryOptimizer, memory_limited_operation
)
from src.utils.tree_utils import (
    prepare_starting_tree, clean_tree_file, validate_tree_file, count_taxa_in_tree
)
from src.utils.script_generators import PAUPScriptGenerator, MrBayesScriptGenerator


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


class TestConfigurationFormatDetector:
    """Test configuration format detection."""
    
    def test_yaml_detection(self):
        """Test YAML format detection."""
        detector = ConfigurationFormatDetector()
        
        yaml_content = """---
analysis_type: ml
model: GTR
gamma: true
"""
        assert detector.detect_format(yaml_content) == 'yaml'
    
    def test_toml_detection(self):
        """Test TOML format detection."""
        detector = ConfigurationFormatDetector()
        
        toml_content = """[analysis]
type = "ml"
model = "GTR"
gamma = true

[[constraints]]
taxa = ["seq1", "seq2"]
"""
        assert detector.detect_format(toml_content) == 'toml'
    
    def test_ini_detection(self):
        """Test INI format detection."""
        detector = ConfigurationFormatDetector()
        
        ini_content = """[analysis]
type = ml
model = GTR
gamma = yes
"""
        assert detector.detect_format(ini_content) == 'ini'
    
    def test_json_detection(self):
        """Test JSON format detection."""
        detector = ConfigurationFormatDetector()
        
        json_content = """{
    "analysis_type": "ml",
    "model": "GTR",
    "gamma": true
}"""
        assert detector.detect_format(json_content) == 'json'
    
    def test_empty_content_defaults_to_ini(self):
        """Test that empty content defaults to INI format."""
        detector = ConfigurationFormatDetector()
        assert detector.detect_format("") == 'ini'
        assert detector.detect_format("   \n  ") == 'ini'
    
    def test_ambiguous_content_defaults_to_ini(self):
        """Test that ambiguous content defaults to INI format."""
        detector = ConfigurationFormatDetector()
        ambiguous_content = "some random text without clear format indicators"
        assert detector.detect_format(ambiguous_content) == 'ini'


class TestAlignmentFormatMapper:
    """Test alignment format mapping."""
    
    def test_nexus_extension_detection(self):
        """Test NEXUS format detection by extension."""
        mapper = AlignmentFormatMapper()
        
        assert mapper.detect_format_by_extension("test.nex") == 'nexus'
        assert mapper.detect_format_by_extension("test.nexus") == 'nexus'
        assert mapper.detect_format_by_extension("test.NEX") == 'nexus'
    
    def test_fasta_extension_detection(self):
        """Test FASTA format detection by extension."""
        mapper = AlignmentFormatMapper()
        
        assert mapper.detect_format_by_extension("test.fa") == 'fasta'
        assert mapper.detect_format_by_extension("test.fas") == 'fasta'
        assert mapper.detect_format_by_extension("test.fasta") == 'fasta'
    
    def test_phylip_extension_detection(self):
        """Test PHYLIP format detection by extension."""
        mapper = AlignmentFormatMapper()
        
        assert mapper.detect_format_by_extension("test.phy") == 'phylip'
        assert mapper.detect_format_by_extension("test.phylip") == 'phylip'
    
    def test_content_based_detection(self):
        """Test format detection based on file content."""
        mapper = AlignmentFormatMapper()
        
        nexus_content = "#NEXUS\nBEGIN DATA;\nsequences here\nEND;"
        assert mapper.detect_format_by_content(nexus_content) == 'nexus'
        
        fasta_content = ">seq1\nATCGATCG\n>seq2\nATCGATCG"
        assert mapper.detect_format_by_content(fasta_content) == 'fasta'
    
    def test_format_detection_fallback(self):
        """Test format detection fallback to FASTA."""
        mapper = AlignmentFormatMapper()
        
        # No extension, no recognizable content
        assert mapper.detect_format("unknown_file", "random content") == 'fasta'
        
        # No extension provided
        assert mapper.detect_format_by_extension("no_extension_file") is None
    
    def test_register_new_extension(self):
        """Test registering a new extension mapping."""
        mapper = AlignmentFormatMapper()
        
        mapper.register_extension(".xyz", "custom")
        assert mapper.detect_format_by_extension("test.xyz") == 'custom'


class TestAnalysisTypeConfiguration:
    """Test analysis type configuration."""
    
    def test_initialization_with_all_types(self):
        """Test initialization with all analysis types enabled."""
        config = AnalysisTypeConfiguration(
            has_ml=True, 
            has_bayesian=True, 
            has_parsimony=True
        )
        
        assert config.has_ml is True
        assert config.has_bayesian is True
        assert config.has_parsimony is True
        assert config.has_all_types is True
    
    def test_initialization_with_single_type(self):
        """Test initialization with single analysis type."""
        config = AnalysisTypeConfiguration(has_ml=True)
        
        assert config.has_ml is True
        assert config.has_bayesian is False
        assert config.has_parsimony is False
        assert config.has_ml_only is True
    
    def test_no_analysis_types_raises_error(self):
        """Test that no analysis types raises validation error."""
        with pytest.raises(ValueError, match="At least one analysis type must be enabled"):
            AnalysisTypeConfiguration()
    
    def test_analysis_types_property(self):
        """Test analysis_types property returns correct list."""
        config = AnalysisTypeConfiguration(has_ml=True, has_bayesian=True)
        
        types = config.analysis_types
        assert 'ml' in types
        assert 'bayesian' in types
        assert 'parsimony' not in types
        assert len(types) == 2
    
    def test_layout_configuration(self):
        """Test layout configuration for different combinations."""
        # All types
        config_all = AnalysisTypeConfiguration(True, True, True)
        assert config_all.get_layout_configuration() == "full_analysis_layout"
        
        # ML + Bayesian
        config_ml_bayes = AnalysisTypeConfiguration(True, True, False)
        assert config_ml_bayes.get_layout_configuration() == "ml_bayesian_layout"
        
        # ML only
        config_ml = AnalysisTypeConfiguration(True, False, False)
        assert config_ml.get_layout_configuration() == "ml_only_layout"
    
    def test_analysis_count(self):
        """Test analysis count property."""
        config1 = AnalysisTypeConfiguration(True, False, False)
        assert config1.analysis_count == 1
        
        config3 = AnalysisTypeConfiguration(True, True, True)
        assert config3.analysis_count == 3
    
    def test_requires_model_configuration(self):
        """Test requires_model_configuration method."""
        # ML and Bayesian require models
        config_ml = AnalysisTypeConfiguration(True, False, False)
        assert config_ml.requires_model_configuration() is True
        
        config_bayes = AnalysisTypeConfiguration(False, True, False)
        assert config_bayes.requires_model_configuration() is True
        
        # Parsimony doesn't require models
        config_parsimony = AnalysisTypeConfiguration(False, False, True)
        assert config_parsimony.requires_model_configuration() is False


class TestConfigurationValueConverter:
    """Test configuration value conversion."""
    
    def test_boolean_conversion(self):
        """Test boolean parameter conversion."""
        converter = ConfigurationValueConverter()
        
        # Test various true values
        assert converter.convert('gamma', 'true') is True
        assert converter.convert('gamma', 'yes') is True
        assert converter.convert('gamma', '1') is True
        assert converter.convert('debug', True) is True
        
        # Test various false values
        assert converter.convert('gamma', 'false') is False
        assert converter.convert('gamma', 'no') is False
        assert converter.convert('gamma', '0') is False
    
    def test_float_conversion(self):
        """Test float parameter conversion."""
        converter = ConfigurationValueConverter()
        
        assert converter.convert('gamma_shape', '0.5') == 0.5
        assert converter.convert('prop_invar', '0.25') == 0.25
        assert converter.convert('bayes_burnin', 0.1) == 0.1
    
    def test_integer_conversion(self):
        """Test integer parameter conversion."""
        converter = ConfigurationValueConverter()
        
        assert converter.convert('nst', '2') == 2
        assert converter.convert('bootstrap_reps', '1000') == 1000
        assert converter.convert('threads', 4) == 4
    
    def test_choice_conversion(self):
        """Test choice parameter conversion and validation."""
        converter = ConfigurationValueConverter()
        
        # Valid choices
        assert converter.convert('data_type', 'dna') == 'dna'
        assert converter.convert('data_type', 'DNA') == 'dna'  # Case insensitive
        assert converter.convert('criterion', 'ml') == 'ml'
        
        # Invalid choice should raise error
        with pytest.raises(ValueError, match="Invalid choice"):
            converter.convert('data_type', 'invalid_type')
    
    def test_normalization_conversion(self):
        """Test normalization parameter conversion."""
        converter = ConfigurationValueConverter()
        
        # Valid normalization levels
        assert converter.convert('normalization', 'none') == 'none'
        assert converter.convert('normalization', 'basic') == 'basic'
        assert converter.convert('normalization', 'full') == 'full'
        assert converter.convert('normalization', 'NONE') == 'none'  # Case insensitive
        
        # Invalid normalization level should raise error
        with pytest.raises(ValueError, match="Invalid normalization level"):
            converter.convert('normalization', 'invalid')
    
    def test_range_validation(self):
        """Test range validation for parameters."""
        converter = ConfigurationValueConverter()
        
        # Valid thread count
        assert converter.convert('threads', '4') == 4
        
        # Invalid thread count (too low)
        with pytest.raises(ValueError, match="below minimum"):
            converter.convert('threads', '0')
        
        # Invalid thread count (too high)
        with pytest.raises(ValueError, match="above maximum"):
            converter.convert('threads', '100')
    
    def test_unknown_parameter_returns_string(self):
        """Test that unknown parameters are returned as strings."""
        converter = ConfigurationValueConverter()
        
        result = converter.convert('unknown_param', 'some_value')
        assert result == 'some_value'
        assert isinstance(result, str)
    
    def test_none_value_handling(self):
        """Test handling of None values."""
        converter = ConfigurationValueConverter()
        
        result = converter.convert('gamma', None)
        assert result is None


class TestThreadCalculator:
    """Test thread calculation utilities."""
    
    @patch('os.cpu_count')
    def test_single_core_strategy(self, mock_cpu_count):
        """Test single core strategy."""
        mock_cpu_count.return_value = 1
        
        calculator = ThreadCalculator()
        threads = calculator.calculate_optimal_threads("auto")
        
        assert threads == 1
    
    @patch('os.cpu_count')
    def test_medium_core_strategy(self, mock_cpu_count):
        """Test medium core strategy."""
        mock_cpu_count.return_value = 4
        
        calculator = ThreadCalculator()
        threads = calculator.calculate_optimal_threads("auto")
        
        assert threads == 3  # 4 - 1
    
    @patch('os.cpu_count')
    def test_high_core_strategy(self, mock_cpu_count):
        """Test high core strategy."""
        mock_cpu_count.return_value = 16
        
        calculator = ThreadCalculator()
        threads = calculator.calculate_optimal_threads("auto")
        
        assert threads == 14  # 16 - 2
    
    def test_explicit_thread_count(self):
        """Test explicit thread count specification."""
        calculator = ThreadCalculator()
        
        threads = calculator.calculate_optimal_threads(8)
        assert threads == 8
        
        threads = calculator.calculate_optimal_threads("6")
        assert threads == 6
    
    def test_invalid_thread_count(self):
        """Test invalid thread count raises error."""
        calculator = ThreadCalculator()
        
        with pytest.raises(ValueError, match="Thread count must be a number"):
            calculator.calculate_optimal_threads("invalid")
        
        with pytest.raises(ValueError, match="Thread count must be at least"):
            calculator.calculate_optimal_threads(0)
        
        with pytest.raises(ValueError, match="Thread count cannot exceed"):
            calculator.calculate_optimal_threads(100)
    
    def test_memory_usage_estimation(self):
        """Test memory usage estimation."""
        calculator = ThreadCalculator()
        
        memory_usage = calculator.estimate_memory_usage(4, 100.0)
        assert memory_usage > 100.0  # Should be more than base data size
        
        # More threads should use more memory
        memory_usage_8 = calculator.estimate_memory_usage(8, 100.0)
        assert memory_usage_8 > memory_usage
    
    def test_thread_memory_validation(self):
        """Test thread-memory combination validation."""
        calculator = ThreadCalculator()
        
        # Should be valid for reasonable combination
        is_valid = calculator.validate_thread_memory_combination(4, 100.0, 2000.0)
        assert is_valid is True
        
        # Should be invalid for excessive thread count
        is_valid = calculator.validate_thread_memory_combination(32, 1000.0, 1000.0)
        assert is_valid is False


class TestAdaptiveThreadCalculator:
    """Test adaptive thread calculation."""
    
    @patch('os.cpu_count')
    def test_analysis_type_adjustment(self, mock_cpu_count):
        """Test thread adjustment for different analysis types."""
        mock_cpu_count.return_value = 8
        
        calculator = AdaptiveThreadCalculator()
        
        ml_threads = calculator.calculate_for_analysis_type('ml', 'auto')
        bayesian_threads = calculator.calculate_for_analysis_type('bayesian', 'auto')
        parsimony_threads = calculator.calculate_for_analysis_type('parsimony', 'auto')
        
        # ML should use most threads, parsimony least
        assert ml_threads >= bayesian_threads >= parsimony_threads
        assert parsimony_threads >= 1
    
    def test_data_size_adjustment_small(self):
        """Test thread adjustment for small datasets."""
        calculator = AdaptiveThreadCalculator()
        
        threads = calculator.calculate_for_data_size(0.5, 'auto')  # 0.5 MB
        assert threads <= 2  # Should limit threads for small datasets
    
    @patch('os.cpu_count')
    def test_data_size_adjustment_large(self, mock_cpu_count):
        """Test thread adjustment for large datasets."""
        mock_cpu_count.return_value = 16
        
        calculator = AdaptiveThreadCalculator()
        
        threads = calculator.calculate_for_data_size(200.0, 'auto')  # 200 MB
        # Should reduce threads due to memory constraints
        base_threads = calculator.calculate_optimal_threads('auto')
        assert threads <= base_threads


class TestMemoryManager:
    """Test memory management utilities."""
    
    def test_memory_monitor_snapshot(self):
        """Test memory monitoring snapshot creation."""
        monitor = MemoryMonitor()
        
        snapshot = monitor.take_snapshot("test_operation")
        
        assert snapshot.description == "test_operation"
        assert snapshot.process_memory_mb >= 0
        assert snapshot.system_memory_mb > 0
        assert 0 <= snapshot.memory_percent <= 100
        assert snapshot.timestamp > 0
    
    def test_memory_monitor_peak_usage(self):
        """Test peak memory usage tracking."""
        monitor = MemoryMonitor()
        
        # Take multiple snapshots
        monitor.take_snapshot("operation_1")
        monitor.take_snapshot("operation_2")
        monitor.take_snapshot("operation_3")
        
        peak_snapshot = monitor.get_peak_memory_usage()
        assert peak_snapshot is not None
        assert peak_snapshot.description in ["operation_1", "operation_2", "operation_3"]
    
    def test_memory_summary(self):
        """Test memory usage summary."""
        monitor = MemoryMonitor()
        
        # Take some snapshots
        monitor.take_snapshot("start")
        monitor.take_snapshot("middle")
        monitor.take_snapshot("end")
        
        summary = monitor.get_memory_summary()
        
        assert "snapshot_count" in summary
        assert "peak_memory_mb" in summary
        assert "avg_memory_mb" in summary
        assert summary["snapshot_count"] == 3
    
    def test_chunked_file_reader(self, temp_dir):
        """Test chunked file reading."""
        # Create a test file
        test_file = temp_dir / "test_large_file.txt"
        test_content = "line1\nline2\nline3\nline4\nline5\n" * 100
        test_file.write_text(test_content)
        
        reader = ChunkedFileReader(chunk_size=100)
        chunks = list(reader.read_chunks(test_file))
        
        assert len(chunks) > 1  # Should be split into multiple chunks
        
        # Verify all content is preserved
        reconstructed = "".join(chunks)
        assert reconstructed == test_content
    
    def test_chunked_line_reader(self, temp_dir):
        """Test chunked line reading."""
        # Create a test file with many lines
        test_file = temp_dir / "test_lines.txt"
        lines = [f"line_{i}" for i in range(50)]
        test_file.write_text("\n".join(lines))
        
        reader = ChunkedFileReader()
        line_chunks = list(reader.read_lines_chunked(test_file, max_lines=10))
        
        assert len(line_chunks) == 5  # 50 lines / 10 per chunk
        assert len(line_chunks[0]) == 10
        assert len(line_chunks[-1]) == 10
    
    def test_memory_efficient_cache(self):
        """Test memory efficient cache."""
        cache = MemoryEfficientCache(max_size=3)
        
        # Add items
        cache.put("key1", "value1")
        cache.put("key2", "value2")
        cache.put("key3", "value3")
        
        assert cache.get("key1") == "value1"
        assert cache.get("key2") == "value2"
        assert cache.get("key3") == "value3"
        assert cache.size() == 3
        
        # Adding fourth item should evict oldest
        cache.put("key4", "value4")
        assert cache.size() == 3
        assert cache.get("key1") is None  # Should be evicted
        assert cache.get("key4") == "value4"
    
    def test_memory_limited_operation(self):
        """Test memory limited operation context manager."""
        with memory_limited_operation(100.0, "test_operation") as monitor:
            assert isinstance(monitor, MemoryMonitor)
            
            # Do some memory allocation
            data = [i for i in range(1000)]
            
            # Monitor should have snapshots
            assert len(monitor.snapshots) >= 2  # start and end


class TestTreeUtils:
    """Test tree utility functions."""
    
    def test_prepare_starting_tree_valid(self, temp_dir):
        """Test preparing a valid starting tree."""
        # Create a valid tree file
        tree_file = temp_dir / "starting_tree.tre"
        tree_content = "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);"
        tree_file.write_text(tree_content)
        
        prepared_tree = prepare_starting_tree(tree_file, temp_dir, "prepared.tre")
        
        assert prepared_tree is not None
        assert prepared_tree.exists()
        assert prepared_tree.name == "prepared.tre"
        assert prepared_tree.read_text() == tree_content
    
    def test_prepare_starting_tree_missing(self, temp_dir):
        """Test preparing a missing starting tree."""
        nonexistent_tree = temp_dir / "nonexistent.tre"
        
        prepared_tree = prepare_starting_tree(nonexistent_tree, temp_dir, "prepared.tre")
        
        assert prepared_tree is None
    
    def test_prepare_starting_tree_none(self, temp_dir):
        """Test preparing with None starting tree."""
        prepared_tree = prepare_starting_tree(None, temp_dir, "prepared.tre")
        
        assert prepared_tree is None
    
    def test_validate_tree_file_valid(self, temp_dir):
        """Test validating a valid tree file."""
        tree_file = temp_dir / "valid_tree.tre"
        tree_content = "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);"
        tree_file.write_text(tree_content)
        
        is_valid = validate_tree_file(tree_file, require_branch_lengths=True)
        assert is_valid is True
        
        is_valid_no_lengths = validate_tree_file(tree_file, require_branch_lengths=False)
        assert is_valid_no_lengths is True
    
    def test_validate_tree_file_no_branch_lengths(self, temp_dir):
        """Test validating tree file without branch lengths."""
        tree_file = temp_dir / "no_lengths_tree.tre"
        tree_content = "(seq1,seq2,(seq3,seq4));"
        tree_file.write_text(tree_content)
        
        is_valid = validate_tree_file(tree_file, require_branch_lengths=False)
        assert is_valid is True
        
        is_valid_with_lengths = validate_tree_file(tree_file, require_branch_lengths=True)
        assert is_valid_with_lengths is False
    
    def test_validate_tree_file_invalid_format(self, temp_dir):
        """Test validating tree file with invalid format."""
        tree_file = temp_dir / "invalid_tree.tre"
        tree_content = "not a valid tree format"
        tree_file.write_text(tree_content)
        
        is_valid = validate_tree_file(tree_file)
        assert is_valid is False
    
    def test_clean_tree_file(self, temp_dir):
        """Test cleaning tree file metadata."""
        tree_file = temp_dir / "messy_tree.tre"
        tree_content = "(seq1:0.1,seq2:0.1);[some metadata after semicolon]"
        tree_file.write_text(tree_content)
        
        clean_tree_file(tree_file, remove_metadata=True)
        
        cleaned_content = tree_file.read_text()
        assert cleaned_content == "(seq1:0.1,seq2:0.1);"
    
    def test_count_taxa_in_tree(self, temp_dir):
        """Test counting taxa in tree file."""
        tree_file = temp_dir / "count_tree.tre"
        tree_content = "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);"
        tree_file.write_text(tree_content)
        
        taxa_count = count_taxa_in_tree(tree_file)
        assert taxa_count == 4  # Should count 4 taxa based on commas


class TestScriptGenerators:
    """Test script generation utilities."""
    
    def test_paup_script_generator(self):
        """Test PAUP* script generation."""
        generator = PAUPScriptGenerator()
        
        script = generator.generate_ml_search_script(
            alignment_file="test.nex",
            threads=2,
            model_settings={'nst': 2, 'rates': 'gamma'}
        )
        
        assert "#NEXUS" in script
        assert "begin paup;" in script
        assert "execute test.nex;" in script
        assert "nthreads=2" in script
        assert "nst=2" in script
        assert "rates=gamma" in script
        assert "quit;" in script
        assert "end;" in script
    
    def test_mrbayes_script_generator(self):
        """Test MrBayes script generation."""
        generator = MrBayesScriptGenerator()
        
        script = generator.generate_bayesian_analysis_script(
            alignment_file="test.nex",
            model_settings={'nst': 2},
            mcmc_settings={'ngen': 1000000, 'nchains': 4}
        )
        
        assert "#NEXUS" in script
        assert "begin mrbayes;" in script
        assert "execute test.nex;" in script
        assert "ngen=1000000" in script
        assert "nchains=4" in script
        assert "mcmc;" in script
        assert "end;" in script


if __name__ == "__main__":
    pytest.main([__file__])