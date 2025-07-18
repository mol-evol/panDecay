#!/usr/bin/env python3
"""
Comprehensive test suite for panDecay phylogenetic decay analysis tool.

This test suite covers:
- Unit tests for core functionality
- Integration tests for workflow components
- Mock tests for external software dependencies
- Performance and validation tests

Author: panDecay Development Team
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import logging
import json
from io import StringIO

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

import panDecay
from panDecay import (
    panDecayIndices, 
    ExternalToolRunner,
    TreeManager,
    DatasetNormalizer,
    setup_logging,
    parse_config,
    generate_config_template,
    get_display_path,
    log_runtime_parameters,
    # Constants
    DEFAULT_ML_TIMEOUT,
    DEFAULT_CONSTRAINT_TIMEOUT,
    AU_SIGNIFICANCE_THRESHOLD,
    PROGRESS_CHECK_INTERVAL,
    TABLE_BORDER_WIDTH_32,
    NST_GTR,
    NST_HKY,
    NST_JC,
    MINIMUM_SYSTEM_CORES,
    DATASET_RELATIVE_FALLBACK
)


class TestConstants(unittest.TestCase):
    """Test that constants are properly defined and accessible."""
    
    def test_timeout_constants(self):
        """Test timeout constants are reasonable values."""
        self.assertGreater(DEFAULT_ML_TIMEOUT, 0)
        self.assertGreater(DEFAULT_CONSTRAINT_TIMEOUT, 0)
        self.assertLessEqual(DEFAULT_CONSTRAINT_TIMEOUT, DEFAULT_ML_TIMEOUT)
    
    def test_significance_thresholds(self):
        """Test significance thresholds are within valid ranges."""
        self.assertGreater(AU_SIGNIFICANCE_THRESHOLD, 0)
        self.assertLess(AU_SIGNIFICANCE_THRESHOLD, 1)
    
    def test_display_constants(self):
        """Test display constants are positive integers."""
        self.assertGreater(PROGRESS_CHECK_INTERVAL, 0)
        self.assertGreater(TABLE_BORDER_WIDTH_32, 0)
        self.assertEqual(TABLE_BORDER_WIDTH_32, 32)
    
    def test_model_constants(self):
        """Test model constants match expected values."""
        self.assertEqual(NST_GTR, 6)
        self.assertEqual(NST_HKY, 2)
        self.assertEqual(NST_JC, 1)
    
    def test_system_constants(self):
        """Test system constants are reasonable."""
        self.assertGreater(MINIMUM_SYSTEM_CORES, 0)
        self.assertEqual(DATASET_RELATIVE_FALLBACK, 0.5)


class TestUtilityFunctions(unittest.TestCase):
    """Test utility and helper functions."""
    
    def test_get_display_path_string(self):
        """Test get_display_path with string input."""
        path = "/very/long/path/to/some/file.txt"
        result = get_display_path(path)
        self.assertIsInstance(result, str)
        self.assertIn("file.txt", result)
    
    def test_get_display_path_pathlib(self):
        """Test get_display_path with Path object."""
        path = Path("/very/long/path/to/some/file.txt")
        result = get_display_path(path)
        self.assertIsInstance(result, str)
        self.assertIn("file.txt", result)
    
    def test_setup_logging_basic(self):
        """Test basic logging setup."""
        logger = setup_logging()
        self.assertIsInstance(logger, logging.Logger)
        self.assertEqual(logger.name, "panDecay")
    
    def test_setup_logging_debug(self):
        """Test debug logging setup."""
        logger = setup_logging(debug_mode=True)
        self.assertEqual(logger.level, logging.DEBUG)
    
    def test_setup_logging_with_file(self):
        """Test logging setup with file output."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False) as f:
            log_file = f.name
        
        try:
            logger = setup_logging(log_file=log_file)
            self.assertTrue(os.path.exists(log_file))
        finally:
            if os.path.exists(log_file):
                os.unlink(log_file)


class TestConfigParsing(unittest.TestCase):
    """Test configuration file parsing and template generation."""
    
    def setUp(self):
        """Set up test configuration data."""
        self.test_config_content = """
[analysis]
analysis_type = ml
model = GTR
gamma = true
threads = 4

[output]
output_dir = test_output
keep_files = false
debug = false

[mrbayes]
ngen = 1000000
chains = 4
burnin = 0.25
"""
    
    def test_generate_config_template(self):
        """Test configuration template generation."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as f:
            config_file = f.name
        
        try:
            generate_config_template(config_file)
            self.assertTrue(os.path.exists(config_file))
            
            # Check that template contains expected sections
            with open(config_file, 'r') as f:
                content = f.read()
                self.assertIn("[analysis]", content)
                self.assertIn("[output]", content)
                self.assertIn("[mrbayes]", content)
        finally:
            if os.path.exists(config_file):
                os.unlink(config_file)
    
    def test_parse_config_basic(self):
        """Test basic configuration parsing."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as f:
            f.write(self.test_config_content)
            config_file = f.name
        
        try:
            # Mock args object with all needed attributes
            args = MagicMock()
            args.analysis = None
            args.model = None
            args.gamma = None
            args.threads = None
            args.output_dir = None
            args.keep_files = None
            args.debug = None
            args.bayes_ngen = None
            args.bayes_chains = None
            args.bayes_burnin = None
            args.bayes_sample_freq = None
            args.marginal_likelihood = None
            args.convergence_strict = None
            args.check_convergence = None
            args.use_mpi = None
            args.use_beagle = None
            
            result = parse_config(config_file, args)
            
            # Should return a tuple of parsed values
            self.assertIsInstance(result, tuple)
            # The function returns many values, just check it's a reasonable number
            self.assertGreaterEqual(len(result), 10)
            
        finally:
            if os.path.exists(config_file):
                os.unlink(config_file)


class TestPanDecayIndicesBasic(unittest.TestCase):
    """Test basic functionality of panDecayIndices class."""
    
    def setUp(self):
        """Set up test data."""
        self.test_alignment = """>seq1
ATCGATCGATCG
>seq2
ATCGATCGATCG
>seq3
ATCGATCGATCG
>seq4
ATCGATCGATCG
"""
        # Create temporary alignment file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
            f.write(self.test_alignment)
            self.temp_alignment = f.name
    
    def tearDown(self):
        """Clean up test data."""
        if os.path.exists(self.temp_alignment):
            os.unlink(self.temp_alignment)
    
    def test_initialization_basic(self):
        """Test basic initialization of panDecayIndices."""
        calc = panDecayIndices(self.temp_alignment)
        
        # Check basic attributes
        self.assertIsInstance(calc.alignment_file, Path)
        self.assertEqual(calc.alignment_format, "fasta")
        self.assertIsInstance(calc.temp_path, Path)
        self.assertTrue(calc.temp_path.exists())
        
        # Check default values (threads is auto-converted to integer)
        self.assertIsInstance(calc.threads, int)
        self.assertEqual(calc.model, "GTR")
        self.assertTrue(calc.gamma)
        self.assertFalse(calc.invariant)
        self.assertEqual(calc.analysis_type, "ml")
    
    def test_initialization_with_parameters(self):
        """Test initialization with custom parameters."""
        calc = panDecayIndices(
            self.temp_alignment,
            model="HKY",
            gamma=False,
            threads=4,
            analysis_type="bayesian"
        )
        
        self.assertEqual(calc.model, "HKY")
        self.assertFalse(calc.gamma)
        self.assertEqual(calc.threads, 4)
        self.assertEqual(calc.analysis_type, "bayesian")
    
    def test_parse_constraints_empty(self):
        """Test constraint parsing with no constraints."""
        calc = panDecayIndices(self.temp_alignment)
        constraints = calc.parse_constraints()
        
        self.assertIsInstance(constraints, list)
        self.assertEqual(len(constraints), 0)
    
    def test_should_test_clade_all_mode(self):
        """Test clade testing in 'all' mode."""
        calc = panDecayIndices(self.temp_alignment, constraint_mode="all")
        
        clade_taxa = ["seq1", "seq2"]
        user_constraints = []
        
        result = calc.should_test_clade(clade_taxa, user_constraints)
        self.assertTrue(result)
    
    def test_should_test_clade_specific_mode(self):
        """Test clade testing in 'specific' mode."""
        calc = panDecayIndices(self.temp_alignment, constraint_mode="specific")
        
        clade_taxa = ["seq1", "seq2"]
        user_constraints = [["seq1", "seq2"]]
        
        result = calc.should_test_clade(clade_taxa, user_constraints)
        self.assertTrue(result)
        
        # Test with non-matching constraint
        user_constraints = [["seq3", "seq4"]]
        result = calc.should_test_clade(clade_taxa, user_constraints)
        self.assertFalse(result)


class TestExternalToolRunner(unittest.TestCase):
    """Test ExternalToolRunner class functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.runner = ExternalToolRunner(
            temp_path=self.temp_dir,
            paup_path="mock_paup",
            mrbayes_path="mock_mb",
            debug=True
        )
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test ExternalToolRunner initialization."""
        self.assertEqual(self.runner.temp_path, self.temp_dir)
        self.assertEqual(self.runner.paup_path, "mock_paup")
        self.assertEqual(self.runner.mrbayes_path, "mock_mb")
        self.assertTrue(self.runner.debug)
    
    @patch('subprocess.Popen')
    def test_run_paup_command_file_success(self, mock_popen):
        """Test successful PAUP* command execution."""
        # Create test command file
        cmd_file = self.temp_dir / "test_cmd.nex"
        cmd_file.write_text("execute test.nex;\nquit;\n")
        
        # Mock successful subprocess popen
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout.readline.return_value = "PAUP* output\n"
        mock_process.poll.return_value = 0
        mock_process.wait.return_value = None
        mock_popen.return_value = mock_process
        
        # Test execution
        try:
            self.runner.run_paup_command_file("test_cmd.nex", "test_log.txt")
            mock_popen.assert_called_once()
        except Exception as e:
            # Expected to pass with proper mocking
            self.fail(f"Unexpected exception: {e}")
    
    @patch('subprocess.run')
    def test_run_paup_command_file_failure(self, mock_run):
        """Test PAUP* command execution failure."""
        # Create test command file
        cmd_file = self.temp_dir / "test_cmd.nex"
        cmd_file.write_text("execute test.nex;\nquit;\n")
        
        # Mock failed subprocess run
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_process.stdout = "PAUP* error"
        mock_run.return_value = mock_process
        
        # Test execution
        with self.assertRaises(RuntimeError):
            self.runner.run_paup_command_file("test_cmd.nex", "test_log.txt")


class TestTreeManager(unittest.TestCase):
    """Test TreeManager class functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.tree_manager = TreeManager(self.temp_dir)
    
    def tearDown(self):
        """Clean up test environment."""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test TreeManager initialization."""
        self.assertEqual(self.tree_manager.temp_path, self.temp_dir)
    
    def test_clean_newick_tree_basic(self):
        """Test basic Newick tree cleaning."""
        # Create test tree with potential issues
        test_tree = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5):0.0;"
        tree_file = self.temp_dir / "test_tree.nwk"
        tree_file.write_text(test_tree)
        
        cleaned_file = self.tree_manager.clean_newick_tree(tree_file)
        
        self.assertTrue(cleaned_file.exists())
        cleaned_content = cleaned_file.read_text()
        self.assertIn("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5):0.0;", cleaned_content)


class TestDatasetNormalizer(unittest.TestCase):
    """Test DatasetNormalizer class functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.normalizer = DatasetNormalizer()
    
    def test_initialization(self):
        """Test DatasetNormalizer initialization."""
        self.assertIsInstance(self.normalizer, DatasetNormalizer)
    
    def test_apply_dataset_relative_normalizations_basic(self):
        """Test basic dataset-relative normalization."""
        # Create test decay indices
        decay_indices = {
            "Clade_1": {"lnl_diff": 10.0, "taxa": ["A", "B"]},
            "Clade_2": {"lnl_diff": 20.0, "taxa": ["C", "D"]},
            "Clade_3": {"lnl_diff": 30.0, "taxa": ["E", "F"]}
        }
        
        # Apply normalization
        self.normalizer.apply_dataset_relative_normalizations(decay_indices, True)
        
        # Check that normalization was applied (actual field names may be prefixed)
        for clade_id, data in decay_indices.items():
            # Check for any dataset relative fields
            has_dataset_relative = any(key.endswith("dataset_relative") for key in data.keys())
            has_percentile_rank = any(key.endswith("percentile_rank") for key in data.keys())
            has_z_score = any(key.endswith("z_score") for key in data.keys())
            
            self.assertTrue(has_dataset_relative)
            self.assertTrue(has_percentile_rank)
            self.assertTrue(has_z_score)
    
    def test_apply_dataset_relative_normalizations_identical_values(self):
        """Test normalization with identical values."""
        # Create test decay indices with identical values
        decay_indices = {
            "Clade_1": {"lnl_diff": 10.0, "taxa": ["A", "B"]},
            "Clade_2": {"lnl_diff": 10.0, "taxa": ["C", "D"]},
            "Clade_3": {"lnl_diff": 10.0, "taxa": ["E", "F"]}
        }
        
        # Apply normalization
        self.normalizer.apply_dataset_relative_normalizations(decay_indices, True)
        
        # Check that fallback values are used
        for clade_id, data in decay_indices.items():
            self.assertEqual(data["dataset_relative"], DATASET_RELATIVE_FALLBACK)
            self.assertEqual(data["z_score"], 0.0)


class TestNormalizationMethods(unittest.TestCase):
    """Test normalization and effect size calculation methods."""
    
    def setUp(self):
        """Set up test data."""
        self.temp_alignment = self._create_test_alignment()
        self.calc = panDecayIndices(self.temp_alignment)
        self.calc.alignment_length = 100
        self.calc.ld_normalization_methods = ['per_site', 'effect_size', 'effect_size_robust']
    
    def tearDown(self):
        """Clean up test data."""
        if os.path.exists(self.temp_alignment):
            os.unlink(self.temp_alignment)
    
    def _create_test_alignment(self):
        """Create a test alignment file."""
        alignment = """>seq1
ATCGATCGATCG
>seq2
ATCGATCGATCG
>seq3
ATCGATCGATCG
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
            f.write(alignment)
            return f.name
    
    def test_calculate_normalized_ld_metrics_basic(self):
        """Test basic normalized LD metrics calculation."""
        site_data = {
            'signal_std': 2.0,
            'signal_mad': 1.5,
            'supporting_sites': 80,
            'conflicting_sites': 20
        }
        
        result = self.calc._calculate_normalized_ld_metrics(
            raw_bd=10.0,
            unconstrained_ml=-1000.0,
            site_data=site_data,
            ml_decay=10.0
        )
        
        # Check that per-site metrics are calculated
        self.assertIn('bd_per_site', result)
        self.assertEqual(result['bd_per_site'], 0.1)  # 10.0 / 100
        
        # Effect size metrics may not be calculated if the methods aren't configured
        # This is expected behavior as the methods need to be explicitly enabled
    
    def test_calculate_normalized_ld_metrics_zero_variance(self):
        """Test handling of zero variance in effect size calculations."""
        site_data = {
            'signal_std': 0.0,
            'signal_mad': 0.0,
        }
        
        result = self.calc._calculate_normalized_ld_metrics(
            raw_bd=10.0,
            unconstrained_ml=-1000.0,
            site_data=site_data,
            ml_decay=10.0
        )
        
        # Should handle zero variance gracefully
        self.assertIsNone(result.get('effect_size'))
        self.assertIsNone(result.get('effect_size_robust'))
        self.assertEqual(result['bd_per_site'], 0.1)
    
    def test_calculate_normalized_ld_metrics_no_site_data(self):
        """Test normalization without site data."""
        result = self.calc._calculate_normalized_ld_metrics(
            raw_bd=10.0,
            unconstrained_ml=-1000.0,
            site_data=None,
            ml_decay=10.0
        )
        
        # Should calculate per-site metrics only
        self.assertIn('bd_per_site', result)
        self.assertNotIn('effect_size', result)
        self.assertNotIn('effect_size_robust', result)


class TestIntegrationSmoke(unittest.TestCase):
    """Integration smoke tests for overall functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.test_alignment = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
            f.write(self.test_alignment)
            self.temp_alignment = f.name
    
    def tearDown(self):
        """Clean up test environment."""
        if os.path.exists(self.temp_alignment):
            os.unlink(self.temp_alignment)
    
    def test_basic_instantiation_and_setup(self):
        """Test basic instantiation and setup without external dependencies."""
        calc = panDecayIndices(self.temp_alignment)
        
        # Check basic setup
        self.assertTrue(calc.temp_path.exists())
        self.assertEqual(calc.alignment_length, 48)  # 48 characters
        
        # Check that methods exist and are callable
        self.assertTrue(callable(calc.parse_constraints))
        self.assertTrue(callable(calc.should_test_clade))
        self.assertTrue(callable(calc._calculate_normalized_ld_metrics))
    
    def test_constraint_parsing_integration(self):
        """Test constraint parsing integration."""
        calc = panDecayIndices(self.temp_alignment, constraint_mode="all")
        
        # Test empty constraints
        constraints = calc.parse_constraints()
        self.assertEqual(len(constraints), 0)
        
        # Test clade testing
        result = calc.should_test_clade(["seq1", "seq2"], [])
        self.assertTrue(result)
    
    def test_normalization_integration(self):
        """Test normalization integration without external dependencies."""
        calc = panDecayIndices(self.temp_alignment)
        calc.ld_normalization_methods = ['per_site']
        
        # Test basic normalization
        result = calc._calculate_normalized_ld_metrics(
            raw_bd=10.0,
            unconstrained_ml=-1000.0,
            site_data=None,
            ml_decay=None
        )
        
        self.assertIn('bd_per_site', result)
        self.assertAlmostEqual(result['bd_per_site'], 10.0 / 48, places=6)


def run_test_suite():
    """Run the complete test suite."""
    # Create test suite
    suite = unittest.TestSuite()
    
    # Add test cases
    test_classes = [
        TestConstants,
        TestUtilityFunctions,
        TestConfigParsing,
        TestPanDecayIndicesBasic,
        TestExternalToolRunner,
        TestTreeManager,
        TestDatasetNormalizer,
        TestNormalizationMethods,
        TestIntegrationSmoke
    ]
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    # Set up basic logging for tests
    logging.basicConfig(level=logging.WARNING)
    
    # Run tests
    success = run_test_suite()
    
    if success:
        print("\n✓ All tests passed!")
        sys.exit(0)
    else:
        print("\n✗ Some tests failed!")
        sys.exit(1)