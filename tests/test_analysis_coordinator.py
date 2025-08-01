#!/usr/bin/env python3
"""
Test Suite for AnalysisCoordinator Component

Tests the analysis coordinator that orchestrates all components to ensure
it behaves identically to the original monolithic implementation.
"""

import unittest
import tempfile
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis.analysis_coordinator import (
    AnalysisCoordinator, AnalysisCoordinatorError, AnalysisConfigurationError,
    AnalysisExecutionError
)


class TestAnalysisCoordinator(unittest.TestCase):
    """Test AnalysisCoordinator component functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create test alignment
        self.test_alignment = """>Homo_sapiens
ATCGATCGATCGATCG
>Pan_troglodytes
ATCGATCGATCGATCG
>Mus_musculus
ATCGATCGATCGATCG
>Rattus_norvegicus
ATCGATCGATCGATCG
"""
        self.alignment_file = self.temp_dir / "test.fasta"
        with open(self.alignment_file, 'w') as f:
            f.write(self.test_alignment)
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization_basic(self):
        """Test basic AnalysisCoordinator initialization"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        self.assertEqual(ac.alignment_file, self.alignment_file)
        self.assertEqual(ac.analysis_mode, "ml")
        self.assertTrue(ac.do_ml)
        self.assertFalse(ac.do_bayesian)
        self.assertFalse(ac.do_parsimony)
        self.assertTrue(ac.debug)
        
        # Check that core components are initialized
        self.assertIsNotNone(ac.file_manager)
        self.assertIsNotNone(ac.external_tools)
        self.assertIsNotNone(ac.data_processor)
        
        # Analysis components should be None until data is loaded
        self.assertIsNone(ac.ml_analyzer)
        self.assertIsNone(ac.bayesian_analyzer)
        self.assertIsNone(ac.parsimony_analyzer)
    
    def test_initialization_all_analyses(self):
        """Test initialization with all analysis types"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="all",
            temp_dir=self.temp_dir
        )
        
        self.assertEqual(ac.analysis_mode, "all")
        self.assertTrue(ac.do_ml)
        self.assertTrue(ac.do_bayesian)
        self.assertTrue(ac.do_parsimony)
    
    def test_initialization_combined_analyses(self):
        """Test initialization with combined analysis modes"""
        # Test ML + Parsimony
        ac1 = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="ml+parsimony",
            temp_dir=self.temp_dir
        )
        
        self.assertTrue(ac1.do_ml)
        self.assertFalse(ac1.do_bayesian)
        self.assertTrue(ac1.do_parsimony)
        
        # Test Bayesian + Parsimony
        ac2 = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="bayesian+parsimony",
            temp_dir=self.temp_dir
        )
        
        self.assertFalse(ac2.do_ml)
        self.assertTrue(ac2.do_bayesian)
        self.assertTrue(ac2.do_parsimony)
    
    def test_initialization_with_parameters(self):
        """Test initialization with various parameters"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            alignment_format="fasta",
            model="GTR+G+I",
            temp_dir=self.temp_dir,
            paup_path="/usr/local/bin/paup",
            threads="4",
            data_type="dna",
            debug=True,
            keep_files=True,
            # Model parameters
            gamma_shape=0.5,
            prop_invar=0.2,
            base_freq="estimate",
            # Bayesian parameters
            bayes_ngen=100000,
            bayes_burnin=0.25,
            bayes_chains=4,
            # Constraint parameters
            constraint_mode="specific",
            test_branches="Homo_sapiens,Pan_troglodytes",
            # Output parameters
            output_style="ascii"
        )
        
        self.assertEqual(ac.model, "GTR+G+I")
        self.assertEqual(ac.params['paup_path'], "/usr/local/bin/paup")
        self.assertEqual(ac.params['threads'], "4")
        self.assertEqual(ac.params['gamma_shape'], 0.5)
        self.assertEqual(ac.params['bayes_ngen'], 100000)
        self.assertEqual(ac.params['constraint_mode'], "specific")
        self.assertEqual(ac.output_style, "ascii")
        self.assertTrue(ac.keep_files)
    
    def test_context_manager(self):
        """Test AnalysisCoordinator as context manager"""
        with AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        ) as ac:
            self.assertIsNotNone(ac.file_manager)
            self.assertEqual(ac.alignment_file, self.alignment_file)
        
        # After exiting context, cleanup should have been called
        # (We can't easily test this without integration, but the structure is correct)
    
    @patch('core.analysis.ml_analyzer.MLAnalyzer.build_ml_tree')
    @patch('core.analysis.ml_analyzer.MLAnalyzer.get_ml_likelihood')
    def test_build_ml_tree_success(self, mock_get_likelihood, mock_build_tree):
        """Test successful ML tree building"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Create mock tree
        mock_tree = Tree(root=Clade())
        mock_build_tree.return_value = mock_tree
        mock_get_likelihood.return_value = -1234.567
        
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        result = ac.build_ml_tree()
        
        self.assertTrue(result)
        self.assertEqual(ac.ml_tree, mock_tree)
        self.assertEqual(ac.ml_likelihood, -1234.567)
        
        # Check that analysis components were initialized
        self.assertIsNotNone(ac.ml_analyzer)
        self.assertIsNotNone(ac.constraint_manager)
        self.assertIsNotNone(ac.result_processor)
    
    @patch('core.analysis.ml_analyzer.MLAnalyzer.build_ml_tree')
    def test_build_ml_tree_failure(self, mock_build_tree):
        """Test ML tree building failure"""
        mock_build_tree.return_value = None
        
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir
        )
        
        with self.assertRaises(AnalysisExecutionError) as context:
            ac.build_ml_tree()
        
        self.assertIn("Failed to build ML tree", str(context.exception))
    
    def test_calculate_decay_indices_without_tree(self):
        """Test decay index calculation without ML tree"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir
        )
        
        with self.assertRaises(AnalysisExecutionError) as context:
            ac.calculate_decay_indices()
        
        self.assertIn("ML tree must be built", str(context.exception))
    
    @patch('core.analysis.analysis_coordinator.AnalysisCoordinator._run_ml_analysis')
    @patch('core.analysis.constraint_manager.ConstraintManager.get_testable_clades_from_tree')
    def test_calculate_decay_indices_ml_only(self, mock_get_clades, mock_run_ml):
        """Test decay index calculation for ML-only analysis"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Set up mocks
        mock_clades = [("Clade_1", {"Homo_sapiens", "Pan_troglodytes"})]
        mock_get_clades.return_value = mock_clades
        
        mock_ml_results = {
            'ml_tree': Tree(root=Clade()),
            'ml_likelihood': -1234.567,
            'constraint_trees': {},
            'constraint_likelihoods': {},
            'au_results': {}
        }
        mock_run_ml.return_value = mock_ml_results
        
        # Initialize coordinator
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="ml",
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Set up ML tree (normally done by build_ml_tree)
        ac.ml_tree = Tree(root=Clade())
        ac.ml_likelihood = -1234.567
        
        # Load alignment data first
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Calculate decay indices
        decay_indices = ac.calculate_decay_indices()
        
        # Verify ML analysis was called
        mock_run_ml.assert_called_once_with(mock_clades)
        
        # Verify result processor was used
        self.assertIsInstance(decay_indices, dict)
    
    @patch('core.analysis.analysis_coordinator.AnalysisCoordinator._run_ml_analysis')
    @patch('core.analysis.analysis_coordinator.AnalysisCoordinator._run_bayesian_analysis')
    @patch('core.analysis.analysis_coordinator.AnalysisCoordinator._run_parsimony_analysis')
    @patch('core.analysis.constraint_manager.ConstraintManager.get_testable_clades_from_tree')
    def test_calculate_decay_indices_all_analyses(self, mock_get_clades, mock_run_pars, 
                                                 mock_run_bayes, mock_run_ml):
        """Test decay index calculation for all analysis types"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Set up mocks
        mock_clades = [("Clade_1", {"Homo_sapiens", "Pan_troglodytes"})]
        mock_get_clades.return_value = mock_clades
        
        mock_tree = Tree(root=Clade())
        
        mock_run_ml.return_value = {
            'ml_tree': mock_tree,
            'ml_likelihood': -1234.567,
            'constraint_trees': {},
            'constraint_likelihoods': {},
            'au_results': {}
        }
        
        mock_run_bayes.return_value = {
            'consensus_tree': mock_tree,
            'marginal_likelihood': -1234.567,
            'constraint_analyses': {},
            'convergence_diagnostics': {}
        }
        
        mock_run_pars.return_value = {
            'parsimony_tree': mock_tree,
            'parsimony_score': 123,
            'constraint_analyses': {},
            'decay_indices': {}
        }
        
        # Initialize coordinator for all analyses
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="all",
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Set up ML tree
        ac.ml_tree = mock_tree
        ac.ml_likelihood = -1234.567
        
        # Load alignment data first
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Calculate decay indices
        decay_indices = ac.calculate_decay_indices()
        
        # Verify all analyses were called
        mock_run_ml.assert_called_once()
        mock_run_bayes.assert_called_once()
        mock_run_pars.assert_called_once()
        
        self.assertIsInstance(decay_indices, dict)
    
    @patch('core.analysis.bootstrap_manager.BootstrapManager.run_ml_bootstrap')
    @patch('core.analysis.bootstrap_manager.BootstrapManager.get_bootstrap_values')
    @patch('core.analysis.bootstrap_manager.BootstrapManager.get_taxa_bootstrap_map')
    def test_run_bootstrap_analysis_success(self, mock_get_taxa_map, mock_get_values, mock_run_bootstrap):
        """Test successful bootstrap analysis"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Set up mocks
        mock_tree = Tree(root=Clade())
        mock_run_bootstrap.return_value = mock_tree
        mock_get_values.return_value = {"node_1": 85.0}
        mock_get_taxa_map.return_value = {frozenset(["Homo_sapiens", "Pan_troglodytes"]): 85.0}
        
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Set up ML tree
        ac.ml_tree = mock_tree
        ac.ml_likelihood = -1234.567
        
        # Load alignment data first
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Run bootstrap analysis
        result = ac.run_bootstrap_analysis(num_replicates=10)
        
        self.assertTrue(result)
        mock_run_bootstrap.assert_called_once()
        
        # Verify bootstrap manager was created
        self.assertIsNotNone(ac.bootstrap_manager)
    
    def test_run_bootstrap_analysis_no_tree(self):
        """Test bootstrap analysis without ML tree"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir
        )
        
        result = ac.run_bootstrap_analysis()
        
        self.assertFalse(result)
    
    def test_write_formatted_results(self):
        """Test formatted results writing"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Initialize components to create result processor
        ac._initialize_components()
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Test writing results (even with empty results)
        output_file = self.temp_dir / "test_results.txt"
        ac.write_formatted_results(output_file)
        
        # File should be created
        self.assertTrue(output_file.exists())
    
    def test_generate_detailed_report(self):
        """Test detailed report generation"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Initialize components
        ac._initialize_components()
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Test report generation
        report_file = self.temp_dir / "test_report.md"
        ac.generate_detailed_report(report_file)
        
        # File should be created
        self.assertTrue(report_file.exists())
        
        # Should contain markdown content
        content = report_file.read_text()
        self.assertIn("# panDecay Analysis Report", content)
    
    def test_annotate_trees(self):
        """Test tree annotation"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            debug=True
        )
        
        # Set up ML tree
        ac.ml_tree = Tree(root=Clade())
        
        # Test tree annotation
        output_dir = self.temp_dir / "trees"
        tree_files = ac.annotate_trees(output_dir, "test_tree")
        
        # Should create output directory and basic tree file
        self.assertTrue(output_dir.exists())
        self.assertIn('basic', tree_files)
        self.assertTrue(tree_files['basic'].exists())
    
    def test_get_analysis_summary(self):
        """Test analysis summary generation"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="ml+parsimony",
            model="GTR+G",
            temp_dir=self.temp_dir,
            debug=True
        )
        
        summary = ac.get_analysis_summary()
        
        # Check basic summary fields
        self.assertIn('alignment_file', summary)
        self.assertIn('analysis_mode', summary)
        self.assertIn('model', summary)
        self.assertIn('analyses_performed', summary)
        
        self.assertEqual(summary['analysis_mode'], 'ml+parsimony')
        self.assertEqual(summary['model'], 'GTR+G')
        self.assertTrue(summary['analyses_performed']['ml'])
        self.assertFalse(summary['analyses_performed']['bayesian'])
        self.assertTrue(summary['analyses_performed']['parsimony'])
        self.assertFalse(summary['tree_built'])
        self.assertEqual(summary['num_decay_indices'], 0)
    
    def test_getter_methods(self):
        """Test getter methods"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir
        )
        
        # Initially None/empty
        self.assertIsNone(ac.get_ml_tree())
        self.assertIsNone(ac.get_ml_likelihood())
        self.assertEqual(ac.get_decay_indices(), {})
        
        # Set some values
        mock_tree = Tree(root=Clade())
        ac.ml_tree = mock_tree
        ac.ml_likelihood = -1234.567
        ac.decay_indices = {"Clade_1": {"taxa": ["A", "B"]}}
        
        # Test getters
        self.assertEqual(ac.get_ml_tree(), mock_tree)
        self.assertEqual(ac.get_ml_likelihood(), -1234.567)
        
        decay_copy = ac.get_decay_indices()
        self.assertEqual(decay_copy, {"Clade_1": {"taxa": ["A", "B"]}})
        
        # Verify it's a copy (modification shouldn't affect original)
        decay_copy["test"] = "value"
        self.assertNotIn("test", ac.decay_indices)
    
    def test_cleanup_intermediate_files(self):
        """Test intermediate file cleanup"""
        # Test with keep_files=False
        ac1 = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            keep_files=False
        )
        
        # Should not raise exception
        ac1.cleanup_intermediate_files()
        
        # Test with keep_files=True
        ac2 = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            temp_dir=self.temp_dir,
            keep_files=True
        )
        
        # Should not raise exception
        ac2.cleanup_intermediate_files()


class TestAnalysisCoordinatorIntegration(unittest.TestCase):
    """Integration tests for AnalysisCoordinator"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create realistic alignment
        self.realistic_alignment = """>Homo_sapiens
ATCGATCGATCGNNNN----ATCGATCG
>Pan_troglodytes
ATCGATCGATCGNNNN----ATCGATCG
>Gorilla_gorilla
ATCGATCGATCGNNNN----ATCGATCG
>Mus_musculus
ATCGATCGATCGNNNN----ATCGATCG
>Rattus_norvegicus
ATCGATCGATCGNNNN----ATCGATCG
"""
        self.alignment_file = self.temp_dir / "realistic.fasta"
        with open(self.alignment_file, 'w') as f:
            f.write(self.realistic_alignment)
    
    def tearDown(self):
        """Clean up integration test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_component_integration(self):
        """Test integration between all components"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="ml",
            model="GTR+G",
            temp_dir=self.temp_dir,
            constraint_mode="all",
            debug=True,
            keep_files=True
        )
        
        # Test that components can be initialized
        ac._initialize_components()
        
        # Load data and initialize analysis components
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Verify all components are initialized
        self.assertIsNotNone(ac.file_manager)
        self.assertIsNotNone(ac.external_tools)
        self.assertIsNotNone(ac.data_processor)
        self.assertIsNotNone(ac.constraint_manager)
        self.assertIsNotNone(ac.result_processor)
        
        # ML analyzer should be initialized for ML analysis
        self.assertIsNotNone(ac.ml_analyzer)
        
        # Bayesian and parsimony analyzers should be None for ML-only
        self.assertIsNone(ac.bayesian_analyzer)
        self.assertIsNone(ac.parsimony_analyzer)
        
        # Test constraint manager integration
        self.assertEqual(len(ac.constraint_manager.available_taxa), 5)
        self.assertEqual(ac.constraint_manager.constraint_mode, "all")
        
        # Test analysis summary
        summary = ac.get_analysis_summary()
        self.assertTrue(summary['analyses_performed']['ml'])
        self.assertEqual(summary['model'], 'GTR+G')
    
    def test_parameter_propagation(self):
        """Test that parameters are properly propagated to components"""
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="all",
            model="HKY+I",
            paup_path="/custom/paup",
            mrbayes_path="/custom/mb",
            threads="8",
            data_type="dna",
            temp_dir=self.temp_dir,
            constraint_mode="specific",
            test_branches="Homo_sapiens,Pan_troglodytes",
            bayes_ngen=50000,
            bayes_chains=2,
            output_style="ascii",
            debug=True
        )
        
        # Check parameter storage
        self.assertEqual(ac.params['paup_path'], "/custom/paup")
        self.assertEqual(ac.params['mrbayes_path'], "/custom/mb")
        self.assertEqual(ac.params['threads'], "8")
        self.assertEqual(ac.params['bayes_ngen'], 50000)
        self.assertEqual(ac.params['bayes_chains'], 2)
        self.assertEqual(ac.params['constraint_mode'], "specific")
        self.assertEqual(ac.params['test_branches'], "Homo_sapiens,Pan_troglodytes")
        
        # Check analysis mode parsing
        self.assertEqual(ac.analysis_mode, "all")
        self.assertTrue(ac.do_ml)
        self.assertTrue(ac.do_bayesian)
        self.assertTrue(ac.do_parsimony)
        
        # Initialize components and check parameter propagation
        ac._initialize_components()
        ac.data_processor.load_alignment(self.alignment_file, "fasta")
        ac._initialize_analysis_components()
        
        # Check external tools configuration
        self.assertEqual(ac.external_tools.paup_path, "/custom/paup")
        self.assertEqual(ac.external_tools.mrbayes_path, "/custom/mb")
        
        # Check data processor configuration
        self.assertEqual(ac.data_processor.data_type, "dna")
        
        # Check constraint manager configuration
        self.assertEqual(ac.constraint_manager.constraint_mode, "specific")
        self.assertEqual(ac.constraint_manager.test_branches, "Homo_sapiens,Pan_troglodytes")
        
        # Check result processor configuration
        self.assertEqual(ac.result_processor.output_style, "ascii")
    
    def test_error_handling_integration(self):
        """Test error handling across component integration"""
        # Test with non-existent alignment file
        with self.assertRaises(Exception):
            ac = AnalysisCoordinator(
                alignment_file="/nonexistent/file.fasta",
                temp_dir=self.temp_dir
            )
            ac.build_ml_tree()
        
        # Test with invalid analysis mode
        ac = AnalysisCoordinator(
            alignment_file=self.alignment_file,
            analysis_mode="invalid_mode",  # Should still work, just no analyses enabled
            temp_dir=self.temp_dir
        )
        
        # All analysis flags should be False for invalid mode
        self.assertFalse(ac.do_ml)
        self.assertFalse(ac.do_bayesian)
        self.assertFalse(ac.do_parsimony)


def run_analysis_coordinator_tests():
    """Run the AnalysisCoordinator test suite"""
    print("=" * 80)
    print("ANALYSIS COORDINATOR COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing orchestration of all analysis components...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestAnalysisCoordinator))
    suite.addTests(loader.loadTestsFromTestCase(TestAnalysisCoordinatorIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("ANALYSIS COORDINATOR TEST RESULTS")
    print("=" * 80)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    if result.failures:
        print("\nFAILURES:")
        for test, traceback in result.failures:
            print(f"- {test}: {traceback}")
    
    if result.errors:
        print("\nERRORS:")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback}")
    
    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\nAnalysisCoordinator Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_analysis_coordinator_tests()
    sys.exit(0 if success else 1)