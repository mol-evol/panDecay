#!/usr/bin/env python3
"""
Test Suite for ResultProcessor Component

Tests the extracted result processing functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
import csv
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis import FileManager, ExternalTools, DataProcessor
from core.analysis.result_processor import (
    ResultProcessor, ResultProcessingError, OutputFormatError,
    TreeAnnotationError, VisualizationError
)


class TestResultProcessor(unittest.TestCase):
    """Test ResultProcessor component functionality"""
    
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
        
        # Initialize components
        self.file_manager = FileManager(temp_dir=self.temp_dir, debug=True)
        self.external_tools = ExternalTools(debug=True)
        self.data_processor = DataProcessor(data_type="dna", debug=True)
        
        # Load alignment
        self.data_processor.load_alignment(self.alignment_file, "fasta")
    
    def tearDown(self):
        """Clean up test fixtures"""
        self.file_manager.cleanup_all()
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def _create_mock_tree(self):
        """Create a mock phylogenetic tree for testing"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Create tree: ((Homo_sapiens,Pan_troglodytes),(Mus_musculus,Rattus_norvegicus))
        homo = Clade(name="Homo_sapiens")
        pan = Clade(name="Pan_troglodytes")
        mus = Clade(name="Mus_musculus")
        rattus = Clade(name="Rattus_norvegicus")
        
        primates = Clade()
        primates.clades = [homo, pan]
        
        rodents = Clade()
        rodents.clades = [mus, rattus]
        
        root = Clade()
        root.clades = [primates, rodents]
        
        return Tree(root=root)
    
    def test_initialization(self):
        """Test ResultProcessor initialization"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            output_style="unicode",
            viz_format="png",
            debug=True
        )
        
        self.assertEqual(rp.output_style, "unicode")
        self.assertEqual(rp.viz_format, "png")
        self.assertTrue(rp.debug)
        self.assertEqual(len(rp.decay_indices), 0)
        self.assertFalse(rp.has_results())
    
    def test_add_ml_results(self):
        """Test adding ML analysis results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_tree()
        ml_likelihood = -1234.567
        constraint_trees = {"Clade_1": mock_tree}
        constraint_likelihoods = {"Clade_1": -1245.678}
        au_results = {"Clade_1": 0.023}
        
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=ml_likelihood,
            constraint_trees=constraint_trees,
            constraint_likelihoods=constraint_likelihoods,
            au_results=au_results
        )
        
        self.assertEqual(rp.ml_results['likelihood'], ml_likelihood)
        self.assertEqual(len(rp.ml_results['constraint_trees']), 1)
        self.assertEqual(len(rp.ml_results['au_results']), 1)
        self.assertTrue(rp.has_results())
    
    def test_add_bayesian_results(self):
        """Test adding Bayesian analysis results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_tree()
        marginal_likelihood = -1234.567
        constraint_analyses = {"Clade_1": {"marginal_likelihood": -1245.678}}
        convergence_diagnostics = {"ess_min": 150, "psrf_max": 1.02}
        
        rp.add_bayesian_results(
            consensus_tree=mock_tree,
            marginal_likelihood=marginal_likelihood,
            constraint_analyses=constraint_analyses,
            convergence_diagnostics=convergence_diagnostics
        )
        
        self.assertEqual(rp.bayesian_results['marginal_likelihood'], marginal_likelihood)
        self.assertEqual(len(rp.bayesian_results['constraint_analyses']), 1)
        self.assertIn('ess_min', rp.bayesian_results['convergence_diagnostics'])
    
    def test_add_parsimony_results(self):
        """Test adding parsimony analysis results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_tree()
        parsimony_score = 123
        constraint_analyses = {"Clade_1": {"score": 125}}
        decay_indices = {"Clade_1": 2.0}
        
        rp.add_parsimony_results(
            parsimony_tree=mock_tree,
            parsimony_score=parsimony_score,
            constraint_analyses=constraint_analyses,
            decay_indices=decay_indices
        )
        
        self.assertEqual(rp.parsimony_results['score'], parsimony_score)
        self.assertEqual(len(rp.parsimony_results['decay_indices']), 1)
    
    def test_add_bootstrap_results(self):
        """Test adding bootstrap analysis results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_tree()
        bootstrap_values = {"node_1": 85.0, "node_2": 92.0}
        taxa_bootstrap_map = {
            frozenset(["Homo_sapiens", "Pan_troglodytes"]): 85.0,
            frozenset(["Mus_musculus", "Rattus_norvegicus"]): 92.0
        }
        
        rp.add_bootstrap_results(
            bootstrap_tree=mock_tree,
            bootstrap_values=bootstrap_values,
            taxa_bootstrap_map=taxa_bootstrap_map
        )
        
        self.assertEqual(len(rp.bootstrap_results['values']), 2)
        self.assertEqual(len(rp.bootstrap_results['taxa_map']), 2)
    
    def test_add_site_analysis_results(self):
        """Test adding site analysis results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        site_data = {
            "Clade_1": {
                "supporting_sites": [1, 3, 5],
                "conflicting_sites": [2, 4],
                "site_likelihoods": [-12.34, -23.45, -34.56]
            }
        }
        site_visualizations = {"Clade_1_hist.png": Path("/tmp/hist.png")}
        
        rp.add_site_analysis_results(
            site_data=site_data,
            site_visualizations=site_visualizations
        )
        
        self.assertEqual(len(rp.site_analysis_results['data']), 1)
        self.assertIn("Clade_1", rp.site_analysis_results['data'])
    
    def test_calculate_decay_indices_ml_only(self):
        """Test decay index calculation for ML-only analysis"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add ML results
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567,
            constraint_likelihoods={"Clade_1": -1245.678},
            au_results={"Clade_1": 0.023}
        )
        
        # Calculate decay indices
        decay_indices = rp.calculate_decay_indices()
        
        self.assertGreater(len(decay_indices), 0)
        
        # Check that at least one clade was processed
        for clade_id, data in decay_indices.items():
            self.assertIn('clade_id', data)
            self.assertIn('taxa', data)
            self.assertIn('analysis_types', data)
            self.assertIn('ml', data['analysis_types'])
    
    def test_support_values_only(self):
        """Test that system provides raw numerical values without categorization"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Verify categorization methods have been removed
        self.assertFalse(hasattr(rp, '_categorize_au_support'))
        self.assertFalse(hasattr(rp, '_categorize_bootstrap_support'))
        self.assertFalse(hasattr(rp, '_categorize_bayesian_decay'))
        self.assertFalse(hasattr(rp, '_categorize_bremer_support'))
        
        # System should only provide raw numerical values for user interpretation
    
    def test_format_support_symbol(self):
        """Test support symbol formatting"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            output_style="unicode"
        )
        
        # Test bootstrap symbols
        symbol_high = rp._format_support_symbol(96.0, "bootstrap")
        symbol_low = rp._format_support_symbol(45.0, "bootstrap")
        
        self.assertEqual(symbol_high, '●')  # Unicode high support
        self.assertEqual(symbol_low, '·')   # Unicode very low support
        
        # Test with ASCII style
        rp.output_style = "ascii"
        symbol_high_ascii = rp._format_support_symbol(96.0, "bootstrap")
        symbol_low_ascii = rp._format_support_symbol(45.0, "bootstrap")
        
        self.assertEqual(symbol_high_ascii, '*')  # ASCII high support
        self.assertEqual(symbol_low_ascii, '.')   # ASCII very low support
    
    def test_format_table_row(self):
        """Test table row formatting"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        columns = ["Clade_1", "85.0", "0.023"]
        widths = [10, 8, 8]
        
        # Test left alignment
        row_left = rp._format_table_row(columns, widths, "left")
        self.assertIn("Clade_1   ", row_left)
        
        # Test right alignment
        row_right = rp._format_table_row(columns, widths, "right")
        self.assertIn("   Clade_1", row_right)
        
        # Test center alignment
        row_center = rp._format_table_row(columns, widths, "center")
        self.assertIn(" Clade_1  ", row_center)
    
    def test_write_formatted_results(self):
        """Test formatted results writing"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add some test data
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567,
            constraint_likelihoods={"Clade_1": -1245.678},
            au_results={"Clade_1": 0.023}
        )
        
        rp.add_bootstrap_results(
            bootstrap_tree=mock_tree,
            bootstrap_values={"node_1": 85.0},
            taxa_bootstrap_map={frozenset(["Homo_sapiens", "Pan_troglodytes"]): 85.0}
        )
        
        # Calculate decay indices
        rp.calculate_decay_indices()
        
        # Write results
        output_file = self.temp_dir / "test_results.txt"
        rp.write_formatted_results(output_file)
        
        # Verify file was created and has content
        self.assertTrue(output_file.exists())
        content = output_file.read_text()
        
        self.assertIn("panDecay Analysis Results", content)
        self.assertIn("Analysis Summary", content)
        self.assertIn("Decay Indices", content)
        self.assertIn("Summary Statistics", content)
    
    def test_generate_detailed_report(self):
        """Test detailed markdown report generation"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add test data
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567,
            constraint_likelihoods={"Clade_1": -1245.678},
            au_results={"Clade_1": 0.023}
        )
        
        rp.calculate_decay_indices()
        
        # Generate report
        report_file = self.temp_dir / "test_report.md"
        rp.generate_detailed_report(report_file)
        
        # Verify file was created and has markdown content
        self.assertTrue(report_file.exists())
        content = report_file.read_text()
        
        self.assertIn("# panDecay Analysis Report", content)
        self.assertIn("## Analysis Overview", content)
        self.assertIn("## Detailed Results", content)
        self.assertIn("## Methodology", content)
        self.assertIn("| Clade |", content)  # Markdown table
    
    def test_get_summary_statistics(self):
        """Test summary statistics generation"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add test data with multiple analysis types
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567,
            au_results={"Clade_1": 0.023, "Clade_2": 0.156}
        )
        
        rp.add_bootstrap_results(
            bootstrap_tree=mock_tree,
            bootstrap_values={"node_1": 85.0, "node_2": 92.0},
            taxa_bootstrap_map={
                frozenset(["Homo_sapiens", "Pan_troglodytes"]): 85.0,
                frozenset(["Mus_musculus", "Rattus_norvegicus"]): 92.0
            }
        )
        
        rp.calculate_decay_indices()
        
        # Get statistics
        stats = rp.get_summary_statistics()
        
        self.assertIn('total_clades', stats)
        self.assertIn('analysis_types', stats)
        self.assertGreater(stats['total_clades'], 0)
        self.assertIn('ml', stats['analysis_types'])
        
        # Check bootstrap statistics if available
        if 'bootstrap' in stats:
            self.assertIn('count', stats['bootstrap'])
            self.assertIn('mean', stats['bootstrap'])
            self.assertIn('high_support', stats['bootstrap'])
        
        # Check AU test statistics if available
        if 'au_test' in stats:
            self.assertIn('count', stats['au_test'])
            self.assertIn('not_rejected', stats['au_test'])
    
    def test_export_results_csv(self):
        """Test CSV export functionality"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add test data
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567,
            constraint_likelihoods={"Clade_1": -1245.678},
            au_results={"Clade_1": 0.023}
        )
        
        rp.calculate_decay_indices()
        
        # Export to CSV
        csv_file = self.temp_dir / "test_results.csv"
        rp.export_results_csv(csv_file)
        
        # Verify CSV file was created and has correct structure
        self.assertTrue(csv_file.exists())
        
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            
            self.assertGreater(len(rows), 0)
            
            # Check that expected columns are present
            if rows:
                first_row = rows[0]
                self.assertIn('clade_id', first_row)
                self.assertIn('taxa_count', first_row)
                self.assertIn('taxa', first_row)
    
    def test_clear_results(self):
        """Test clearing all results"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add some results
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567
        )
        
        self.assertTrue(rp.has_results())
        
        # Clear results
        rp.clear_results()
        
        self.assertFalse(rp.has_results())
        self.assertEqual(len(rp.decay_indices), 0)
        self.assertEqual(len(rp.ml_results), 0)
    
    def test_has_results(self):
        """Test result availability checking"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Initially no results
        self.assertFalse(rp.has_results())
        
        # Add some results
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567
        )
        
        self.assertTrue(rp.has_results())
        
        # Calculate decay indices
        rp.calculate_decay_indices()
        self.assertTrue(rp.has_results())
    
    def test_get_decay_indices(self):
        """Test decay indices retrieval"""
        rp = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Add results and calculate indices
        mock_tree = self._create_mock_tree()
        rp.add_ml_results(
            ml_tree=mock_tree,
            ml_likelihood=-1234.567
        )
        
        rp.calculate_decay_indices()
        
        # Get decay indices
        indices = rp.get_decay_indices()
        
        self.assertIsInstance(indices, dict)
        # Verify it's a copy (modifications shouldn't affect original)
        if indices:
            first_key = next(iter(indices.keys()))
            indices[first_key]['test_modification'] = 'test'
            
            original_indices = rp.get_decay_indices()
            self.assertNotIn('test_modification', original_indices[first_key])


class TestResultProcessorIntegration(unittest.TestCase):
    """Integration tests for ResultProcessor component"""
    
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
    
    def _create_realistic_tree(self):
        """Create a realistic phylogenetic tree"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        homo = Clade(name="Homo_sapiens", branch_length=0.01)
        pan = Clade(name="Pan_troglodytes", branch_length=0.02)
        gorilla = Clade(name="Gorilla_gorilla", branch_length=0.03)
        mus = Clade(name="Mus_musculus", branch_length=0.15)
        rattus = Clade(name="Rattus_norvegicus", branch_length=0.16)
        
        # Create nested clades
        hominins = Clade(branch_length=0.001)
        hominins.clades = [homo, pan]
        
        hominids = Clade(branch_length=0.005)
        hominids.clades = [hominins, gorilla]
        
        rodents = Clade(branch_length=0.01)
        rodents.clades = [mus, rattus]
        
        root = Clade(branch_length=0.0)
        root.clades = [hominids, rodents]
        
        return Tree(root=root)
    
    def test_complete_analysis_workflow(self):
        """Test complete result processing workflow"""
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load alignment
            dp.load_alignment(self.alignment_file, "fasta")
            
            # Create result processor
            rp = ResultProcessor(
                file_manager=fm,
                external_tools=et,
                data_processor=dp,
                output_style="unicode",
                viz_format="png",
                debug=True
            )
            
            # Create realistic results
            realistic_tree = self._create_realistic_tree()
            
            # Add ML results
            rp.add_ml_results(
                ml_tree=realistic_tree,
                ml_likelihood=-2345.678,
                constraint_likelihoods={
                    "Clade_1": -2350.123,
                    "Clade_2": -2355.456,
                    "Clade_3": -2348.789
                },
                au_results={
                    "Clade_1": 0.012,
                    "Clade_2": 0.089,
                    "Clade_3": 0.156
                }
            )
            
            # Add Bayesian results
            rp.add_bayesian_results(
                consensus_tree=realistic_tree,
                marginal_likelihood=-2345.678,
                constraint_analyses={
                    "Clade_1": {"marginal_likelihood": -2350.123},
                    "Clade_2": {"marginal_likelihood": -2355.456}
                },
                convergence_diagnostics={
                    "ess_min": 180,
                    "psrf_max": 1.05,
                    "asdsf_final": 0.008
                }
            )
            
            # Add bootstrap results
            rp.add_bootstrap_results(
                bootstrap_tree=realistic_tree,
                bootstrap_values={
                    "node_1": 98.0,
                    "node_2": 89.0,
                    "node_3": 76.0,
                    "node_4": 94.0
                },
                taxa_bootstrap_map={
                    frozenset(["Homo_sapiens", "Pan_troglodytes"]): 98.0,
                    frozenset(["Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla"]): 89.0,
                    frozenset(["Mus_musculus", "Rattus_norvegicus"]): 94.0
                }
            )
            
            # Add parsimony results
            rp.add_parsimony_results(
                parsimony_tree=realistic_tree,
                parsimony_score=156,
                decay_indices={
                    "Clade_1": 3.0,
                    "Clade_2": 1.5,
                    "Clade_3": 4.2
                }
            )
            
            # Calculate decay indices
            decay_indices = rp.calculate_decay_indices()
            
            # Verify comprehensive results
            self.assertGreater(len(decay_indices), 0)
            
            # Test formatted output
            results_file = fm.get_temp_path("complete_results.txt")
            rp.write_formatted_results(results_file)
            self.assertTrue(results_file.exists())
            
            # Test markdown report
            report_file = fm.get_temp_path("complete_report.md")
            rp.generate_detailed_report(report_file)
            self.assertTrue(report_file.exists())
            
            # Test CSV export
            csv_file = fm.get_temp_path("complete_results.csv")
            rp.export_results_csv(csv_file)
            self.assertTrue(csv_file.exists())
            
            # Test summary statistics
            stats = rp.get_summary_statistics()
            self.assertIn('total_clades', stats)
            self.assertIn('analysis_types', stats)
            
            # Verify multiple analysis types
            expected_types = {'ml', 'bayesian', 'parsimony'}
            actual_types = set(stats['analysis_types'])
            self.assertTrue(expected_types.issubset(actual_types))
            
            # Test that files have meaningful content
            results_content = results_file.read_text()
            self.assertIn("panDecay Analysis Results", results_content)
            self.assertIn("Bootstrap", results_content)
            
            report_content = report_file.read_text()
            self.assertIn("# panDecay Analysis Report", report_content)
            self.assertIn("Maximum Likelihood", report_content)
            
            # Verify CSV structure
            with open(csv_file, 'r') as f:
                csv_reader = csv.DictReader(f)
                csv_rows = list(csv_reader)
                
                if csv_rows:
                    self.assertIn('clade_id', csv_rows[0])
                    self.assertIn('taxa', csv_rows[0])
    
    def test_output_format_comparison(self):
        """Test different output formats produce consistent results"""
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            dp.load_alignment(self.alignment_file, "fasta")
            
            # Test both unicode and ascii output styles
            for style in ["unicode", "ascii"]:
                rp = ResultProcessor(
                    file_manager=fm,
                    external_tools=et,
                    data_processor=dp,
                    output_style=style,
                    debug=True
                )
                
                # Add test results
                realistic_tree = self._create_realistic_tree()
                rp.add_ml_results(
                    ml_tree=realistic_tree,
                    ml_likelihood=-2345.678,
                    au_results={"Clade_1": 0.023}
                )
                
                rp.add_bootstrap_results(
                    bootstrap_tree=realistic_tree,
                    bootstrap_values={"node_1": 85.0},
                    taxa_bootstrap_map={frozenset(["Homo_sapiens", "Pan_troglodytes"]): 85.0}
                )
                
                # Calculate and output
                rp.calculate_decay_indices()
                
                # Test symbol formatting
                symbol = rp._format_support_symbol(85.0, "bootstrap")
                
                if style == "unicode":
                    self.assertIn(symbol, ['●', '◐', '○', '·'])
                else:  # ascii
                    self.assertIn(symbol, ['*', '+', 'o', '.'])
                
                # Test formatted output
                output_file = fm.get_temp_path(f"results_{style}.txt")
                rp.write_formatted_results(output_file)
                self.assertTrue(output_file.exists())
                
                # Both should have same basic structure
                content = output_file.read_text()
                self.assertIn("panDecay Analysis Results", content)
                self.assertIn("Decay Indices", content)


def run_result_processor_tests():
    """Run the ResultProcessor test suite"""
    print("=" * 80)
    print("RESULT PROCESSOR COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted result processing functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestResultProcessor))
    suite.addTests(loader.loadTestsFromTestCase(TestResultProcessorIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("RESULT PROCESSOR TEST RESULTS")
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
    print(f"\nResultProcessor Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_result_processor_tests()
    sys.exit(0 if success else 1)