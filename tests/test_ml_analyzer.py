#!/usr/bin/env python3
"""
Test Suite for MLAnalyzer Component

Tests the extracted ML analysis functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock, mock_open

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis import FileManager, ExternalTools, DataProcessor
from core.analysis.ml_analyzer import (
    MLAnalyzer, MLAnalysisError, TreeBuildError, ConstraintTreeError, LikelihoodParsingError
)


class TestMLAnalyzer(unittest.TestCase):
    """Test MLAnalyzer component functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create test alignment
        self.test_alignment = """>seq1
ATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCG
>seq4
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
        
        # Create NEXUS file
        self.nexus_file = self.file_manager.get_temp_path("alignment.nex")
        self.data_processor.convert_to_nexus(self.nexus_file)
        
        self.model_commands = "lset nst=6 rmatrix=estimate basefreq=estimate rates=gamma;"
    
    def tearDown(self):
        """Clean up test fixtures"""
        self.file_manager.cleanup_all()
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test MLAnalyzer initialization"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands,
            debug=True
        )
        
        self.assertEqual(ml_analyzer.model_commands, self.model_commands)
        self.assertTrue(ml_analyzer.debug)
        self.assertIsNone(ml_analyzer.ml_tree)
        self.assertIsNone(ml_analyzer.ml_likelihood)
        self.assertEqual(ml_analyzer.constraint_trees, {})
    
    def test_initialization_no_alignment(self):
        """Test initialization fails without loaded alignment"""
        empty_processor = DataProcessor()
        
        with self.assertRaises(MLAnalysisError) as context:
            MLAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=empty_processor,
                model_commands=self.model_commands
            )
        
        self.assertIn("must have alignment loaded", str(context.exception))
    
    @patch('Bio.Phylo.read')
    def test_build_ml_tree_success(self, mock_phylo_read):
        """Test successful ML tree building"""
        # Mock successful PAUP* execution
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = "Log-likelihood = 1234.56"
            mock_paup.return_value = mock_result
            
            # Mock tree file creation
            tree_file = self.file_manager.get_temp_path("ml_tree.tre")
            with open(tree_file, 'w') as f:
                f.write("(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
            
            # Mock score file creation
            score_file = self.file_manager.get_temp_path("ml_score.txt")
            with open(score_file, 'w') as f:
                f.write("Tree    -ln L\n1    1234.56\n")
            
            # Mock phylo read
            mock_tree = MagicMock()
            mock_phylo_read.return_value = mock_tree
            
            # Initialize and test
            ml_analyzer = MLAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                model_commands=self.model_commands
            )
            
            result_tree = ml_analyzer.build_ml_tree()
            
            # Verify results
            self.assertEqual(result_tree, mock_tree)
            self.assertEqual(ml_analyzer.ml_tree, mock_tree)
            self.assertEqual(ml_analyzer.ml_likelihood, 1234.56)
            
            # Verify PAUP* was called correctly
            mock_paup.assert_called_once()
            args = mock_paup.call_args
            self.assertIn(self.model_commands, args[0][0].read_text())
    
    def test_build_ml_tree_with_starting_tree(self):
        """Test ML tree building with starting tree"""
        # Create starting tree file
        starting_tree = self.temp_dir / "starting.tre"
        with open(starting_tree, 'w') as f:
            f.write("(seq1:0.1,(seq2:0.1,(seq3:0.1,seq4:0.1):0.1):0.1);")
        
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = ""
            mock_paup.return_value = mock_result
            
            # Mock tree and score files
            tree_file = self.file_manager.get_temp_path("ml_tree.tre")
            score_file = self.file_manager.get_temp_path("ml_score.txt")
            with open(tree_file, 'w') as f:
                f.write("(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
            with open(score_file, 'w') as f:
                f.write("Tree    -ln L\n1    1000.00\n")
            
            with patch('Bio.Phylo.read') as mock_phylo:
                mock_phylo.return_value = MagicMock()
                
                ml_analyzer = MLAnalyzer(
                    file_manager=self.file_manager,
                    external_tools=self.external_tools,
                    data_processor=self.data_processor,
                    model_commands=self.model_commands,
                    starting_tree=starting_tree
                )
                
                ml_analyzer.build_ml_tree()
                
                # Verify starting tree commands were included
                script_content = mock_paup.call_args[0][0].read_text()
                self.assertIn("gettrees file=start_tree.tre", script_content)
                self.assertIn("hsearch start=current", script_content)
    
    def test_build_ml_tree_missing_tree_file(self):
        """Test ML tree building failure when tree file missing"""
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = ""
            mock_paup.return_value = mock_result
            
            # Don't create tree file to simulate failure
            
            ml_analyzer = MLAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                model_commands=self.model_commands
            )
            
            with self.assertRaises(TreeBuildError) as context:
                ml_analyzer.build_ml_tree()
            
            self.assertIn("not found or is empty", str(context.exception))
    
    def test_generate_constraint_tree_success(self):
        """Test successful constraint tree generation"""
        # Initialize ML analyzer with ML tree
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        ml_analyzer.ml_likelihood = 1000.0  # Set baseline likelihood
        
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_paup.return_value = mock_result
            
            # Create constraint tree and score files
            tree_file = self.file_manager.get_temp_path("constraint_tree_1.tre")
            score_file = self.file_manager.get_temp_path("constraint_score_1.txt")
            
            with open(tree_file, 'w') as f:
                f.write("((seq1:0.1,seq2:0.1):0.1,(seq3:0.1,seq4:0.1):0.1);")
            with open(score_file, 'w') as f:
                f.write("Tree    -ln L\n1    1010.50\n")
            
            # Test constraint generation
            tree_filename, likelihood = ml_analyzer.generate_constraint_tree(
                clade_taxa=["seq1", "seq2"],
                tree_idx=1
            )
            
            # Verify results
            self.assertEqual(tree_filename, "constraint_tree_1.tre")
            self.assertEqual(likelihood, 1010.50)
            
            # Check constraint was stored
            self.assertIn("Clade_1", ml_analyzer.constraint_trees)
            constraint_info = ml_analyzer.constraint_trees["Clade_1"]
            self.assertEqual(constraint_info['taxa'], ["seq1", "seq2"])
            self.assertEqual(constraint_info['constrained_lnl'], 1010.50)
            self.assertEqual(constraint_info['lnl_diff'], 10.50)  # 1010.50 - 1000.0
    
    def test_generate_constraint_tree_all_taxa(self):
        """Test constraint tree generation with all taxa (should skip)"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Try to constrain all taxa (invalid)
        all_taxa = self.data_processor.get_taxon_names()
        tree_filename, likelihood = ml_analyzer.generate_constraint_tree(
            clade_taxa=all_taxa,
            tree_idx=1
        )
        
        # Should return None, None
        self.assertIsNone(tree_filename)
        self.assertIsNone(likelihood)
    
    def test_generate_constraint_tree_empty_taxa(self):
        """Test constraint tree generation with empty taxa list"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        tree_filename, likelihood = ml_analyzer.generate_constraint_tree(
            clade_taxa=[],
            tree_idx=1
        )
        
        # Should return None, None
        self.assertIsNone(tree_filename)
        self.assertIsNone(likelihood)
    
    def test_parse_likelihood_from_score_file(self):
        """Test likelihood parsing from score file"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Create test score file
        score_content = """Tree    -ln L    Diff -ln L
1    12345.67    0.00
2    12350.00    4.33
"""
        score_file = self.file_manager.get_temp_path("test_scores.txt")
        with open(score_file, 'w') as f:
            f.write(score_content)
        
        likelihood = ml_analyzer._parse_likelihood_from_score_file(score_file)
        self.assertEqual(likelihood, 12345.67)
    
    def test_parse_likelihood_from_score_file_missing(self):
        """Test likelihood parsing from missing score file"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        missing_file = self.file_manager.get_temp_path("missing_scores.txt")
        likelihood = ml_analyzer._parse_likelihood_from_score_file(missing_file)
        self.assertIsNone(likelihood)
    
    def test_parse_likelihood_from_log(self):
        """Test likelihood parsing from log content"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        log_content = """
PAUP* output...
Tree score = 1234.567
More output...
-ln L = 5678.901
Final output...
"""
        
        likelihood = ml_analyzer._parse_likelihood_from_log(log_content)
        self.assertEqual(likelihood, 5678.901)  # Should get the last match
    
    def test_parse_likelihood_from_log_no_match(self):
        """Test likelihood parsing from log with no matches"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        log_content = "No likelihood values here"
        likelihood = ml_analyzer._parse_likelihood_from_log(log_content)
        self.assertIsNone(likelihood)
    
    def test_clean_newick_tree(self):
        """Test Newick tree cleaning"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Create tree file with metadata
        tree_content = "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1); [Some metadata here]"
        tree_file = self.file_manager.get_temp_path("messy_tree.tre")
        with open(tree_file, 'w') as f:
            f.write(tree_content)
        
        cleaned_path = ml_analyzer._clean_newick_tree(tree_file, delete_cleaned=False)
        
        # Verify cleaned content
        cleaned_content = cleaned_path.read_text()
        self.assertEqual(cleaned_content, "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
        self.assertNotIn("metadata", cleaned_content.lower())
    
    def test_get_analysis_summary(self):
        """Test analysis summary generation"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Set some results
        ml_analyzer.ml_tree = MagicMock()
        ml_analyzer.ml_likelihood = 1000.0
        ml_analyzer.constraint_trees = {
            'Clade_1': {'taxa': ['seq1', 'seq2'], 'constrained_lnl': 1010.0},
            'Clade_2': {'taxa': ['seq3', 'seq4'], 'constrained_lnl': None}
        }
        
        summary = ml_analyzer.get_analysis_summary()
        
        self.assertTrue(summary['ml_tree_available'])
        self.assertEqual(summary['ml_likelihood'], 1000.0)
        self.assertEqual(summary['num_constraint_trees'], 2)
        self.assertEqual(summary['constraint_trees_with_likelihood'], 1)
        self.assertEqual(summary['constraint_clades'], ['Clade_1', 'Clade_2'])
    
    def test_calculate_likelihood_differences(self):
        """Test likelihood difference calculation"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Set up test data
        ml_analyzer.ml_likelihood = 1000.0
        ml_analyzer.constraint_trees = {
            'Clade_1': {'constrained_lnl': 1010.0},
            'Clade_2': {'constrained_lnl': 1005.5},
            'Clade_3': {'constrained_lnl': None}  # Should be skipped
        }
        
        differences = ml_analyzer.calculate_likelihood_differences()
        
        self.assertEqual(differences['Clade_1'], 10.0)
        self.assertEqual(differences['Clade_2'], 5.5)
        self.assertNotIn('Clade_3', differences)
    
    def test_validate_ml_analysis(self):
        """Test ML analysis validation"""
        ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands=self.model_commands
        )
        
        # Initially invalid
        self.assertFalse(ml_analyzer.validate_ml_analysis())
        
        # Add tree but no likelihood
        ml_analyzer.ml_tree = MagicMock()
        self.assertFalse(ml_analyzer.validate_ml_analysis())
        
        # Add likelihood
        ml_analyzer.ml_likelihood = 1000.0
        self.assertTrue(ml_analyzer.validate_ml_analysis())


class TestMLAnalyzerIntegration(unittest.TestCase):
    """Integration tests for MLAnalyzer component"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create realistic alignment
        self.realistic_alignment = """>Homo_sapiens
ATCGATCGATCGNNNN----ATCGATCG
>Pan_troglodytes
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
    
    def test_complete_ml_workflow(self):
        """Test complete ML analysis workflow"""
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load and process alignment
            dp.load_alignment(self.alignment_file, "fasta")
            nexus_file = fm.get_temp_path("alignment.nex")
            dp.convert_to_nexus(nexus_file)
            
            # Mock PAUP* execution for ML tree building
            with patch.object(et, 'run_paup_command_file') as mock_paup:
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_result.stdout = "Tree likelihood = 2000.0"
                mock_paup.return_value = mock_result
                
                # Create mock tree and score files
                tree_file = fm.get_temp_path("ml_tree.tre")
                score_file = fm.get_temp_path("ml_score.txt")
                
                with open(tree_file, 'w') as f:
                    f.write("(Homo_sapiens:0.1,(Pan_troglodytes:0.1,(Mus_musculus:0.1,Rattus_norvegicus:0.1):0.1):0.1);")
                with open(score_file, 'w') as f:
                    f.write("Tree    -ln L\n1    2000.0\n")
                
                # Mock phylo.read
                with patch('Bio.Phylo.read') as mock_phylo:
                    mock_tree = MagicMock()
                    mock_phylo.return_value = mock_tree
                    
                    # Initialize ML analyzer
                    ml_analyzer = MLAnalyzer(
                        file_manager=fm,
                        external_tools=et,
                        data_processor=dp,
                        model_commands="lset nst=6 rmatrix=estimate basefreq=estimate rates=gamma;",
                        debug=True
                    )
                    
                    # Build ML tree
                    result_tree = ml_analyzer.build_ml_tree()
                    
                    # Verify results
                    self.assertEqual(result_tree, mock_tree)
                    self.assertEqual(ml_analyzer.get_ml_likelihood(), 2000.0)
                    self.assertTrue(ml_analyzer.validate_ml_analysis())
                    
                    # Test constraint tree generation
                    clade_taxa = ["Homo_sapiens", "Pan_troglodytes"]
                    
                    # Create constraint files
                    constraint_tree = fm.get_temp_path("constraint_tree_1.tre") 
                    constraint_score = fm.get_temp_path("constraint_score_1.txt")
                    
                    with open(constraint_tree, 'w') as f:
                        f.write("((Homo_sapiens:0.1,Pan_troglodytes:0.1):0.1,(Mus_musculus:0.1,Rattus_norvegicus:0.1):0.1);")
                    with open(constraint_score, 'w') as f:
                        f.write("Tree    -ln L\n1    2015.5\n")
                    
                    tree_filename, likelihood = ml_analyzer.generate_constraint_tree(
                        clade_taxa=clade_taxa,
                        tree_idx=1
                    )
                    
                    # Verify constraint results
                    self.assertEqual(tree_filename, "constraint_tree_1.tre")
                    self.assertEqual(likelihood, 2015.5)
                    
                    # Check likelihood differences
                    differences = ml_analyzer.calculate_likelihood_differences()
                    self.assertEqual(differences['Clade_1'], 15.5)  # 2015.5 - 2000.0
                    
                    # Get analysis summary
                    summary = ml_analyzer.get_analysis_summary()
                    self.assertTrue(summary['ml_tree_available'])
                    self.assertEqual(summary['ml_likelihood'], 2000.0)
                    self.assertEqual(summary['num_constraint_trees'], 1)


def run_ml_analyzer_tests():
    """Run the MLAnalyzer test suite"""
    print("=" * 80)
    print("ML ANALYZER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted ML analysis functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestMLAnalyzer))
    suite.addTests(loader.loadTestsFromTestCase(TestMLAnalyzerIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\\n" + "=" * 80)
    print("ML ANALYZER TEST RESULTS")
    print("=" * 80)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    
    if result.failures:
        print("\\nFAILURES:")
        for test, traceback in result.failures:
            print(f"- {test}: {traceback}")
    
    if result.errors:
        print("\\nERRORS:")
        for test, traceback in result.errors:
            print(f"- {test}: {traceback}")
    
    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\\nMLAnalyzer Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_ml_analyzer_tests()
    sys.exit(0 if success else 1)