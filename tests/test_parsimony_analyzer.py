#!/usr/bin/env python3
"""
Test Suite for ParsimonyAnalyzer Component

Tests the extracted parsimony analysis functionality to ensure
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
from core.analysis.parsimony_analyzer import (
    ParsimonyAnalyzer, ParsimonyAnalysisError, TreeBuildError, 
    ConstraintAnalysisError, ScoreParsingError
)


class TestParsimonyAnalyzer(unittest.TestCase):
    """Test ParsimonyAnalyzer component functionality"""
    
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
    
    def tearDown(self):
        """Clean up test fixtures"""
        self.file_manager.cleanup_all()
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test ParsimonyAnalyzer initialization"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            search_replicates=5,
            swap_method="tbr",
            debug=True
        )
        
        self.assertEqual(pa.search_replicates, 5)
        self.assertEqual(pa.swap_method, "tbr")
        self.assertTrue(pa.multrees)
        self.assertTrue(pa.debug)
        self.assertIsNone(pa.parsimony_tree)
        self.assertIsNone(pa.parsimony_score)
        self.assertEqual(pa.constraint_analyses, {})
        self.assertEqual(pa.decay_indices, {})
    
    def test_initialization_invalid_swap_method(self):
        """Test initialization with invalid swap method"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            swap_method="invalid"
        )
        
        # Should default to tbr
        self.assertEqual(pa.swap_method, "tbr")
    
    def test_initialization_no_alignment(self):
        """Test initialization fails without loaded alignment"""
        empty_processor = DataProcessor()
        
        with self.assertRaises(ParsimonyAnalysisError) as context:
            ParsimonyAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=empty_processor
            )
        
        self.assertIn("must have alignment loaded", str(context.exception))
    
    @patch('Bio.Phylo.read')
    def test_build_parsimony_tree_success(self, mock_phylo_read):
        """Test successful parsimony tree building"""
        # Mock successful PAUP* execution
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = "Tree length = 42"
            mock_paup.return_value = mock_result
            
            # Mock tree file creation
            tree_file = self.file_manager.get_temp_path("parsimony_trees.tre")
            with open(tree_file, 'w') as f:
                f.write("(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
            
            # Mock score file creation
            score_file = self.file_manager.get_temp_path("parsimony_score.txt")
            with open(score_file, 'w') as f:
                f.write("Tree    Length\n1    42\n")
            
            # Mock phylo read
            mock_tree = MagicMock()
            mock_phylo_read.return_value = mock_tree
            
            # Initialize and test
            pa = ParsimonyAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                search_replicates=5
            )
            
            result_tree = pa.build_parsimony_tree()
            
            # Verify results
            self.assertEqual(result_tree, mock_tree)
            self.assertEqual(pa.parsimony_tree, mock_tree)
            self.assertEqual(pa.parsimony_score, 42.0)
            
            # Verify PAUP* was called correctly
            mock_paup.assert_called_once()
            args = mock_paup.call_args
            script_content = args[0][0].read_text()
            self.assertIn("set criterion=parsimony;", script_content)
            self.assertIn("nreps=5", script_content)
            self.assertIn("swap=tbr", script_content)
    
    def test_build_parsimony_tree_missing_tree_file(self):
        """Test parsimony tree building failure when tree file missing"""
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = ""
            mock_paup.return_value = mock_result
            
            # Don't create tree file to simulate failure
            
            pa = ParsimonyAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor
            )
            
            with self.assertRaises(TreeBuildError) as context:
                pa.build_parsimony_tree()
            
            self.assertIn("not found or is empty", str(context.exception))
    
    def test_extract_score_from_text_pattern1(self):
        """Test score extraction pattern 1: 'Length = 123'"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "Tree analysis complete\nLength = 45\nOther output"
        score = pa._extract_score_from_text(text)
        self.assertEqual(score, 45.0)
    
    def test_extract_score_from_text_pattern2(self):
        """Test score extraction pattern 2: 'Tree length = 123'"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "PAUP output\nTree length = 67.5\nMore output"
        score = pa._extract_score_from_text(text)
        self.assertEqual(score, 67.5)
    
    def test_extract_score_from_text_pattern3(self):
        """Test score extraction pattern 3: 'Score = 123'"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "Analysis results\nScore = 123\nEnd"
        score = pa._extract_score_from_text(text)
        self.assertEqual(score, 123.0)
    
    def test_extract_score_from_text_pattern4(self):
        """Test score extraction pattern 4: 'Length: 123'"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "Results summary\nLength: 89\nComplete"
        score = pa._extract_score_from_text(text)
        self.assertEqual(score, 89.0)
    
    def test_extract_score_from_text_pattern5(self):
        """Test score extraction pattern 5: keyword lines with numbers"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "PAUP analysis\nTotal tree length is 156 steps\nDone"
        score = pa._extract_score_from_text(text)
        self.assertEqual(score, 156.0)
    
    def test_extract_score_from_text_no_match(self):
        """Test score extraction when no patterns match"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        text = "No scores here, just random text"
        score = pa._extract_score_from_text(text)
        self.assertIsNone(score)
    
    def test_parse_parsimony_score_from_file(self):
        """Test parsimony score parsing from score file"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Create test score file
        score_content = """Tree    Length
1    234
2    235
"""
        score_file = self.file_manager.get_temp_path("test_scores.txt")
        with open(score_file, 'w') as f:
            f.write(score_content)
        
        score = pa._parse_parsimony_score(score_file)
        self.assertEqual(score, 234.0)
    
    def test_parse_parsimony_score_from_stdout(self):
        """Test parsimony score parsing from PAUP* stdout"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        stdout_content = "PAUP analysis\nLength = 78\nComplete"
        score = pa._parse_parsimony_score(None, stdout_content)
        self.assertEqual(score, 78.0)
    
    def test_parse_parsimony_score_missing_file(self):
        """Test parsimony score parsing with missing file"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        missing_file = self.file_manager.get_temp_path("missing_scores.txt")
        score = pa._parse_parsimony_score(missing_file)
        self.assertIsNone(score)
    
    def test_generate_constraint_analysis_success(self):
        """Test successful constraint analysis generation"""
        # Initialize parsimony analyzer with parsimony score
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        pa.parsimony_score = 50.0  # Set baseline score
        
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = "Length = 55"
            mock_paup.return_value = mock_result
            
            # Create constraint tree and score files
            tree_file = self.file_manager.get_temp_path("pars_constraint_clade1.tre")
            score_file = self.file_manager.get_temp_path("pars_constraint_score_clade1.txt")
            
            with open(tree_file, 'w') as f:
                f.write("((seq1:0.1,seq2:0.1):0.1,(seq3:0.1,seq4:0.1):0.1);")
            with open(score_file, 'w') as f:
                f.write("Tree    Length\n1    55\n")
            
            # Test constraint generation
            tree_filename, tree_length = pa.generate_constraint_analysis(
                clade_taxa=["seq1", "seq2"],
                constraint_id="clade1"
            )
            
            # Verify results
            self.assertEqual(tree_filename, "pars_constraint_clade1.tre")
            self.assertEqual(tree_length, 55.0)
            
            # Check constraint was stored
            self.assertIn("clade1", pa.constraint_analyses)
            constraint_info = pa.constraint_analyses["clade1"]
            self.assertEqual(constraint_info['taxa'], ["seq1", "seq2"])
            self.assertEqual(constraint_info['constrained_score'], 55.0)
            self.assertEqual(constraint_info['decay_value'], 5.0)  # 55 - 50
            
            # Check decay index was calculated
            self.assertEqual(pa.get_decay_indices()["clade1"], 5.0)
    
    def test_generate_constraint_analysis_all_taxa(self):
        """Test constraint analysis with all taxa (should skip)"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Try to constrain all taxa (invalid)
        all_taxa = self.data_processor.get_taxon_names()
        tree_filename, tree_length = pa.generate_constraint_analysis(
            clade_taxa=all_taxa,
            constraint_id="all_taxa_test"
        )
        
        # Should return None, None
        self.assertIsNone(tree_filename)
        self.assertIsNone(tree_length)
    
    def test_generate_constraint_analysis_empty_taxa(self):
        """Test constraint analysis with empty taxa list"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        tree_filename, tree_length = pa.generate_constraint_analysis(
            clade_taxa=[],
            constraint_id="empty_test"
        )
        
        # Should return None, None
        self.assertIsNone(tree_filename)
        self.assertIsNone(tree_length)
    
    def test_calculate_bremer_support(self):
        """Test Bremer support calculation"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Set up test data
        pa.parsimony_score = 40.0
        pa.constraint_analyses = {
            'clade1': {'constrained_score': 45.0},
            'clade2': {'constrained_score': 42.0},
            'clade3': {'constrained_score': None}  # Should return None
        }
        
        # Test direct calculation
        bremer1 = pa.calculate_bremer_support('clade1')
        self.assertEqual(bremer1, 5.0)
        
        bremer2 = pa.calculate_bremer_support('clade2')
        self.assertEqual(bremer2, 2.0)
        
        bremer3 = pa.calculate_bremer_support('clade3')
        self.assertIsNone(bremer3)
        
        # Test non-existent constraint
        bremer_none = pa.calculate_bremer_support('nonexistent')
        self.assertIsNone(bremer_none)
    
    def test_clean_newick_tree(self):
        """Test Newick tree cleaning"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Create tree file with metadata
        tree_content = "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1); [Some metadata here]"
        tree_file = self.file_manager.get_temp_path("messy_tree.tre")
        with open(tree_file, 'w') as f:
            f.write(tree_content)
        
        cleaned_path = pa._clean_newick_tree(tree_file, delete_cleaned=False)
        
        # Verify cleaned content
        cleaned_content = cleaned_path.read_text()
        self.assertEqual(cleaned_content, "(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
        self.assertNotIn("metadata", cleaned_content.lower())
    
    def test_validate_parsimony_analysis(self):
        """Test parsimony analysis validation"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Initially invalid
        self.assertFalse(pa.validate_parsimony_analysis())
        
        # Add tree but no score
        pa.parsimony_tree = MagicMock()
        self.assertFalse(pa.validate_parsimony_analysis())
        
        # Add score
        pa.parsimony_score = 50.0
        self.assertTrue(pa.validate_parsimony_analysis())
    
    def test_get_analysis_summary(self):
        """Test analysis summary generation"""
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            search_replicates=5,
            swap_method="spr",
            multrees=False
        )
        
        # Set some results
        pa.parsimony_tree = MagicMock()
        pa.parsimony_score = 40.0
        pa.constraint_analyses = {
            'clade1': {'constrained_score': 45.0},
            'clade2': {'constrained_score': 42.0}
        }
        pa.decay_indices = {
            'clade1': 5.0,
            'clade2': 2.0
        }
        
        summary = pa.get_analysis_summary()
        
        self.assertTrue(summary['parsimony_tree_available'])
        self.assertEqual(summary['parsimony_score'], 40.0)
        self.assertEqual(summary['num_constraint_analyses'], 2)
        self.assertEqual(summary['num_decay_indices'], 2)
        self.assertEqual(summary['search_replicates'], 5)
        self.assertEqual(summary['swap_method'], "spr")
        self.assertFalse(summary['multrees'])
        self.assertEqual(summary['constraint_ids'], ['clade1', 'clade2'])
        self.assertEqual(summary['available_decay_values'], ['clade1', 'clade2'])
        
        # Check decay statistics
        decay_stats = summary['decay_stats']
        self.assertEqual(decay_stats['min'], 2.0)
        self.assertEqual(decay_stats['max'], 5.0)
        self.assertEqual(decay_stats['mean'], 3.5)
        self.assertEqual(decay_stats['count'], 2)
    
    def test_swap_method_validation(self):
        """Test different swap method settings"""
        # Test valid swap methods
        for method in ["tbr", "spr", "nni"]:
            pa = ParsimonyAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                swap_method=method
            )
            self.assertEqual(pa.swap_method, method)
        
        # Test case insensitive
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            swap_method="TBR"
        )
        self.assertEqual(pa.swap_method, "tbr")
    
    def test_multrees_setting(self):
        """Test multrees parameter setting"""
        # Test multrees enabled
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            multrees=True
        )
        self.assertTrue(pa.multrees)
        
        # Test multrees disabled
        pa = ParsimonyAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            multrees=False
        )
        self.assertFalse(pa.multrees)


class TestParsimonyAnalyzerIntegration(unittest.TestCase):
    """Integration tests for ParsimonyAnalyzer component"""
    
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
    
    @patch('Bio.Phylo.read')
    def test_complete_parsimony_workflow(self, mock_phylo_read):
        """Test complete parsimony analysis workflow"""
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load and process alignment
            dp.load_alignment(self.alignment_file, "fasta")
            nexus_file = fm.get_temp_path("alignment.nex")
            dp.convert_to_nexus(nexus_file)
            
            # Mock PAUP* execution for parsimony tree building
            with patch.object(et, 'run_paup_command_file') as mock_paup:
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_result.stdout = "Tree length = 100"
                mock_paup.return_value = mock_result
                
                # Create mock tree and score files
                tree_file = fm.get_temp_path("parsimony_trees.tre")
                score_file = fm.get_temp_path("parsimony_score.txt")
                
                with open(tree_file, 'w') as f:
                    f.write("(Homo_sapiens,(Pan_troglodytes,(Mus_musculus,Rattus_norvegicus)));")
                with open(score_file, 'w') as f:
                    f.write("Tree    Length\n1    100\n")
                
                # Mock phylo.read
                mock_tree = MagicMock()
                mock_phylo_read.return_value = mock_tree
                
                # Initialize parsimony analyzer
                pa = ParsimonyAnalyzer(
                    file_manager=fm,
                    external_tools=et,
                    data_processor=dp,
                    search_replicates=3,
                    swap_method="tbr",
                    debug=True
                )
                
                # Build parsimony tree
                result_tree = pa.build_parsimony_tree()
                
                # Verify results
                self.assertEqual(result_tree, mock_tree)
                self.assertEqual(pa.get_parsimony_score(), 100.0)
                self.assertTrue(pa.validate_parsimony_analysis())
                
                # Verify PAUP* was called correctly
                mock_paup.assert_called_once()
                
                # Test constraint analysis
                clade_taxa = ["Homo_sapiens", "Pan_troglodytes"]
                
                # Create constraint files
                constraint_tree = fm.get_temp_path("pars_constraint_clade1.tre")
                constraint_score = fm.get_temp_path("pars_constraint_score_clade1.txt")
                
                with open(constraint_tree, 'w') as f:
                    f.write("((Homo_sapiens,Pan_troglodytes),(Mus_musculus,Rattus_norvegicus));")
                with open(constraint_score, 'w') as f:
                    f.write("Tree    Length\n1    105\n")
                
                tree_filename, tree_length = pa.generate_constraint_analysis(
                    clade_taxa=clade_taxa,
                    constraint_id="clade1"
                )
                
                # Verify constraint results
                self.assertEqual(tree_filename, "pars_constraint_clade1.tre")
                self.assertEqual(tree_length, 105.0)
                
                # Check Bremer support
                bremer_support = pa.calculate_bremer_support("clade1")
                self.assertEqual(bremer_support, 5.0)  # 105 - 100
                
                # Get analysis summary
                summary = pa.get_analysis_summary()
                self.assertTrue(summary['parsimony_tree_available'])
                self.assertEqual(summary['parsimony_score'], 100.0)
                self.assertEqual(summary['num_constraint_analyses'], 1)
                self.assertEqual(summary['decay_stats']['min'], 5.0)
                self.assertEqual(summary['decay_stats']['max'], 5.0)


def run_parsimony_analyzer_tests():
    """Run the ParsimonyAnalyzer test suite"""
    print("=" * 80)
    print("PARSIMONY ANALYZER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted parsimony analysis functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestParsimonyAnalyzer))
    suite.addTests(loader.loadTestsFromTestCase(TestParsimonyAnalyzerIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("PARSIMONY ANALYZER TEST RESULTS")
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
    print(f"\nParsimonyAnalyzer Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_parsimony_analyzer_tests()
    sys.exit(0 if success else 1)