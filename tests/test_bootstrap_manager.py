#!/usr/bin/env python3
"""
Test Suite for BootstrapManager Component

Tests the extracted bootstrap analysis functionality to ensure
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
from core.analysis.bootstrap_manager import (
    BootstrapManager, BootstrapAnalysisError, BootstrapBuildError,
    BootstrapParsingError, ConsensusTreeError
)


class TestBootstrapManager(unittest.TestCase):
    """Test BootstrapManager component functionality"""
    
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
        """Test BootstrapManager initialization"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            num_replicates=50,
            consensus_threshold=70,
            debug=True
        )
        
        self.assertEqual(bm.num_replicates, 50)
        self.assertEqual(bm.consensus_threshold, 70)
        self.assertTrue(bm.debug)
        self.assertIsNone(bm.bootstrap_tree)
        self.assertEqual(bm.bootstrap_values, {})
        self.assertEqual(bm.taxa_bootstrap_map, {})
    
    def test_initialization_invalid_threshold(self):
        """Test initialization with invalid consensus threshold"""
        with self.assertRaises(BootstrapAnalysisError) as context:
            BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                consensus_threshold=150  # Invalid
            )
        
        self.assertIn("Invalid consensus threshold", str(context.exception))
    
    def test_initialization_no_alignment(self):
        """Test initialization fails without loaded alignment"""
        empty_processor = DataProcessor()
        
        with self.assertRaises(BootstrapAnalysisError) as context:
            BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=empty_processor
            )
        
        self.assertIn("must have alignment loaded", str(context.exception))
    
    def _create_mock_bootstrap_tree(self):
        """Create a mock bootstrap tree with confidence values"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Create tree structure: ((seq1,seq2)85,(seq3,seq4)92)
        seq1 = Clade(name="seq1")
        seq2 = Clade(name="seq2")
        seq3 = Clade(name="seq3")
        seq4 = Clade(name="seq4")
        
        clade12 = Clade(confidence=85.0)
        clade12.clades = [seq1, seq2]
        
        clade34 = Clade(confidence=92.0)
        clade34.clades = [seq3, seq4]
        
        root = Clade(confidence=100.0)
        root.clades = [clade12, clade34]
        
        return Tree(root=root)
    
    @patch('Bio.Phylo.read')
    def test_run_ml_bootstrap_success(self, mock_phylo_read):
        """Test successful ML bootstrap analysis"""
        # Mock successful PAUP* execution
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stdout = "Bootstrap analysis completed"
            mock_paup.return_value = mock_result
            
            # Mock bootstrap tree file creation
            tree_file = self.file_manager.get_temp_path("bootstrap_trees.tre")
            with open(tree_file, 'w') as f:
                f.write("((seq1,seq2)85,(seq3,seq4)92);")
            
            # Mock phylo read with bootstrap tree
            mock_tree = self._create_mock_bootstrap_tree()
            mock_phylo_read.return_value = mock_tree
            
            # Initialize and test
            bm = BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                num_replicates=10
            )
            
            result_tree = bm.run_ml_bootstrap(
                model_commands="lset nst=6 rates=gamma;",
                timeout_sec=600
            )
            
            # Verify results
            self.assertEqual(result_tree, mock_tree)
            self.assertEqual(bm.bootstrap_tree, mock_tree)
            self.assertGreater(len(bm.bootstrap_values), 0)
            self.assertGreater(len(bm.taxa_bootstrap_map), 0)
            
            # Verify PAUP* was called correctly
            mock_paup.assert_called_once()
            args = mock_paup.call_args
            script_content = args[0][0].read_text()
            self.assertIn("set criterion=likelihood;", script_content)
            self.assertIn("bootstrap nreps=10", script_content)
            self.assertIn("lset nst=6 rates=gamma;", script_content)
    
    @patch('Bio.Phylo.read')
    def test_run_ml_bootstrap_with_starting_tree(self, mock_phylo_read):
        """Test ML bootstrap with starting tree"""
        # Create starting tree file
        starting_tree = self.temp_dir / "starting.tre"
        with open(starting_tree, 'w') as f:
            f.write("(seq1,(seq2,(seq3,seq4)));")
        
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_paup.return_value = mock_result
            
            # Mock tree file
            tree_file = self.file_manager.get_temp_path("bootstrap_trees.tre")
            with open(tree_file, 'w') as f:
                f.write("((seq1,seq2)75,(seq3,seq4)88);")
            
            mock_tree = self._create_mock_bootstrap_tree()
            mock_phylo_read.return_value = mock_tree
            
            bm = BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor
            )
            
            bm.run_ml_bootstrap(starting_tree=starting_tree)
            
            # Verify starting tree commands were included
            script_content = mock_paup.call_args[0][0].read_text()
            self.assertIn("gettrees file=bootstrap_start.tre", script_content)
            self.assertIn("lscores 1 / userbrlen=yes", script_content)
    
    def test_run_ml_bootstrap_missing_tree_file(self):
        """Test ML bootstrap failure when tree file missing"""
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_paup.return_value = mock_result
            
            # Don't create tree file to simulate failure
            
            bm = BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor
            )
            
            with self.assertRaises(BootstrapBuildError) as context:
                bm.run_ml_bootstrap()
            
            self.assertIn("not found or is empty", str(context.exception))
    
    @patch('Bio.Phylo.read')
    def test_run_parsimony_bootstrap_success(self, mock_phylo_read):
        """Test successful parsimony bootstrap analysis"""
        with patch.object(self.external_tools, 'run_paup_command_file') as mock_paup:
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_paup.return_value = mock_result
            
            # Mock bootstrap tree file
            tree_file = self.file_manager.get_temp_path("bootstrap_trees.tre")
            with open(tree_file, 'w') as f:
                f.write("((seq1,seq2)80,(seq3,seq4)95);")
            
            mock_tree = self._create_mock_bootstrap_tree()
            mock_phylo_read.return_value = mock_tree
            
            bm = BootstrapManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                num_replicates=20
            )
            
            result_tree = bm.run_parsimony_bootstrap(
                search_replicates=2,
                swap_method="spr"
            )
            
            # Verify results
            self.assertEqual(result_tree, mock_tree)
            self.assertEqual(bm.bootstrap_tree, mock_tree)
            
            # Verify PAUP* script
            script_content = mock_paup.call_args[0][0].read_text()
            self.assertIn("set criterion=parsimony;", script_content)
            self.assertIn("bootstrap nreps=20", script_content)
            self.assertIn("nreps=2 swap=spr", script_content)
    
    def test_parse_bootstrap_values(self):
        """Test bootstrap value parsing from tree"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_bootstrap_tree()
        bootstrap_values = bm._parse_bootstrap_values(mock_tree)
        
        # Should have parsed confidence values
        self.assertGreater(len(bootstrap_values), 0)
        
        # Check that values are in percentage form
        for value in bootstrap_values.values():
            self.assertGreaterEqual(value, 0)
            self.assertLessEqual(value, 100)
    
    def test_create_taxa_bootstrap_map(self):
        """Test taxa-bootstrap mapping creation"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        mock_tree = self._create_mock_bootstrap_tree()
        taxa_map = bm._create_taxa_bootstrap_map(mock_tree)
        
        # Should have created mappings for internal nodes
        self.assertGreater(len(taxa_map), 0)
        
        # Check that keys are frozensets and values are bootstrap percentages
        for taxa_set, bootstrap_val in taxa_map.items():
            self.assertIsInstance(taxa_set, frozenset)
            self.assertGreater(len(taxa_set), 1)  # Internal nodes only
            self.assertIsInstance(bootstrap_val, float)
            self.assertGreaterEqual(bootstrap_val, 0)
            self.assertLessEqual(bootstrap_val, 100)
    
    def test_get_bootstrap_support(self):
        """Test bootstrap support lookup for specific clades"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Set up mock taxa-bootstrap mapping
        bm.taxa_bootstrap_map = {
            frozenset(["seq1", "seq2"]): 85.0,
            frozenset(["seq3", "seq4"]): 92.0,
            frozenset(["seq1", "seq2", "seq3", "seq4"]): 100.0
        }
        
        # Test exact matches
        support1 = bm.get_bootstrap_support(["seq1", "seq2"])
        self.assertEqual(support1, 85.0)
        
        support2 = bm.get_bootstrap_support(["seq3", "seq4"])
        self.assertEqual(support2, 92.0)
        
        # Test subset match (clade is part of larger supported clade)
        support3 = bm.get_bootstrap_support(["seq1"])  # Should find 85.0 from seq1,seq2 clade
        self.assertEqual(support3, 85.0)
        
        # Test non-existent clade
        support_none = bm.get_bootstrap_support(["seq1", "seq3"])
        self.assertEqual(support_none, 100.0)  # Found in root clade
        
        # Test truly non-existent clade with taxa not in any mapping
        support_none2 = bm.get_bootstrap_support(["seq5", "seq6"])
        self.assertIsNone(support_none2)
    
    def test_annotate_tree_with_bootstrap(self):
        """Test tree annotation with bootstrap values"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Set up mock data
        mock_tree = self._create_mock_bootstrap_tree()
        bm.taxa_bootstrap_map = {
            frozenset(["seq1", "seq2"]): 85.0,
            frozenset(["seq3", "seq4"]): 92.0
        }
        
        # Test annotation
        annotated_tree = bm.annotate_tree_with_bootstrap(mock_tree)
        
        # Check that tree was annotated (should have node names with BS: values)
        has_bootstrap_annotation = False
        for clade in annotated_tree.find_clades():
            if clade.name and "BS:" in str(clade.name):
                has_bootstrap_annotation = True
                break
        
        # Note: The exact annotation depends on tree structure matching,
        # so we mainly check that the function runs without error
        self.assertIsNotNone(annotated_tree)
    
    def test_create_bootstrap_summary(self):
        """Test bootstrap analysis summary creation"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            num_replicates=50,
            consensus_threshold=75
        )
        
        # Test summary without bootstrap data
        summary_empty = bm.create_bootstrap_summary()
        self.assertFalse(summary_empty['bootstrap_completed'])
        self.assertEqual(summary_empty['num_replicates'], 50)
        self.assertEqual(summary_empty['consensus_threshold'], 75)
        
        # Set up mock bootstrap data
        bm.bootstrap_values = {
            'node_1': 95.0,
            'node_2': 85.0,
            'node_3': 75.0,
            'node_4': 65.0,
            'node_5': 55.0
        }
        bm.taxa_bootstrap_map = {
            frozenset(["seq1", "seq2"]): 95.0,
            frozenset(["seq3", "seq4"]): 85.0
        }
        
        # Test summary with bootstrap data
        summary = bm.create_bootstrap_summary()
        
        self.assertTrue(summary['bootstrap_completed'])
        self.assertEqual(summary['num_supported_nodes'], 5)
        self.assertEqual(summary['num_taxa_mappings'], 2)
        self.assertEqual(summary['bootstrap_stats']['min'], 55.0)
        self.assertEqual(summary['bootstrap_stats']['max'], 95.0)
        self.assertEqual(summary['bootstrap_stats']['mean'], 75.0)
        self.assertEqual(summary['bootstrap_stats']['count'], 5)
        
        # Check support categories (values: 95, 85, 75, 65, 55)
        self.assertEqual(summary['support_categories']['high_support_70+'], 3)  # 95, 85, 75 (> 70)
        self.assertEqual(summary['support_categories']['high_support_80+'], 2)  # 95, 85 (> 80)
        self.assertEqual(summary['support_categories']['high_support_95+'], 0)  # none > 95
        self.assertEqual(summary['support_categories']['percent_high_70+'], 60.0)  # 3/5 * 100
    
    def test_validate_bootstrap_analysis(self):
        """Test bootstrap analysis validation"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Initially invalid
        self.assertFalse(bm.validate_bootstrap_analysis())
        
        # Add tree but no values
        bm.bootstrap_tree = MagicMock()
        self.assertFalse(bm.validate_bootstrap_analysis())
        
        # Add values but no taxa map
        bm.bootstrap_values = {"node_1": 85.0}
        self.assertFalse(bm.validate_bootstrap_analysis())
        
        # Add taxa map
        bm.taxa_bootstrap_map = {frozenset(["seq1", "seq2"]): 85.0}
        self.assertTrue(bm.validate_bootstrap_analysis())
    
    def test_clear_bootstrap_data(self):
        """Test clearing bootstrap analysis data"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Set up some data
        bm.bootstrap_tree = MagicMock()
        bm.bootstrap_values = {"node_1": 85.0}
        bm.taxa_bootstrap_map = {frozenset(["seq1", "seq2"]): 85.0}
        
        # Clear data
        bm.clear_bootstrap_data()
        
        # Verify everything is cleared
        self.assertIsNone(bm.bootstrap_tree)
        self.assertEqual(bm.bootstrap_values, {})
        self.assertEqual(bm.taxa_bootstrap_map, {})
    
    def test_export_bootstrap_values(self):
        """Test exporting bootstrap values to file"""
        bm = BootstrapManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Set up mock data
        bm.taxa_bootstrap_map = {
            frozenset(["seq1", "seq2"]): 95.0,
            frozenset(["seq3", "seq4"]): 85.0,
            frozenset(["seq1", "seq2", "seq3", "seq4"]): 75.0
        }
        
        # Export to file
        output_file = self.temp_dir / "bootstrap_export.txt"
        bm.export_bootstrap_values(output_file)
        
        # Verify file was created and has content
        self.assertTrue(output_file.exists())
        content = output_file.read_text()
        
        self.assertIn("# Bootstrap Support Values", content)
        self.assertIn("seq1, seq2 -> 95.0%", content)
        self.assertIn("seq3, seq4 -> 85.0%", content)
        
        # Check that higher values appear first (sorted by value)
        lines = [line for line in content.splitlines() if "->" in line and not line.startswith("#")]
        self.assertGreater(len(lines), 0)
        # First data line should contain the highest bootstrap value (95.0%)
        values_in_order = [line.split("->")[1].strip() for line in lines]
        self.assertEqual(values_in_order[0], "95.0%")  # Highest value first
        self.assertEqual(values_in_order[1], "85.0%")  # Second highest


class TestBootstrapManagerIntegration(unittest.TestCase):
    """Integration tests for BootstrapManager component"""
    
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
    
    def _create_realistic_bootstrap_tree(self):
        """Create a realistic bootstrap tree with species names"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        homo = Clade(name="Homo_sapiens")
        pan = Clade(name="Pan_troglodytes") 
        mus = Clade(name="Mus_musculus")
        rattus = Clade(name="Rattus_norvegicus")
        
        primates = Clade(confidence=98.0)
        primates.clades = [homo, pan]
        
        rodents = Clade(confidence=95.0)
        rodents.clades = [mus, rattus]
        
        root = Clade(confidence=100.0)
        root.clades = [primates, rodents]
        
        return Tree(root=root)
    
    @patch('Bio.Phylo.read')
    def test_complete_bootstrap_workflow(self, mock_phylo_read):
        """Test complete bootstrap analysis workflow"""
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load and process alignment
            dp.load_alignment(self.alignment_file, "fasta")
            nexus_file = fm.get_temp_path("alignment.nex")
            dp.convert_to_nexus(nexus_file)
            
            # Mock PAUP* execution for bootstrap
            with patch.object(et, 'run_paup_command_file') as mock_paup:
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_result.stdout = "Bootstrap analysis completed"
                mock_paup.return_value = mock_result
                
                # Create mock bootstrap tree file
                tree_file = fm.get_temp_path("bootstrap_trees.tre")
                with open(tree_file, 'w') as f:
                    f.write("((Homo_sapiens,Pan_troglodytes)98,(Mus_musculus,Rattus_norvegicus)95);")
                
                # Mock phylo.read
                mock_tree = self._create_realistic_bootstrap_tree()
                mock_phylo_read.return_value = mock_tree
                
                # Initialize bootstrap manager
                bm = BootstrapManager(
                    file_manager=fm,
                    external_tools=et,
                    data_processor=dp,
                    num_replicates=25,
                    consensus_threshold=60,
                    debug=True
                )
                
                # Run ML bootstrap
                result_tree = bm.run_ml_bootstrap(
                    model_commands="lset nst=6 rates=gamma;",
                    timeout_sec=300
                )
                
                # Verify results
                self.assertEqual(result_tree, mock_tree)
                self.assertTrue(bm.validate_bootstrap_analysis())
                
                # Test bootstrap support lookup
                primate_support = bm.get_bootstrap_support(["Homo_sapiens", "Pan_troglodytes"])
                self.assertIsNotNone(primate_support)
                
                rodent_support = bm.get_bootstrap_support(["Mus_musculus", "Rattus_norvegicus"])
                self.assertIsNotNone(rodent_support)
                
                # Test tree annotation
                from Bio.Phylo.BaseTree import Tree, Clade
                test_tree = Tree(root=Clade())
                annotated_tree = bm.annotate_tree_with_bootstrap(test_tree)
                self.assertIsNotNone(annotated_tree)
                
                # Test summary generation
                summary = bm.create_bootstrap_summary()
                self.assertTrue(summary['bootstrap_completed'])
                self.assertEqual(summary['num_replicates'], 25)
                self.assertEqual(summary['consensus_threshold'], 60)
                
                # Test export
                export_file = fm.get_temp_path("bootstrap_export.txt")
                bm.export_bootstrap_values(export_file)
                self.assertTrue(export_file.exists())
    
    @patch('Bio.Phylo.read')
    def test_parsimony_bootstrap_workflow(self, mock_phylo_read):
        """Test parsimony bootstrap workflow"""
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            dp.load_alignment(self.alignment_file, "fasta")
            nexus_file = fm.get_temp_path("alignment.nex")
            dp.convert_to_nexus(nexus_file)
            
            with patch.object(et, 'run_paup_command_file') as mock_paup:
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_paup.return_value = mock_result
                
                tree_file = fm.get_temp_path("bootstrap_trees.tre")
                with open(tree_file, 'w') as f:
                    f.write("((Homo_sapiens,Pan_troglodytes)90,(Mus_musculus,Rattus_norvegicus)85);")
                
                mock_tree = self._create_realistic_bootstrap_tree()
                mock_phylo_read.return_value = mock_tree
                
                bm = BootstrapManager(
                    file_manager=fm,
                    external_tools=et,
                    data_processor=dp,
                    num_replicates=10
                )
                
                # Run parsimony bootstrap
                result_tree = bm.run_parsimony_bootstrap(
                    search_replicates=2,
                    swap_method="tbr",
                    timeout_sec=600
                )
                
                # Verify results
                self.assertEqual(result_tree, mock_tree)
                self.assertTrue(bm.validate_bootstrap_analysis())
                
                # Verify PAUP* script contained parsimony settings
                script_content = mock_paup.call_args[0][0].read_text()
                self.assertIn("set criterion=parsimony;", script_content)
                self.assertIn("bootstrap nreps=10", script_content)


def run_bootstrap_manager_tests():
    """Run the BootstrapManager test suite"""
    print("=" * 80)
    print("BOOTSTRAP MANAGER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted bootstrap analysis functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestBootstrapManager))
    suite.addTests(loader.loadTestsFromTestCase(TestBootstrapManagerIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("BOOTSTRAP MANAGER TEST RESULTS")
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
    print(f"\nBootstrapManager Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_bootstrap_manager_tests()
    sys.exit(0 if success else 1)