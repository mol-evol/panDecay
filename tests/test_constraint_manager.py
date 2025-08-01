#!/usr/bin/env python3
"""
Test Suite for ConstraintManager Component

Tests the extracted constraint management functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis import FileManager, ExternalTools, DataProcessor
from core.analysis.constraint_manager import (
    ConstraintManager, ConstraintError, ConstraintParsingError,
    ConstraintValidationError, ConstraintTestingError
)


class TestConstraintManager(unittest.TestCase):
    """Test ConstraintManager component functionality"""
    
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
>Canis_lupus
ATCGATCGATCGATCG
>Felis_catus
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
    
    def test_initialization_default(self):
        """Test ConstraintManager initialization with defaults"""
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        self.assertEqual(cm.constraint_mode, "all")
        self.assertIsNone(cm.test_branches)
        self.assertIsNone(cm.constraint_file)
        self.assertEqual(cm.config_constraints, {})
        self.assertEqual(cm.user_constraints, {})
        self.assertEqual(len(cm.available_taxa), 6)
    
    def test_initialization_invalid_mode(self):
        """Test initialization with invalid constraint mode"""
        with self.assertRaises(ConstraintError) as context:
            ConstraintManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                constraint_mode="invalid"
            )
        
        self.assertIn("Invalid constraint mode", str(context.exception))
    
    def test_parse_command_line_constraints(self):
        """Test parsing constraints from command line"""
        test_branches = "Homo_sapiens,Pan_troglodytes;Mus_musculus,Rattus_norvegicus"
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            test_branches=test_branches
        )
        
        self.assertEqual(len(cm.user_constraints), 2)
        self.assertIn("cmd_constraint_1", cm.user_constraints)
        self.assertIn("cmd_constraint_2", cm.user_constraints)
        
        constraint1 = cm.user_constraints["cmd_constraint_1"]
        self.assertEqual(constraint1, {"Homo_sapiens", "Pan_troglodytes"})
        
        constraint2 = cm.user_constraints["cmd_constraint_2"]
        self.assertEqual(constraint2, {"Mus_musculus", "Rattus_norvegicus"})
    
    def test_parse_constraint_file(self):
        """Test parsing constraints from file"""
        # Create constraint file
        constraint_file = self.temp_dir / "constraints.txt"
        with open(constraint_file, 'w') as f:
            f.write("# Primate clade\n")
            f.write("Homo_sapiens,Pan_troglodytes\n")
            f.write("\n")  # Empty line
            f.write("Mus_musculus,Rattus_norvegicus\n")
            f.write("# Comment line\n")
            f.write("Canis_lupus,Felis_catus\n")
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_file=constraint_file
        )
        
        self.assertEqual(len(cm.user_constraints), 3)
        self.assertIn("file_constraint_2", cm.user_constraints)  # Line 2
        self.assertIn("file_constraint_4", cm.user_constraints)  # Line 4
        self.assertIn("file_constraint_6", cm.user_constraints)  # Line 6
    
    def test_parse_config_constraints(self):
        """Test parsing constraints from config"""
        config_constraints = {
            "primates": "Homo_sapiens,Pan_troglodytes",
            "rodents": "Mus_musculus,Rattus_norvegicus",
            "carnivores": "Canis_lupus,Felis_catus"
        }
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            config_constraints=config_constraints
        )
        
        self.assertEqual(len(cm.user_constraints), 3)
        self.assertIn("primates", cm.user_constraints)
        self.assertIn("rodents", cm.user_constraints)
        self.assertIn("carnivores", cm.user_constraints)
        
        self.assertEqual(cm.user_constraints["primates"], {"Homo_sapiens", "Pan_troglodytes"})
    
    def test_constraint_validation_missing_taxa(self):
        """Test constraint validation with missing taxa"""
        # Create constraint with non-existent taxon
        test_branches = "Homo_sapiens,Non_existent_taxon"
        
        with self.assertLogs(level='WARNING') as log:
            cm = ConstraintManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                test_branches=test_branches
            )
        
        # Should have no valid constraints
        self.assertEqual(len(cm.user_constraints), 0)
        self.assertTrue(any("Skipping invalid constraint" in message for message in log.output))
    
    def test_constraint_validation_too_few_taxa(self):
        """Test constraint validation with too few taxa"""
        test_branches = "Homo_sapiens"  # Only one taxon
        
        with self.assertLogs(level='WARNING') as log:
            cm = ConstraintManager(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=self.data_processor,
                test_branches=test_branches
            )
        
        # Should warn about too few taxa
        self.assertTrue(any("fewer than 2 taxa" in message for message in log.output))
    
    def test_should_test_clade_all_mode(self):
        """Test clade testing in 'all' mode"""
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_mode="all"
        )
        
        # Should test all clades
        clade1 = {"Homo_sapiens", "Pan_troglodytes"}
        clade2 = {"Mus_musculus", "Rattus_norvegicus"}
        
        self.assertTrue(cm.should_test_clade(clade1))
        self.assertTrue(cm.should_test_clade(clade2))
    
    def test_should_test_clade_specific_mode(self):
        """Test clade testing in 'specific' mode"""
        test_branches = "Homo_sapiens,Pan_troglodytes"
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_mode="specific",
            test_branches=test_branches
        )
        
        # Should only test specified clades
        clade1 = {"Homo_sapiens", "Pan_troglodytes"}  # Specified
        clade2 = {"Mus_musculus", "Rattus_norvegicus"}  # Not specified
        
        self.assertTrue(cm.should_test_clade(clade1))
        self.assertFalse(cm.should_test_clade(clade2))
    
    def test_should_test_clade_exclude_mode(self):
        """Test clade testing in 'exclude' mode"""
        test_branches = "Homo_sapiens,Pan_troglodytes"
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_mode="exclude",
            test_branches=test_branches
        )
        
        # Should test all except specified clades
        clade1 = {"Homo_sapiens", "Pan_troglodytes"}  # Excluded
        clade2 = {"Mus_musculus", "Rattus_norvegicus"}  # Not excluded
        
        self.assertFalse(cm.should_test_clade(clade1))
        self.assertTrue(cm.should_test_clade(clade2))
    
    def test_get_testable_clades_from_tree(self):
        """Test extracting testable clades from tree"""
        from Bio.Phylo.BaseTree import Tree, Clade
        
        # Create test tree: ((Homo_sapiens,Pan_troglodytes),(Mus_musculus,Rattus_norvegicus))
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
        
        tree = Tree(root=root)
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_mode="all"
        )
        
        testable_clades = cm.get_testable_clades_from_tree(tree)
        
        # Should find 3 internal nodes (root + 2 clades)
        self.assertEqual(len(testable_clades), 3)
        
        # Check that we found the expected clades
        clade_taxa_sets = [taxa_set for _, taxa_set in testable_clades]
        
        expected_primates = {"Homo_sapiens", "Pan_troglodytes"}
        expected_rodents = {"Mus_musculus", "Rattus_norvegicus"}
        expected_root = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus", "Rattus_norvegicus"}
        
        self.assertIn(expected_primates, clade_taxa_sets)
        self.assertIn(expected_rodents, clade_taxa_sets)
        self.assertIn(expected_root, clade_taxa_sets)
    
    def test_filename_generation(self):
        """Test constraint filename generation"""
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        clade_id = "test_clade"
        
        tree_path = cm.get_constraint_tree_path(clade_id)
        score_path = cm.get_constraint_score_path(clade_id)
        search_path = cm.get_constraint_search_path(clade_id)
        log_path = cm.get_constraint_log_path(clade_id)
        
        self.assertTrue(str(tree_path).endswith("constraint_tree_test_clade.tre"))
        self.assertTrue(str(score_path).endswith("constraint_score_test_clade.txt"))
        self.assertTrue(str(search_path).endswith("constraint_search_test_clade.nex"))
        self.assertTrue(str(log_path).endswith("paup_constraint_test_clade.log"))
    
    def test_format_taxa_for_paup(self):
        """Test PAUP* taxa formatting"""
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        taxa_set = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus"}
        formatted = cm.format_taxa_for_paup(taxa_set)
        
        # Should be sorted and formatted for PAUP*
        expected = "((Homo_sapiens, Mus_musculus, Pan_troglodytes))"
        self.assertEqual(formatted, expected)
    
    def test_format_taxa_for_mrbayes(self):
        """Test MrBayes taxa formatting"""
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        taxa_set = {"Homo_sapiens", "Pan_troglodytes", "Mus_musculus"}
        formatted = cm.format_taxa_for_mrbayes(taxa_set)
        
        # Should be sorted and space-separated
        expected = "Homo_sapiens Mus_musculus Pan_troglodytes"
        self.assertEqual(formatted, expected)
    
    def test_create_constraint_summary(self):
        """Test constraint summary creation"""
        config_constraints = {
            "primates": "Homo_sapiens,Pan_troglodytes",
            "rodents": "Mus_musculus,Rattus_norvegicus"
        }
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            constraint_mode="specific",
            config_constraints=config_constraints
        )
        
        summary = cm.create_constraint_summary()
        
        self.assertEqual(summary['constraint_mode'], 'specific')
        self.assertEqual(summary['num_user_constraints'], 2)
        self.assertIn('config_file', summary['constraint_sources'])
        self.assertEqual(summary['available_taxa_count'], 6)
        
        # Check constraint details
        self.assertIn('primates', summary['constraints'])
        self.assertEqual(summary['constraints']['primates']['taxa_count'], 2)
        self.assertEqual(summary['constraints']['primates']['taxa'], ['Homo_sapiens', 'Pan_troglodytes'])
    
    def test_validate_all_constraints(self):
        """Test constraint validation"""
        # Valid constraints
        config_constraints = {
            "primates": "Homo_sapiens,Pan_troglodytes",
            "rodents": "Mus_musculus,Rattus_norvegicus"
        }
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            config_constraints=config_constraints
        )
        
        self.assertTrue(cm.validate_all_constraints())
        
        # Test with no constraints
        cm_empty = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        self.assertTrue(cm_empty.validate_all_constraints())
    
    def test_utility_methods(self):
        """Test utility methods"""
        config_constraints = {
            "primates": "Homo_sapiens,Pan_troglodytes"
        }
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            config_constraints=config_constraints
        )
        
        # Test getter methods
        self.assertEqual(cm.get_constraint_mode(), "all")
        self.assertTrue(cm.has_user_constraints())
        
        user_constraints = cm.get_user_constraints()
        self.assertEqual(len(user_constraints), 1)
        self.assertIn("primates", user_constraints)
        
        # Test clear
        cm.clear_constraints()
        self.assertFalse(cm.has_user_constraints())
        self.assertEqual(len(cm.get_user_constraints()), 0)
    
    def test_export_constraints(self):
        """Test constraint export to file"""
        config_constraints = {
            "primates": "Homo_sapiens,Pan_troglodytes",
            "rodents": "Mus_musculus,Rattus_norvegicus"
        }
        
        cm = ConstraintManager(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            config_constraints=config_constraints
        )
        
        # Export constraints
        export_file = self.temp_dir / "exported_constraints.txt"
        cm.export_constraints(export_file)
        
        # Verify file was created and has content
        self.assertTrue(export_file.exists())
        content = export_file.read_text()
        
        self.assertIn("# Phylogenetic Constraints", content)
        self.assertIn("# Mode: all", content)
        self.assertIn("# Total constraints: 2", content)
        self.assertIn("Homo_sapiens, Pan_troglodytes", content)
        self.assertIn("Mus_musculus, Rattus_norvegicus", content)


class TestConstraintManagerIntegration(unittest.TestCase):
    """Integration tests for ConstraintManager component"""
    
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
>Pongo_pygmaeus
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
    
    def test_complete_constraint_workflow(self):
        """Test complete constraint management workflow"""
        # Create constraint file
        constraint_file = self.temp_dir / "phylo_constraints.txt"
        with open(constraint_file, 'w') as f:
            f.write("# Great apes\n")
            f.write("Homo_sapiens,Pan_troglodytes,Gorilla_gorilla\n")
            f.write("# Hominidae\n")
            f.write("Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Pongo_pygmaeus\n")
        
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load alignment
            dp.load_alignment(self.alignment_file, "fasta")
            
            # Create constraint manager with multiple sources
            cm = ConstraintManager(
                file_manager=fm,
                external_tools=et,
                data_processor=dp,
                constraint_mode="specific",
                test_branches="Mus_musculus,Rattus_norvegicus",
                constraint_file=constraint_file,
                config_constraints={
                    "all_primates": "Homo_sapiens,Pan_troglodytes,Gorilla_gorilla,Pongo_pygmaeus"
                },
                debug=True
            )
            
            # Verify constraint parsing
            self.assertGreater(len(cm.user_constraints), 0)
            self.assertTrue(cm.validate_all_constraints())
            
            # Test clade testing logic
            great_apes = {"Homo_sapiens", "Pan_troglodytes", "Gorilla_gorilla"}
            rodents = {"Mus_musculus", "Rattus_norvegicus"}
            random_clade = {"Homo_sapiens", "Mus_musculus"}
            
            # Should test constraints that match user-defined ones
            self.assertTrue(cm.should_test_clade(rodents))  # From command line
            
            # Test filename generation
            tree_path = cm.get_constraint_tree_path("test_clade")
            self.assertTrue(str(tree_path).endswith(".tre"))
            
            # Test formatting
            paup_format = cm.format_taxa_for_paup(great_apes)
            self.assertIn("((", paup_format)
            self.assertIn("))", paup_format)
            
            mrbayes_format = cm.format_taxa_for_mrbayes(great_apes)
            self.assertIn("Homo_sapiens", mrbayes_format)
            self.assertIn(" ", mrbayes_format)  # Space-separated
            
            # Test summary
            summary = cm.create_constraint_summary()
            self.assertEqual(summary['constraint_mode'], 'specific')
            self.assertGreater(summary['num_user_constraints'], 0)
            
            # Test export
            export_file = fm.get_temp_path("exported_constraints.txt")
            cm.export_constraints(export_file)
            self.assertTrue(export_file.exists())
    
    def test_constraint_mode_workflows(self):
        """Test different constraint mode workflows"""
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            dp.load_alignment(self.alignment_file, "fasta")
            
            # Test 'all' mode
            cm_all = ConstraintManager(
                file_manager=fm,
                external_tools=et,
                data_processor=dp,
                constraint_mode="all"
            )
            
            # Should test any clade
            test_clade = {"Homo_sapiens", "Pan_troglodytes"}
            self.assertTrue(cm_all.should_test_clade(test_clade))
            
            # Test 'specific' mode with constraints
            cm_specific = ConstraintManager(
                file_manager=fm,
                external_tools=et,
                data_processor=dp,
                constraint_mode="specific",
                test_branches="Homo_sapiens,Pan_troglodytes"
            )
            
            # Should only test specified clades
            specified_clade = {"Homo_sapiens", "Pan_troglodytes"}
            other_clade = {"Mus_musculus", "Rattus_norvegicus"}
            
            self.assertTrue(cm_specific.should_test_clade(specified_clade))
            self.assertFalse(cm_specific.should_test_clade(other_clade))
            
            # Test 'exclude' mode
            cm_exclude = ConstraintManager(
                file_manager=fm,
                external_tools=et,
                data_processor=dp,
                constraint_mode="exclude",
                test_branches="Homo_sapiens,Pan_troglodytes"
            )
            
            # Should test all except specified clades
            self.assertFalse(cm_exclude.should_test_clade(specified_clade))
            self.assertTrue(cm_exclude.should_test_clade(other_clade))


def run_constraint_manager_tests():
    """Run the ConstraintManager test suite"""
    print("=" * 80)
    print("CONSTRAINT MANAGER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted constraint management functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestConstraintManager))
    suite.addTests(loader.loadTestsFromTestCase(TestConstraintManagerIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("CONSTRAINT MANAGER TEST RESULTS")
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
    print(f"\nConstraintManager Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_constraint_manager_tests()
    sys.exit(0 if success else 1)