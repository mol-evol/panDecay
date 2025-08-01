#!/usr/bin/env python3
"""
Test Suite for DataProcessor Component

Tests the extracted data processing functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
import os
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis.data_processor import (
    DataProcessor, DataProcessingError, AlignmentLoadError, AlignmentValidationError
)


class TestDataProcessor(unittest.TestCase):
    """Test DataProcessor component functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create test DNA alignment
        self.dna_content = """>seq1
ATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCG
>seq4
ATCGATCGATCGATCG
"""
        self.dna_file = self.temp_dir / "test_dna.fasta"
        with open(self.dna_file, 'w') as f:
            f.write(self.dna_content)
        
        # Create test protein alignment
        self.protein_content = """>seq1
ACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWY
>seq3
ACDEFGHIKLMNPQRSTVWY
>seq4
ACDEFGHIKLMNPQRSTVWY
"""
        self.protein_file = self.temp_dir / "test_protein.fasta"
        with open(self.protein_file, 'w') as f:
            f.write(self.protein_content)
        
        # Create test discrete alignment
        self.discrete_content = """>seq1
01010101
>seq2
10101010
>seq3
01010101
>seq4
10101010
"""
        self.discrete_file = self.temp_dir / "test_discrete.fasta"
        with open(self.discrete_file, 'w') as f:
            f.write(self.discrete_content)
        
        # Create invalid discrete alignment
        self.invalid_discrete_content = """>seq1
ATCGATCG
>seq2
ATCGATCG
"""
        self.invalid_discrete_file = self.temp_dir / "invalid_discrete.fasta"
        with open(self.invalid_discrete_file, 'w') as f:
            f.write(self.invalid_discrete_content)
        
        # Create PHYLIP format alignment
        self.phylip_content = """4 16
seq1        ATCGATCGATCGATCG
seq2        ATCGATCGATCGATCG  
seq3        ATCGATCGATCGATCG
seq4        ATCGATCGATCGATCG
"""
        self.phylip_file = self.temp_dir / "test.phy"
        with open(self.phylip_file, 'w') as f:
            f.write(self.phylip_content)
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test DataProcessor initialization"""
        # Test default initialization
        dp = DataProcessor()
        self.assertEqual(dp.data_type, "dna")
        self.assertFalse(dp.debug)
        self.assertIsNone(dp.alignment)
        self.assertEqual(dp.alignment_stats, {})
        
        # Test custom initialization
        dp_custom = DataProcessor(data_type="protein", debug=True)
        self.assertEqual(dp_custom.data_type, "protein")
        self.assertTrue(dp_custom.debug)
        
        # Test invalid data type defaults to DNA
        dp_invalid = DataProcessor(data_type="invalid")
        self.assertEqual(dp_invalid.data_type, "dna")
    
    def test_load_alignment_dna(self):
        """Test loading DNA alignment"""
        dp = DataProcessor(data_type="dna")
        alignment = dp.load_alignment(self.dna_file, "fasta")
        
        self.assertIsNotNone(alignment)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.get_alignment_length(), 16)
        
        # Check alignment stats
        stats = dp.get_alignment_stats()
        self.assertEqual(stats["num_sequences"], 4)
        self.assertEqual(stats["alignment_length"], 16)
        self.assertEqual(stats["format"], "fasta")
        self.assertEqual(stats["data_type"], "dna")
    
    def test_load_alignment_protein(self):
        """Test loading protein alignment"""
        dp = DataProcessor(data_type="protein")
        alignment = dp.load_alignment(self.protein_file, "fasta")
        
        self.assertIsNotNone(alignment)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.get_alignment_length(), 20)
    
    def test_load_alignment_discrete(self):
        """Test loading discrete alignment"""
        dp = DataProcessor(data_type="discrete")
        alignment = dp.load_alignment(self.discrete_file, "fasta")
        
        self.assertIsNotNone(alignment)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.get_alignment_length(), 8)
    
    def test_load_alignment_phylip(self):
        """Test loading PHYLIP format alignment"""
        dp = DataProcessor(data_type="dna")
        alignment = dp.load_alignment(self.phylip_file, "phylip")
        
        self.assertIsNotNone(alignment)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(alignment.get_alignment_length(), 16)
    
    def test_load_alignment_missing_file(self):
        """Test loading non-existent alignment file"""
        dp = DataProcessor()
        missing_file = self.temp_dir / "missing.fasta"
        
        with self.assertRaises(AlignmentLoadError) as context:
            dp.load_alignment(missing_file)
        
        self.assertIn("not found", str(context.exception))
    
    def test_validate_discrete_data_valid(self):
        """Test validation of valid discrete data"""
        dp = DataProcessor(data_type="discrete")
        dp.load_alignment(self.discrete_file, "fasta")
        
        self.assertTrue(dp.validate_discrete_data())
    
    def test_validate_discrete_data_invalid(self):
        """Test validation of invalid discrete data"""
        dp = DataProcessor(data_type="discrete")
        dp.load_alignment(self.invalid_discrete_file, "fasta")
        
        self.assertFalse(dp.validate_discrete_data())
    
    def test_format_taxon_for_paup(self):
        """Test taxon name formatting for PAUP*"""
        dp = DataProcessor()
        
        test_cases = [
            ("Normal_Name", "Normal_Name"),
            ("Name with spaces", "'Name with spaces'"),
            ("Name-with-dashes", "Name-with-dashes"),
            ("Name(with)parens", "'Name(with)parens'"),
            ("Name,with,commas", "'Name,with,commas'"),
            ("Name'with'quotes", "'Name_with_quotes'"),  # Single quotes replaced
            ("Name:with:colons", "'Name:with:colons'")
        ]
        
        for input_name, expected_output in test_cases:
            formatted = dp.format_taxon_for_paup(input_name)
            self.assertEqual(formatted, expected_output)
    
    def test_format_taxon_non_string(self):
        """Test taxon formatting with non-string input"""
        dp = DataProcessor()
        
        formatted = dp.format_taxon_for_paup(123)
        self.assertEqual(formatted, "123")
    
    def test_convert_to_nexus_dna(self):
        """Test NEXUS conversion for DNA data"""
        dp = DataProcessor(data_type="dna")
        dp.load_alignment(self.dna_file, "fasta")
        
        nexus_file = self.temp_dir / "test_dna.nex"
        result_path = dp.convert_to_nexus(nexus_file)
        
        self.assertEqual(result_path, nexus_file)
        self.assertTrue(nexus_file.exists())
        self.assertGreater(nexus_file.stat().st_size, 0)
        
        # Check content
        content = nexus_file.read_text()
        self.assertIn("#NEXUS", content)
        self.assertIn("BEGIN DATA", content)
        self.assertIn("DATATYPE=DNA", content)
        self.assertIn("NTAX=4", content)
        self.assertIn("NCHAR=16", content)
        self.assertIn("seq1", content)
    
    def test_convert_to_nexus_protein(self):
        """Test NEXUS conversion for protein data"""
        dp = DataProcessor(data_type="protein")
        dp.load_alignment(self.protein_file, "fasta")
        
        nexus_file = self.temp_dir / "test_protein.nex"
        dp.convert_to_nexus(nexus_file)
        
        content = nexus_file.read_text()
        self.assertIn("DATATYPE=PROTEIN", content)
    
    def test_convert_to_nexus_discrete(self):
        """Test NEXUS conversion for discrete data"""
        dp = DataProcessor(data_type="discrete")
        dp.load_alignment(self.discrete_file, "fasta")
        
        nexus_file = self.temp_dir / "test_discrete.nex"
        dp.convert_to_nexus(nexus_file)
        
        content = nexus_file.read_text()
        self.assertIn("DATATYPE=STANDARD", content)
        self.assertIn('SYMBOLS="01"', content)
        self.assertIn("BEGIN ASSUMPTIONS", content)
    
    def test_convert_to_nexus_no_assumptions(self):
        """Test NEXUS conversion without assumptions block"""
        dp = DataProcessor(data_type="discrete")
        dp.load_alignment(self.discrete_file, "fasta")
        
        nexus_file = self.temp_dir / "test_discrete_no_assumptions.nex"
        dp.convert_to_nexus(nexus_file, include_assumptions=False)
        
        content = nexus_file.read_text()
        self.assertNotIn("BEGIN ASSUMPTIONS", content)
    
    def test_convert_to_nexus_no_alignment(self):
        """Test NEXUS conversion without loaded alignment"""
        dp = DataProcessor()
        
        nexus_file = self.temp_dir / "no_alignment.nex"
        
        with self.assertRaises(DataProcessingError) as context:
            dp.convert_to_nexus(nexus_file)
        
        self.assertIn("No alignment loaded", str(context.exception))
    
    def test_get_taxon_names(self):
        """Test getting taxon names"""
        dp = DataProcessor()
        dp.load_alignment(self.dna_file, "fasta")
        
        # Test raw names
        names = dp.get_taxon_names()
        expected = ["seq1", "seq2", "seq3", "seq4"]
        self.assertEqual(names, expected)
        
        # Test formatted names
        formatted_names = dp.get_taxon_names(formatted=True)
        self.assertEqual(formatted_names, expected)  # These names don't need formatting
    
    def test_get_taxon_names_with_spaces(self):
        """Test getting formatted taxon names with spaces"""
        # Create alignment with spaces in names (using | to avoid BioPython parsing issues)
        content_with_spaces = """>seq|1
ATCG
>seq|2  
ATCG
"""
        file_with_spaces = self.temp_dir / "spaces.fasta"
        with open(file_with_spaces, 'w') as f:
            f.write(content_with_spaces)
        
        dp = DataProcessor()
        dp.load_alignment(file_with_spaces, "fasta")
        
        # Test raw names
        names = dp.get_taxon_names()
        self.assertEqual(names, ["seq|1", "seq|2"])
        
        # Test formatted names (pipes don't need quoting)
        formatted_names = dp.get_taxon_names(formatted=True)
        self.assertEqual(formatted_names, ["seq|1", "seq|2"])
    
    def test_get_sequence_by_name(self):
        """Test getting specific sequence by name"""
        dp = DataProcessor()
        dp.load_alignment(self.dna_file, "fasta")
        
        # Test existing sequence
        seq = dp.get_sequence_by_name("seq1")
        self.assertIsNotNone(seq)
        self.assertEqual(seq.id, "seq1")
        
        # Test non-existing sequence
        seq_none = dp.get_sequence_by_name("nonexistent")
        self.assertIsNone(seq_none)
        
        # Test without loaded alignment
        dp_empty = DataProcessor()
        seq_empty = dp_empty.get_sequence_by_name("seq1")
        self.assertIsNone(seq_empty)
    
    def test_filter_sequences(self):
        """Test filtering sequences by taxon names"""
        dp = DataProcessor()
        dp.load_alignment(self.dna_file, "fasta")
        
        # Test filtering subset
        filtered = dp.filter_sequences(["seq1", "seq3"])
        self.assertEqual(len(filtered), 2)
        self.assertEqual([record.id for record in filtered], ["seq1", "seq3"])
        
        # Test filtering with non-existent taxa
        with self.assertRaises(DataProcessingError) as context:
            dp.filter_sequences(["nonexistent"])
        
        self.assertIn("No sequences found", str(context.exception))
        
        # Test without loaded alignment
        dp_empty = DataProcessor()
        with self.assertRaises(DataProcessingError):
            dp_empty.filter_sequences(["seq1"])
    
    def test_get_alignment_stats_dna(self):
        """Test getting alignment statistics for DNA"""
        dp = DataProcessor(data_type="dna")
        dp.load_alignment(self.dna_file, "fasta")
        
        stats = dp.get_alignment_stats()
        
        self.assertEqual(stats["num_sequences"], 4)
        self.assertEqual(stats["alignment_length"], 16)
        self.assertEqual(stats["data_type"], "dna")
        self.assertIn("composition", stats)
        
        # Check DNA composition
        composition = stats["composition"]
        self.assertIn("A", composition)
        self.assertIn("T", composition)
        self.assertIn("C", composition)
        self.assertIn("G", composition)
    
    def test_get_alignment_stats_protein(self):
        """Test getting alignment statistics for protein"""
        dp = DataProcessor(data_type="protein")
        dp.load_alignment(self.protein_file, "fasta")
        
        stats = dp.get_alignment_stats()
        composition = stats["composition"]
        
        self.assertIn("hydrophobic", composition)
        self.assertIn("polar", composition)
        self.assertIn("charged", composition)
    
    def test_get_alignment_stats_discrete(self):
        """Test getting alignment statistics for discrete"""
        dp = DataProcessor(data_type="discrete")
        dp.load_alignment(self.discrete_file, "fasta")
        
        stats = dp.get_alignment_stats()
        composition = stats["composition"]
        
        self.assertIn("0", composition)
        self.assertIn("1", composition)
        self.assertEqual(composition["0"], 50.0)  # Half 0s, half 1s
        self.assertEqual(composition["1"], 50.0)
    
    def test_validate_clade_taxa(self):
        """Test validating clade taxa against alignment"""
        dp = DataProcessor()
        dp.load_alignment(self.dna_file, "fasta")
        
        # Test valid taxa
        valid_taxa = dp.validate_clade_taxa(["seq1", "seq2"])
        self.assertEqual(valid_taxa, ["seq1", "seq2"])
        
        # Test mixed valid/invalid taxa
        mixed_taxa = dp.validate_clade_taxa(["seq1", "nonexistent", "seq3"])
        self.assertEqual(mixed_taxa, ["seq1", "seq3"])
        
        # Test all invalid taxa
        with self.assertRaises(DataProcessingError) as context:
            dp.validate_clade_taxa(["nonexistent1", "nonexistent2"])
        
        self.assertIn("No valid taxa found", str(context.exception))
        
        # Test without loaded alignment
        dp_empty = DataProcessor()
        with self.assertRaises(DataProcessingError):
            dp_empty.validate_clade_taxa(["seq1"])
    
    def test_is_alignment_loaded(self):
        """Test checking if alignment is loaded"""
        dp = DataProcessor()
        
        # Initially no alignment
        self.assertFalse(dp.is_alignment_loaded())
        
        # After loading
        dp.load_alignment(self.dna_file, "fasta")
        self.assertTrue(dp.is_alignment_loaded())
    
    def test_data_validation_dna(self):
        """Test DNA data validation"""
        dp = DataProcessor(data_type="dna", debug=True)
        dp.load_alignment(self.dna_file, "fasta")
        
        # Should validate successfully (no exceptions)
        self.assertTrue(dp._validate_dna_data())
    
    def test_data_validation_protein(self):
        """Test protein data validation"""
        dp = DataProcessor(data_type="protein", debug=True)
        dp.load_alignment(self.protein_file, "fasta")
        
        # Should validate successfully
        self.assertTrue(dp._validate_protein_data())
    
    def test_edge_cases(self):
        """Test various edge cases"""
        dp = DataProcessor()
        
        # Test operations without loaded alignment
        self.assertEqual(dp.get_taxon_names(), [])
        stats = dp.get_alignment_stats()
        self.assertIn("data_type", stats)
        self.assertEqual(stats["data_type"], "dna")
        
        # Test empty clade validation
        dp.load_alignment(self.dna_file, "fasta")
        with self.assertRaises(DataProcessingError):
            dp.validate_clade_taxa([])


class TestDataProcessorIntegration(unittest.TestCase):
    """Integration tests for DataProcessor component"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create realistic alignment with various characters
        self.realistic_dna = """>Homo_sapiens
ATCGATCGATCGNNNN----ATCGATCG
>Pan_troglodytes  
ATCGATCGATCGNNNN----ATCGATCG
>Mus_musculus
ATCGATCGATCGNNNN----ATCGATCG
>Rattus_norvegicus
ATCGATCGATCGNNNN----ATCGATCG
"""
        self.realistic_file = self.temp_dir / "realistic.fasta"
        with open(self.realistic_file, 'w') as f:
            f.write(self.realistic_dna)
    
    def tearDown(self):
        """Clean up integration test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_complete_workflow(self):
        """Test complete data processing workflow"""
        dp = DataProcessor(data_type="dna", debug=True)
        
        # Load alignment
        alignment = dp.load_alignment(self.realistic_file, "fasta")
        self.assertIsNotNone(alignment)
        
        # Get taxon names (with underscores instead of spaces)
        taxa = dp.get_taxon_names()
        self.assertIn("Mus_musculus", taxa)  # Has underscore
        
        formatted_taxa = dp.get_taxon_names(formatted=True)
        self.assertIn("Mus_musculus", formatted_taxa)  # Underscores don't need quoting
        
        # Convert to NEXUS
        nexus_file = self.temp_dir / "realistic.nex"
        dp.convert_to_nexus(nexus_file)
        
        # Verify NEXUS content
        content = nexus_file.read_text()
        self.assertIn("Mus_musculus", content)  # Underscores preserved
        self.assertIn("NTAX=4", content)
        
        # Get stats
        stats = dp.get_alignment_stats()
        self.assertEqual(stats["num_sequences"], 4)
        self.assertIn("composition", stats)
        
        # Validate clade
        valid_clade = dp.validate_clade_taxa(["Homo_sapiens", "Mus_musculus"])
        self.assertEqual(len(valid_clade), 2)


def run_data_processor_tests():
    """Run the DataProcessor test suite"""
    print("=" * 80)
    print("DATA PROCESSOR COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted data processing functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestDataProcessor))
    suite.addTests(loader.loadTestsFromTestCase(TestDataProcessorIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\\n" + "=" * 80)
    print("DATA PROCESSOR TEST RESULTS")
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
    print(f"\\nDataProcessor Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_data_processor_tests()
    sys.exit(0 if success else 1)