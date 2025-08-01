#!/usr/bin/env python3
"""
Component Integration Test Suite

Tests that the extracted components (FileManager, ExternalTools, DataProcessor)
work together correctly and maintain the same interface as the original.
"""

import unittest
import tempfile
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis import FileManager, ExternalTools, DataProcessor


class TestComponentIntegration(unittest.TestCase):
    """Test integration between extracted components"""
    
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
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_filemanager_dataprocessor_integration(self):
        """Test FileManager and DataProcessor working together"""
        
        # Initialize components
        with FileManager(debug=True) as fm:
            dp = DataProcessor(data_type="dna")
            
            # Copy alignment to temporary directory
            temp_alignment = fm.copy_to_temp(self.alignment_file, "analysis.fasta")
            self.assertTrue(temp_alignment.exists())
            
            # Load alignment using DataProcessor
            alignment = dp.load_alignment(temp_alignment, "fasta")
            self.assertIsNotNone(alignment)
            
            # Convert to NEXUS in temporary directory
            nexus_file = fm.get_temp_path("analysis.nex")
            dp.convert_to_nexus(nexus_file)
            
            # Verify NEXUS file exists and has correct content
            self.assertTrue(nexus_file.exists())
            content = nexus_file.read_text()
            self.assertIn("#NEXUS", content)
            self.assertIn("NTAX=4", content)
            
            # Get alignment stats
            stats = dp.get_alignment_stats()
            self.assertEqual(stats["num_sequences"], 4)
            self.assertEqual(stats["data_type"], "dna")
    
    @patch('subprocess.Popen')
    def test_filemanager_externaltools_integration(self, mock_popen):
        """Test FileManager and ExternalTools working together"""
        
        # Mock successful PAUP* execution
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("PAUP output", "")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        
        # Initialize components
        with FileManager(debug=True) as fm:
            et = ExternalTools(debug=True)
            
            # Create command file in temporary directory
            cmd_content = """#NEXUS
begin paup;
    execute analysis.nex;
    set criterion=likelihood;
    hsearch;
    savetrees file=tree.tre;
    quit;
end;
"""
            cmd_file = fm.write_file("commands.nex", cmd_content)
            log_file = fm.get_temp_path("analysis.log")
            
            # Run PAUP* command using ExternalTools
            result = et.run_paup_command_file(
                cmd_file,
                log_file,
                fm.get_temp_path()
            )
            
            # Verify execution
            self.assertEqual(result.returncode, 0)
            self.assertEqual(result.stdout, "PAUP output")
            self.assertTrue(log_file.exists())
    
    def test_three_component_workflow(self):
        """Test complete workflow using all three components"""
        
        # Mock external tool execution
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = ("Analysis complete", "")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            # Initialize all components
            with FileManager(debug=True) as fm:
                dp = DataProcessor(data_type="dna", debug=True)
                et = ExternalTools(debug=True)
                
                # Step 1: Copy and process alignment data
                temp_alignment = fm.copy_to_temp(self.alignment_file)
                alignment = dp.load_alignment(temp_alignment, "fasta")
                
                # Step 2: Convert to NEXUS format
                nexus_file = fm.get_temp_path("analysis.nex")
                dp.convert_to_nexus(nexus_file)
                
                # Step 3: Create PAUP* command file
                taxon_names = dp.get_taxon_names(formatted=True)
                cmd_content = f"""#NEXUS
begin paup;
    execute {nexus_file.name};
    set criterion=likelihood;
    hsearch;
    quit;
end;
"""
                cmd_file = fm.write_file("analysis_commands.nex", cmd_content)
                
                # Step 4: Run analysis
                log_file = fm.get_temp_path("analysis.log")
                result = et.run_paup_command_file(
                    cmd_file,
                    log_file,
                    fm.get_temp_path()
                )
                
                # Step 5: Verify results
                self.assertEqual(result.returncode, 0)
                self.assertTrue(nexus_file.exists())
                self.assertTrue(cmd_file.exists())
                self.assertTrue(log_file.exists())
                
                # Check file contents
                nexus_content = nexus_file.read_text()
                self.assertIn("seq1", nexus_content)
                self.assertIn("DATATYPE=DNA", nexus_content)
                
                cmd_content_read = cmd_file.read_text()
                self.assertIn("execute analysis.nex", cmd_content_read)
                
                # Verify alignment stats
                stats = dp.get_alignment_stats()
                self.assertEqual(stats["num_sequences"], 4)
                self.assertEqual(stats["alignment_length"], 16)
    
    def test_component_error_handling(self):
        """Test error handling across components"""
        
        with FileManager() as fm:
            dp = DataProcessor()
            et = ExternalTools()
            
            # Test DataProcessor error propagation
            with self.assertRaises(Exception):
                dp.load_alignment(fm.get_temp_path("nonexistent.fasta"))
            
            # Test ExternalTools error propagation
            with self.assertRaises(Exception):
                et.run_paup_command_file(
                    fm.get_temp_path("nonexistent.nex"),
                    fm.get_temp_path("log.txt"),
                    fm.get_temp_path()
                )
    
    def test_component_memory_cleanup(self):
        """Test that components properly clean up resources"""
        
        # Create components in separate scope
        temp_paths = []
        
        # Test FileManager cleanup
        with FileManager() as fm:
            temp_path = fm.get_temp_path()
            temp_paths.append(temp_path)
            
            # Create files
            test_file = fm.write_file("test.txt", "content")
            self.assertTrue(test_file.exists())
        
        # After context exit, temp directory should be cleaned up
        # (unless debug mode is enabled)
        
        # Test DataProcessor memory cleanup  
        dp = DataProcessor()
        dp.load_alignment(self.alignment_file, "fasta")
        
        # Alignment should be loaded
        self.assertTrue(dp.is_alignment_loaded())
        
        # After going out of scope, memory should be managed by Python GC
        del dp


class TestBackwardCompatibility(unittest.TestCase):
    """Test that extracted components maintain backward compatibility"""
    
    def setUp(self):
        """Set up compatibility test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create test alignment
        self.test_alignment = """>Homo_sapiens
ATCGATCGATCGATCG
>Pan_troglodytes
ATCGATCGATCGATCG
"""
        self.alignment_file = self.temp_dir / "test.fasta"
        with open(self.alignment_file, 'w') as f:
            f.write(self.test_alignment)
    
    def tearDown(self):
        """Clean up compatibility test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_original_interface_preserved(self):
        """Test that original method signatures are preserved"""
        
        # FileManager interface
        fm = FileManager()
        self.assertTrue(hasattr(fm, 'get_temp_path'))
        self.assertTrue(hasattr(fm, 'write_file'))
        self.assertTrue(hasattr(fm, 'copy_to_temp'))
        self.assertTrue(hasattr(fm, 'cleanup_all'))
        
        # DataProcessor interface
        dp = DataProcessor()
        self.assertTrue(hasattr(dp, 'load_alignment'))
        self.assertTrue(hasattr(dp, 'convert_to_nexus'))
        self.assertTrue(hasattr(dp, 'format_taxon_for_paup'))
        self.assertTrue(hasattr(dp, 'validate_discrete_data'))
        
        # ExternalTools interface
        et = ExternalTools()
        self.assertTrue(hasattr(et, 'run_paup_command_file'))
        self.assertTrue(hasattr(et, 'run_mrbayes'))
        self.assertTrue(hasattr(et, 'validate_tool'))
    
    def test_parameter_compatibility(self):
        """Test that method parameters match original expectations"""
        
        # Test DataProcessor taxon formatting (original behavior)
        dp = DataProcessor()
        
        test_cases = [
            ("Normal_Name", "Normal_Name"),
            ("Name with spaces", "'Name with spaces'"),
            ("Name'with'quotes", "'Name_with_quotes'")
        ]
        
        for input_name, expected in test_cases:
            result = dp.format_taxon_for_paup(input_name)
            self.assertEqual(result, expected)


def run_integration_tests():
    """Run the component integration test suite"""
    print("=" * 80)
    print("COMPONENT INTEGRATION TEST SUITE")
    print("=" * 80)
    print("Testing integration between extracted components...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestComponentIntegration))
    suite.addTests(loader.loadTestsFromTestCase(TestBackwardCompatibility))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\\n" + "=" * 80)
    print("INTEGRATION TEST RESULTS")
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
    print(f"\\nIntegration Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_integration_tests()
    sys.exit(0 if success else 1)