#!/usr/bin/env python3
"""
Test Suite for FileManager Component

Tests the extracted file management functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
import os
import shutil
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis.file_manager import FileManager


class TestFileManager(unittest.TestCase):
    """Test FileManager component functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.test_content = "Test file content"
        
    def tearDown(self):
        """Clean up any remaining test directories"""
        # Clean up any debug_runs directories created during testing
        debug_runs_path = Path.cwd() / "debug_runs"
        if debug_runs_path.exists():
            try:
                shutil.rmtree(debug_runs_path)
            except OSError:
                # If directory is not empty, try to clean it up more thoroughly
                for root, dirs, files in os.walk(debug_runs_path, topdown=False):
                    for file in files:
                        try:
                            os.unlink(os.path.join(root, file))
                        except OSError:
                            pass
                    for dir in dirs:
                        try:
                            os.rmdir(os.path.join(root, dir))
                        except OSError:
                            pass
                try:
                    os.rmdir(debug_runs_path)
                except OSError:
                    pass  # If it still can't be removed, leave it
    
    def test_auto_cleanup_mode(self):
        """Test automatic cleanup mode (default)"""
        with FileManager() as fm:
            temp_path = fm.get_temp_path()
            self.assertTrue(temp_path.exists())
            self.assertTrue(temp_path.is_dir())
            
            # Write a test file
            test_file = fm.write_file("test.txt", self.test_content)
            self.assertTrue(test_file.exists())
        
        # After context exit, directory should be cleaned up
        self.assertFalse(temp_path.exists())
    
    def test_debug_mode(self):
        """Test debug mode (keeps files)"""
        fm = FileManager(debug=True)
        temp_path = fm.get_temp_path()
        
        self.assertTrue(temp_path.exists())
        self.assertTrue("debug_runs" in str(temp_path))
        
        # Write test file
        test_file = fm.write_file("debug_test.txt", self.test_content)
        self.assertTrue(test_file.exists())
        
        # Manual cleanup
        fm.cleanup_all()
        # In debug mode, files should still exist
        self.assertTrue(temp_path.exists())
    
    def test_keep_files_mode(self):
        """Test keep_files mode"""
        fm = FileManager(keep_files=True)
        temp_path = fm.get_temp_path()
        
        test_file = fm.write_file("keep_test.txt", self.test_content)
        self.assertTrue(test_file.exists())
        
        fm.cleanup_all()
        # Files should be kept
        self.assertTrue(temp_path.exists())
        self.assertTrue(test_file.exists())
    
    def test_custom_temp_dir(self):
        """Test custom temporary directory"""
        custom_dir = Path(tempfile.mkdtemp()) / "custom_analysis"
        
        try:
            fm = FileManager(temp_dir=custom_dir)
            self.assertEqual(fm.get_temp_path(), custom_dir)
            self.assertTrue(custom_dir.exists())
            
            test_file = fm.write_file("custom_test.txt", self.test_content)
            self.assertTrue(test_file.exists())
            
        finally:
            if custom_dir.exists():
                shutil.rmtree(custom_dir.parent)
    
    def test_file_operations(self):
        """Test basic file operations"""
        with FileManager() as fm:
            # Test write_file
            test_file = fm.write_file("test.txt", self.test_content)
            self.assertTrue(fm.file_exists("test.txt"))
            self.assertEqual(fm.get_file_size("test.txt"), len(self.test_content))
            
            # Test file content
            with open(test_file, 'r') as f:
                content = f.read()
            self.assertEqual(content, self.test_content)
    
    def test_copy_to_temp(self):
        """Test copying files to temporary directory"""
        # Create source file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(self.test_content)
            source_path = f.name
        
        try:
            with FileManager() as fm:
                # Test copy with default name
                copied_file = fm.copy_to_temp(source_path)
                self.assertTrue(copied_file.exists())
                self.assertEqual(copied_file.name, Path(source_path).name)
                
                # Test copy with custom name
                custom_copied = fm.copy_to_temp(source_path, "custom_name.txt")
                self.assertTrue(fm.file_exists("custom_name.txt"))
                
                # Verify content
                with open(custom_copied, 'r') as f:
                    content = f.read()
                self.assertEqual(content, self.test_content)
        
        finally:
            os.unlink(source_path)
    
    def test_copy_nonexistent_file(self):
        """Test copying non-existent file raises error"""
        with FileManager() as fm:
            with self.assertRaises(FileNotFoundError):
                fm.copy_to_temp("/nonexistent/file.txt")
    
    def test_directory_operations(self):
        """Test directory creation and management"""
        with FileManager() as fm:
            # Test ensure_directory_exists
            subdir = fm.ensure_directory_exists("subdir")
            self.assertTrue(subdir.exists())
            self.assertTrue(subdir.is_dir())
            
            # Test nested directory
            nested_dir = fm.ensure_directory_exists("level1/level2/level3")
            self.assertTrue(nested_dir.exists())
    
    def test_file_listing(self):
        """Test file listing functionality"""
        with FileManager() as fm:
            # Create multiple files
            fm.write_file("test1.txt", "content1")
            fm.write_file("test2.txt", "content2")
            fm.write_file("data.csv", "csv content")
            
            # Test list all files
            all_files = fm.list_files()
            self.assertEqual(len(all_files), 3)
            
            # Test pattern matching
            txt_files = fm.list_files("*.txt")
            self.assertEqual(len(txt_files), 2)
            
            csv_files = fm.list_files("*.csv")
            self.assertEqual(len(csv_files), 1)
    
    def test_cleanup_marked_files(self):
        """Test cleanup of marked files"""
        with FileManager() as fm:
            # Create test files
            test_file1 = fm.write_file("keep.txt", "keep this")
            test_file2 = fm.write_file("cleanup.txt", "cleanup this")
            
            # Mark one for cleanup
            fm.mark_for_cleanup(test_file2)
            
            # Run cleanup
            fm.cleanup_intermediate_files()
            
            # Check results
            self.assertTrue(test_file1.exists())  # Should be kept
            self.assertFalse(test_file2.exists())  # Should be cleaned up
    
    def test_cleanup_patterns(self):
        """Test cleanup with glob patterns"""
        with FileManager() as fm:
            # Create test files
            fm.write_file("keep.txt", "keep")
            fm.write_file("temp1.tmp", "temp1")
            fm.write_file("temp2.tmp", "temp2")
            fm.write_file("log.log", "log")
            
            # Cleanup with patterns
            fm.cleanup_intermediate_files(patterns=["*.tmp", "*.log"])
            
            # Check results
            self.assertTrue(fm.file_exists("keep.txt"))
            self.assertFalse(fm.file_exists("temp1.tmp"))
            self.assertFalse(fm.file_exists("temp2.tmp"))
            self.assertFalse(fm.file_exists("log.log"))
    
    def test_cleanup_in_debug_mode(self):
        """Test that cleanup is skipped in debug mode"""
        fm = FileManager(debug=True)
        
        try:
            test_file = fm.write_file("debug_file.txt", "content")
            fm.mark_for_cleanup(test_file)
            
            # Cleanup should be skipped
            fm.cleanup_intermediate_files()
            
            # File should still exist
            self.assertTrue(test_file.exists())
        
        finally:
            fm.cleanup_all()
    
    def test_get_file_size_nonexistent(self):
        """Test get_file_size with non-existent file"""
        with FileManager() as fm:
            with self.assertRaises(FileNotFoundError):
                fm.get_file_size("nonexistent.txt")
    
    def test_work_dir_prefix(self):
        """Test custom work directory prefix"""
        fm = FileManager(debug=True, work_dir_prefix="custom_prefix")
        
        try:
            temp_path = fm.get_temp_path()
            self.assertIn("custom_prefix", str(temp_path))
        finally:
            fm.cleanup_all()


def run_file_manager_tests():
    """Run the FileManager test suite"""
    print("=" * 80)
    print("FILE MANAGER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted file management functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestFileManager)
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("FILE MANAGER TEST RESULTS")
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
    print(f"\nFileManager Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_file_manager_tests()
    sys.exit(0 if success else 1)