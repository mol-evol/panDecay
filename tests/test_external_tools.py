#!/usr/bin/env python3
"""
Test Suite for ExternalTools Component

Tests the extracted external tool execution functionality with mock testing
to ensure it behaves identically to the original implementation.
"""

import unittest
import tempfile
import os
import subprocess
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock, call
import time

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis.external_tools import (
    ExternalTools, ToolExecutionResult, ExternalToolError, 
    ToolTimeoutError, ToolExecutionError
)


class TestExternalTools(unittest.TestCase):
    """Test ExternalTools component functionality"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.tools = ExternalTools(debug=True)
        
        # Create test command file
        self.test_cmd_file = self.temp_dir / "test_commands.nex"
        with open(self.test_cmd_file, 'w') as f:
            f.write("begin paup; quit; end;")
        
        # Create test nexus file
        self.test_nexus_file = self.temp_dir / "test.nex"
        with open(self.test_nexus_file, 'w') as f:
            f.write("#NEXUS\\nbegin mrbayes; quit; end;")
        
        self.log_file = self.temp_dir / "test.log"
    
    def tearDown(self):
        """Clean up test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test ExternalTools initialization"""
        tools = ExternalTools(
            paup_path="/usr/bin/paup",
            mrbayes_path="/usr/bin/mb",
            mpirun_path="/usr/bin/mpirun",
            debug=True
        )
        
        self.assertEqual(tools.paup_path, "/usr/bin/paup")
        self.assertEqual(tools.mrbayes_path, "/usr/bin/mb")
        self.assertEqual(tools.mpirun_path, "/usr/bin/mpirun")
        self.assertTrue(tools.debug)
        self.assertEqual(tools._tool_validation_cache, {})
    
    def test_tool_execution_result(self):
        """Test ToolExecutionResult object"""
        result = ToolExecutionResult(
            args=["paup", "-n", "test.nex"],
            returncode=0,
            stdout="PAUP output",
            stderr="",
            execution_time=1.5
        )
        
        self.assertEqual(result.args, ["paup", "-n", "test.nex"])
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "PAUP output")
        self.assertEqual(result.stderr, "")
        self.assertEqual(result.execution_time, 1.5)
    
    @patch('subprocess.run')
    def test_validate_tool_success(self, mock_subprocess):
        """Test successful tool validation"""
        # Mock successful tool validation
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_subprocess.return_value = mock_result
        
        result = self.tools.validate_tool("/usr/bin/paup", "PAUP*")
        
        self.assertTrue(result)
        mock_subprocess.assert_called_once_with(
            ["/usr/bin/paup", "-h"],
            capture_output=True,
            text=True,
            timeout=10
        )
        
        # Test caching
        result2 = self.tools.validate_tool("/usr/bin/paup", "PAUP*")
        self.assertTrue(result2)
        # Should not call subprocess again due to caching
        self.assertEqual(mock_subprocess.call_count, 1)
    
    @patch('subprocess.run')
    def test_validate_tool_failure(self, mock_subprocess):
        """Test tool validation failure"""
        # Mock tool not found
        mock_subprocess.side_effect = FileNotFoundError("Tool not found")
        
        result = self.tools.validate_tool("/nonexistent/tool", "NonexistentTool")
        
        self.assertFalse(result)
        self.assertIn("/nonexistent/tool", self.tools._tool_validation_cache)
        self.assertFalse(self.tools._tool_validation_cache["/nonexistent/tool"])
    
    @patch('subprocess.Popen')
    def test_run_paup_command_file_success(self, mock_popen):
        """Test successful PAUP* command execution"""
        # Mock successful PAUP* execution
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("PAUP output", "")
        mock_process.returncode = 0
        mock_popen.return_value = mock_process
        
        result = self.tools.run_paup_command_file(
            self.test_cmd_file,
            self.log_file,
            self.temp_dir
        )
        
        self.assertIsInstance(result, ToolExecutionResult)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "PAUP output")
        self.assertEqual(result.stderr, "")
        self.assertTrue(self.log_file.exists())
        
        # Verify Popen was called correctly
        mock_popen.assert_called_once_with(
            ["paup", "-n", "test_commands.nex"],
            cwd=str(self.temp_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )
    
    @patch('subprocess.Popen')
    def test_run_paup_command_file_failure(self, mock_popen):
        """Test PAUP* command execution failure"""
        # Mock failed PAUP* execution
        mock_process = MagicMock()
        mock_process.communicate.return_value = ("", "PAUP error")
        mock_process.returncode = 1
        mock_popen.return_value = mock_process
        
        with self.assertRaises(ToolExecutionError) as context:
            self.tools.run_paup_command_file(
                self.test_cmd_file,
                self.log_file,
                self.temp_dir
            )
        
        self.assertIn("PAUP* failed with exit code 1", str(context.exception))
    
    @patch('subprocess.Popen')
    def test_run_paup_command_file_timeout(self, mock_popen):
        """Test PAUP* command execution timeout"""
        # Mock timeout
        mock_process = MagicMock()
        # Set up the timeout to occur on the first call with timeout, 
        # then succeed on the second call (after kill)
        mock_process.communicate.side_effect = [
            subprocess.TimeoutExpired("paup", 5),  # First call times out
            ("", "")  # Second call after kill succeeds
        ]
        mock_process.returncode = None
        mock_popen.return_value = mock_process
        
        with self.assertRaises(ToolTimeoutError) as context:
            self.tools.run_paup_command_file(
                self.test_cmd_file,
                self.log_file,
                self.temp_dir,
                timeout_sec=5
            )
        
        self.assertIn("timed out after 5 seconds", str(context.exception))
        mock_process.kill.assert_called_once()
    
    def test_run_paup_command_file_missing_file(self):
        """Test PAUP* execution with missing command file"""
        missing_file = self.temp_dir / "missing.nex"
        
        with self.assertRaises(FileNotFoundError) as context:
            self.tools.run_paup_command_file(
                missing_file,
                self.log_file,
                self.temp_dir
            )
        
        self.assertIn("command file not found", str(context.exception))
    
    @patch('subprocess.run')
    def test_run_mrbayes_success(self, mock_subprocess):
        """Test successful MrBayes execution"""
        # Mock successful MrBayes execution
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "MrBayes completed"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result
        
        result = self.tools.run_mrbayes(
            self.test_nexus_file,
            self.temp_dir
        )
        
        self.assertIsInstance(result, ToolExecutionResult)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "MrBayes completed")
        self.assertEqual(result.stderr, "")
        
        mock_subprocess.assert_called_once_with(
            ["mb", "test.nex"],
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
            timeout=None
        )
    
    @patch('subprocess.run')
    def test_run_mrbayes_with_mpi(self, mock_subprocess):
        """Test MrBayes execution with MPI"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "MrBayes MPI completed"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result
        
        result = self.tools.run_mrbayes(
            self.test_nexus_file,
            self.temp_dir,
            use_mpi=True,
            mpi_processors=4
        )
        
        self.assertEqual(result.returncode, 0)
        
        mock_subprocess.assert_called_once_with(
            ["mpirun", "-np", "4", "mb", "test.nex"],
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
            timeout=None
        )
    
    @patch('subprocess.run')
    def test_run_mrbayes_failure(self, mock_subprocess):
        """Test MrBayes execution failure"""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = "Error in MrBayes"
        mock_result.stderr = "MrBayes stderr"
        mock_subprocess.return_value = mock_result
        
        with self.assertRaises(ToolExecutionError) as context:
            self.tools.run_mrbayes(
                self.test_nexus_file,
                self.temp_dir
            )
        
        self.assertIn("MrBayes failed with exit code 1", str(context.exception))
    
    @patch('subprocess.run')
    def test_run_mrbayes_timeout(self, mock_subprocess):
        """Test MrBayes execution timeout"""
        mock_subprocess.side_effect = subprocess.TimeoutExpired("mb", 10)
        
        with self.assertRaises(ToolTimeoutError) as context:
            self.tools.run_mrbayes(
                self.test_nexus_file,
                self.temp_dir,
                timeout_sec=10
            )
        
        self.assertIn("timed out after 10 seconds", str(context.exception))
    
    def test_run_mrbayes_missing_file(self):
        """Test MrBayes execution with missing nexus file"""
        missing_file = self.temp_dir / "missing.nex"
        
        with self.assertRaises(FileNotFoundError) as context:
            self.tools.run_mrbayes(
                missing_file,
                self.temp_dir
            )
        
        self.assertIn("nexus file not found", str(context.exception))
    
    @patch('subprocess.run')
    def test_run_generic_command_success(self, mock_subprocess):
        """Test successful generic command execution"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "Command output"
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result
        
        result = self.tools.run_generic_command(
            ["echo", "test"],
            self.temp_dir
        )
        
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "Command output")
        
        mock_subprocess.assert_called_once_with(
            ["echo", "test"],
            cwd=str(self.temp_dir),
            capture_output=True,
            text=True,
            timeout=None
        )
    
    @patch('subprocess.run')
    def test_run_generic_command_timeout(self, mock_subprocess):
        """Test generic command timeout"""
        mock_subprocess.side_effect = subprocess.TimeoutExpired(["sleep", "10"], 5)
        
        with self.assertRaises(ToolTimeoutError):
            self.tools.run_generic_command(
                ["sleep", "10"],
                self.temp_dir,
                timeout_sec=5
            )
    
    @patch('subprocess.run')
    def test_get_tool_version_success(self, mock_subprocess):
        """Test successful tool version retrieval"""
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "PAUP* Version 4.0\\nBuild info..."
        mock_result.stderr = ""
        mock_subprocess.return_value = mock_result
        
        version = self.tools.get_tool_version("paup", "--version")
        
        self.assertEqual(version, "PAUP* Version 4.0")
        
        mock_subprocess.assert_called_once_with(
            ["paup", "--version"],
            capture_output=True,
            text=True,
            timeout=10
        )
    
    @patch('subprocess.run')
    def test_get_tool_version_stderr(self, mock_subprocess):
        """Test tool version from stderr"""
        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stdout = ""
        mock_result.stderr = "Tool Version 1.0\\nUsage..."
        mock_subprocess.return_value = mock_result
        
        version = self.tools.get_tool_version("tool", "-v")
        
        self.assertEqual(version, "Tool Version 1.0")
    
    @patch('subprocess.run')
    def test_get_tool_version_failure(self, mock_subprocess):
        """Test tool version retrieval failure"""
        mock_subprocess.side_effect = FileNotFoundError("Tool not found")
        
        version = self.tools.get_tool_version("nonexistent")
        
        self.assertIsNone(version)
    
    def test_exception_hierarchy(self):
        """Test exception class hierarchy"""
        # Test that specific exceptions inherit from base
        self.assertTrue(issubclass(ToolTimeoutError, ExternalToolError))
        self.assertTrue(issubclass(ToolExecutionError, ExternalToolError))
        
        # Test exception creation
        timeout_error = ToolTimeoutError("Timeout message")
        self.assertIn("Timeout message", str(timeout_error))
        
        exec_error = ToolExecutionError("Execution failed", 1, "stdout", "stderr")
        self.assertIn("Execution failed", str(exec_error))


class TestExternalToolsIntegration(unittest.TestCase):
    """Integration tests for ExternalTools component"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        self.tools = ExternalTools(debug=False)
    
    def tearDown(self):
        """Clean up integration test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_realistic_command_workflow(self):
        """Test a realistic command execution workflow"""
        # Create a realistic command file
        cmd_file = self.temp_dir / "realistic.nex"
        with open(cmd_file, 'w') as f:
            f.write("""
#NEXUS
begin paup;
    set warnroot=no;
    log file=output.log;
    quit;
end;
""")
        
        log_file = self.temp_dir / "execution.log"
        
        # This would actually try to run PAUP*, so we mock it
        with patch('subprocess.Popen') as mock_popen:
            mock_process = MagicMock()
            mock_process.communicate.return_value = ("Execution complete", "")
            mock_process.returncode = 0
            mock_popen.return_value = mock_process
            
            result = self.tools.run_paup_command_file(
                cmd_file,
                log_file,
                self.temp_dir,
                timeout_sec=30
            )
            
            self.assertEqual(result.returncode, 0)
            self.assertTrue(log_file.exists())


def run_external_tools_tests():
    """Run the ExternalTools test suite"""
    print("=" * 80)
    print("EXTERNAL TOOLS COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted external tool execution functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestExternalTools))
    suite.addTests(loader.loadTestsFromTestCase(TestExternalToolsIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\\n" + "=" * 80)
    print("EXTERNAL TOOLS TEST RESULTS")
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
    print(f"\\nExternalTools Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_external_tools_tests()
    sys.exit(0 if success else 1)