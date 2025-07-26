"""
Comprehensive tests for exception handling.

Tests the custom exception hierarchy and error handling
throughout the panDecay codebase.
"""

import pytest
from pathlib import Path

from src.exceptions.analysis_exceptions import (
    PanDecayError, AnalysisError, ValidationError, ConfigurationError
)
from src.exceptions.tool_exceptions import (
    ExternalToolError, ToolNotFoundError, ToolExecutionError, ToolTimeoutError
)
from src.exceptions.io_exceptions import (
    FileOperationError, AlignmentParsingError, TreeParsingError
)


class TestPanDecayError:
    """Test the base PanDecayError exception."""
    
    def test_basic_error_creation(self):
        """Test creating a basic PanDecayError."""
        error = PanDecayError("Test error message")
        
        assert str(error) == "Test error message"
        assert error.context == {}
    
    def test_error_with_context(self):
        """Test creating PanDecayError with context information."""
        context = {
            'operation': 'test_operation',
            'file_path': '/path/to/file.txt',
            'timestamp': '2024-01-01T12:00:00'
        }
        
        error = PanDecayError("Test error with context", context=context)
        
        assert str(error) == "Test error with context"
        assert error.context == context
        assert error.context['operation'] == 'test_operation'
        assert error.context['file_path'] == '/path/to/file.txt'
    
    def test_error_inheritance(self):
        """Test that PanDecayError properly inherits from Exception."""
        error = PanDecayError("Test error")
        
        assert isinstance(error, Exception)
        assert isinstance(error, PanDecayError)
    
    def test_error_context_modification(self):
        """Test that error context can be modified after creation."""
        error = PanDecayError("Test error")
        
        error.context['new_key'] = 'new_value'
        assert error.context['new_key'] == 'new_value'
    
    def test_error_repr(self):
        """Test string representation of PanDecayError."""
        context = {'key': 'value'}
        error = PanDecayError("Test message", context=context)
        
        repr_str = repr(error)
        assert "PanDecayError" in repr_str
        assert "Test message" in repr_str


class TestAnalysisError:
    """Test the AnalysisError exception."""
    
    def test_analysis_error_creation(self):
        """Test creating an AnalysisError."""
        error = AnalysisError("Analysis failed", analysis_type="ml")
        
        assert str(error) == "Analysis failed"
        assert error.analysis_type == "ml"
        assert error.context['analysis_type'] == "ml"
    
    def test_analysis_error_without_type(self):
        """Test creating AnalysisError without analysis type."""
        error = AnalysisError("Analysis failed")
        
        assert str(error) == "Analysis failed"
        assert error.analysis_type is None
        assert 'analysis_type' not in error.context
    
    def test_analysis_error_with_context(self):
        """Test AnalysisError with additional context."""
        context = {'constraint_id': 'clade_1', 'iteration': 5}
        error = AnalysisError("Analysis failed", analysis_type="bayesian", context=context)
        
        assert error.analysis_type == "bayesian"
        assert error.context['analysis_type'] == "bayesian"
        assert error.context['constraint_id'] == 'clade_1'
        assert error.context['iteration'] == 5
    
    def test_analysis_error_inheritance(self):
        """Test AnalysisError inheritance."""
        error = AnalysisError("Test error")
        
        assert isinstance(error, PanDecayError)
        assert isinstance(error, AnalysisError)


class TestValidationError:
    """Test the ValidationError exception."""
    
    def test_validation_error_creation(self):
        """Test creating a ValidationError."""
        error = ValidationError("Invalid input", field="threads", value=0)
        
        assert str(error) == "Invalid input"
        assert error.field == "threads"
        assert error.value == 0
        assert error.context['field'] == "threads"
        assert error.context['value'] == 0
    
    def test_validation_error_without_field_and_value(self):
        """Test ValidationError without field and value."""
        error = ValidationError("General validation error")
        
        assert str(error) == "General validation error"
        assert error.field is None
        assert error.value is None
        assert 'field' not in error.context
        assert 'value' not in error.context
    
    def test_validation_error_with_complex_value(self):
        """Test ValidationError with complex value types."""
        complex_value = {'nested': {'key': 'value'}, 'list': [1, 2, 3]}
        error = ValidationError("Complex validation error", field="config", value=complex_value)
        
        assert error.value == complex_value
        assert error.context['value'] == complex_value
    
    def test_validation_error_inheritance(self):
        """Test ValidationError inheritance."""
        error = ValidationError("Test error")
        
        assert isinstance(error, PanDecayError)
        assert isinstance(error, ValidationError)


class TestConfigurationError:
    """Test the ConfigurationError exception."""
    
    def test_configuration_error_creation(self):
        """Test creating a ConfigurationError."""
        error = ConfigurationError("Config error", config_file="config.ini", section="analysis")
        
        assert str(error) == "Config error"
        assert error.config_file == "config.ini"
        assert error.section == "analysis"
        assert error.context['config_file'] == "config.ini"
        assert error.context['section'] == "analysis"
    
    def test_configuration_error_minimal(self):
        """Test ConfigurationError with minimal parameters."""
        error = ConfigurationError("Config error")
        
        assert str(error) == "Config error"
        assert error.config_file is None
        assert error.section is None
    
    def test_configuration_error_inheritance(self):
        """Test ConfigurationError inheritance."""
        error = ConfigurationError("Test error")
        
        assert isinstance(error, PanDecayError)
        assert isinstance(error, ConfigurationError)


class TestExternalToolError:
    """Test the ExternalToolError exception."""
    
    def test_external_tool_error_creation(self):
        """Test creating an ExternalToolError."""
        error = ExternalToolError("Tool failed", tool_name="paup", return_code=1)
        
        assert str(error) == "Tool failed"
        assert error.tool_name == "paup"
        assert error.return_code == 1
        assert error.context['tool_name'] == "paup"
        assert error.context['return_code'] == 1
    
    def test_external_tool_error_with_output(self):
        """Test ExternalToolError with stdout and stderr."""
        context = {
            'stdout': 'Tool output message',
            'stderr': 'Tool error message',
            'command': 'paup -n script.nex'
        }
        error = ExternalToolError("Tool execution failed", tool_name="paup", context=context)
        
        assert error.context['stdout'] == 'Tool output message'
        assert error.context['stderr'] == 'Tool error message'
        assert error.context['command'] == 'paup -n script.nex'
    
    def test_external_tool_error_minimal(self):
        """Test ExternalToolError with minimal parameters."""
        error = ExternalToolError("Tool error")
        
        assert str(error) == "Tool error"
        assert error.tool_name is None
        assert error.return_code is None
    
    def test_external_tool_error_inheritance(self):
        """Test ExternalToolError inheritance."""
        error = ExternalToolError("Test error")
        
        assert isinstance(error, PanDecayError)
        assert isinstance(error, ExternalToolError)


class TestToolNotFoundError:
    """Test the ToolNotFoundError exception."""
    
    def test_tool_not_found_error_creation(self):
        """Test creating a ToolNotFoundError."""
        error = ToolNotFoundError("PAUP* not found", tool_name="paup", search_paths=["/usr/bin", "/usr/local/bin"])
        
        assert str(error) == "PAUP* not found"
        assert error.tool_name == "paup"
        assert error.search_paths == ["/usr/bin", "/usr/local/bin"]
        assert error.context['search_paths'] == ["/usr/bin", "/usr/local/bin"]
    
    def test_tool_not_found_error_inheritance(self):
        """Test ToolNotFoundError inheritance."""
        error = ToolNotFoundError("Test error")
        
        assert isinstance(error, ExternalToolError)
        assert isinstance(error, ToolNotFoundError)


class TestToolExecutionError:
    """Test the ToolExecutionError exception."""
    
    def test_tool_execution_error_creation(self):
        """Test creating a ToolExecutionError."""
        error = ToolExecutionError("Execution failed", tool_name="mrbayes", command="mb script.nex")
        
        assert str(error) == "Execution failed"
        assert error.tool_name == "mrbayes"
        assert error.command == "mb script.nex"
        assert error.context['command'] == "mb script.nex"
    
    def test_tool_execution_error_inheritance(self):
        """Test ToolExecutionError inheritance."""
        error = ToolExecutionError("Test error")
        
        assert isinstance(error, ExternalToolError)
        assert isinstance(error, ToolExecutionError)


class TestToolTimeoutError:
    """Test the ToolTimeoutError exception."""
    
    def test_tool_timeout_error_creation(self):
        """Test creating a ToolTimeoutError."""
        error = ToolTimeoutError("Tool timed out", tool_name="paup", timeout_seconds=3600)
        
        assert str(error) == "Tool timed out"
        assert error.tool_name == "paup"
        assert error.timeout_seconds == 3600
        assert error.context['timeout_seconds'] == 3600
    
    def test_tool_timeout_error_inheritance(self):
        """Test ToolTimeoutError inheritance."""
        error = ToolTimeoutError("Test error")
        
        assert isinstance(error, ExternalToolError)
        assert isinstance(error, ToolTimeoutError)


class TestFileOperationError:
    """Test the FileOperationError exception."""
    
    def test_file_operation_error_creation(self):
        """Test creating a FileOperationError."""
        file_path = Path("/path/to/file.txt")
        error = FileOperationError("File operation failed", file_path=file_path, operation="read")
        
        assert str(error) == "File operation failed"
        assert error.file_path == file_path
        assert error.operation == "read"
        assert error.context['file_path'] == str(file_path)
        assert error.context['operation'] == "read"
    
    def test_file_operation_error_minimal(self):
        """Test FileOperationError with minimal parameters."""
        error = FileOperationError("File error")
        
        assert str(error) == "File error"
        assert error.file_path is None
        assert error.operation is None
    
    def test_file_operation_error_inheritance(self):
        """Test FileOperationError inheritance."""
        error = FileOperationError("Test error")
        
        assert isinstance(error, PanDecayError)
        assert isinstance(error, FileOperationError)


class TestAlignmentParsingError:
    """Test the AlignmentParsingError exception."""
    
    def test_alignment_parsing_error_creation(self):
        """Test creating an AlignmentParsingError."""
        file_path = Path("/path/to/alignment.fas")
        error = AlignmentParsingError("Parsing failed", file_path=file_path, format_type="fasta", line_number=10)
        
        assert str(error) == "Parsing failed"
        assert error.file_path == file_path
        assert error.format_type == "fasta"
        assert error.line_number == 10
        assert error.context['format_type'] == "fasta"
        assert error.context['line_number'] == 10
    
    def test_alignment_parsing_error_inheritance(self):
        """Test AlignmentParsingError inheritance."""
        error = AlignmentParsingError("Test error")
        
        assert isinstance(error, FileOperationError)
        assert isinstance(error, AlignmentParsingError)


class TestTreeParsingError:
    """Test the TreeParsingError exception."""
    
    def test_tree_parsing_error_creation(self):
        """Test creating a TreeParsingError."""
        file_path = Path("/path/to/tree.nwk")
        error = TreeParsingError("Tree parsing failed", file_path=file_path, format_type="newick", position=25)
        
        assert str(error) == "Tree parsing failed"
        assert error.file_path == file_path
        assert error.format_type == "newick"
        assert error.position == 25
        assert error.context['format_type'] == "newick"
        assert error.context['position'] == 25
    
    def test_tree_parsing_error_inheritance(self):
        """Test TreeParsingError inheritance."""
        error = TreeParsingError("Test error")
        
        assert isinstance(error, FileOperationError)
        assert isinstance(error, TreeParsingError)


class TestExceptionHierarchy:
    """Test the overall exception hierarchy."""
    
    def test_hierarchy_structure(self):
        """Test that the exception hierarchy is properly structured."""
        # Test base exception
        base_error = PanDecayError("Base error")
        assert isinstance(base_error, Exception)
        
        # Test analysis exceptions
        analysis_error = AnalysisError("Analysis error")
        assert isinstance(analysis_error, PanDecayError)
        
        validation_error = ValidationError("Validation error")
        assert isinstance(validation_error, PanDecayError)
        
        config_error = ConfigurationError("Config error")
        assert isinstance(config_error, PanDecayError)
        
        # Test tool exceptions
        tool_error = ExternalToolError("Tool error")
        assert isinstance(tool_error, PanDecayError)
        
        not_found_error = ToolNotFoundError("Not found")
        assert isinstance(not_found_error, ExternalToolError)
        assert isinstance(not_found_error, PanDecayError)
        
        execution_error = ToolExecutionError("Execution error")
        assert isinstance(execution_error, ExternalToolError)
        assert isinstance(execution_error, PanDecayError)
        
        timeout_error = ToolTimeoutError("Timeout error")
        assert isinstance(timeout_error, ExternalToolError)
        assert isinstance(timeout_error, PanDecayError)
        
        # Test IO exceptions
        file_error = FileOperationError("File error")
        assert isinstance(file_error, PanDecayError)
        
        alignment_error = AlignmentParsingError("Alignment error")
        assert isinstance(alignment_error, FileOperationError)
        assert isinstance(alignment_error, PanDecayError)
        
        tree_error = TreeParsingError("Tree error")
        assert isinstance(tree_error, FileOperationError)
        assert isinstance(tree_error, PanDecayError)
    
    def test_exception_catching_hierarchy(self):
        """Test that exceptions can be caught at different levels."""
        # Test catching specific exception
        try:
            raise AnalysisError("Specific analysis error")
        except AnalysisError as e:
            assert isinstance(e, AnalysisError)
        
        # Test catching at base level
        try:
            raise ValidationError("Validation error")
        except PanDecayError as e:
            assert isinstance(e, ValidationError)
            assert isinstance(e, PanDecayError)
        
        # Test catching tool exceptions
        try:
            raise ToolNotFoundError("Tool not found")
        except ExternalToolError as e:
            assert isinstance(e, ToolNotFoundError)
        except PanDecayError as e:
            pytest.fail("Should have been caught as ExternalToolError")


class TestExceptionUsagePatterns:
    """Test common exception usage patterns."""
    
    def test_chaining_exceptions(self):
        """Test exception chaining with from clause."""
        original_error = ValueError("Original error")
        
        try:
            raise AnalysisError("Analysis failed") from original_error
        except AnalysisError as e:
            assert e.__cause__ is original_error
            assert str(e.__cause__) == "Original error"
    
    def test_context_preservation_across_raises(self):
        """Test that context is preserved when re-raising."""
        context = {'operation': 'test', 'data': 'important'}
        original_error = AnalysisError("Original error", context=context)
        
        try:
            try:
                raise original_error
            except AnalysisError as e:
                # Re-raise with additional context
                e.context['re_raised'] = True
                raise
        except AnalysisError as final_error:
            assert final_error.context['operation'] == 'test'
            assert final_error.context['data'] == 'important'
            assert final_error.context['re_raised'] is True
    
    def test_exception_with_nested_context(self):
        """Test exceptions with deeply nested context."""
        nested_context = {
            'analysis': {
                'type': 'ml',
                'constraints': [
                    {'id': 'clade_1', 'taxa': ['seq1', 'seq2']},
                    {'id': 'clade_2', 'taxa': ['seq3', 'seq4']}
                ]
            },
            'system': {
                'memory_mb': 1024,
                'threads': 4
            }
        }
        
        error = AnalysisError("Complex analysis error", context=nested_context)
        
        assert error.context['analysis']['type'] == 'ml'
        assert len(error.context['analysis']['constraints']) == 2
        assert error.context['system']['memory_mb'] == 1024
    
    def test_exception_message_formatting(self):
        """Test consistent exception message formatting."""
        # Test that error messages can include context information
        file_path = Path("/path/to/file.txt")
        error = FileOperationError(
            f"Failed to read file: {file_path}",
            file_path=file_path,
            operation="read"
        )
        
        assert str(file_path) in str(error)
        assert error.context['file_path'] == str(file_path)


if __name__ == "__main__":
    pytest.main([__file__])