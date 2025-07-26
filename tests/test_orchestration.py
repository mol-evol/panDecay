"""
Comprehensive tests for analysis orchestration.

Tests the orchestration system that coordinates multiple analysis engines
and manages the overall analysis workflow.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from dataclasses import dataclass

from src.orchestration.analysis_orchestrator import (
    AnalysisOrchestrator, OrchestrationResult
)
from src.analysis.engines.base_engine import AnalysisData, AnalysisResult
from src.config.constants import AnalysisModes, AnalysisTimeouts
from src.exceptions.analysis_exceptions import ValidationError


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def sample_alignment_file(temp_dir):
    """Create a sample alignment file for testing."""
    alignment_content = """#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX=4 NCHAR=10;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    seq1 ATCGATCGAT
    seq2 ATCGATCGAA
    seq3 ATCGATCGTT
    seq4 ATCGATCGCC
  ;
END;
"""
    alignment_file = temp_dir / "test_alignment.nex"
    alignment_file.write_text(alignment_content)
    return alignment_file


@pytest.fixture
def sample_analysis_data(temp_dir, sample_alignment_file):
    """Create sample analysis data for testing."""
    return AnalysisData(
        alignment_file=sample_alignment_file,
        tree_file=None,
        temp_dir=temp_dir,
        model_settings={
            'nst': 2,
            'rates': 'gamma',
            'threads': 1
        },
        constraints=[
            {'id': 'test_constraint_1', 'taxa': ['seq1', 'seq2']},
            {'id': 'test_constraint_2', 'taxa': ['seq3', 'seq4']}
        ],
        timeouts=AnalysisTimeouts()
    )


@pytest.fixture
def orchestrator(temp_dir):
    """Create an analysis orchestrator for testing."""
    return AnalysisOrchestrator(temp_dir, debug=True)


class TestAnalysisOrchestrator:
    """Test the analysis orchestrator functionality."""
    
    def test_initialization(self, orchestrator, temp_dir):
        """Test orchestrator initialization."""
        assert orchestrator.temp_dir == temp_dir
        assert orchestrator.debug is True
        
        # Should have ML engine registered by default
        supported_modes = orchestrator.get_supported_analysis_modes()
        assert AnalysisModes.ML in supported_modes
        assert AnalysisModes.BAYESIAN in supported_modes
        assert AnalysisModes.PARSIMONY in supported_modes
    
    def test_engine_registration(self, orchestrator):
        """Test registering new analysis engines."""
        # Create a mock engine class
        MockEngine = Mock()
        MockEngine.__name__ = "MockEngine"
        
        initial_count = len(orchestrator.get_supported_analysis_modes())
        orchestrator.register_engine("mock", MockEngine)
        
        new_count = len(orchestrator.get_supported_analysis_modes())
        assert new_count == initial_count + 1
        assert "mock" in orchestrator.get_supported_analysis_modes()
    
    def test_create_engines(self, orchestrator):
        """Test creating analysis engines."""
        orchestrator.create_engines([AnalysisModes.ML], threads=2)
        
        assert AnalysisModes.ML in orchestrator._engines
        assert len(orchestrator._engines) == 1
    
    def test_create_engines_invalid_mode(self, orchestrator):
        """Test creating engines with invalid mode logs warning."""
        with patch.object(orchestrator.logger, 'warning') as mock_warning:
            orchestrator.create_engines(["invalid_mode"])
            mock_warning.assert_called_once()
    
    def test_validate_analysis_request_valid(self, orchestrator, sample_analysis_data):
        """Test validation of valid analysis request."""
        # Should not raise any exceptions
        orchestrator.validate_analysis_request([AnalysisModes.ML], sample_analysis_data)
    
    def test_validate_analysis_request_unsupported_mode(self, orchestrator, sample_analysis_data):
        """Test validation fails with unsupported analysis mode."""
        with pytest.raises(ValidationError, match="Unsupported analysis modes"):
            orchestrator.validate_analysis_request(["invalid_mode"], sample_analysis_data)
    
    def test_validate_analysis_request_no_modes(self, orchestrator, sample_analysis_data):
        """Test validation fails with no analysis modes."""
        with pytest.raises(ValidationError, match="At least one analysis mode must be specified"):
            orchestrator.validate_analysis_request([], sample_analysis_data)
    
    def test_validate_analysis_request_missing_alignment(self, orchestrator, temp_dir):
        """Test validation fails with missing alignment file."""
        invalid_data = AnalysisData(
            alignment_file=temp_dir / "nonexistent.nex",
            tree_file=None,
            temp_dir=temp_dir,
            model_settings={},
            constraints=[],
            timeouts=AnalysisTimeouts()
        )
        
        with pytest.raises(ValidationError, match="Alignment file not found"):
            orchestrator.validate_analysis_request([AnalysisModes.ML], invalid_data)
    
    def test_successful_single_analysis(self, orchestrator, sample_analysis_data):
        """Test successful single analysis orchestration."""
        # Mock the ML engine
        mock_engine = Mock()
        mock_result = AnalysisResult(
            analysis_type="ml",
            success=True,
            branches_tested=2,
            decay_indices={'branch1': 0.5, 'branch2': 0.8},
            execution_time=5.0
        )
        mock_engine.analyze.return_value = mock_result
        mock_engine.validate_inputs.return_value = True
        mock_engine.cleanup.return_value = None
        
        # Replace the engine in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_engine
        
        result = orchestrator.run_analysis([AnalysisModes.ML], sample_analysis_data)
        
        assert isinstance(result, OrchestrationResult)
        assert result.success is True
        assert result.total_execution_time > 0
        assert AnalysisModes.ML in result.analysis_results
        assert result.analysis_results[AnalysisModes.ML].success is True
    
    def test_successful_multiple_analysis(self, orchestrator, sample_analysis_data):
        """Test successful multiple analysis orchestration."""
        # Mock ML engine
        mock_ml_engine = Mock()
        mock_ml_result = AnalysisResult(
            analysis_type="ml",
            success=True,
            branches_tested=2,
            decay_indices={'branch1': 0.5},
            execution_time=3.0
        )
        mock_ml_engine.analyze.return_value = mock_ml_result
        mock_ml_engine.validate_inputs.return_value = True
        mock_ml_engine.cleanup.return_value = None
        
        # Mock Bayesian engine
        mock_bayes_engine = Mock()
        mock_bayes_result = AnalysisResult(
            analysis_type="bayesian",
            success=True,
            branches_tested=2,
            decay_indices={'branch1': 0.7},
            execution_time=7.0
        )
        mock_bayes_engine.analyze.return_value = mock_bayes_result
        mock_bayes_engine.validate_inputs.return_value = True
        mock_bayes_engine.cleanup.return_value = None
        
        # Replace engines in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_ml_engine
        orchestrator._engines[AnalysisModes.BAYESIAN] = mock_bayes_engine
        
        result = orchestrator.run_analysis(
            [AnalysisModes.ML, AnalysisModes.BAYESIAN], 
            sample_analysis_data
        )
        
        assert result.success is True
        assert len(result.analysis_results) == 2
        assert AnalysisModes.ML in result.analysis_results
        assert AnalysisModes.BAYESIAN in result.analysis_results
        assert result.metadata['successful_analyses'] == [AnalysisModes.ML, AnalysisModes.BAYESIAN]
    
    def test_partial_failure_analysis(self, orchestrator, sample_analysis_data):
        """Test orchestration with partial failures."""
        # Mock ML engine (success)
        mock_ml_engine = Mock()
        mock_ml_result = AnalysisResult(
            analysis_type="ml",
            success=True,
            branches_tested=2,
            decay_indices={'branch1': 0.5},
            execution_time=3.0
        )
        mock_ml_engine.analyze.return_value = mock_ml_result
        mock_ml_engine.validate_inputs.return_value = True
        mock_ml_engine.cleanup.return_value = None
        
        # Mock Bayesian engine (failure)
        mock_bayes_engine = Mock()
        mock_bayes_result = AnalysisResult(
            analysis_type="bayesian",
            success=False,
            branches_tested=0,
            decay_indices={},
            execution_time=1.0,
            error_message="Bayesian analysis failed"
        )
        mock_bayes_engine.analyze.return_value = mock_bayes_result
        mock_bayes_engine.validate_inputs.return_value = True
        mock_bayes_engine.cleanup.return_value = None
        
        # Replace engines in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_ml_engine
        orchestrator._engines[AnalysisModes.BAYESIAN] = mock_bayes_engine
        
        result = orchestrator.run_analysis(
            [AnalysisModes.ML, AnalysisModes.BAYESIAN], 
            sample_analysis_data
        )
        
        # Should still be successful overall if at least one analysis succeeds
        assert result.success is True
        assert len(result.analysis_results) == 2
        assert result.analysis_results[AnalysisModes.ML].success is True
        assert result.analysis_results[AnalysisModes.BAYESIAN].success is False
        assert result.metadata['successful_analyses'] == [AnalysisModes.ML]
    
    def test_complete_failure_analysis(self, orchestrator, sample_analysis_data):
        """Test orchestration with complete failure."""
        # Mock ML engine (failure)
        mock_ml_engine = Mock()
        mock_ml_result = AnalysisResult(
            analysis_type="ml",
            success=False,
            branches_tested=0,
            decay_indices={},
            execution_time=1.0,
            error_message="ML analysis failed"
        )
        mock_ml_engine.analyze.return_value = mock_ml_result
        mock_ml_engine.validate_inputs.return_value = True
        mock_ml_engine.cleanup.return_value = None
        
        # Replace engine in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_ml_engine
        
        result = orchestrator.run_analysis([AnalysisModes.ML], sample_analysis_data)
        
        assert result.success is False
        assert result.error_message == "All analyses failed"
        assert len(result.analysis_results) == 1
        assert result.analysis_results[AnalysisModes.ML].success is False
    
    def test_exception_during_analysis(self, orchestrator, sample_analysis_data):
        """Test handling of exceptions during analysis."""
        # Mock ML engine that raises exception
        mock_ml_engine = Mock()
        mock_ml_engine.analyze.side_effect = Exception("Unexpected error")
        mock_ml_engine.validate_inputs.return_value = True
        mock_ml_engine.cleanup.return_value = None
        
        # Replace engine in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_ml_engine
        
        result = orchestrator.run_analysis([AnalysisModes.ML], sample_analysis_data)
        
        assert result.success is False
        assert "Unexpected error" in result.analysis_results[AnalysisModes.ML].error_message
    
    def test_validation_error_during_orchestration(self, orchestrator):
        """Test handling of validation errors during orchestration."""
        invalid_data = AnalysisData(
            alignment_file=Path("nonexistent.nex"),
            tree_file=None,
            temp_dir=orchestrator.temp_dir,
            model_settings={},
            constraints=[],
            timeouts=AnalysisTimeouts()
        )
        
        result = orchestrator.run_analysis([AnalysisModes.ML], invalid_data)
        
        assert result.success is False
        assert "Orchestrated analysis failed" in result.error_message
        assert len(result.analysis_results) == 0
    
    def test_cleanup_engines_called(self, orchestrator, sample_analysis_data):
        """Test that engine cleanup is called after analysis."""
        # Mock ML engine
        mock_ml_engine = Mock()
        mock_ml_result = AnalysisResult(
            analysis_type="ml",
            success=True,
            branches_tested=2,
            decay_indices={'branch1': 0.5},
            execution_time=3.0
        )
        mock_ml_engine.analyze.return_value = mock_ml_result
        mock_ml_engine.validate_inputs.return_value = True
        mock_ml_engine.cleanup.return_value = None
        
        # Replace engine in orchestrator
        orchestrator._engines[AnalysisModes.ML] = mock_ml_engine
        
        result = orchestrator.run_analysis([AnalysisModes.ML], sample_analysis_data)
        
        # Verify cleanup was called
        mock_ml_engine.cleanup.assert_called_once()
    
    def test_parse_analysis_mode_string(self, orchestrator):
        """Test parsing analysis mode strings."""
        # Test single mode
        modes = orchestrator.parse_analysis_mode_string("ml")
        assert modes == ["ml"]
        
        # Test compound mode
        modes = orchestrator.parse_analysis_mode_string("ml+bayesian")
        assert modes == ["ml", "bayesian"]
        
        # Test all mode
        modes = orchestrator.parse_analysis_mode_string("all")
        supported_modes = orchestrator.get_supported_analysis_modes()
        assert set(modes) == set(supported_modes)
    
    def test_create_analysis_data(self, orchestrator, sample_alignment_file, temp_dir):
        """Test creating analysis data object."""
        data = orchestrator.create_analysis_data(
            alignment_file=sample_alignment_file,
            temp_dir=temp_dir,
            model_settings={'nst': 2},
            constraints=[{'id': 'test', 'taxa': ['seq1', 'seq2']}],
            timeouts=None,
            tree_file=None
        )
        
        assert isinstance(data, AnalysisData)
        assert data.alignment_file == sample_alignment_file
        assert data.temp_dir == temp_dir
        assert data.model_settings == {'nst': 2}
        assert len(data.constraints) == 1
        assert isinstance(data.timeouts, AnalysisTimeouts)


class TestOrchestrationResult:
    """Test the OrchestrationResult dataclass."""
    
    def test_successful_result_creation(self):
        """Test creating a successful orchestration result."""
        analysis_results = {
            'ml': AnalysisResult(
                analysis_type="ml",
                success=True,
                branches_tested=2,
                decay_indices={'branch1': 0.5},
                execution_time=5.0
            )
        }
        
        result = OrchestrationResult(
            success=True,
            total_execution_time=10.5,
            analysis_results=analysis_results,
            metadata={'test_key': 'test_value'}
        )
        
        assert result.success is True
        assert result.total_execution_time == 10.5
        assert len(result.analysis_results) == 1
        assert 'ml' in result.analysis_results
        assert result.metadata == {'test_key': 'test_value'}
        assert result.error_message is None
    
    def test_failed_result_creation(self):
        """Test creating a failed orchestration result."""
        result = OrchestrationResult(
            success=False,
            total_execution_time=1.0,
            analysis_results={},
            error_message="Orchestration failed"
        )
        
        assert result.success is False
        assert result.error_message == "Orchestration failed"
        assert len(result.analysis_results) == 0
        assert result.metadata is None


class TestOrchestrationIntegration:
    """Test integration scenarios for orchestration."""
    
    def test_full_workflow_simulation(self, orchestrator, sample_analysis_data):
        """Test a complete analysis workflow simulation."""
        # Create mock engines for all types
        engines = {}
        for analysis_type in [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY]:
            mock_engine = Mock()
            mock_result = AnalysisResult(
                analysis_type=analysis_type,
                success=True,
                branches_tested=2,
                decay_indices={f'branch_{analysis_type}': 0.5 + ord(analysis_type[0]) * 0.1},
                execution_time=3.0 + ord(analysis_type[0])
            )
            mock_engine.analyze.return_value = mock_result
            mock_engine.validate_inputs.return_value = True
            mock_engine.cleanup.return_value = None
            engines[analysis_type] = mock_engine
        
        # Replace engines in orchestrator
        orchestrator._engines.update(engines)
        
        # Run all analysis types
        result = orchestrator.run_analysis(
            [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY],
            sample_analysis_data
        )
        
        assert result.success is True
        assert len(result.analysis_results) == 3
        assert result.metadata['successful_analyses'] == [
            AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY
        ]
        
        # Verify all engines were cleaned up
        for engine in engines.values():
            engine.cleanup.assert_called_once()
    
    def test_engine_creation_with_parameters(self, orchestrator):
        """Test engine creation with custom parameters."""
        orchestrator.create_engines(
            [AnalysisModes.ML],
            paup_path="/custom/paup",
            threads=8,
            custom_param="test_value"
        )
        
        assert AnalysisModes.ML in orchestrator._engines
        # Note: We can't easily test the parameter passing without mocking 
        # the engine constructors, but this tests the interface
    
    def test_concurrent_orchestration_safety(self, temp_dir):
        """Test that multiple orchestrators can run safely."""
        orchestrator1 = AnalysisOrchestrator(temp_dir / "orch1", debug=True)
        orchestrator2 = AnalysisOrchestrator(temp_dir / "orch2", debug=True)
        
        # Both should have independent engine registries
        orchestrator1.register_engine("test1", Mock)
        orchestrator2.register_engine("test2", Mock)
        
        modes1 = orchestrator1.get_supported_analysis_modes()
        modes2 = orchestrator2.get_supported_analysis_modes()
        
        assert "test1" in modes1
        assert "test1" not in modes2
        assert "test2" in modes2
        assert "test2" not in modes1


if __name__ == "__main__":
    pytest.main([__file__])