"""
Comprehensive tests for analysis engines.

Tests the modular analysis engine architecture including ML, Bayesian,
and Parsimony engines with proper mocking of external tools.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from dataclasses import dataclass

from src.analysis.engines.base_engine import AnalysisEngine, AnalysisData, AnalysisResult
from src.analysis.engines.ml_engine import MLAnalysisEngine
from src.analysis.engines.bayesian_engine import BayesianAnalysisEngine
from src.analysis.engines.parsimony_engine import ParsimonyAnalysisEngine
from src.config.constants import AnalysisTimeouts, FileNames
from src.exceptions.analysis_exceptions import AnalysisError, ValidationError


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


class TestBaseAnalysisEngine:
    """Test the base analysis engine functionality."""
    
    def test_abstract_methods_raise_not_implemented(self, temp_dir):
        """Test that abstract methods raise NotImplementedError."""
        
        class TestEngine(AnalysisEngine):
            def get_analysis_type(self):
                return "test"
            
            def validate_inputs(self, data):
                return True
        
        engine = TestEngine(temp_dir, debug=False)
        
        with pytest.raises(NotImplementedError):
            engine.analyze(Mock())
    
    def test_prepare_analysis_creates_directory(self, temp_dir, sample_analysis_data):
        """Test that prepare_analysis creates analysis directory."""
        
        class TestEngine(AnalysisEngine):
            def get_analysis_type(self):
                return "test"
            
            def validate_inputs(self, data):
                return True
            
            def analyze(self, data):
                return Mock()
        
        engine = TestEngine(temp_dir, debug=False)
        engine.prepare_analysis(sample_analysis_data)
        
        assert engine.analysis_temp_dir.exists()
        assert engine.analysis_temp_dir.is_dir()
    
    def test_cleanup_removes_temp_directory(self, temp_dir, sample_analysis_data):
        """Test that cleanup removes temporary directory."""
        
        class TestEngine(AnalysisEngine):
            def get_analysis_type(self):
                return "test"
            
            def validate_inputs(self, data):
                return True
            
            def analyze(self, data):
                return Mock()
        
        engine = TestEngine(temp_dir, debug=False)
        engine.prepare_analysis(sample_analysis_data)
        
        temp_analysis_dir = engine.analysis_temp_dir
        assert temp_analysis_dir.exists()
        
        engine.cleanup()
        assert not temp_analysis_dir.exists()
    
    def test_debug_mode_preserves_temp_directory(self, temp_dir, sample_analysis_data):
        """Test that debug mode preserves temporary directory."""
        
        class TestEngine(AnalysisEngine):
            def get_analysis_type(self):
                return "test"
            
            def validate_inputs(self, data):
                return True
            
            def analyze(self, data):
                return Mock()
        
        engine = TestEngine(temp_dir, debug=True)
        engine.prepare_analysis(sample_analysis_data)
        
        temp_analysis_dir = engine.analysis_temp_dir
        assert temp_analysis_dir.exists()
        
        engine.cleanup()
        assert temp_analysis_dir.exists()  # Should still exist in debug mode


class TestMLAnalysisEngine:
    """Test the Maximum Likelihood analysis engine."""
    
    @pytest.fixture
    def ml_engine(self, temp_dir):
        """Create ML analysis engine for testing."""
        return MLAnalysisEngine(temp_dir, paup_path="mock_paup", debug=True)
    
    def test_initialization(self, ml_engine):
        """Test ML engine initialization."""
        assert ml_engine.get_analysis_type() == "ml"
        assert ml_engine.paup_path == "mock_paup"
        assert ml_engine.debug is True
    
    def test_validate_inputs_missing_alignment(self, ml_engine, temp_dir):
        """Test validation fails with missing alignment file."""
        data = AnalysisData(
            alignment_file=temp_dir / "nonexistent.nex",
            tree_file=None,
            temp_dir=temp_dir,
            model_settings={},
            constraints=[],
            timeouts=AnalysisTimeouts()
        )
        
        with pytest.raises(AnalysisError, match="Alignment file not found"):
            ml_engine.validate_inputs(data)
    
    @patch('src.analysis.engines.ml_engine.ExternalToolManager')
    def test_validate_inputs_paup_unavailable(self, mock_tool_manager, ml_engine, sample_analysis_data):
        """Test validation fails when PAUP* is unavailable."""
        mock_manager = Mock()
        mock_manager.__enter__ = Mock(return_value=mock_manager)
        mock_manager.__exit__ = Mock(return_value=None)
        mock_manager.check_tool_availability.side_effect = Exception("PAUP* not found")
        mock_tool_manager.return_value = mock_manager
        
        with pytest.raises(AnalysisError, match="PAUP\\* not available"):
            ml_engine.validate_inputs(sample_analysis_data)
    
    @patch('src.analysis.engines.ml_engine.ExternalToolManager')
    def test_successful_analysis_workflow(self, mock_tool_manager, ml_engine, sample_analysis_data):
        """Test successful ML analysis workflow."""
        # Mock ExternalToolManager
        mock_manager = Mock()
        mock_manager.__enter__ = Mock(return_value=mock_manager)
        mock_manager.__exit__ = Mock(return_value=None)
        mock_manager.check_tool_availability.return_value = None
        mock_manager.execute_script_file.return_value = None
        mock_tool_manager.return_value = mock_manager
        
        # Mock likelihood parsing
        with patch.object(ml_engine, '_parse_likelihood_score', return_value=-1000.5):
            with patch.object(ml_engine, '_test_constraints', return_value={}):
                with patch.object(ml_engine, '_calculate_decay_indices', return_value={}):
                    
                    result = ml_engine.analyze(sample_analysis_data)
                    
                    assert isinstance(result, AnalysisResult)
                    assert result.analysis_type == "ml"
                    assert result.success is True
                    assert result.execution_time > 0
    
    def test_parse_likelihood_score_valid_file(self, ml_engine, temp_dir):
        """Test parsing likelihood score from valid score file."""
        score_file = temp_dir / "test_score.txt"
        score_content = "Tree\t-lnL\tLength\n1\t1000.5\t0.123\n"
        score_file.write_text(score_content)
        
        likelihood = ml_engine._parse_likelihood_score(score_file)
        assert likelihood == 1000.5
    
    def test_parse_likelihood_score_missing_file(self, ml_engine, temp_dir):
        """Test parsing likelihood score from missing file raises error."""
        score_file = temp_dir / "nonexistent.txt"
        
        with pytest.raises(AnalysisError, match="Score file not found"):
            ml_engine._parse_likelihood_score(score_file)
    
    def test_parse_likelihood_score_invalid_format(self, ml_engine, temp_dir):
        """Test parsing likelihood score from invalid format raises error."""
        score_file = temp_dir / "test_score.txt"
        score_content = "Invalid format\n"
        score_file.write_text(score_content)
        
        with pytest.raises(AnalysisError, match="Invalid score file format"):
            ml_engine._parse_likelihood_score(score_file)


class TestBayesianAnalysisEngine:
    """Test the Bayesian analysis engine."""
    
    @pytest.fixture
    def bayesian_engine(self, temp_dir):
        """Create Bayesian analysis engine for testing."""
        return BayesianAnalysisEngine(temp_dir, mrbayes_path="mock_mb", debug=True)
    
    def test_initialization(self, bayesian_engine):
        """Test Bayesian engine initialization."""
        assert bayesian_engine.get_analysis_type() == "bayesian"
        assert bayesian_engine.mrbayes_path == "mock_mb"
        assert bayesian_engine.debug is True
    
    def test_validate_inputs_missing_alignment(self, bayesian_engine, temp_dir):
        """Test validation fails with missing alignment file."""
        data = AnalysisData(
            alignment_file=temp_dir / "nonexistent.nex",
            tree_file=None,
            temp_dir=temp_dir,
            model_settings={},
            constraints=[],
            timeouts=AnalysisTimeouts()
        )
        
        with pytest.raises(AnalysisError, match="Alignment file not found"):
            bayesian_engine.validate_inputs(data)
    
    @patch('src.analysis.engines.bayesian_engine.ExternalToolManager')
    def test_validate_inputs_mrbayes_unavailable(self, mock_tool_manager, bayesian_engine, sample_analysis_data):
        """Test validation fails when MrBayes is unavailable."""
        mock_manager = Mock()
        mock_manager.__enter__ = Mock(return_value=mock_manager)
        mock_manager.__exit__ = Mock(return_value=None)
        mock_manager.check_tool_availability.side_effect = Exception("MrBayes not found")
        mock_tool_manager.return_value = mock_manager
        
        with pytest.raises(AnalysisError, match="MrBayes not available"):
            bayesian_engine.validate_inputs(sample_analysis_data)
    
    def test_parse_marginal_likelihood_valid_log(self, bayesian_engine, temp_dir):
        """Test parsing marginal likelihood from valid log file."""
        log_file = temp_dir / "test.log"
        log_content = """
MrBayes log file

Stepping-stone sampler marginal likelihood = -1500.25

Analysis completed.
"""
        log_file.write_text(log_content)
        
        marginal_likelihood = bayesian_engine._parse_marginal_likelihood(log_file)
        assert marginal_likelihood == -1500.25
    
    def test_parse_marginal_likelihood_missing_file(self, bayesian_engine, temp_dir):
        """Test parsing marginal likelihood from missing file returns None."""
        log_file = temp_dir / "nonexistent.log"
        
        marginal_likelihood = bayesian_engine._parse_marginal_likelihood(log_file)
        assert marginal_likelihood is None


class TestParsimonyAnalysisEngine:
    """Test the Parsimony analysis engine."""
    
    @pytest.fixture
    def parsimony_engine(self, temp_dir):
        """Create Parsimony analysis engine for testing."""
        return ParsimonyAnalysisEngine(temp_dir, paup_path="mock_paup", debug=True)
    
    def test_initialization(self, parsimony_engine):
        """Test Parsimony engine initialization."""
        assert parsimony_engine.get_analysis_type() == "parsimony"
        assert parsimony_engine.paup_path == "mock_paup"
        assert parsimony_engine.debug is True
    
    def test_parse_parsimony_score_valid_file(self, parsimony_engine, temp_dir):
        """Test parsing parsimony score from valid score file."""
        score_file = temp_dir / "test_score.txt"
        score_content = "Tree\tLength\tSteps\n1\t100\t150\n"
        score_file.write_text(score_content)
        
        score = parsimony_engine._parse_parsimony_score(score_file)
        assert score == 150
    
    def test_parse_parsimony_score_missing_file(self, parsimony_engine, temp_dir):
        """Test parsing parsimony score from missing file raises error."""
        score_file = temp_dir / "nonexistent.txt"
        
        with pytest.raises(AnalysisError, match="Score file not found"):
            parsimony_engine._parse_parsimony_score(score_file)


class TestAnalysisResult:
    """Test the AnalysisResult dataclass."""
    
    def test_analysis_result_creation(self):
        """Test creating an AnalysisResult instance."""
        result = AnalysisResult(
            analysis_type="test",
            success=True,
            branches_tested=5,
            decay_indices={'branch1': 0.5, 'branch2': 0.8},
            execution_time=10.5,
            metadata={'test_key': 'test_value'}
        )
        
        assert result.analysis_type == "test"
        assert result.success is True
        assert result.branches_tested == 5
        assert result.decay_indices == {'branch1': 0.5, 'branch2': 0.8}
        assert result.execution_time == 10.5
        assert result.metadata == {'test_key': 'test_value'}
        assert result.error_message is None
    
    def test_analysis_result_failure(self):
        """Test creating a failed AnalysisResult instance."""
        result = AnalysisResult(
            analysis_type="test",
            success=False,
            branches_tested=0,
            decay_indices={},
            execution_time=1.0,
            error_message="Test error"
        )
        
        assert result.success is False
        assert result.error_message == "Test error"
        assert result.branches_tested == 0
        assert result.decay_indices == {}


class TestAnalysisData:
    """Test the AnalysisData dataclass."""
    
    def test_analysis_data_creation(self, temp_dir, sample_alignment_file):
        """Test creating an AnalysisData instance."""
        data = AnalysisData(
            alignment_file=sample_alignment_file,
            tree_file=None,
            temp_dir=temp_dir,
            model_settings={'nst': 2},
            constraints=[{'id': 'test', 'taxa': ['seq1', 'seq2']}],
            timeouts=AnalysisTimeouts()
        )
        
        assert data.alignment_file == sample_alignment_file
        assert data.tree_file is None
        assert data.temp_dir == temp_dir
        assert data.model_settings == {'nst': 2}
        assert len(data.constraints) == 1
        assert isinstance(data.timeouts, AnalysisTimeouts)


if __name__ == "__main__":
    pytest.main([__file__])