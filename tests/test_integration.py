"""
Integration tests for panDecay.

Tests the complete system working together with realistic scenarios,
including end-to-end workflows and cross-module interactions.
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from src.orchestration.analysis_orchestrator import AnalysisOrchestrator
from src.analysis.engines.base_engine import AnalysisData, AnalysisResult
from src.config.constants import AnalysisModes, AnalysisTimeouts, FileNames
from src.utils.format_detectors import AlignmentFormatMapper, AnalysisTypeConfiguration
from src.utils.config_converters import ConfigurationValueConverter
from src.utils.thread_calculator import ThreadCalculator
from src.utils.memory_manager import MemoryOptimizer
from src.utils.tree_utils import prepare_starting_tree, validate_tree_file
from src.exceptions.analysis_exceptions import AnalysisError, ValidationError


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def sample_alignment_content():
    """Create sample alignment content in NEXUS format."""
    return """#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX=6 NCHAR=20;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    Human     ATCGATCGATCGATCGATCG
    Chimp     ATCGATCGATCGATCGATCG
    Gorilla   ATCGATCGATCGATCGATCC
    Orangutan ATCGATCGATCGATCGATAA
    Gibbon    ATCGATCGATCGATCGTTTT
    Macaque   ATCGATCGATCGATCGAAAA
  ;
END;
"""


@pytest.fixture
def sample_tree_content():
    """Create sample tree content in Newick format."""
    return "(Human:0.1,Chimp:0.1,(Gorilla:0.15,(Orangutan:0.2,(Gibbon:0.25,Macaque:0.3):0.1):0.05):0.05);"


@pytest.fixture
def sample_constraints():
    """Create sample constraint definitions."""
    return [
        {
            'id': 'primates',
            'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon']
        },
        {
            'id': 'great_apes',
            'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan']
        },
        {
            'id': 'hominids',
            'taxa': ['Human', 'Chimp', 'Gorilla']
        }
    ]


@pytest.fixture
def setup_test_files(temp_dir, sample_alignment_content, sample_tree_content):
    """Set up test files for integration testing."""
    alignment_file = temp_dir / "test_alignment.nex"
    alignment_file.write_text(sample_alignment_content)
    
    tree_file = temp_dir / "starting_tree.tre"
    tree_file.write_text(sample_tree_content)
    
    return {
        'alignment_file': alignment_file,
        'tree_file': tree_file,
        'temp_dir': temp_dir
    }


class TestEndToEndWorkflow:
    """Test complete end-to-end analysis workflows."""
    
    def test_single_analysis_workflow(self, setup_test_files, sample_constraints):
        """Test a complete single analysis workflow."""
        files = setup_test_files
        
        # Create orchestrator
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        # Create analysis data
        analysis_data = orchestrator.create_analysis_data(
            alignment_file=files['alignment_file'],
            temp_dir=files['temp_dir'],
            model_settings={'nst': 2, 'rates': 'gamma', 'threads': 1},
            constraints=sample_constraints,
            timeouts=AnalysisTimeouts(),
            tree_file=files['tree_file']
        )
        
        # Mock the ML engine to simulate successful analysis
        with patch('src.analysis.engines.ml_engine.ExternalToolManager') as mock_tool_manager:
            mock_manager = Mock()
            mock_manager.__enter__ = Mock(return_value=mock_manager)
            mock_manager.__exit__ = Mock(return_value=None)
            mock_manager.check_tool_availability.return_value = None
            mock_manager.execute_script_file.return_value = None
            mock_tool_manager.return_value = mock_manager
            
            # Mock likelihood parsing
            orchestrator.create_engines([AnalysisModes.ML])
            ml_engine = orchestrator._engines[AnalysisModes.ML]
            
            with patch.object(ml_engine, '_parse_likelihood_score', return_value=-1000.5):
                with patch.object(ml_engine, '_test_constraints') as mock_constraints:
                    mock_constraints.return_value = {
                        'primates': {
                            'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon'],
                            'constrained_likelihood': -1010.5,
                            'likelihood_difference': -10.0
                        }
                    }
                    
                    with patch.object(ml_engine, '_calculate_decay_indices') as mock_decay:
                        mock_decay.return_value = {
                            'primates': {
                                'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon'],
                                'likelihood_difference': -10.0,
                                'au_pvalue': 0.02,
                                'significant_au': True
                            }
                        }
                        
                        # Run the analysis
                        result = orchestrator.run_analysis([AnalysisModes.ML], analysis_data)
                        
                        # Verify results
                        assert result.success is True
                        assert AnalysisModes.ML in result.analysis_results
                        ml_result = result.analysis_results[AnalysisModes.ML]
                        assert ml_result.success is True
                        assert ml_result.branches_tested > 0
                        assert 'primates' in ml_result.decay_indices
    
    def test_multi_analysis_workflow(self, setup_test_files, sample_constraints):
        """Test workflow with multiple analysis types."""
        files = setup_test_files
        
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        # Create analysis data
        analysis_data = orchestrator.create_analysis_data(
            alignment_file=files['alignment_file'],
            temp_dir=files['temp_dir'],
            model_settings={'nst': 2, 'rates': 'gamma', 'threads': 1},
            constraints=sample_constraints
        )
        
        # Mock all analysis engines
        mock_results = {}
        for analysis_type in [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY]:
            mock_engine = Mock()
            mock_result = AnalysisResult(
                analysis_type=analysis_type,
                success=True,
                branches_tested=len(sample_constraints),
                decay_indices={
                    constraint['id']: {
                        'taxa': constraint['taxa'],
                        f'{analysis_type}_metric': 0.5 + hash(analysis_type) % 10 * 0.1
                    }
                    for constraint in sample_constraints
                },
                execution_time=5.0,
                metadata={f'{analysis_type}_specific': True}
            )
            mock_engine.analyze.return_value = mock_result
            mock_engine.validate_inputs.return_value = True
            mock_engine.cleanup.return_value = None
            
            orchestrator._engines[analysis_type] = mock_engine
            mock_results[analysis_type] = mock_result
        
        # Run multi-analysis
        result = orchestrator.run_analysis(
            [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY],
            analysis_data
        )
        
        # Verify results
        assert result.success is True
        assert len(result.analysis_results) == 3
        assert result.metadata['successful_analyses'] == [
            AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY
        ]
        
        # Check that all constraint IDs appear in all analyses
        for analysis_type in [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY]:
            analysis_result = result.analysis_results[analysis_type]
            for constraint in sample_constraints:
                assert constraint['id'] in analysis_result.decay_indices


class TestSystemIntegration:
    """Test integration between different system components."""
    
    def test_format_detection_integration(self, temp_dir):
        """Test integration of format detection with file processing."""
        # Create files in different formats
        fasta_content = ">seq1\nATCGATCG\n>seq2\nATCGATCG\n"
        nexus_content = "#NEXUS\nBEGIN DATA;\nsequences here\nEND;\n"
        
        fasta_file = temp_dir / "test.fas"
        nexus_file = temp_dir / "test.nex"
        
        fasta_file.write_text(fasta_content)
        nexus_file.write_text(nexus_content)
        
        mapper = AlignmentFormatMapper()
        
        # Test format detection integration
        fasta_format = mapper.detect_format("test.fas", fasta_content)
        nexus_format = mapper.detect_format("test.nex", nexus_content)
        
        assert fasta_format == 'fasta'
        assert nexus_format == 'nexus'
        
        # Test fallback behavior
        unknown_format = mapper.detect_format("unknown.xyz", "random content")
        assert unknown_format == 'fasta'  # Should fall back to fasta
    
    def test_configuration_processing_integration(self):
        """Test integration of configuration processing."""
        converter = ConfigurationValueConverter()
        
        # Test processing of realistic configuration
        config_values = {
            'gamma': 'true',
            'nst': '2',
            'threads': '4',
            'gamma_shape': '0.5',
            'data_type': 'DNA',
            'normalization': 'basic'
        }
        
        processed_config = {}
        for param, value in config_values.items():
            processed_config[param] = converter.convert(param, value)
        
        # Verify correct type conversions
        assert processed_config['gamma'] is True
        assert processed_config['nst'] == 2
        assert processed_config['threads'] == 4
        assert processed_config['gamma_shape'] == 0.5
        assert processed_config['data_type'] == 'dna'
        assert processed_config['normalization'] == 'basic'
    
    def test_thread_calculation_integration(self):
        """Test integration of thread calculation with system detection."""
        calculator = ThreadCalculator()
        
        # Test auto calculation
        auto_threads = calculator.calculate_optimal_threads("auto")
        assert auto_threads >= 1
        
        # Test explicit specification
        explicit_threads = calculator.calculate_optimal_threads(2)
        assert explicit_threads == 2
        
        # Test memory estimation integration
        memory_usage = calculator.estimate_memory_usage(auto_threads, 50.0)
        assert memory_usage > 50.0  # Should be more than base data size
        
        # Test validation integration
        is_valid = calculator.validate_thread_memory_combination(auto_threads, 50.0, 2000.0)
        assert is_valid is True
    
    def test_memory_optimization_integration(self, temp_dir):
        """Test integration of memory optimization components."""
        optimizer = MemoryOptimizer()
        
        # Test large alignment optimization recommendations
        recommendations = optimizer.optimize_for_large_alignment(150.0)  # 150 MB
        
        assert recommendations['use_chunked_reading'] is True
        assert recommendations['enable_garbage_collection'] is True
        assert recommendations['reduce_cache_size'] is True
        assert 'suggested_chunk_size' in recommendations
        
        # Test memory estimation integration
        estimates = optimizer.estimate_analysis_memory_requirements(
            alignment_size_mb=100.0,
            num_sequences=50,
            analysis_types=['ml', 'bayesian']
        )
        
        assert estimates['base_memory'] == 100.0
        assert estimates['total_estimated'] > 100.0
        assert estimates['analysis_memory'] > 0
    
    def test_tree_utilities_integration(self, temp_dir, sample_tree_content):
        """Test integration of tree utility functions."""
        # Create a tree file
        tree_file = temp_dir / "integration_tree.tre"
        tree_file.write_text(sample_tree_content)
        
        # Test tree preparation integration
        prepared_tree = prepare_starting_tree(tree_file, temp_dir, "prepared.tre")
        assert prepared_tree is not None
        assert prepared_tree.exists()
        
        # Test tree validation integration
        is_valid = validate_tree_file(prepared_tree, require_branch_lengths=True)
        assert is_valid is True
        
        # Test tree cleaning integration
        from src.utils.tree_utils import clean_tree_file
        messy_tree = temp_dir / "messy.tre"
        messy_content = sample_tree_content + "[metadata]"
        messy_tree.write_text(messy_content)
        
        clean_tree_file(messy_tree, remove_metadata=True)
        cleaned_content = messy_tree.read_text()
        assert "[metadata]" not in cleaned_content


class TestErrorHandlingIntegration:
    """Test integration of error handling across components."""
    
    def test_analysis_error_propagation(self, setup_test_files):
        """Test that analysis errors propagate correctly through the system."""
        files = setup_test_files
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        # Create analysis data with missing alignment
        invalid_data = AnalysisData(
            alignment_file=files['temp_dir'] / "nonexistent.nex",
            tree_file=None,
            temp_dir=files['temp_dir'],
            model_settings={},
            constraints=[],
            timeouts=AnalysisTimeouts()
        )
        
        # Should handle validation error gracefully
        result = orchestrator.run_analysis([AnalysisModes.ML], invalid_data)
        
        assert result.success is False
        assert "failed" in result.error_message.lower()
    
    def test_tool_error_integration(self, setup_test_files, sample_constraints):
        """Test integration of external tool error handling."""
        files = setup_test_files
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        analysis_data = orchestrator.create_analysis_data(
            alignment_file=files['alignment_file'],
            temp_dir=files['temp_dir'],
            model_settings={'nst': 2},
            constraints=sample_constraints
        )
        
        # Mock external tool failure
        with patch('src.analysis.engines.ml_engine.ExternalToolManager') as mock_tool_manager:
            mock_manager = Mock()
            mock_manager.__enter__ = Mock(return_value=mock_manager)
            mock_manager.__exit__ = Mock(return_value=None)
            mock_manager.check_tool_availability.side_effect = Exception("PAUP* not found")
            mock_tool_manager.return_value = mock_manager
            
            orchestrator.create_engines([AnalysisModes.ML])
            result = orchestrator.run_analysis([AnalysisModes.ML], analysis_data)
            
            assert result.success is False
            assert "failed" in result.error_message.lower()
    
    def test_configuration_error_integration(self):
        """Test integration of configuration error handling."""
        converter = ConfigurationValueConverter()
        
        # Test that invalid values raise appropriate errors
        with pytest.raises(ValueError):
            converter.convert('threads', 'invalid_number')
        
        with pytest.raises(ValueError):
            converter.convert('data_type', 'invalid_type')
        
        with pytest.raises(ValueError):
            converter.convert('normalization', 'invalid_level')


class TestPerformanceIntegration:
    """Test performance-related integration scenarios."""
    
    def test_large_dataset_simulation(self, temp_dir):
        """Test system behavior with simulated large datasets."""
        # Create a simulated large alignment
        large_alignment_content = """#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX=20 NCHAR=1000;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
"""
        
        # Add sequences
        for i in range(20):
            sequence = "ATCG" * 250  # 1000 characters
            large_alignment_content += f"    seq{i:02d}  {sequence}\n"
        
        large_alignment_content += "  ;\nEND;\n"
        
        alignment_file = temp_dir / "large_alignment.nex"
        alignment_file.write_text(large_alignment_content)
        
        # Test memory optimization recommendations
        optimizer = MemoryOptimizer()
        file_size_mb = len(large_alignment_content) / (1024 * 1024)
        
        recommendations = optimizer.optimize_for_large_alignment(file_size_mb)
        memory_estimates = optimizer.estimate_analysis_memory_requirements(
            alignment_size_mb=file_size_mb,
            num_sequences=20,
            analysis_types=['ml']
        )
        
        assert 'total_estimated' in memory_estimates
        assert memory_estimates['total_estimated'] > file_size_mb
        
        # Test thread calculation for large dataset
        calculator = ThreadCalculator()
        if hasattr(calculator, 'calculate_for_data_size'):
            from src.utils.thread_calculator import AdaptiveThreadCalculator
            adaptive_calculator = AdaptiveThreadCalculator()
            optimal_threads = adaptive_calculator.calculate_for_data_size(file_size_mb, 'auto')
            assert optimal_threads >= 1
    
    def test_resource_management_integration(self, temp_dir):
        """Test integration of resource management across components."""
        from src.utils.memory_manager import memory_limited_operation
        
        # Test memory-limited operation context manager
        with memory_limited_operation(50.0, "integration_test") as monitor:
            # Simulate some memory-intensive operations
            data_structures = []
            for i in range(100):
                data_structures.append([j for j in range(100)])
            
            # Monitor should track memory usage
            assert len(monitor.snapshots) >= 2
            
            # Get memory summary
            summary = monitor.get_memory_summary()
            assert 'peak_memory_mb' in summary
            assert 'avg_memory_mb' in summary


class TestConfigurationIntegration:
    """Test integration of configuration management."""
    
    def test_analysis_type_configuration_integration(self):
        """Test integration of analysis type configuration."""
        # Test various configurations
        configs = [
            AnalysisTypeConfiguration(has_ml=True),
            AnalysisTypeConfiguration(has_ml=True, has_bayesian=True),
            AnalysisTypeConfiguration(has_ml=True, has_bayesian=True, has_parsimony=True)
        ]
        
        for config in configs:
            # Test layout configuration
            layout = config.get_layout_configuration()
            assert layout is not None
            assert isinstance(layout, str)
            
            # Test analysis types
            types = config.analysis_types
            assert len(types) >= 1
            assert all(t in [AnalysisModes.ML, AnalysisModes.BAYESIAN, AnalysisModes.PARSIMONY] for t in types)
            
            # Test requirements
            if config.has_ml or config.has_bayesian:
                assert config.requires_model_configuration() is True
            else:
                assert config.requires_model_configuration() is False
    
    def test_orchestrator_mode_parsing_integration(self, temp_dir):
        """Test integration of orchestrator mode parsing."""
        orchestrator = AnalysisOrchestrator(temp_dir)
        
        # Test various mode string formats
        test_cases = [
            ("ml", ["ml"]),
            ("ml+bayesian", ["ml", "bayesian"]),
            ("bayesian+parsimony", ["bayesian", "parsimony"]),
            ("all", orchestrator.get_supported_analysis_modes())
        ]
        
        for mode_string, expected_modes in test_cases:
            parsed_modes = orchestrator.parse_analysis_mode_string(mode_string)
            assert set(parsed_modes) == set(expected_modes)


class TestRealWorldScenarios:
    """Test realistic usage scenarios."""
    
    def test_typical_ml_analysis_scenario(self, setup_test_files, sample_constraints):
        """Test a typical ML analysis scenario."""
        files = setup_test_files
        
        # Step 1: Set up orchestrator
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        # Step 2: Process configuration
        converter = ConfigurationValueConverter()
        config = {
            'analysis_type': converter.convert('criterion', 'ml'),
            'model': 'GTR',
            'gamma': converter.convert('gamma', 'true'),
            'threads': converter.convert('threads', '2')
        }
        
        # Step 3: Calculate optimal resources
        calculator = ThreadCalculator()
        optimal_threads = calculator.calculate_optimal_threads(config['threads'])
        
        # Step 4: Create analysis data
        analysis_data = orchestrator.create_analysis_data(
            alignment_file=files['alignment_file'],
            temp_dir=files['temp_dir'],
            model_settings={
                'nst': 6,  # GTR model
                'rates': 'gamma',
                'threads': optimal_threads
            },
            constraints=sample_constraints
        )
        
        # Step 5: Mock and run analysis
        with patch('src.analysis.engines.ml_engine.ExternalToolManager') as mock_tool_manager:
            mock_manager = Mock()
            mock_manager.__enter__ = Mock(return_value=mock_manager)
            mock_manager.__exit__ = Mock(return_value=None)
            mock_manager.check_tool_availability.return_value = None
            mock_manager.execute_script_file.return_value = None
            mock_tool_manager.return_value = mock_manager
            
            orchestrator.create_engines([AnalysisModes.ML])
            ml_engine = orchestrator._engines[AnalysisModes.ML]
            
            with patch.object(ml_engine, '_parse_likelihood_score', return_value=-1500.0):
                with patch.object(ml_engine, '_test_constraints', return_value={}):
                    with patch.object(ml_engine, '_calculate_decay_indices', return_value={}):
                        result = orchestrator.run_analysis([AnalysisModes.ML], analysis_data)
                        
                        assert result.success is True
                        assert result.total_execution_time > 0
    
    def test_multi_constraint_analysis_scenario(self, setup_test_files):
        """Test analysis with multiple complex constraints."""
        files = setup_test_files
        
        # Create complex constraint set
        complex_constraints = [
            {'id': 'mammals', 'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon', 'Macaque']},
            {'id': 'primates', 'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan', 'Gibbon']},
            {'id': 'great_apes', 'taxa': ['Human', 'Chimp', 'Gorilla', 'Orangutan']},
            {'id': 'african_apes', 'taxa': ['Human', 'Chimp', 'Gorilla']},
            {'id': 'hominini', 'taxa': ['Human', 'Chimp']},
        ]
        
        orchestrator = AnalysisOrchestrator(files['temp_dir'], debug=True)
        
        analysis_data = orchestrator.create_analysis_data(
            alignment_file=files['alignment_file'],
            temp_dir=files['temp_dir'],
            model_settings={'nst': 2, 'rates': 'gamma'},
            constraints=complex_constraints
        )
        
        # Verify constraint processing
        assert len(analysis_data.constraints) == 5
        assert all('taxa' in constraint for constraint in analysis_data.constraints)
        assert all('id' in constraint for constraint in analysis_data.constraints)
        
        # Test validation passes
        try:
            orchestrator.validate_analysis_request([AnalysisModes.ML], analysis_data)
        except ValidationError:
            pytest.fail("Validation should pass for valid constraints")


if __name__ == "__main__":
    pytest.main([__file__])