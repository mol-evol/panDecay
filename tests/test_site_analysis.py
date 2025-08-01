"""
Test suite for site analysis functionality.

Tests the integration of SiteAnalyzer, MLAnalyzer, ResultProcessor, and AnalysisCoordinator
for site-specific likelihood analysis.
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import Mock, MagicMock

# Import the components
from src.core.analysis.site_analyzer import SiteAnalyzer, SiteData, CladesSiteAnalysis
from src.core.analysis.ml_analyzer import MLAnalyzer
from src.core.analysis.result_processor import ResultProcessor
from src.core.analysis.file_manager import FileManager
from src.core.analysis.external_tools import ExternalTools
from src.core.analysis.data_processor import DataProcessor


class TestSiteAnalyzer:
    """Test the core SiteAnalyzer component."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create mock dependencies
        self.file_manager = Mock(spec=FileManager)
        self.file_manager.get_temp_path.return_value = self.temp_dir / "test.txt"
        self.file_manager.write_file.return_value = self.temp_dir / "test.nex"
        
        self.external_tools = Mock(spec=ExternalTools)
        self.external_tools.run_paup.return_value = Mock(returncode=0, stderr="")
        
        self.data_processor = Mock(spec=DataProcessor)
        self.data_processor.is_alignment_loaded.return_value = True
        self.data_processor.get_taxon_names.return_value = ["taxon1", "taxon2", "taxon3"]
    
    def test_site_analyzer_initialization(self):
        """Test SiteAnalyzer initialization."""
        analyzer = SiteAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            paup_model_commands="lset nst=6 rates=gamma;",
            debug=True
        )
        
        assert analyzer is not None
        assert analyzer.debug is True
        assert analyzer.site_likelihood_threshold == 0.1  # Default
        assert len(analyzer.site_analyses) == 0
    
    def test_site_data_structure(self):
        """Test SiteData dataclass."""
        site_data = SiteData(
            site_position=1,
            ml_likelihood=-12.34,
            constraint_likelihood=-15.67,
            delta_lnl=3.33,
            site_type="supporting"
        )
        
        assert site_data.site_position == 1
        assert site_data.ml_likelihood == -12.34
        assert site_data.constraint_likelihood == -15.67
        assert site_data.delta_lnl == 3.33
        assert site_data.site_type == "supporting"
    
    def test_clades_site_analysis_structure(self):
        """Test CladesSiteAnalysis dataclass."""
        site_data = {
            1: SiteData(1, -10.0, -12.0, 2.0, "supporting"),
            2: SiteData(2, -8.0, -7.0, -1.0, "conflicting")
        }
        
        analysis = CladesSiteAnalysis(
            clade_id="Clade_1",
            site_data=site_data,
            supporting_sites=1,
            conflicting_sites=1,
            neutral_sites=0,
            support_ratio=0.5,
            sum_supporting_delta=2.0,
            sum_conflicting_delta=-1.0,
            weighted_support_ratio=0.67
        )
        
        assert analysis.clade_id == "Clade_1"
        assert len(analysis.site_data) == 2
        assert analysis.supporting_sites == 1
        assert analysis.conflicting_sites == 1


class TestMLAnalyzerSiteIntegration:
    """Test MLAnalyzer site analysis integration."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create mock dependencies
        self.file_manager = Mock(spec=FileManager)
        self.external_tools = Mock(spec=ExternalTools)
        self.data_processor = Mock(spec=DataProcessor)
        self.data_processor.is_alignment_loaded.return_value = True
        
        self.ml_analyzer = MLAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            model_commands="lset nst=6;",
            debug=True
        )
    
    def test_site_analysis_enablement(self):
        """Test enabling site analysis in MLAnalyzer."""
        assert not self.ml_analyzer.has_site_analysis_enabled()
        
        self.ml_analyzer.enable_site_analysis()
        
        assert self.ml_analyzer.has_site_analysis_enabled()
        assert self.ml_analyzer.site_analyzer is not None
    
    def test_site_analysis_summary(self):
        """Test site analysis summary generation."""
        # Enable site analysis
        self.ml_analyzer.enable_site_analysis()
        
        # Create mock site analysis results
        mock_site_analysis = Mock()
        mock_site_analysis.site_data = {1: Mock(), 2: Mock(), 3: Mock()}
        mock_site_analysis.supporting_sites = 2
        mock_site_analysis.conflicting_sites = 1
        mock_site_analysis.neutral_sites = 0
        
        self.ml_analyzer.site_analysis_results["Clade_1"] = mock_site_analysis
        
        summary = self.ml_analyzer.get_site_analysis_summary()
        
        assert summary['site_analysis_enabled'] is True
        assert summary['clades_analyzed'] == 1
        assert summary['total_sites_analyzed'] == 3
        assert summary['total_supporting_sites'] == 2
        assert summary['total_conflicting_sites'] == 1


class TestResultProcessorSiteIntegration:
    """Test ResultProcessor site analysis integration."""
    
    def setup_method(self):
        """Set up test fixtures."""
        self.file_manager = Mock(spec=FileManager)
        self.external_tools = Mock(spec=ExternalTools)
        self.data_processor = Mock(spec=DataProcessor)
        
        self.result_processor = ResultProcessor(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            debug=True
        )
    
    def test_site_analysis_integration(self):
        """Test integration of site analysis data."""
        # Create mock site analysis data
        mock_site_analysis = Mock()
        mock_site_analysis.site_data = {1: Mock(), 2: Mock()}
        mock_site_analysis.supporting_sites = 1
        mock_site_analysis.conflicting_sites = 1
        mock_site_analysis.neutral_sites = 0
        mock_site_analysis.support_ratio = 0.5
        mock_site_analysis.weighted_support_ratio = 0.6
        mock_site_analysis.sum_supporting_delta = 2.0
        mock_site_analysis.sum_conflicting_delta = -1.0
        
        site_data = {"Clade_1": mock_site_analysis}
        
        # Create a decay index entry first
        self.result_processor.decay_indices["Clade_1"] = {
            'clade_id': 'Clade_1',
            'taxa': ['taxon1', 'taxon2'],
            'analysis_types': ['ml']
        }
        
        # Add site analysis results
        self.result_processor.add_site_analysis_results(site_data)
        
        # Check integration
        decay_entry = self.result_processor.decay_indices["Clade_1"]
        assert 'supporting_sites' in decay_entry
        assert 'conflicting_sites' in decay_entry
        assert 'support_ratio' in decay_entry
        assert decay_entry['supporting_sites'] == 1
        assert decay_entry['conflicting_sites'] == 1


class TestSiteAnalysisWorkflow:
    """Test the complete site analysis workflow."""
    
    def test_workflow_integration(self):
        """Test that all components work together for site analysis."""
        # This is a high-level integration test
        # In practice, this would use real test data and verify the complete pipeline
        
        # Create minimal test setup
        file_manager = Mock()
        external_tools = Mock()
        data_processor = Mock()
        data_processor.is_alignment_loaded.return_value = True
        
        # Create site analyzer
        site_analyzer = SiteAnalyzer(
            file_manager=file_manager,
            external_tools=external_tools,
            data_processor=data_processor,
            paup_model_commands="lset nst=6;"
        )
        
        # Create ML analyzer with site analysis
        ml_analyzer = MLAnalyzer(
            file_manager=file_manager,
            external_tools=external_tools,
            data_processor=data_processor,
            model_commands="lset nst=6;"
        )
        
        # Create result processor
        result_processor = ResultProcessor(
            file_manager=file_manager,
            external_tools=external_tools,
            data_processor=data_processor
        )
        
        # Test basic initialization and compatibility
        assert site_analyzer is not None
        assert ml_analyzer is not None
        assert result_processor is not None
        
        # Test site analysis can be enabled
        ml_analyzer.enable_site_analysis()
        assert ml_analyzer.has_site_analysis_enabled()


if __name__ == "__main__":
    # Run basic smoke tests
    print("Running site analysis tests...")
    
    # Test basic imports
    try:
        from src.core.analysis.site_analyzer import SiteAnalyzer, SiteData, CladesSiteAnalysis
        print("✓ SiteAnalyzer imports successful")
    except ImportError as e:
        print(f"✗ SiteAnalyzer import failed: {e}")
    
    # Test data structures
    try:
        site_data = SiteData(1, -10.0, -12.0, 2.0, "supporting")
        print("✓ SiteData structure works")
    except Exception as e:
        print(f"✗ SiteData failed: {e}")
    
    # Test CladesSiteAnalysis structure
    try:
        analysis = CladesSiteAnalysis(
            clade_id="test",
            site_data={},
            supporting_sites=0,
            conflicting_sites=0,
            neutral_sites=0,
            support_ratio=0.0,
            sum_supporting_delta=0.0,
            sum_conflicting_delta=0.0,
            weighted_support_ratio=0.0
        )
        print("✓ CladesSiteAnalysis structure works")
    except Exception as e:
        print(f"✗ CladesSiteAnalysis failed: {e}")
    
    print("Site analysis test suite completed.")