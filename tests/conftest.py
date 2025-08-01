"""
Pytest configuration and shared fixtures.

This module contains shared test fixtures and configuration
for the panDecay test suite.
"""

import pytest
import tempfile
import shutil
import logging
from pathlib import Path
from unittest.mock import Mock
from typing import Dict, Any

# Configure logging for tests
logging.basicConfig(level=logging.DEBUG)


@pytest.fixture(scope="session")
def test_data_dir():
    """Get the test data directory."""
    return Path(__file__).parent / "data"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for each test."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def sample_sequences():
    """Provide sample DNA sequences for testing."""
    return {
        'seq1': 'ATCGATCGATCGATCG',
        'seq2': 'ATCGATCGATCGATCG', 
        'seq3': 'ATCGATCGATCGATCC',
        'seq4': 'ATCGATCGATCGATAA',
        'seq5': 'ATCGATCGATCGTTTT',
        'seq6': 'ATCGATCGATCGAAAA'
    }


@pytest.fixture
def sample_nexus_alignment(sample_sequences):
    """Create a sample NEXUS alignment."""
    ntax = len(sample_sequences)
    nchar = len(next(iter(sample_sequences.values())))
    
    content = f"""#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX={ntax} NCHAR={nchar};
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
"""
    
    for seq_name, sequence in sample_sequences.items():
        content += f"    {seq_name:<10} {sequence}\n"
    
    content += "  ;\nEND;\n"
    return content


@pytest.fixture
def sample_fasta_alignment(sample_sequences):
    """Create a sample FASTA alignment."""
    content = ""
    for seq_name, sequence in sample_sequences.items():
        content += f">{seq_name}\n{sequence}\n"
    return content


@pytest.fixture
def sample_tree():
    """Provide a sample phylogenetic tree."""
    return "(seq1:0.1,seq2:0.1,(seq3:0.15,(seq4:0.2,(seq5:0.25,seq6:0.3):0.1):0.05):0.05);"


@pytest.fixture
def sample_constraints():
    """Provide sample constraint definitions."""
    return [
        {'id': 'group1', 'taxa': ['seq1', 'seq2', 'seq3']},
        {'id': 'group2', 'taxa': ['seq4', 'seq5', 'seq6']},
        {'id': 'pair1', 'taxa': ['seq1', 'seq2']},
        {'id': 'pair2', 'taxa': ['seq5', 'seq6']}
    ]


@pytest.fixture
def mock_external_tool():
    """Create a mock external tool for testing."""
    mock_tool = Mock()
    mock_tool.check_availability.return_value = True
    mock_tool.execute.return_value = (0, "Success", "")
    mock_tool.execute_script.return_value = (0, "Success", "")
    return mock_tool


@pytest.fixture
def sample_model_settings():
    """Provide sample model settings."""
    return {
        'nst': 2,
        'rates': 'gamma',
        'gamma_shape': 0.5,
        'prop_invar': 0.25,
        'threads': 1
    }


@pytest.fixture
def sample_analysis_config():
    """Provide sample analysis configuration."""
    return {
        'analysis_type': 'ml',
        'model': 'GTR',
        'gamma': True,
        'invariable': True,
        'threads': 2,
        'bootstrap': False,
        'bootstrap_reps': 1000
    }


@pytest.fixture(autouse=True)
def reset_singletons():
    """Reset any singleton instances between tests."""
    # This fixture automatically runs before each test
    # Add any singleton reset logic here if needed
    yield
    # Cleanup after test


@pytest.fixture
def caplog_debug(caplog):
    """Capture debug logs during tests."""
    with caplog.at_level(logging.DEBUG):
        yield caplog


def create_test_file(temp_dir: Path, filename: str, content: str) -> Path:
    """Helper function to create test files."""
    file_path = temp_dir / filename
    file_path.write_text(content)
    return file_path


def pytest_configure(config):
    """Configure pytest with custom markers and settings."""
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "external: mark test as requiring external tools"
    )
    config.addinivalue_line(
        "markers", "memory: mark test as memory-related"
    )
    config.addinivalue_line(
        "markers", "performance: mark test as performance-related"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers based on test names."""
    for item in items:
        # Add unit marker to test_* files
        if "test_" in item.nodeid and not "integration" in item.nodeid:
            item.add_marker(pytest.mark.unit)
        
        # Add integration marker to integration tests
        if "integration" in item.nodeid:
            item.add_marker(pytest.mark.integration)
        
        # Add slow marker to tests that might be slow
        if any(keyword in item.name.lower() for keyword in ['large', 'performance', 'benchmark']):
            item.add_marker(pytest.mark.slow)
        
        # Add external marker to tests requiring external tools
        if any(keyword in item.name.lower() for keyword in ['paup', 'mrbayes', 'external']):
            item.add_marker(pytest.mark.external)
        
        # Add memory marker to memory-related tests
        if any(keyword in item.name.lower() for keyword in ['memory', 'cache', 'optimization']):
            item.add_marker(pytest.mark.memory)


@pytest.fixture
def mock_analysis_result():
    """Create a mock analysis result for testing."""
    return {
        'analysis_type': "test",
        'success': True,
        'branches_tested': 3,
        'decay_indices': {
            'group1': {'taxa': ['seq1', 'seq2', 'seq3'], 'decay_value': 0.5},
            'group2': {'taxa': ['seq4', 'seq5', 'seq6'], 'decay_value': 0.8}
        },
        'execution_time': 5.0,
        'metadata': {'test_run': True}
    }


class TestHelpers:
    """Helper class with utility methods for tests."""
    
    @staticmethod
    def create_alignment_file(temp_dir: Path, format_type: str = "nexus", 
                            sequences: Dict[str, str] = None) -> Path:
        """Create an alignment file in the specified format."""
        if sequences is None:
            sequences = {
                'seq1': 'ATCGATCG',
                'seq2': 'ATCGATCG',
                'seq3': 'ATCGATCC'
            }
        
        if format_type.lower() == "nexus":
            content = f"""#NEXUS

BEGIN DATA;
  DIMENSIONS NTAX={len(sequences)} NCHAR={len(next(iter(sequences.values())))};
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
"""
            for name, seq in sequences.items():
                content += f"    {name:<10} {seq}\n"
            content += "  ;\nEND;\n"
            filename = "test_alignment.nex"
            
        elif format_type.lower() == "fasta":
            content = ""
            for name, seq in sequences.items():
                content += f">{name}\n{seq}\n"
            filename = "test_alignment.fas"
            
        else:
            raise ValueError(f"Unsupported format: {format_type}")
        
        file_path = temp_dir / filename
        file_path.write_text(content)
        return file_path
    
    @staticmethod
    def create_tree_file(temp_dir: Path, tree_string: str = None) -> Path:
        """Create a tree file."""
        if tree_string is None:
            tree_string = "(seq1:0.1,seq2:0.1,seq3:0.15);"
        
        tree_file = temp_dir / "test_tree.tre"
        tree_file.write_text(tree_string)
        return tree_file


@pytest.fixture
def test_helpers():
    """Provide test helper utilities."""
    return TestHelpers