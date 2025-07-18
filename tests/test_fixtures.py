#!/usr/bin/env python3
"""
Test fixtures and mock data for panDecay test suite.

This module provides reusable test fixtures, mock data, and helper functions
for testing panDecay functionality.
"""

import os
import tempfile
import shutil
from pathlib import Path
from unittest.mock import MagicMock
import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_path = Path(tempfile.mkdtemp())
    yield temp_path
    # Cleanup
    if temp_path.exists():
        shutil.rmtree(temp_path)


@pytest.fixture
def simple_alignment():
    """Create a simple test alignment file."""
    alignment_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
        f.write(alignment_content)
        temp_file = f.name
    
    yield temp_file
    
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)


@pytest.fixture
def complex_alignment():
    """Create a more complex test alignment with variations."""
    alignment_content = """>Homo_sapiens
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Pan_troglodytes
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Gorilla_gorilla
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Pongo_abelii
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Macaca_mulatta
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Papio_anubis
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Chlorocebus_sabaeus
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Saimiri_boliviensis
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fas', delete=False) as f:
        f.write(alignment_content)
        temp_file = f.name
    
    yield temp_file
    
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)


@pytest.fixture
def test_config():
    """Create a test configuration file."""
    config_content = """
[analysis]
analysis_type = ml
model = GTR
gamma = true
threads = 4

[output]
output_dir = test_output
keep_files = false
debug = false
output_style = ascii

[mrbayes]
ngen = 1000000
chains = 4
burnin = 0.25
sample_freq = 1000

[constraints]
constraint1 = seq1,seq2
constraint2 = seq3,seq4

[paup]
ml_timeout = 3600
constraint_timeout = 600
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.ini', delete=False) as f:
        f.write(config_content)
        temp_file = f.name
    
    yield temp_file
    
    # Cleanup
    if os.path.exists(temp_file):
        os.unlink(temp_file)


@pytest.fixture
def mock_paup_output():
    """Mock PAUP* output data."""
    return {
        'ml_score': "-1234.567890",
        'constraint_score': "-1244.567890",
        'au_results': {
            1: {'lnL': -1234.567890, 'AU_pvalue': 0.95},
            2: {'lnL': -1244.567890, 'AU_pvalue': 0.02}
        },
        'site_analysis': {
            'site_1': {'ml_lnl': -2.5, 'constraint_lnl': -3.0},
            'site_2': {'ml_lnl': -1.8, 'constraint_lnl': -2.2}
        }
    }


@pytest.fixture
def mock_mrbayes_output():
    """Mock MrBayes output data."""
    return {
        'unconstrained_ml': -1230.456789,
        'constrained_ml': -1240.456789,
        'convergence': {
            'min_ess': 250.0,
            'max_psrf': 1.008,
            'asdsf': 0.005,
            'converged': True,
            'warnings': []
        }
    }


@pytest.fixture
def mock_tree_string():
    """Mock Newick tree string."""
    return "(((Homo_sapiens:0.01,Pan_troglodytes:0.01):0.02,Gorilla_gorilla:0.03):0.04,Pongo_abelii:0.05);"


@pytest.fixture
def mock_decay_indices():
    """Mock decay indices data structure."""
    return {
        'Clade_1': {
            'taxa': ['Homo_sapiens', 'Pan_troglodytes'],
            'lnl_diff': 10.0,
            'AU_pvalue': 0.02,
            'significant_AU': True,
            'constrained_lnl': -1244.567890
        },
        'Clade_2': {
            'taxa': ['Gorilla_gorilla', 'Pongo_abelii'],
            'lnl_diff': 5.0,
            'AU_pvalue': 0.15,
            'significant_AU': False,
            'constrained_lnl': -1239.567890
        }
    }


@pytest.fixture
def mock_site_data():
    """Mock site-specific analysis data."""
    return {
        'signal_std': 2.5,
        'signal_mad': 1.8,
        'supporting_sites': 85,
        'conflicting_sites': 15,
        'neutral_sites': 0,
        'support_ratio': 0.85,
        'sum_supporting_delta': 42.5,
        'sum_conflicting_delta': -7.5,
        'weighted_support_ratio': 0.92
    }


@pytest.fixture
def mock_external_software():
    """Mock external software dependencies."""
    mock_paup = MagicMock()
    mock_paup.return_value = MagicMock()
    mock_paup.return_value.returncode = 0
    
    mock_mrbayes = MagicMock()
    mock_mrbayes.return_value = MagicMock()
    mock_mrbayes.return_value.returncode = 0
    
    return {
        'paup': mock_paup,
        'mrbayes': mock_mrbayes
    }


class MockPAUPRunner:
    """Mock PAUP* runner for testing."""
    
    def __init__(self, temp_path, success=True):
        self.temp_path = temp_path
        self.success = success
    
    def run_paup_command_file(self, cmd_file, log_file, timeout=None):
        """Mock PAUP* command execution."""
        if not self.success:
            raise RuntimeError("Mock PAUP* execution failed")
        
        # Create mock output files
        log_path = self.temp_path / log_file
        log_path.write_text("Mock PAUP* output")
        
        return True


class MockMrBayesRunner:
    """Mock MrBayes runner for testing."""
    
    def __init__(self, temp_path, success=True):
        self.temp_path = temp_path
        self.success = success
    
    def run_mrbayes(self, nexus_file, output_prefix):
        """Mock MrBayes execution."""
        if not self.success:
            return None
        
        # Create mock output files
        lstat_file = self.temp_path / f"{output_prefix}.nex.lstat"
        lstat_file.write_text("all\t-1230.456789\t-1230.456789\t0\n")
        
        return -1230.456789


@pytest.fixture
def mock_subprocess_success():
    """Mock successful subprocess execution."""
    mock_result = MagicMock()
    mock_result.returncode = 0
    mock_result.stdout = "Mock successful output"
    mock_result.stderr = ""
    return mock_result


@pytest.fixture
def mock_subprocess_failure():
    """Mock failed subprocess execution."""
    mock_result = MagicMock()
    mock_result.returncode = 1
    mock_result.stdout = ""
    mock_result.stderr = "Mock error output"
    return mock_result


# Test data constants
TEST_ALIGNMENT_DNA = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq3
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""

TEST_ALIGNMENT_PROTEIN = """>seq1
MKTVRQERLKSIVRILYFIKQFRCDATVNKKQSKAQLVVPNGTKGCVSVSDDVYSDMVRKQMKAGVDDDTDVDVDVDVDVDVDVDVDVDVDVDVDVDVDVD
>seq2
MKTVRQERLKSIVRILYFIKQFRCDATVNKKQSKAQLVVPNGTKGCVSVSDDVYSDMVRKQMKAGVDDDTDVDVDVDVDVDVDVDVDVDVDVDVDVDVDVD
>seq3
MKTVRQERLKSIVRILYFIKQFRCDATVNKKQSKAQLVVPNGTKGCVSVSDDVYSDMVRKQMKAGVDDDTDVDVDVDVDVDVDVDVDVDVDVDVDVDVDVD
>seq4
MKTVRQERLKSIVRILYFIKQFRCDATVNKKQSKAQLVVPNGTKGCVSVSDDVYSDMVRKQMKAGVDDDTDVDVDVDVDVDVDVDVDVDVDVDVDVDVDVD
"""

TEST_ALIGNMENT_MORPHOLOGY = """>seq1
01001101010110101001011010101010
>seq2
01001101010110101001011010101010
>seq3
01001101010110101001011010101010
>seq4
01001101010110101001011010101010
"""

TEST_TREE_NEWICK = "(((seq1:0.01,seq2:0.01):0.02,seq3:0.03):0.04,seq4:0.05);"

TEST_CONSTRAINTS = [
    ["seq1", "seq2"],
    ["seq3", "seq4"],
    ["seq1", "seq2", "seq3"]
]