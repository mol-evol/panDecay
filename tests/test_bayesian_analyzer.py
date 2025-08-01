#!/usr/bin/env python3
"""
Test Suite for BayesianAnalyzer Component

Tests the extracted Bayesian analysis functionality to ensure
it behaves identically to the original implementation.
"""

import unittest
import tempfile
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock, mock_open

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

from core.analysis import FileManager, ExternalTools, DataProcessor
from core.analysis.bayesian_analyzer import (
    BayesianAnalyzer, BayesianAnalysisError, MrBayesExecutionError, 
    ConvergenceError, MarginalLikelihoodError
)


class TestBayesianAnalyzer(unittest.TestCase):
    """Test BayesianAnalyzer component functionality"""
    
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
        
        # Initialize components
        self.file_manager = FileManager(temp_dir=self.temp_dir, debug=True)
        self.external_tools = ExternalTools(debug=True)
        self.data_processor = DataProcessor(data_type="dna", debug=True)
        
        # Load alignment
        self.data_processor.load_alignment(self.alignment_file, "fasta")
        
        # Create NEXUS file
        self.nexus_file = self.file_manager.get_temp_path("alignment.nex")
        self.data_processor.convert_to_nexus(self.nexus_file)
    
    def tearDown(self):
        """Clean up test fixtures"""
        self.file_manager.cleanup_all()
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test BayesianAnalyzer initialization"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="GTR",
            debug=True
        )
        
        self.assertEqual(ba.bayes_model, "GTR")
        self.assertTrue(ba.debug)
        self.assertIsNone(ba.consensus_tree)
        self.assertIsNone(ba.marginal_lnl)
        self.assertEqual(ba.constraint_analyses, {})
        self.assertTrue(ba.check_convergence)
    
    def test_initialization_no_alignment(self):
        """Test initialization fails without loaded alignment"""
        empty_processor = DataProcessor()
        
        with self.assertRaises(BayesianAnalysisError) as context:
            BayesianAnalyzer(
                file_manager=self.file_manager,
                external_tools=self.external_tools,
                data_processor=empty_processor
            )
        
        self.assertIn("must have alignment loaded", str(context.exception))
    
    def test_get_dna_model_commands_gtr(self):
        """Test DNA model commands for GTR"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="GTR"
        )
        
        commands = ba._get_dna_model_commands()
        
        self.assertIn("  lset nst=6 rates=gamma;", commands)
        self.assertIn("  prset revmatpr=dirichlet(1,1,1,1,1,1);", commands)
        self.assertIn("  prset statefreqpr=dirichlet(1,1,1,1);", commands)
        self.assertIn("  prset shapepr=exponential(2.0);", commands)
    
    def test_get_dna_model_commands_hky(self):
        """Test DNA model commands for HKY"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="HKY"
        )
        
        commands = ba._get_dna_model_commands()
        
        self.assertIn("  lset nst=2 rates=gamma;", commands)
        self.assertIn("  prset revmatpr=dirichlet(1,1);", commands)
    
    def test_get_dna_model_commands_jc(self):
        """Test DNA model commands for JC"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="JC"
        )
        
        commands = ba._get_dna_model_commands()
        
        self.assertIn("  lset nst=1 rates=equal;", commands)
    
    def test_get_protein_model_commands(self):
        """Test protein model commands"""
        # Create protein data processor
        protein_processor = DataProcessor(data_type="protein", debug=True)
        protein_content = """>seq1
ACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWY"""
        
        protein_file = self.temp_dir / "protein.fasta"
        with open(protein_file, 'w') as f:
            f.write(protein_content)
        
        protein_processor.load_alignment(protein_file, "fasta")
        
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=protein_processor,
            bayes_model="WAG"
        )
        
        commands = ba._get_protein_model_commands()
        
        self.assertIn("  prset aamodelpr=fixed(wag);", commands)
        self.assertIn("  lset rates=gamma;", commands)
        self.assertIn("  prset shapepr=exponential(2.0);", commands)
    
    def test_get_constraint_commands(self):
        """Test constraint command generation"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        commands = ba._get_constraint_commands(["seq1", "seq2"], "test_constraint")
        
        self.assertIn("  constraint test_constraint = (seq1 seq2);", commands)
        self.assertIn("  prset topologypr=constraints(test_constraint);", commands)
    
    def test_get_mcmc_commands(self):
        """Test MCMC parameter commands"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_ngen=100000,
            bayes_chains=4,
            bayes_sample_freq=100,
            bayes_burnin=0.25
        )
        
        commands = ba._get_mcmc_commands("test_output")
        
        self.assertIn("  mcmcp ngen=100000 nchains=4 nruns=2;", commands)
        self.assertIn("  mcmcp samplefreq=100 printfreq=1000;", commands)
        self.assertIn("  mcmcp burnin=25000;", commands)  # 25% of 100000
        self.assertIn("  mcmcp filename=test_output;", commands)
    
    def test_get_stepping_stone_commands(self):
        """Test stepping-stone sampling commands"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            ss_alpha=0.4,
            ss_nsteps=50
        )
        
        commands = ba._get_stepping_stone_commands()
        
        self.assertIn("  ss alpha=0.4 nsteps=50;", commands)
    
    def test_get_beagle_commands(self):
        """Test BEAGLE acceleration commands"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            use_beagle=True,
            beagle_device="auto",
            beagle_precision="double",
            beagle_scaling="dynamic"
        )
        
        commands = ba._get_beagle_commands()
        
        self.assertIn("  set usebeagle=yes;", commands)
        self.assertIn("  set beagledevice=auto;", commands)
        self.assertIn("  set beagleprecision=double;", commands)
        self.assertIn("  set beaglescaling=dynamic;", commands)
    
    def test_generate_mrbayes_nexus_basic(self):
        """Test basic MrBayes NEXUS generation"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="GTR"
        )
        
        nexus_file = ba._generate_mrbayes_nexus(output_prefix="test")
        
        self.assertTrue(nexus_file.exists())
        content = nexus_file.read_text()
        
        self.assertIn("#NEXUS", content)
        self.assertIn("begin mrbayes;", content)
        self.assertIn("execute alignment.nex;", content)
        self.assertIn("lset nst=6 rates=gamma;", content)
        self.assertIn("mcmc;", content)
        self.assertIn("sump;", content)
        self.assertIn("sumt;", content)
        self.assertIn("end;", content)
    
    def test_generate_mrbayes_nexus_with_constraint(self):
        """Test MrBayes NEXUS generation with constraints"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        nexus_file = ba._generate_mrbayes_nexus(
            clade_taxa=["seq1", "seq2"],
            constraint_id="test_constraint",
            output_prefix="constrained"
        )
        
        content = nexus_file.read_text()
        
        self.assertIn("constraint test_constraint = (seq1 seq2);", content)
        self.assertIn("prset topologypr=constraints(test_constraint);", content)
    
    def test_generate_mrbayes_nexus_with_beagle(self):
        """Test MrBayes NEXUS generation with BEAGLE"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            use_beagle=True
        )
        
        nexus_file = ba._generate_mrbayes_nexus(output_prefix="beagle_test")
        content = nexus_file.read_text()
        
        self.assertIn("set usebeagle=yes;", content)
    
    def test_generate_mrbayes_nexus_stepping_stone(self):
        """Test MrBayes NEXUS generation with stepping-stone"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            marginal_likelihood="ss"
        )
        
        nexus_file = ba._generate_mrbayes_nexus(output_prefix="ss_test")
        content = nexus_file.read_text()
        
        self.assertIn("ss alpha=", content)
    
    @patch('Bio.Phylo.read')
    def test_parse_consensus_tree_success(self, mock_phylo_read):
        """Test successful consensus tree parsing"""
        # Mock tree file and Phylo.read
        mock_tree = MagicMock()
        mock_phylo_read.return_value = mock_tree
        
        # Create mock consensus tree file
        con_file = self.file_manager.get_temp_path("test.con.tre")
        with open(con_file, 'w') as f:
            f.write("(seq1:0.1,seq2:0.1,(seq3:0.1,seq4:0.1):0.1);")
        
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        tree = ba._parse_consensus_tree("test")
        
        self.assertEqual(tree, mock_tree)
        mock_phylo_read.assert_called_once()
    
    def test_parse_consensus_tree_missing_file(self):
        """Test consensus tree parsing with missing file"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        tree = ba._parse_consensus_tree("nonexistent")
        
        self.assertIsNone(tree)
    
    def test_parse_stepping_stone_likelihood(self):
        """Test stepping-stone likelihood parsing"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Create mock .ss file
        ss_content = """
MrBayes output...
Stepping-stone estimate of marginal likelihood = -1234.567
More output...
"""
        ss_file = self.file_manager.get_temp_path("test.ss")
        with open(ss_file, 'w') as f:
            f.write(ss_content)
        
        nexus_file = self.file_manager.get_temp_path("test.nex")
        
        lnl = ba._parse_stepping_stone_likelihood(nexus_file, "test")
        
        self.assertEqual(lnl, -1234.567)
    
    def test_parse_stepping_stone_likelihood_missing(self):
        """Test stepping-stone likelihood parsing with missing file"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        nexus_file = self.file_manager.get_temp_path("test.nex")
        lnl = ba._parse_stepping_stone_likelihood(nexus_file, "nonexistent")
        
        self.assertIsNone(lnl)
    
    def test_parse_harmonic_mean_likelihood(self):
        """Test harmonic mean likelihood parsing"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Create mock .lstat file
        lstat_content = """
Gen	LnL	LnPr	
1000	-2345.67	-12.34	
2000	-2344.89	-11.98	
...
Harmonic mean of likelihoods = -2350.123
"""
        lstat_file = self.file_manager.get_temp_path("test.run1.lstat")
        with open(lstat_file, 'w') as f:
            f.write(lstat_content)
        
        lnl = ba._parse_harmonic_mean_likelihood("test")
        
        self.assertEqual(lnl, -2350.123)
    
    def test_parse_harmonic_mean_likelihood_missing(self):
        """Test harmonic mean likelihood parsing with missing file"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        lnl = ba._parse_harmonic_mean_likelihood("nonexistent")
        
        self.assertIsNone(lnl)
    
    def test_parse_convergence_diagnostics(self):
        """Test convergence diagnostics parsing"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Create mock .pstat file
        pstat_content = """
# Parameter statistics
Gen	LnL	LnPr	TL	pi(A)	ESS	PSRF
1000	-2345.67	-12.34	0.5	0.25	150.5	1.002
2000	-2344.89	-11.98	0.52	0.24	200.1	1.001
"""
        pstat_file = self.file_manager.get_temp_path("test.run1.pstat")
        with open(pstat_file, 'w') as f:
            f.write(pstat_content)
        
        # Create mock .mcmc file
        mcmc_content = """
MrBayes output...
Generation 1000: Average standard deviation of split frequencies: 0.012345
Generation 2000: Average standard deviation of split frequencies: 0.005678
"""
        mcmc_file = self.file_manager.get_temp_path("test.mcmc")
        with open(mcmc_file, 'w') as f:
            f.write(mcmc_content)
        
        nexus_file = self.file_manager.get_temp_path("test.nex")
        
        convergence_data = ba._parse_convergence_diagnostics(nexus_file, "test")
        
        self.assertEqual(convergence_data['min_ess'], 150.5)
        self.assertEqual(convergence_data['max_psrf'], 1.002)
        self.assertEqual(convergence_data['asdsf'], 0.005678)
    
    def test_check_convergence_success(self):
        """Test convergence checking with good values"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            min_ess=200,
            max_psrf=1.01,
            max_asdsf=0.01
        )
        
        ba.convergence_data = {
            'min_ess': 250.0,
            'max_psrf': 1.005,
            'asdsf': 0.005
        }
        
        # Should not raise exception
        ba._check_convergence("test")
        
        self.assertTrue(ba.convergence_data['converged'])
    
    def test_check_convergence_failure_strict(self):
        """Test convergence checking failure in strict mode"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            min_ess=200,
            max_psrf=1.01,
            max_asdsf=0.01,
            convergence_strict=True
        )
        
        ba.convergence_data = {
            'min_ess': 150.0,  # Too low
            'max_psrf': 1.005,
            'asdsf': 0.005
        }
        
        with self.assertRaises(ConvergenceError):
            ba._check_convergence("test")
    
    def test_check_convergence_failure_lenient(self):
        """Test convergence checking failure in lenient mode"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            min_ess=200,
            convergence_strict=False
        )
        
        ba.convergence_data = {
            'min_ess': 150.0,  # Too low
            'max_psrf': 1.005,
            'asdsf': 0.005
        }
        
        # Should not raise exception but log warning
        ba._check_convergence("test")
        
        self.assertFalse(ba.convergence_data['converged'])
    
    def test_calculate_bayes_factor(self):
        """Test Bayes factor calculation"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        ba.marginal_lnl = -1000.0
        
        # Test normal case
        bf = ba._calculate_bayes_factor(-1005.0)
        self.assertAlmostEqual(bf, 148.413, places=2)  # exp(5)
        
        # Test very large BF (should cap at infinity)
        bf_large = ba._calculate_bayes_factor(-1200.0)
        self.assertEqual(bf_large, float('inf'))
        
        # Test very small BF (should cap at 0)
        bf_small = ba._calculate_bayes_factor(-900.0)
        self.assertEqual(bf_small, 0.0)
    
    def test_calculate_bayes_factor_no_unconstrained(self):
        """Test Bayes factor calculation without unconstrained likelihood"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # marginal_lnl is None by default
        bf = ba._calculate_bayes_factor(-1005.0)
        self.assertIsNone(bf)
    
    def test_generate_constraint_analysis_all_taxa(self):
        """Test constraint analysis with all taxa (should skip)"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        all_taxa = self.data_processor.get_taxon_names()
        tree_filename, likelihood = ba.generate_constraint_analysis(
            clade_taxa=all_taxa,
            constraint_id="all_taxa_test",
            timeout_sec=60
        )
        
        self.assertIsNone(tree_filename)
        self.assertIsNone(likelihood)
    
    def test_generate_constraint_analysis_empty_taxa(self):
        """Test constraint analysis with empty taxa list"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        tree_filename, likelihood = ba.generate_constraint_analysis(
            clade_taxa=[],
            constraint_id="empty_test",
            timeout_sec=60
        )
        
        self.assertIsNone(tree_filename)
        self.assertIsNone(likelihood)
    
    def test_validate_bayesian_analysis(self):
        """Test Bayesian analysis validation"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        # Initially invalid
        self.assertFalse(ba.validate_bayesian_analysis())
        
        # Add consensus tree but no likelihood
        ba.consensus_tree = MagicMock()
        self.assertFalse(ba.validate_bayesian_analysis())
        
        # Add likelihood - but still won't validate because convergence is checked by default
        ba.marginal_lnl = -1000.0
        # This will fail because check_convergence=True by default but convergence_data is empty
        self.assertFalse(ba.validate_bayesian_analysis())
        
        # Test convergence failure with convergence checking disabled
        ba.check_convergence = False
        self.assertTrue(ba.validate_bayesian_analysis())
        
        # Test convergence failure with convergence checking enabled
        ba.check_convergence = True
        ba.convergence_data = {'converged': False}
        self.assertFalse(ba.validate_bayesian_analysis())
        
        # Fix convergence
        ba.convergence_data = {'converged': True}
        self.assertTrue(ba.validate_bayesian_analysis())
    
    def test_get_analysis_summary(self):
        """Test analysis summary generation"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor,
            bayes_model="GTR",
            bayes_ngen=100000,
            bayes_chains=4
        )
        
        # Set some results
        ba.consensus_tree = MagicMock()
        ba.marginal_lnl = -1000.0
        ba.convergence_data = {
            'converged': True,
            'min_ess': 250.0,
            'max_psrf': 1.005,
            'asdsf': 0.008
        }
        ba.constraint_analyses = {
            'clade1': {'marginal_lnl': -1005.0},
            'clade2': {'marginal_lnl': -1010.0}
        }
        
        summary = ba.get_analysis_summary()
        
        self.assertTrue(summary['consensus_tree_available'])
        self.assertEqual(summary['marginal_likelihood'], -1000.0)
        self.assertTrue(summary['convergence_checked'])
        self.assertTrue(summary['converged'])
        self.assertEqual(summary['num_constraint_analyses'], 2)
        self.assertEqual(summary['bayes_model'], "GTR")
        self.assertEqual(summary['mcmc_generations'], 100000)
        self.assertEqual(summary['mcmc_chains'], 4)
        self.assertEqual(summary['min_ess'], 250.0)
        self.assertEqual(summary['max_psrf'], 1.005)
        self.assertEqual(summary['asdsf'], 0.008)
        self.assertEqual(summary['constraint_ids'], ['clade1', 'clade2'])
    
    def test_calculate_bayes_factors(self):
        """Test Bayes factor calculation for all constraints"""
        ba = BayesianAnalyzer(
            file_manager=self.file_manager,
            external_tools=self.external_tools,
            data_processor=self.data_processor
        )
        
        ba.marginal_lnl = -1000.0
        ba.constraint_analyses = {
            'clade1': {'marginal_lnl': -1005.0},
            'clade2': {'marginal_lnl': -1010.0},
            'clade3': {'marginal_lnl': None}  # Should be skipped
        }
        
        bayes_factors = ba.calculate_bayes_factors()
        
        self.assertAlmostEqual(bayes_factors['clade1'], 148.413, places=2)
        self.assertAlmostEqual(bayes_factors['clade2'], 22026.466, places=0)
        self.assertNotIn('clade3', bayes_factors)


class TestBayesianAnalyzerIntegration(unittest.TestCase):
    """Integration tests for BayesianAnalyzer component"""
    
    def setUp(self):
        """Set up integration test fixtures"""
        self.temp_dir = Path(tempfile.mkdtemp())
        
        # Create realistic alignment
        self.realistic_alignment = """>Homo_sapiens
ATCGATCGATCGNNNN----ATCGATCG
>Pan_troglodytes
ATCGATCGATCGNNNN----ATCGATCG
>Mus_musculus
ATCGATCGATCGNNNN----ATCGATCG
>Rattus_norvegicus
ATCGATCGATCGNNNN----ATCGATCG
"""
        self.alignment_file = self.temp_dir / "realistic.fasta"
        with open(self.alignment_file, 'w') as f:
            f.write(self.realistic_alignment)
    
    def tearDown(self):
        """Clean up integration test fixtures"""
        import shutil
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)
    
    @patch('Bio.Phylo.read')
    def test_complete_bayesian_workflow(self, mock_phylo_read):
        """Test complete Bayesian analysis workflow"""
        # Initialize all components
        with FileManager(temp_dir=self.temp_dir, debug=True) as fm:
            et = ExternalTools(debug=True)
            dp = DataProcessor(data_type="dna", debug=True)
            
            # Load and process alignment
            dp.load_alignment(self.alignment_file, "fasta")
            nexus_file = fm.get_temp_path("alignment.nex")
            dp.convert_to_nexus(nexus_file)
            
            # Mock MrBayes execution
            with patch.object(et, 'run_mrbayes') as mock_mrbayes:
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_result.stdout = "MrBayes analysis completed"
                mock_mrbayes.return_value = mock_result
                
                # Create mock output files
                con_tree_file = fm.get_temp_path("bayesian.con.tre")
                ss_file = fm.get_temp_path("bayesian.ss")
                pstat_file = fm.get_temp_path("bayesian.run1.pstat")
                mcmc_file = fm.get_temp_path("bayesian.mcmc")
                
                with open(con_tree_file, 'w') as f:
                    f.write("(Homo_sapiens:0.1,(Pan_troglodytes:0.1,(Mus_musculus:0.1,Rattus_norvegicus:0.1):0.1):0.1);")
                
                with open(ss_file, 'w') as f:
                    f.write("Stepping-stone estimate of marginal likelihood = -2000.567")
                
                with open(pstat_file, 'w') as f:
                    f.write("Gen\tLnL\tLnPr\tTL\tpi(A)\tESS\tPSRF\n")
                    f.write("1000\t-2345.67\t-12.34\t0.5\t0.25\t250.5\t1.002\n")
                
                with open(mcmc_file, 'w') as f:
                    f.write("Average standard deviation of split frequencies: 0.005678")
                
                # Mock phylo.read
                mock_tree = MagicMock()
                mock_phylo_read.return_value = mock_tree
                
                # Initialize Bayesian analyzer
                ba = BayesianAnalyzer(
                    file_manager=fm,
                    external_tools=et,
                    data_processor=dp,
                    bayes_model="GTR",
                    bayes_ngen=10000,  # Small for testing
                    marginal_likelihood="ss",
                    debug=True
                )
                
                # Run Bayesian analysis
                result_tree = ba.run_bayesian_analysis(timeout_sec=60)
                
                # Verify results
                self.assertEqual(result_tree, mock_tree)
                self.assertEqual(ba.get_marginal_likelihood(), -2000.567)
                self.assertTrue(ba.validate_bayesian_analysis())
                
                # Verify MrBayes was called correctly
                mock_mrbayes.assert_called_once()
                
                # Test constraint analysis
                clade_taxa = ["Homo_sapiens", "Pan_troglodytes"]
                
                # Create constraint files
                constraint_con_file = fm.get_temp_path("constraint_clade1.con.tre")
                constraint_ss_file = fm.get_temp_path("constraint_clade1.ss")
                
                with open(constraint_con_file, 'w') as f:
                    f.write("((Homo_sapiens:0.1,Pan_troglodytes:0.1):0.1,(Mus_musculus:0.1,Rattus_norvegicus:0.1):0.1);")
                with open(constraint_ss_file, 'w') as f:
                    f.write("Stepping-stone estimate of marginal likelihood = -2010.234")
                
                tree_filename, likelihood = ba.generate_constraint_analysis(
                    clade_taxa=clade_taxa,
                    constraint_id="clade1",
                    timeout_sec=60
                )
                
                # Verify constraint results
                self.assertEqual(tree_filename, "constraint_clade1.con.tre")
                self.assertEqual(likelihood, -2010.234)
                
                # Check Bayes factors
                bayes_factors = ba.calculate_bayes_factors()
                self.assertIn('clade1', bayes_factors)
                self.assertGreater(bayes_factors['clade1'], 0)
                
                # Get analysis summary
                summary = ba.get_analysis_summary()
                self.assertTrue(summary['consensus_tree_available'])
                self.assertEqual(summary['marginal_likelihood'], -2000.567)
                self.assertEqual(summary['num_constraint_analyses'], 1)


def run_bayesian_analyzer_tests():
    """Run the BayesianAnalyzer test suite"""
    print("=" * 80)
    print("BAYESIAN ANALYZER COMPONENT TEST SUITE")
    print("=" * 80)
    print("Testing extracted Bayesian analysis functionality...")
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestBayesianAnalyzer))
    suite.addTests(loader.loadTestsFromTestCase(TestBayesianAnalyzerIntegration))
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2, stream=sys.stdout)
    result = runner.run(suite)
    
    # Summary
    print("\n" + "=" * 80)
    print("BAYESIAN ANALYZER TEST RESULTS")
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
    print(f"\nBayesianAnalyzer Test Status: {'✅ PASS' if success else '❌ FAIL'}")
    
    return success


if __name__ == "__main__":
    success = run_bayesian_analyzer_tests()
    sys.exit(0 if success else 1)