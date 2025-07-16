#!/usr/bin/env python3

import os
import sys
import argparse
import configparser
import json
import numpy as np
from Bio import Phylo, AlignIO, SeqIO
import tempfile
import shutil
import subprocess
import logging
import re
import multiprocessing
import time
import datetime
import glob
from pathlib import Path

VERSION = "1.1.0"
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

# --- Constants for Filenames ---
NEXUS_ALIGNMENT_FN = "alignment.nex"
ML_TREE_FN = "ml_tree.tre"
ML_SCORE_FN = "ml_score.txt"
ML_SEARCH_NEX_FN = "ml_search.nex"
ML_LOG_FN = "paup_ml.log"

AU_TEST_NEX_FN = "au_test.nex"
AU_TEST_SCORE_FN = "au_test_results.txt"
AU_LOG_FN = "paup_au.log"


class panDecayIndices:
    """
    Implements phylogenetic decay indices (Bremer support) using multiple approaches.
    Calculates support by comparing optimal trees with constrained trees using:
    - ML (Maximum Likelihood) with AU test
    - Bayesian analysis with marginal likelihood comparisons
    - Parsimony analysis with step differences
    """

    def __init__(self, alignment_file, alignment_format="fasta", model="GTR+G",
                 temp_dir: Path = None, paup_path="paup", threads="auto",
                 starting_tree: Path = None, data_type="dna",
                 debug=False, keep_files=False, gamma_shape=None, prop_invar=None,
                 base_freq=None, rates=None, protein_model=None, nst=None,
                 parsmodel=None, paup_block=None, analysis_mode="ml",
                 bayesian_software=None, mrbayes_path="mb",
                 bayes_model=None, bayes_ngen=1000000, bayes_burnin=0.25,
                 bayes_chains=4, bayes_sample_freq=1000, marginal_likelihood="ss",
                 ss_alpha=0.4, ss_nsteps=50, use_mpi=False, mpi_processors=None,
                 mpirun_path="mpirun", use_beagle=False, beagle_device="auto",
                 beagle_precision="double", beagle_scaling="dynamic",
                 constraint_mode="all", test_branches=None, constraint_file=None,
                 config_constraints=None, check_convergence=True, min_ess=200,
                 max_psrf=1.01, max_asdsf=0.01, convergence_strict=False,
                 mrbayes_parse_timeout=30.0, output_style="unicode",
                 normalize_bd=True, bd_normalization_methods=None):

        self.alignment_file = Path(alignment_file)
        self.alignment_format = alignment_format
        self.model_str = model # Keep original model string for reference
        self.paup_path = paup_path
        self.starting_tree = starting_tree # Already a Path or None from main
        self.debug = debug
        self.keep_files = keep_files or debug
        self.gamma_shape_arg = gamma_shape
        self.prop_invar_arg = prop_invar
        self.base_freq_arg = base_freq
        self.rates_arg = rates
        self.protein_model_arg = protein_model
        self.nst_arg = nst
        self.parsmodel_arg = parsmodel # For discrete data, used in _convert_model_to_paup
        self.user_paup_block = paup_block # Raw user block content
        self._files_to_cleanup = []
        
        # Parse analysis mode to set boolean flags
        self.analysis_mode = analysis_mode
        self.do_ml = "ml" in analysis_mode or analysis_mode == "all"
        self.do_bayesian = "bayesian" in analysis_mode or analysis_mode == "all"
        self.do_parsimony = "parsimony" in analysis_mode or analysis_mode == "all"
        
        # Bayesian analysis parameters
        self.bayesian_software = bayesian_software
        self.mrbayes_path = mrbayes_path
        self.bayes_model = bayes_model or model  # Use ML model if not specified
        self.bayes_ngen = bayes_ngen
        self.bayes_burnin = bayes_burnin
        self.bayes_chains = bayes_chains
        self.bayes_sample_freq = bayes_sample_freq
        self.marginal_likelihood = marginal_likelihood
        self.ss_alpha = ss_alpha
        self.ss_nsteps = ss_nsteps
        
        # MPI and BEAGLE parameters
        self.use_mpi = use_mpi
        self.mpi_processors = mpi_processors
        self.mpirun_path = mpirun_path
        self.use_beagle = use_beagle
        self.beagle_device = beagle_device
        self.beagle_precision = beagle_precision
        self.beagle_scaling = beagle_scaling
        
        # Constraint selection parameters
        self.constraint_mode = constraint_mode
        self.test_branches = test_branches
        self.constraint_file = constraint_file
        self.config_constraints = config_constraints or {}
        
        # Convergence checking parameters
        self.check_convergence = check_convergence
        self.min_ess = min_ess
        self.max_psrf = max_psrf
        self.max_asdsf = max_asdsf
        self.convergence_strict = convergence_strict
        
        # Output and parsing parameters
        self.mrbayes_parse_timeout = mrbayes_parse_timeout
        self.output_style = output_style
        
        # Normalization parameters
        self.normalize_bd = normalize_bd
        self.bd_normalization_methods = bd_normalization_methods or ["per_site", "relative"]

        self.data_type = data_type.lower()
        if self.data_type not in ["dna", "protein", "discrete"]:
            logger.warning(f"Unknown data type: {data_type}, defaulting to DNA")
            self.data_type = "dna"

        if threads == "auto":
            total_cores = multiprocessing.cpu_count()
            if total_cores > 2:
                self.threads = total_cores - 2 # Leave 2 cores for OS/other apps
            elif total_cores > 1:
                self.threads = total_cores - 1 # Leave 1 core
            else:
                self.threads = 1 # Use 1 core if only 1 is available
            logger.info(f"Using 'auto' threads: PAUP* will be configured for {self.threads} thread(s) (leaving some for system).")
        elif str(threads).lower() == "all": # Add an explicit "all" option if you really want it
            self.threads = multiprocessing.cpu_count()
            logger.warning(f"PAUP* configured to use ALL {self.threads} threads. System may become unresponsive.")
        else:
            try:
                self.threads = int(threads)
                if self.threads < 1:
                    logger.warning(f"Thread count {self.threads} is invalid, defaulting to 1.")
                    self.threads = 1
                elif self.threads > multiprocessing.cpu_count():
                    logger.warning(f"Requested {self.threads} threads, but only {multiprocessing.cpu_count()} cores available. Using {multiprocessing.cpu_count()}.")
                    self.threads = multiprocessing.cpu_count()
            except ValueError:
                logger.warning(f"Invalid thread count '{threads}', defaulting to 1.")
                self.threads = 1

        logger.info(f"PAUP* will be configured to use up to {self.threads} thread(s).")

        # --- Temporary Directory Setup ---
        self._temp_dir_obj = None  # For TemporaryDirectory lifecycle
        if self.debug or self.keep_files or temp_dir:
            if temp_dir: # User-provided temp_dir (already a Path object)
                self.temp_path = temp_dir
                self.temp_path.mkdir(parents=True, exist_ok=True)
            else: # Debug/keep_files, create a timestamped dir
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                debug_runs_path = Path.cwd() / "debug_runs"
                debug_runs_path.mkdir(parents=True, exist_ok=True)
                self.work_dir_name = f"mldecay_{timestamp}"
                self.temp_path = debug_runs_path / self.work_dir_name
                self.temp_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Using temporary directory: {self.temp_path}")
        else: # Auto-cleanup
            self._temp_dir_obj = tempfile.TemporaryDirectory(prefix="mldecay_")
            self.temp_path = Path(self._temp_dir_obj.name)
            self.work_dir_name = self.temp_path.name
            logger.info(f"Using temporary directory (auto-cleanup): {self.temp_path}")

        # --- PAUP* Model Settings ---
        self.parsmodel = False # Default, will be set by _convert_model_to_paup if discrete
        if self.user_paup_block is None:
            self.paup_model_cmds = self._convert_model_to_paup(
                self.model_str, gamma_shape=self.gamma_shape_arg, prop_invar=self.prop_invar_arg,
                base_freq=self.base_freq_arg, rates=self.rates_arg,
                protein_model=self.protein_model_arg, nst=self.nst_arg,
                parsmodel_user_intent=self.parsmodel_arg # Pass user intent
            )
        else:
            logger.info("Using user-provided PAUP block for model specification.")
            self.paup_model_cmds = self.user_paup_block # This is the content of the block

        # --- Alignment Handling ---
        try:
            self.alignment = AlignIO.read(str(self.alignment_file), self.alignment_format)
            self.alignment_length = self.alignment.get_alignment_length()
            logger.info(f"Loaded alignment: {len(self.alignment)} sequences, {self.alignment_length} sites.")
        except Exception as e:
            logger.error(f"Failed to load alignment '{self.alignment_file}': {e}")
            if self._temp_dir_obj: self._temp_dir_obj.cleanup() # Manual cleanup if init fails early
            raise

        if self.data_type == "discrete":
            if not self._validate_discrete_data():
                logger.warning("Discrete data validation failed based on content, proceeding but results may be unreliable.")

        if self.keep_files or self.debug: # Copy original alignment for debugging
             shutil.copy(str(self.alignment_file), self.temp_path / f"original_alignment.{self.alignment_format}")

        self.nexus_file_path = self.temp_path / NEXUS_ALIGNMENT_FN
        
        # If input is already NEXUS, copy it directly instead of converting
        if self.alignment_format.lower() == "nexus":
            logger.info(f"Input is already NEXUS format, copying directly to temp directory")
            shutil.copy(str(self.alignment_file), str(self.nexus_file_path))
        else:
            self._convert_to_nexus() # Writes to self.nexus_file_path
            
        # Validate that NEXUS file exists and has content
        if not self.nexus_file_path.exists():
            raise FileNotFoundError(f"NEXUS file was not created at {self.nexus_file_path}")
        if self.nexus_file_path.stat().st_size == 0:
            raise ValueError(f"NEXUS file at {self.nexus_file_path} is empty")

        self.ml_tree = None
        self.ml_likelihood = None
        self.decay_indices = {}

    def __del__(self):
        """Cleans up temporary files if TemporaryDirectory object was used."""
        if hasattr(self, '_temp_dir_obj') and self._temp_dir_obj:
            logger.debug(f"Attempting to cleanup temp_dir_obj for {self.temp_path}")
            self._temp_dir_obj.cleanup()
            logger.info(f"Auto-cleaned temporary directory: {self.temp_path}")
        elif hasattr(self, 'temp_path') and self.temp_path.exists() and not self.keep_files:
            logger.info(f"Manually cleaning up temporary directory: {self.temp_path}")
            shutil.rmtree(self.temp_path)
        elif hasattr(self, 'temp_path') and self.keep_files:
            logger.info(f"Keeping temporary directory: {self.temp_path}")

    def _convert_model_to_paup(self, model_str, gamma_shape, prop_invar, base_freq, rates, protein_model, nst, parsmodel_user_intent):
        """Converts model string and params to PAUP* 'lset' command part (without 'lset' itself)."""
        cmd_parts = []
        has_gamma = "+G" in model_str.upper()
        has_invar = "+I" in model_str.upper()
        base_model_name = model_str.split("+")[0].upper()
        
        if self.debug:
            logger.debug(f"Model conversion debug - model_str: {model_str}, base_model_name: {base_model_name}, data_type: {self.data_type}")

        if self.data_type == "dna":
            if nst is not None: cmd_parts.append(f"nst={nst}")
            elif base_model_name == "GTR": cmd_parts.append("nst=6")
            elif base_model_name in ["HKY", "K2P", "K80", "TN93"]: cmd_parts.append("nst=2")
            elif base_model_name in ["JC", "JC69", "F81"]: cmd_parts.append("nst=1")
            else:
                logger.warning(f"Unknown DNA model: {base_model_name}, defaulting to GTR (nst=6).")
                cmd_parts.append("nst=6")

            current_nst = next((p.split('=')[1] for p in cmd_parts if "nst=" in p), None)
            if current_nst == '6' or (base_model_name == "GTR" and nst is None):
                cmd_parts.append("rmatrix=estimate")
            elif current_nst == '2' or (base_model_name in ["HKY", "K2P"] and nst is None):
                cmd_parts.append("tratio=estimate")

            if base_freq: cmd_parts.append(f"basefreq={base_freq}")
            elif base_model_name in ["JC", "K2P", "JC69", "K80"] : cmd_parts.append("basefreq=equal")
            else: cmd_parts.append("basefreq=estimate") # GTR, HKY, F81, TN93 default to estimate

        elif self.data_type == "protein":
            valid_protein_models = ["JTT", "WAG", "LG", "DAYHOFF", "MTREV", "CPREV", "BLOSUM62", "HIVB", "HIVW"]
            if protein_model: cmd_parts.append(f"protein={protein_model.lower()}")
            elif base_model_name.upper() in valid_protein_models: cmd_parts.append(f"protein={base_model_name.lower()}")
            else:
                logger.warning(f"Unknown protein model: {base_model_name}, defaulting to JTT.")
                cmd_parts.append("protein=jtt")

        elif self.data_type == "discrete": # Typically Mk model
            cmd_parts.append("nst=1") # For standard Mk
            if base_freq: cmd_parts.append(f"basefreq={base_freq}")
            else: cmd_parts.append("basefreq=equal") # Default for Mk

            if parsmodel_user_intent is None: # If user didn't specify, default to True for discrete
                self.parsmodel = True
            else:
                self.parsmodel = bool(parsmodel_user_intent)


        # Common rate variation and invariable sites for all data types
        if rates: cmd_parts.append(f"rates={rates}")
        elif has_gamma: cmd_parts.append("rates=gamma")
        else: cmd_parts.append("rates=equal")

        current_rates = next((p.split('=')[1] for p in cmd_parts if "rates=" in p), "equal")
        if gamma_shape is not None and (current_rates == "gamma" or has_gamma):
            cmd_parts.append(f"shape={gamma_shape}")
        elif current_rates == "gamma" or has_gamma:
            cmd_parts.append("shape=estimate")

        if prop_invar is not None:
            cmd_parts.append(f"pinvar={prop_invar}")
        elif has_invar:
            cmd_parts.append("pinvar=estimate")
        else: # No +I and no explicit prop_invar
            cmd_parts.append("pinvar=0")

        return "lset " + " ".join(cmd_parts) + ";"

    def _validate_discrete_data(self):
        """Validate that discrete data contains only 0, 1, -, ? characters."""
        if self.data_type == "discrete":
            valid_chars = set("01-?")
            for record in self.alignment:
                seq_chars = set(str(record.seq).upper()) # Convert to upper for case-insensitivity if needed
                invalid_chars = seq_chars - valid_chars
                if invalid_chars:
                    logger.warning(f"Sequence {record.id} contains invalid discrete characters: {invalid_chars}. Expected only 0, 1, -, ?.")
                    return False
        return True

    def _format_taxon_for_paup(self, taxon_name):
        """Format a taxon name for PAUP* (handles spaces, special chars by quoting)."""
        if not isinstance(taxon_name, str): taxon_name = str(taxon_name)
        # PAUP* needs quotes if name contains whitespace or NEXUS special chars: ( ) [ ] { } / \ , ; = * ` " ' < >
        if re.search(r'[\s\(\)\[\]\{\}/\\,;=\*`"\'<>]', taxon_name) or ':' in taxon_name: # Colon also problematic
            return f"'{taxon_name.replace(chr(39), '_')}'" # chr(39) is single quote

        return taxon_name

    def _convert_to_nexus(self):
        """Converts alignment to NEXUS, writes to self.nexus_file_path."""
        try:
            with open(self.nexus_file_path, 'w') as f:
                f.write("#NEXUS\n\n")
                f.write("BEGIN DATA;\n")
                dt = "DNA"
                if self.data_type == "protein": dt = "PROTEIN"
                elif self.data_type == "discrete": dt = "STANDARD"

                f.write(f"  DIMENSIONS NTAX={len(self.alignment)} NCHAR={self.alignment.get_alignment_length()};\n")
                format_line = f"  FORMAT DATATYPE={dt} MISSING=? GAP=- INTERLEAVE=NO"
                if self.data_type == "discrete":
                    format_line += " SYMBOLS=\"01\"" # Assuming binary discrete data
                f.write(format_line + ";\n")
                f.write("  MATRIX\n")
                for record in self.alignment:
                    f.write(f"  {self._format_taxon_for_paup(record.id)} {record.seq}\n")
                f.write("  ;\nEND;\n")

                if self.data_type == "discrete":
                    f.write("\nBEGIN ASSUMPTIONS;\n")
                    f.write("  OPTIONS DEFTYPE=UNORD POLYTCOUNT=MINSTEPS;\n") # Common for Mk
                    f.write("END;\n")
            logger.info(f"Converted alignment to NEXUS: {self.nexus_file_path}")
        except Exception as e:
            logger.error(f"Failed to convert alignment to NEXUS: {e}")
            raise

    def _get_paup_model_setup_cmds(self):
        """Returns the model setup command string(s) for PAUP* script."""
        if self.user_paup_block is None:
            # self.paup_model_cmds is like "lset nst=6 ...;"
            # Remove "lset " for combining with nthreads, keep ";"
            model_params_only = self.paup_model_cmds.replace("lset ", "", 1)
            base_cmds = [
                f"lset nthreads={self.threads} {model_params_only}", # model_params_only includes the trailing ";"
                "set criterion=likelihood;"
            ]
            if self.data_type == "discrete":
                base_cmds.append("options deftype=unord polytcount=minsteps;")
                if self.parsmodel: # self.parsmodel is set by _convert_model_to_paup
                    base_cmds.append("set parsmodel=yes;")
            return "\n".join(f"    {cmd}" for cmd in base_cmds)
        else:
            # self.paup_model_cmds is the user's raw block content
            # Assume it sets threads, model, criterion, etc.
            return self.paup_model_cmds # Return as is, for direct insertion

    def _run_paup_command_file(self, paup_cmd_filename_str: str, log_filename_str: str, timeout_sec: int = None):
        """Runs a PAUP* .nex command file located in self.temp_path."""
        paup_cmd_file = self.temp_path / paup_cmd_filename_str
        # The main log file will capture both stdout and stderr from PAUP*
        combined_log_file_path = self.temp_path / log_filename_str

        if not paup_cmd_file.exists():
            logger.error(f"PAUP* command file not found: {paup_cmd_file}")
            raise FileNotFoundError(f"PAUP* command file not found: {paup_cmd_file}")

        logger.debug(f"Running PAUP* command file: {paup_cmd_filename_str} (Log: {log_filename_str})")

        # stdout_content and stderr_content will be filled for logging/debugging if needed
        stdout_capture = ""
        stderr_capture = ""

        try:
            # Open the log file once for both stdout and stderr
            with open(combined_log_file_path, 'w') as f_log:
                process = subprocess.Popen(
                    [self.paup_path, "-n", paup_cmd_filename_str],
                    cwd=str(self.temp_path),
                    stdout=subprocess.PIPE, # Capture stdout
                    stderr=subprocess.PIPE, # Capture stderr
                    text=True,
                    universal_newlines=True # For text=True
                )

                # Read stdout and stderr in a non-blocking way or use communicate
                # communicate() is simpler and safer for handling potential deadlocks
                try:
                    stdout_capture, stderr_capture = process.communicate(timeout=timeout_sec)
                except subprocess.TimeoutExpired:
                    process.kill() # Ensure process is killed on timeout
                    stdout_capture, stderr_capture = process.communicate() # Try to get any remaining output
                    logger.error(f"PAUP* command {paup_cmd_filename_str} timed out after {timeout_sec}s.")
                    f_log.write(f"--- PAUP* Execution Timed Out ({timeout_sec}s) ---\n")
                    if stdout_capture: f_log.write("--- STDOUT (partial) ---\n" + stdout_capture)
                    if stderr_capture: f_log.write("\n--- STDERR (partial) ---\n" + stderr_capture)
                    raise # Re-raise the TimeoutExpired exception

                # Write captured output to the log file
                f_log.write("--- STDOUT ---\n")
                f_log.write(stdout_capture if stdout_capture else "No stdout captured.\n")
                if stderr_capture:
                    f_log.write("\n--- STDERR ---\n")
                    f_log.write(stderr_capture)

                retcode = process.returncode
                if retcode != 0:
                    logger.error(f"PAUP* execution failed for {paup_cmd_filename_str}. Exit code: {retcode}")
                    # The log file already contains stdout/stderr
                    logger.error(f"PAUP* stdout/stderr saved to {combined_log_file_path}. Stderr sample: {stderr_capture[:500]}...")
                    # Raise an equivalent of CalledProcessError
                    raise subprocess.CalledProcessError(retcode, process.args, output=stdout_capture, stderr=stderr_capture)

            if self.debug:
                logger.debug(f"PAUP* output saved to: {combined_log_file_path}")
                logger.debug(f"PAUP* stdout sample (from capture):\n{stdout_capture[:500]}...")
                if stderr_capture: logger.debug(f"PAUP* stderr sample (from capture):\n{stderr_capture[:500]}...")

            # Return a simple object that mimics CompletedProcess for the parts we use
            # Or adjust callers to expect (stdout_str, stderr_str, retcode) tuple
            class MockCompletedProcess:
                def __init__(self, args, returncode, stdout, stderr):
                    self.args = args
                    self.returncode = returncode
                    self.stdout = stdout
                    self.stderr = stderr

            return MockCompletedProcess(process.args, retcode, stdout_capture, stderr_capture)

        except subprocess.CalledProcessError: # Already logged, just re-raise
            raise
        except subprocess.TimeoutExpired: # Already logged, just re-raise
            raise
        except Exception as e:
            # Fallback for other errors during Popen or communicate
            logger.error(f"Unexpected error running PAUP* for {paup_cmd_filename_str}: {e}")
            # Attempt to write to log if f_log was opened
            if 'f_log' in locals() and not f_log.closed:
                 f_log.write(f"\n--- Script Error during PAUP* execution ---\n{str(e)}\n")
            raise

    def _parse_likelihood_from_score_file(self, score_file_path: Path):
        if not score_file_path.exists():
            logger.warning(f"Score file not found: {score_file_path}")
            return None
        try:
            content = score_file_path.read_text()
            if self.debug: logger.debug(f"Score file ({score_file_path}) content:\n{content}")

            lines = content.splitlines()
            header_idx, lnl_col_idx = -1, -1

            for i, line_text in enumerate(lines):
                norm_line = ' '.join(line_text.strip().lower().split()) # Normalize
                if "tree" in norm_line and ("-lnl" in norm_line or "loglk" in norm_line or "likelihood" in norm_line):
                    header_idx = i
                    headers = norm_line.split()
                    for col_name in ["-lnl", "loglk", "likelihood", "-loglk"]:
                        if col_name in headers:
                            lnl_col_idx = headers.index(col_name)
                            break
                    if lnl_col_idx != -1: break

            if header_idx == -1 or lnl_col_idx == -1:
                logger.warning(f"Could not find valid header or likelihood column in {score_file_path}.")
                return None
            logger.debug(f"Found LNL column at index {lnl_col_idx} in header: {lines[header_idx].strip()}")

            for i in range(header_idx + 1, len(lines)):
                data_line_text = lines[i].strip()
                if not data_line_text: continue # Skip empty

                parts = data_line_text.split()
                if len(parts) > lnl_col_idx:
                    try:
                        val_str = parts[lnl_col_idx]
                        if '*' in val_str : # Handle cases like '**********' or if PAUP adds flags
                            logger.warning(f"Likelihood value problematic (e.g., '******') in {score_file_path}, line: '{data_line_text}'")
                            continue # Try next line if multiple scores
                        likelihood = float(val_str)
                        logger.info(f"Parsed log-likelihood from {score_file_path}: {likelihood}")
                        return likelihood
                    except ValueError:
                        logger.warning(f"Could not convert LNL value to float: '{parts[lnl_col_idx]}' from line '{data_line_text}' in {score_file_path}")
                else: logger.warning(f"Insufficient columns in data line: '{data_line_text}' in {score_file_path}")
            logger.warning(f"No parsable data lines found after header in {score_file_path}")
            return None
        except Exception as e:
            logger.warning(f"Error reading/parsing score file {score_file_path}: {e}")
            return None

    def build_ml_tree(self):
        logger.info("Building maximum likelihood tree...")
        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", self._get_paup_model_setup_cmds()]

        if self.user_paup_block is None: # Standard model processing, add search commands
            if self.starting_tree and self.starting_tree.exists():
                start_tree_fn_temp = "start_tree.tre" # Relative to temp_path
                shutil.copy(str(self.starting_tree), str(self.temp_path / start_tree_fn_temp))
                script_cmds.extend([
                    f"gettrees file={start_tree_fn_temp};",
                    "lscores 1 / userbrlen=yes;", "hsearch start=current;"
                ])
            elif self.starting_tree: # Path provided but not found
                 logger.warning(f"Starting tree file not found: {self.starting_tree}. Performing standard search.")
                 script_cmds.append("hsearch start=stepwise addseq=random nreps=10;")
            else: # No starting tree
                script_cmds.append("hsearch start=stepwise addseq=random nreps=10;")

            script_cmds.extend([
                f"savetrees file={ML_TREE_FN} format=newick brlens=yes replace=yes;",
                f"lscores 1 / scorefile={ML_SCORE_FN} replace=yes;"
            ])
        else: # User-provided PAUP block, assume it handles search & save. Add defensively if not detected.
            block_lower = self.user_paup_block.lower()
            if "savetrees" not in block_lower:
                script_cmds.append(f"savetrees file={ML_TREE_FN} format=newick brlens=yes replace=yes;")
            if "lscores" not in block_lower and "lscore" not in block_lower : # Check for lscore too
                script_cmds.append(f"lscores 1 / scorefile={ML_SCORE_FN} replace=yes;")

        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        ml_search_cmd_path = self.temp_path / ML_SEARCH_NEX_FN
        ml_search_cmd_path.write_text(paup_script_content)
        if self.debug: 
            logger.debug(f"ML search PAUP* script ({ml_search_cmd_path}):\n{paup_script_content}")
        else:
            # Always log the PAUP* commands being executed for troubleshooting
            logger.debug(f"Executing PAUP* with model: {self.model_str}, threads: {self.threads}")

        try:
            paup_result = self._run_paup_command_file(ML_SEARCH_NEX_FN, ML_LOG_FN, timeout_sec=3600) # 1hr timeout

            self.ml_likelihood = self._parse_likelihood_from_score_file(self.temp_path / ML_SCORE_FN)
            if self.ml_likelihood is None and paup_result.stdout: # Fallback to log
                logger.info(f"Fallback: Parsing ML likelihood from PAUP* log {ML_LOG_FN}")
                patterns = [r'-ln\s*L\s*=\s*([0-9.]+)', r'likelihood\s*=\s*([0-9.]+)', r'score\s*=\s*([0-9.]+)']
                for p in patterns:
                    m = re.findall(p, paup_result.stdout, re.IGNORECASE)
                    if m: self.ml_likelihood = float(m[-1]); break
                if self.ml_likelihood: logger.info(f"Extracted ML likelihood from log: {self.ml_likelihood}")
                else: logger.warning("Could not extract ML likelihood from PAUP* log.")

            ml_tree_path = self.temp_path / ML_TREE_FN
            if ml_tree_path.exists() and ml_tree_path.stat().st_size > 0:
                # Clean the tree file if it has metadata after semicolon
                cleaned_tree_path = self._clean_newick_tree(ml_tree_path)
                self.ml_tree = Phylo.read(str(cleaned_tree_path), "newick")
                logger.info(f"Successfully built ML tree. Log-likelihood: {self.ml_likelihood if self.ml_likelihood is not None else 'N/A'}")
                if self.ml_likelihood is None:
                    logger.error("ML tree built, but likelihood could not be determined. Analysis may be compromised.")
                    # Decide if this is a fatal error for downstream steps
            else:
                logger.error(f"ML tree file {ml_tree_path} not found or is empty after PAUP* run.")
                raise FileNotFoundError(f"ML tree file missing or empty: {ml_tree_path}")
        except Exception as e:
            logger.error(f"ML tree construction failed: {e}")
            raise # Re-raise to be handled by the main try-except block

    def _clean_newick_tree(self, tree_path, delete_cleaned=True):
        """
        Clean Newick tree files that may have metadata after the semicolon.

        Args:
            tree_path: Path to the tree file
            delete_cleaned: Whether to delete the cleaned file after use (if caller manages reading)

        Returns:
            Path to a cleaned tree file or the original path if no cleaning was needed
        """
        try:
            content = Path(tree_path).read_text()

            # Check if there's any text after a semicolon (including whitespace)
            semicolon_match = re.search(r';(.+)', content, re.DOTALL)
            if semicolon_match:
                # Get everything up to the first semicolon
                clean_content = content.split(';')[0] + ';'

                # Write the cleaned tree to a new file
                cleaned_path = Path(str(tree_path) + '.cleaned')
                cleaned_path.write_text(clean_content)

                # Mark the file for later deletion if requested
                if delete_cleaned:
                    self._files_to_cleanup.append(cleaned_path)

                if self.debug:
                    logger.debug(f"Original tree content: '{content}'")
                    logger.debug(f"Cleaned tree content: '{clean_content}'")

                logger.debug(f"Cleaned tree file {tree_path} - removed metadata after semicolon")
                return cleaned_path

            return tree_path  # No cleaning needed
        except Exception as e:
            logger.warning(f"Error cleaning Newick tree {tree_path}: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Traceback for tree cleaning error: {traceback.format_exc()}")
            return tree_path  # Return original path if cleaning fails

    def run_bootstrap_analysis(self, num_replicates=100):
        """
        Run bootstrap analysis with PAUP* to calculate support values.

        Args:
            num_replicates: Number of bootstrap replicates to perform

        Returns:
            The bootstrap consensus tree with support values, or None if analysis failed
        """
        # Define bootstrap constants
        BOOTSTRAP_NEX_FN = "bootstrap_search.nex"
        BOOTSTRAP_LOG_FN = "paup_bootstrap.log"
        BOOTSTRAP_TREE_FN = "bootstrap_trees.tre"

        logger.info(f"Running bootstrap analysis with {num_replicates} replicates...")

        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", self._get_paup_model_setup_cmds()]

        # Add bootstrap commands
        script_cmds.extend([
            f"set criterion=likelihood;",
            f"hsearch;",  # Find the ML tree first
            f"bootstrap nreps={num_replicates} search=heuristic keepall=no conlevel=50 / start=stepwise addseq=random nreps=1;",
            # The bootstrap command creates a consensus tree with support values
            # We'll extract the ML tree topology with bootstrap values
            f"describetrees 1 / brlens=yes;",  # Show the tree with bootstrap values
            f"savetrees from=1 to=1 file={BOOTSTRAP_TREE_FN} format=newick brlens=yes replace=yes supportValues=nodeLabels;"
        ])

        # Create and execute PAUP script
        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        bootstrap_cmd_path = self.temp_path / BOOTSTRAP_NEX_FN
        bootstrap_cmd_path.write_text(paup_script_content)

        if self.debug: logger.debug(f"Bootstrap PAUP* script ({bootstrap_cmd_path}):\n{paup_script_content}")

        try:
            # Run the bootstrap analysis - timeout based on number of replicates
            self._run_paup_command_file(BOOTSTRAP_NEX_FN, BOOTSTRAP_LOG_FN,
                                      timeout_sec=max(3600, 60 * num_replicates))

            # Get the bootstrap tree
            bootstrap_tree_path = self.temp_path / BOOTSTRAP_TREE_FN

            if bootstrap_tree_path.exists() and bootstrap_tree_path.stat().st_size > 0:
                # Log the bootstrap tree file content for debugging
                if self.debug:
                    bootstrap_content = bootstrap_tree_path.read_text()
                    logger.debug(f"Bootstrap tree file content:\n{bootstrap_content}")

                # Clean the tree file if it has metadata after semicolon
                cleaned_tree_path = self._clean_newick_tree(bootstrap_tree_path)

                # Log the cleaned bootstrap tree file for debugging
                if self.debug:
                    cleaned_content = cleaned_tree_path.read_text() if Path(cleaned_tree_path).exists() else "Cleaning failed"
                    logger.debug(f"Cleaned bootstrap tree file content:\n{cleaned_content}")

                try:
                    # Parse bootstrap values from tree file
                    bootstrap_tree = Phylo.read(str(cleaned_tree_path), "newick")
                    self.bootstrap_tree = bootstrap_tree

                    # Verify that bootstrap values are present
                    has_bootstrap_values = False
                    for node in bootstrap_tree.get_nonterminals():
                        if node.confidence is not None:
                            has_bootstrap_values = True
                            break

                    if has_bootstrap_values:
                        logger.info(f"Bootstrap analysis complete with {num_replicates} replicates and bootstrap values")
                    else:
                        logger.warning(f"Bootstrap tree found, but no bootstrap values detected. Check PAUP* output format.")

                    return bootstrap_tree
                except Exception as parse_error:
                    logger.error(f"Error parsing bootstrap tree: {parse_error}")
                    if self.debug:
                        import traceback
                        logger.debug(f"Traceback for bootstrap parse error: {traceback.format_exc()}")
                    return None
            else:
                logger.error(f"Bootstrap tree file not found or empty: {bootstrap_tree_path}")
                return None
        except Exception as e:
            logger.error(f"Bootstrap analysis failed: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Traceback: {traceback.format_exc()}")
            return None

    def _generate_and_score_constraint_tree(self, clade_taxa: list, tree_idx: int):
        # Returns (relative_tree_filename_str_or_None, likelihood_float_or_None)
        formatted_clade_taxa = [self._format_taxon_for_paup(t) for t in clade_taxa]
        if not formatted_clade_taxa : # Should not happen if called correctly
            logger.warning(f"Constraint {tree_idx}: No taxa provided for clade. Skipping.")
            return None, None

        # All taxa in alignment (already formatted by _format_taxon_for_paup in _convert_to_nexus if that logic was used)
        # For safety, re-format here if needed or assume names are simple.
        # Here, get raw IDs then format.
        all_raw_taxa_ids = [rec.id for rec in self.alignment]
        if len(clade_taxa) == len(all_raw_taxa_ids):
             logger.warning(f"Constraint {tree_idx}: Clade contains all taxa. Skipping as no outgroup possible for MONOPHYLY constraint.")
             return None, None


        clade_spec = "((" + ", ".join(formatted_clade_taxa) + "));"

        constr_tree_fn = f"constraint_tree_{tree_idx}.tre"
        constr_score_fn = f"constraint_score_{tree_idx}.txt"
        constr_cmd_fn = f"constraint_search_{tree_idx}.nex"
        constr_log_fn = f"paup_constraint_{tree_idx}.log"

        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", self._get_paup_model_setup_cmds()]
        script_cmds.extend([
            f"constraints clade_constraint (MONOPHYLY) = {clade_spec}",
            "set maxtrees=100 increase=auto;" # Sensible default for constrained search
        ])

        if self.user_paup_block is None: # Standard search
            script_cmds.extend([
                "hsearch start=stepwise addseq=random nreps=1;", # Initial unconstrained to get a tree in memory
                "hsearch start=1 enforce=yes converse=yes constraints=clade_constraint;",
                f"savetrees file={constr_tree_fn} format=newick brlens=yes replace=yes;",
                f"lscores 1 / scorefile={constr_score_fn} replace=yes;"
            ])
        else: # User PAUP block
            block_lower = self.user_paup_block.lower()
            if not any(cmd in block_lower for cmd in ["hsearch", "bandb", "alltrees"]): # If no search specified
                script_cmds.append("hsearch start=stepwise addseq=random nreps=1;")
            # Add enforce to the existing search or a new one. This is tricky.
            # Simplest: add a new constrained search. User might need to adjust their block.
            script_cmds.append("hsearch start=1 enforce=yes converse=yes constraints=clade_constraint;")
            if "savetrees" not in block_lower:
                script_cmds.append(f"savetrees file={constr_tree_fn} format=newick brlens=yes replace=yes;")
            if "lscores" not in block_lower and "lscore" not in block_lower:
                script_cmds.append(f"lscores 1 / scorefile={constr_score_fn} replace=yes;")

        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        cmd_file_path = self.temp_path / constr_cmd_fn
        cmd_file_path.write_text(paup_script_content)
        if self.debug: logger.debug(f"Constraint search {tree_idx} script ({cmd_file_path}):\n{paup_script_content}")

        try:
            self._run_paup_command_file(constr_cmd_fn, constr_log_fn, timeout_sec=600)

            score_file_path = self.temp_path / constr_score_fn
            constrained_lnl = self._parse_likelihood_from_score_file(score_file_path)

            tree_file_path = self.temp_path / constr_tree_fn
            if tree_file_path.exists() and tree_file_path.stat().st_size > 0:
                return constr_tree_fn, constrained_lnl # Return relative filename
            else:
                logger.error(f"Constraint tree file {tree_file_path} (idx {tree_idx}) not found or empty.")
                # Try to get LNL from log if score file failed and tree missing
                if constrained_lnl is None:
                    log_content = (self.temp_path / constr_log_fn).read_text()
                    patterns = [r'-ln\s*L\s*=\s*([0-9.]+)', r'likelihood\s*=\s*([0-9.]+)', r'score\s*=\s*([0-9.]+)']
                    for p in patterns:
                        m = re.findall(p, log_content, re.IGNORECASE)
                        if m: constrained_lnl = float(m[-1]); break
                    if constrained_lnl: logger.info(f"Constraint {tree_idx}: LNL from log: {constrained_lnl} (tree file missing)")
                return None, constrained_lnl
        except Exception as e:
            logger.error(f"Constraint tree generation/scoring failed for index {tree_idx}: {e}")
            return None, None
    
    # ===== Bayesian Analysis Methods =====
    
    def _generate_mrbayes_nexus(self, constraint_tree_file=None, clade_taxa=None, constraint_id=None):
        """
        Generate MrBayes NEXUS file with optional constraint for non-monophyly.
        
        Args:
            constraint_tree_file: Path to save the constraint tree (for file-based constraints)
            clade_taxa: List of taxa to constrain as non-monophyletic
            constraint_id: Identifier for the constraint
            
        Returns:
            String containing MrBayes block
        """
        blocks = []
        
        # MrBayes block
        blocks.append("begin mrbayes;")
        
        # If constraint is specified, add it
        if clade_taxa and constraint_id:
            # Format taxa for MrBayes
            formatted_taxa = [self._format_taxon_for_paup(t) for t in clade_taxa]
            taxa_string = " ".join(formatted_taxa)
            
            # Create a negative constraint to force NON-monophyly
            # Based on MrBayes manual section 6.8.1: negative constraints 'ban' trees
            # where listed taxa form a monophyletic group
            blocks.append(f"    constraint broken_{constraint_id} negative = {taxa_string};")
            blocks.append(f"    prset topologypr = constraints(broken_{constraint_id});")
        
        # Model settings based on data type
        if self.data_type == "dna":
            if self.debug:
                logger.debug(f"MrBayes model debug - bayes_model: {self.bayes_model}, data_type: {self.data_type}")
            
            # Determine nst parameter
            if "GTR" in self.bayes_model.upper():
                nst_val = "6"
            elif "HKY" in self.bayes_model.upper():
                nst_val = "2"
            elif "JC" in self.bayes_model.upper():
                nst_val = "1"
            else:
                nst_val = "6"  # Default to GTR
                
            # Determine rates parameter
            if "+G" in self.bayes_model and "+I" in self.bayes_model:
                rates_val = "invgamma"
            elif "+G" in self.bayes_model:
                rates_val = "gamma"
            elif "+I" in self.bayes_model:
                rates_val = "propinv"
            else:
                rates_val = "equal"
                
            # Combine into single lset command
            blocks.append(f"    lset nst={nst_val} rates={rates_val};")
                
        elif self.data_type == "protein":
            protein_models = {
                "JTT": "jones", "WAG": "wag", "LG": "lg", 
                "DAYHOFF": "dayhoff", "CPREV": "cprev", "MTREV": "mtrev"
            }
            model_name = "wag"  # default
            for pm, mb_name in protein_models.items():
                if pm in self.bayes_model.upper():
                    model_name = mb_name
                    break
            blocks.append(f"    prset aamodelpr=fixed({model_name});")
            
            if "+G" in self.bayes_model:
                blocks.append("    lset rates=gamma;")
        
        # BEAGLE settings if enabled
        if self.use_beagle:
            beagle_cmd = f"    set usebeagle=yes beagledevice={self.beagle_device} "
            beagle_cmd += f"beagleprecision={self.beagle_precision} "
            beagle_cmd += f"beaglescaling={self.beagle_scaling};"
            blocks.append(beagle_cmd)
        
        # MCMC settings
        blocks.append(f"    mcmc ngen={self.bayes_ngen} samplefreq={self.bayes_sample_freq} "
                     f"nchains={self.bayes_chains} savebrlens=yes printfreq=1000 diagnfreq=5000;")
        
        # Summary commands first
        burnin_samples = int(self.bayes_ngen / self.bayes_sample_freq * self.bayes_burnin)
        blocks.append(f"    sump burnin={burnin_samples};")
        blocks.append(f"    sumt burnin={burnin_samples};")
        
        # Add stepping-stone sampling if requested
        if self.marginal_likelihood == "ss":
            # Stepping-stone sampling parameters
            # alpha: shape parameter for Beta distribution (default 0.4)
            # nsteps: number of steps between prior and posterior (default 50)
            blocks.append(f"    ss alpha={self.ss_alpha} nsteps={self.ss_nsteps} "
                         f"burnin={burnin_samples};")
        
        blocks.append("end;")
        blocks.append("")  # Empty line
        blocks.append("quit;")  # Ensure MrBayes exits
        
        return "\n".join(blocks)
    
    def _run_mrbayes(self, nexus_file, output_prefix):
        """
        Execute MrBayes and return the marginal likelihood.
        
        Args:
            nexus_file: Path to NEXUS file with data and MrBayes block
            output_prefix: Prefix for output files
            
        Returns:
            Marginal likelihood value or None if failed
        """
        try:
            # MrBayes needs absolute paths and proper quoting for paths with spaces
            # We'll use a relative path instead to avoid issues with spaces
            nexus_filename = nexus_file.name
            
            # Build MrBayes command
            if self.use_mpi:
                # For MPI, determine number of processors
                n_procs = self.mpi_processors
                if n_procs is None:
                    # Default: one processor per chain
                    n_procs = self.bayes_chains
                cmd = [self.mpirun_path, "-np", str(n_procs), self.mrbayes_path, nexus_filename]
            else:
                # Standard MrBayes command
                cmd = [self.mrbayes_path, nexus_filename]
            
            logger.debug(f"Running MrBayes: {' '.join(cmd)} in directory {self.temp_path}")
            logger.debug(f"Working directory: {self.temp_path}")
            logger.debug(f"NEXUS file: {nexus_file}")
            
            # Run MrBayes
            # More realistic timeout: assume ~1000 generations/second, multiply by safety factor
            timeout_seconds = max(7200, (self.bayes_ngen / 500) * self.bayes_chains)
            logger.debug(f"MrBayes timeout set to {timeout_seconds} seconds")
            
            result = subprocess.run(cmd, cwd=str(self.temp_path), 
                                  capture_output=True, text=True, 
                                  timeout=timeout_seconds)
            
            if result.returncode != 0:
                logger.error(f"MrBayes failed with return code {result.returncode}")
                logger.error(f"MrBayes stdout: {result.stdout[:500]}")  # First 500 chars
                logger.error(f"MrBayes stderr: {result.stderr[:500]}")  # First 500 chars
                # Check for specific error patterns
                if "Error" in result.stdout or "Could not" in result.stdout:
                    error_lines = [line for line in result.stdout.split('\n') if 'Error' in line or 'Could not' in line]
                    for line in error_lines[:5]:
                        logger.error(f"MrBayes error: {line}")
                return None
            
            # Log successful completion at debug level to avoid redundancy
            logger.debug(f"MrBayes completed successfully for {output_prefix}")
            
            # Parse marginal likelihood from output
            ml_value = None
            
            # Use stepping-stone if requested
            if self.marginal_likelihood == "ss":
                logger.debug(f"Looking for stepping-stone output for {output_prefix}")
                ml_value = self._parse_mrbayes_stepping_stone(nexus_file, output_prefix)
                
                # Fall back to harmonic mean if stepping-stone failed
                if ml_value is None:
                    logger.warning("Stepping-stone parsing failed, falling back to harmonic mean")
                    lstat_file_path = self.temp_path / f"{nexus_file.name}.lstat"
                    ml_value = self._parse_mrbayes_marginal_likelihood(lstat_file_path, output_prefix)
            else:
                # Use harmonic mean from .lstat file
                lstat_file_path = self.temp_path / f"{nexus_file.name}.lstat"
                logger.debug(f"Looking for MrBayes lstat file: {lstat_file_path}")
                ml_value = self._parse_mrbayes_marginal_likelihood(lstat_file_path, output_prefix)
            
            # Parse posterior probabilities from consensus tree (only for unconstrained analysis)
            if output_prefix == "unc":
                con_tree_path = self.temp_path / f"{nexus_file.name}.con.tre"
                if con_tree_path.exists():
                    logger.debug(f"Parsing posterior probabilities from {con_tree_path}")
                    self.posterior_probs = self._parse_mrbayes_posterior_probs(con_tree_path)
                else:
                    logger.warning(f"Consensus tree not found: {con_tree_path}")
            
            # Check convergence diagnostics
            convergence_data = self._parse_mrbayes_convergence_diagnostics(nexus_file, output_prefix)
            if not self._check_mrbayes_convergence(convergence_data, output_prefix):
                # If strict mode and convergence failed, return None
                if self.convergence_strict:
                    return None
            
            # Store convergence data for reporting
            if not hasattr(self, 'convergence_diagnostics'):
                self.convergence_diagnostics = {}
            self.convergence_diagnostics[output_prefix] = convergence_data
            
            return ml_value
            
        except subprocess.TimeoutExpired:
            logger.error(f"MrBayes timed out for {nexus_file}")
            return None
        except Exception as e:
            logger.error(f"Error running MrBayes: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return None
    
    def run_bayesian_decay_analysis(self):
        """
        Run Bayesian decay analysis for all clades identified in the ML tree.
        
        Returns:
            Dictionary mapping clade_id to Bayesian decay metrics
        """
        if not self.ml_tree:
            logger.error("No ML tree available. Build ML tree first to identify clades.")
            return {}
            
        if not self.bayesian_software:
            logger.error("No Bayesian software specified.")
            return {}
            
        logger.info(f"Running Bayesian decay analysis using {self.bayesian_software}")
        logger.info(f"MCMC settings: {self.bayes_ngen} generations, {self.bayes_chains} chains, sampling every {self.bayes_sample_freq}")
        
        # Log parallel processing settings
        if self.use_mpi:
            n_procs = self.mpi_processors if self.mpi_processors else self.bayes_chains
            logger.info(f"MPI enabled: Using {n_procs} processors with {self.mpirun_path}")
        if self.use_beagle:
            logger.info(f"BEAGLE enabled: device={self.beagle_device}, precision={self.beagle_precision}, scaling={self.beagle_scaling}")
        
        # Warn if this will take a long time
        estimated_time = (self.bayes_ngen * self.bayes_chains) / 10000  # Very rough estimate
        if estimated_time > 60:
            logger.warning(f"Bayesian analysis may take a long time (~{estimated_time:.0f} seconds per run)")
        
        # First, run unconstrained Bayesian analysis
        box_content = [
            "Unconstrained analysis (baseline)",
            "Running MrBayes without constraints",
            f"{self.bayes_ngen:,} generations, {self.bayes_chains} chains"
        ]
        logger.info(self._format_progress_box("Bayesian Analysis", box_content))
        
        # Create NEXUS file with MrBayes block
        nexus_content = self.nexus_file_path.read_text()
        
        # Remove PAUP-specific 'options' commands that MrBayes doesn't understand
        nexus_lines = nexus_content.split('\n')
        filtered_lines = []
        for line in nexus_lines:
            # Skip lines that contain 'options' command (case-insensitive)
            if line.strip().lower().startswith('options'):
                logger.debug(f"Filtering out PAUP-specific line for MrBayes: {line.strip()}")
                continue
            filtered_lines.append(line)
        
        filtered_nexus = '\n'.join(filtered_lines)
        mrbayes_block = self._generate_mrbayes_nexus()
        combined_nexus = filtered_nexus + "\n" + mrbayes_block
        
        unconstrained_nexus = self.temp_path / "unc.nex"
        unconstrained_nexus.write_text(combined_nexus)
        
        # Run unconstrained analysis
        unconstrained_ml = self._run_mrbayes(unconstrained_nexus, "unc")
        
        if unconstrained_ml is None:
            logger.error("Unconstrained Bayesian analysis failed")
            return {}
            
        logger.info(f"Unconstrained marginal likelihood: {unconstrained_ml}")
        
        # Now run constrained analyses for each clade
        bayesian_results = {}
        
        # Get all clades from ML tree (same as in ML analysis)
        internal_clades = [cl for cl in self.ml_tree.get_nonterminals() if cl and cl.clades]
        
        # Parse user constraints if constraint mode is not "all"
        user_constraints = []
        if self.constraint_mode != "all":
            user_constraints = self.parse_constraints()
            if not user_constraints and self.constraint_mode == "specific":
                logger.warning("Constraint mode is 'specific' but no constraints were provided. No branches will be tested.")
                return {}
            logger.info(f"Parsed {len(user_constraints)} user-defined constraints for Bayesian analysis")
        
        # Count testable branches (same logic as ML analysis)
        testable_branches = []
        for i, clade_obj in enumerate(internal_clades):
            clade_log_idx = i + 1
            clade_taxa = [leaf.name for leaf in clade_obj.get_terminals()]
            total_taxa_count = len(self.ml_tree.get_terminals())
            
            if len(clade_taxa) <= 1 or len(clade_taxa) >= total_taxa_count - 1:
                continue
            if not self.should_test_clade(clade_taxa, user_constraints):
                continue
            testable_branches.append((i, clade_obj, clade_log_idx, clade_taxa))
        
        logger.info(f"Testing {len(testable_branches)} branches for Bayesian decay...")
        
        for branch_num, (i, clade_obj, clade_log_idx, clade_taxa) in enumerate(testable_branches, 1):
            clade_id = f"Clade_{clade_log_idx}"
            
            # Display progress box
            # Show first few taxa names for context
            taxa_sample = ", ".join(clade_taxa[:2])
            if len(clade_taxa) > 2:
                taxa_sample += "..."
            
            box_content = [
                f"Testing clade {clade_log_idx} ({branch_num} of {len(testable_branches)}) • {len(clade_taxa)} taxa",
                f"Testing constraint on: {taxa_sample}",
                "Running MrBayes with negative constraint",
                f"{self.bayes_ngen:,} generations, {self.bayes_chains} chains"
            ]
            logger.info(self._format_progress_box("Bayesian Analysis", box_content))
            
            # Create constrained NEXUS file
            mrbayes_block = self._generate_mrbayes_nexus(
                clade_taxa=clade_taxa, 
                constraint_id=clade_id
            )
            combined_nexus = filtered_nexus + "\n" + mrbayes_block
            
            constrained_nexus = self.temp_path / f"c_{clade_log_idx}.nex"
            constrained_nexus.write_text(combined_nexus)
            
            # Debug: save first constraint file for inspection
            if clade_log_idx == 3 and self.debug:
                debug_copy = self.temp_path.parent / "debug_mrbayes_constraint.nex"
                debug_copy.write_text(combined_nexus)
                logger.info(f"Debug: Saved constraint file to {debug_copy}")
            
            # Run constrained analysis
            constrained_ml = self._run_mrbayes(constrained_nexus, f"c_{clade_log_idx}")
            
            if constrained_ml is not None:
                # Calculate Bayesian decay (marginal likelihood difference)
                bayes_decay = unconstrained_ml - constrained_ml
                
                # Calculate site data for effect size calculations if normalization is enabled
                site_data = None
                if self.normalize_bd and any(method.startswith('effect_size') or method == 'signal_to_noise' 
                                           for method in self.bd_normalization_methods):
                    # We need to perform site analysis for effect size calculations
                    # Use the ML tree files if available
                    if hasattr(self, 'ml_tree_file') and self.ml_tree_file:
                        # Check if we have a constraint tree file for this clade
                        constraint_tree_file = self.temp_path / f"constrained_tree_{clade_log_idx}.tre"
                        if constraint_tree_file.exists():
                            tree_files = [str(self.ml_tree_file), str(constraint_tree_file)]
                            site_data = self._calculate_site_likelihoods(tree_files, clade_id)
                            if site_data:
                                logger.debug(f"Site analysis completed for {clade_id} (Bayesian normalization)")
                
                # Calculate normalized BD metrics if enabled
                normalized_bd_metrics = {}
                if self.normalize_bd:
                    normalized_bd_metrics = self._calculate_normalized_bd_metrics(
                        bayes_decay, unconstrained_ml, site_data=site_data
                    )
                
                bayesian_results[clade_id] = {
                    'unconstrained_ml': unconstrained_ml,
                    'constrained_ml': constrained_ml,
                    'bayes_decay': bayes_decay,
                    'taxa': clade_taxa,
                    **normalized_bd_metrics  # Add normalized metrics
                }
                
                # Update progress box with results
                result_lines = [
                    f"MrBayes completed successfully",
                    f"Marginal likelihood (stepping-stone): {constrained_ml:.3f}",
                    f"Bayes Decay: {bayes_decay:.2f}"
                ]
                
                # Add normalized metrics to progress display if available
                if self.normalize_bd and normalized_bd_metrics:
                    if 'bd_per_site' in normalized_bd_metrics:
                        bd_ps = normalized_bd_metrics['bd_per_site']
                        result_lines.append(f"BD/site: {bd_ps:.6f}")
                    if 'bd_relative' in normalized_bd_metrics:
                        bd_rel = normalized_bd_metrics['bd_relative']
                        result_lines.append(f"BD%: {bd_rel*100:.3f}%")
                logger.info(self._format_progress_box("Bayesian Analysis Results", result_lines))
                if bayes_decay < 0:
                    logger.warning(f"⚠️  {clade_id} has negative Bayes Decay ({bayes_decay:.4f}), suggesting potential convergence or estimation issues")
            else:
                logger.warning(f"Constrained analysis failed for {clade_id}")
        
        # Check for negative Bayes Decay values and issue summary warning
        if bayesian_results:
            negative_clades = [cid for cid, data in bayesian_results.items() if data['bayes_decay'] < 0]
            if negative_clades:
                logger.warning(f"\n⚠️  WARNING: {len(negative_clades)}/{len(bayesian_results)} clades have negative Bayes Decay values!")
                logger.warning("This suggests potential issues with MCMC convergence or marginal likelihood estimation.")
                logger.warning("Consider:")
                logger.warning("  1. Increasing MCMC generations (--bayes-ngen 5000000 or higher)")
                logger.warning("  2. Using more chains (--bayes-chains 8)")
                logger.warning("  3. Checking MCMC convergence diagnostics in MrBayes output")
                logger.warning("  4. Verifying chain convergence (check .stat files for ESS values)")
                logger.warning(f"Affected clades: {', '.join(negative_clades)}\n")
                
        return bayesian_results
    
    def _parse_mrbayes_marginal_likelihood(self, lstat_file, output_prefix):
        """
        Parse marginal likelihood from MrBayes .lstat file.
        
        Args:
            lstat_file: Path to MrBayes .lstat file
            output_prefix: Prefix to identify the run
            
        Returns:
            Marginal likelihood value or None
        """
        if not lstat_file.exists():
            logger.warning(f"MrBayes lstat file not found: {lstat_file}")
            return None
            
        try:
            lstat_content = lstat_file.read_text()
            
            # Parse the .lstat file format
            # Format: run  arithmetic_mean  harmonic_mean  values_discarded
            # We want the harmonic mean from the "all" row
            lines = lstat_content.strip().split('\n')
            
            for line in lines:
                if line.startswith('all'):
                    parts = line.split()
                    if len(parts) >= 3:
                        # harmonic_mean is the third column
                        harmonic_mean = float(parts[2])
                        logger.info(f"Parsed harmonic mean marginal likelihood for {output_prefix}: {harmonic_mean}")
                        return harmonic_mean
            
            # If no "all" row, try to get from individual runs
            for line in lines:
                if line[0].isdigit():  # Run number
                    parts = line.split()
                    if len(parts) >= 3:
                        harmonic_mean = float(parts[2])
                        logger.info(f"Parsed harmonic mean marginal likelihood for {output_prefix}: {harmonic_mean}")
                        return harmonic_mean
                
            logger.warning(f"Could not find marginal likelihood in {lstat_file}")
            return None
            
        except Exception as e:
            logger.error(f"Error parsing MrBayes lstat file: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return None

    def _parse_mrbayes_stepping_stone(self, nexus_file_path, output_prefix):
        """
        Parse stepping-stone marginal likelihood from MrBayes .ss output file.
        
        Args:
            nexus_file_path: Path to the NEXUS file (base name for output files)
            output_prefix: Prefix to identify the run
            
        Returns:
            Marginal likelihood value or None
        """
        # MrBayes creates .ss file with stepping-stone results
        ss_file_path = self.temp_path / f"{nexus_file_path.name}.ss"
        
        if not ss_file_path.exists():
            logger.warning(f"MrBayes stepping-stone file not found: {ss_file_path}")
            return None
            
        try:
            ss_content = ss_file_path.read_text()
            
            # Parse the stepping-stone output
            # Look for line with "Marginal likelihood (ln)"
            for line in ss_content.splitlines():
                if "Marginal likelihood (ln)" in line:
                    # Extract the value - format varies slightly between MrBayes versions
                    # Common formats:
                    # "Marginal likelihood (ln) = -1234.56"
                    # "Marginal likelihood (ln):     -1234.56"
                    parts = line.split("=")
                    if len(parts) < 2:
                        parts = line.split(":")
                    
                    if len(parts) >= 2:
                        try:
                            ml_value = float(parts[-1].strip())
                            logger.info(f"Parsed stepping-stone marginal likelihood for {output_prefix}: {ml_value}")
                            return ml_value
                        except ValueError:
                            logger.warning(f"Could not parse ML value from line: {line}")
            
            # Alternative: Calculate marginal likelihood from step contributions
            # The .ss file contains contributions for each step
            # We need to sum the contributions for each run
            
            # Parse the step contributions table
            run1_sum = 0.0
            run2_sum = 0.0
            found_data = False
            
            for line in ss_content.splitlines():
                # Skip header lines
                if line.startswith('[') or line.startswith('Step'):
                    continue
                    
                parts = line.split()
                if len(parts) >= 5:  # Step, Power, run1, run2, aSplit0
                    try:
                        step = int(parts[0])
                        run1_contrib = float(parts[2])
                        run2_contrib = float(parts[3])
                        
                        run1_sum += run1_contrib
                        run2_sum += run2_contrib
                        found_data = True
                    except (ValueError, IndexError):
                        continue
            
            if found_data:
                # Average the marginal likelihoods from both runs
                ml_value = (run1_sum + run2_sum) / 2.0
                logger.debug(f"Calculated stepping-stone marginal likelihood for {output_prefix}: {ml_value}")
                logger.debug(f"Run1 ML: {run1_sum}, Run2 ML: {run2_sum}")
                return ml_value
            
            logger.warning(f"Could not calculate stepping-stone marginal likelihood from {ss_file_path}")
            if self.debug:
                logger.debug(f"Stepping-stone file content:\n{ss_content[:500]}")
            return None
            
        except Exception as e:
            logger.error(f"Error parsing MrBayes stepping-stone file: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return None

    def _parse_mrbayes_posterior_probs(self, con_tree_path):
        """
        Parse posterior probabilities from MrBayes consensus tree file.
        
        Args:
            con_tree_path: Path to .con.tre file
            
        Returns:
            Dictionary mapping clade (as frozenset of taxa) to posterior probability
        """
        try:
            logger.debug(f"Parsing MrBayes consensus tree from: {con_tree_path}")
            
            if not con_tree_path.exists():
                logger.warning(f"MrBayes consensus tree file not found: {con_tree_path}")
                return {}
            
            posterior_probs = {}
            
            # Read the consensus tree file
            tree_content = con_tree_path.read_text()
            logger.debug(f"Consensus tree file has {len(tree_content)} characters")
            
            # Find the tree line (starts with "tree con_50")
            tree_line = None
            for line in tree_content.splitlines():
                if line.strip().startswith("tree con_"):
                    tree_line = line
                    break
            
            if not tree_line:
                logger.warning("Could not find consensus tree line in .con.tre file")
                return {}
            
            # Extract the Newick string
            # Format: tree con_50 = [&U] (taxon1:0.1,taxon2:0.1)[1.00]:0.0;
            parts = tree_line.split("=", 1)
            if len(parts) < 2:
                logger.warning("Invalid consensus tree format")
                return {}
            
            newick_str = parts[1].strip()
            
            # Debug: log the original newick string
            if self.debug or True:  # Always log for debugging
                logger.debug(f"Original newick string: {newick_str[:300]}...")
            
            # Remove [&U] or other tree attributes at the beginning
            if newick_str.startswith("["):
                end_bracket = newick_str.find("]")
                if end_bracket != -1:
                    newick_str = newick_str[end_bracket+1:].strip()
            
            logger.debug(f"Newick after removing tree attributes: {newick_str[:200]}...")
            
            # MrBayes puts posterior probabilities in square brackets after clades
            # We need to convert these to BioPython confidence values
            # First, let's parse the tree without the posterior probabilities
            
            # Create a version without posterior probabilities for BioPython
            import re
            # Pattern to match posterior probabilities
            # MrBayes extended format: [&prob=1.00000000e+00,...]
            prob_pattern = r'\[&prob=([0-9.eE+-]+)[,\]]'
            
            # Extract posterior probabilities before removing them
            prob_matches = list(re.finditer(prob_pattern, newick_str))
            logger.info(f"Found {len(prob_matches)} posterior probability annotations in consensus tree")
            
            # Also check for other possible formats
            if len(prob_matches) == 0:
                # Check if posterior probs might be in a different format
                alt_patterns = [
                    r'\)(\d+\.\d+):',  # )0.95:
                    r'\)(\d+):',       # )95:
                    r'\{(\d+\.\d+)\}', # {0.95}
                ]
                for pattern in alt_patterns:
                    alt_matches = list(re.finditer(pattern, newick_str))
                    if alt_matches:
                        logger.info(f"Found {len(alt_matches)} matches with alternative pattern: {pattern}")
            
            # Remove the posterior probabilities for BioPython parsing
            clean_newick = re.sub(r'\[[0-9.]+\]', '', newick_str)
            logger.debug(f"Clean newick for BioPython: {clean_newick[:200]}...")
            
            # Parse using BioPython
            from io import StringIO
            tree_io = StringIO(clean_newick)
            
            try:
                tree = Phylo.read(tree_io, "newick")
                logger.debug(f"Successfully parsed tree with {len(list(tree.get_nonterminals()))} internal nodes")
                
                # Now we need to match the posterior probabilities to the clades
                # This is tricky because we need to traverse the tree in the same order
                # as the Newick string
                
                # Use regex to extract clades and their posterior probabilities
                posterior_probs = self._extract_mrbayes_posterior_probs(newick_str)
                
                logger.debug(f"Extracted posterior probabilities for {len(posterior_probs)} clades")
                
            except Exception as e:
                logger.warning(f"Could not parse consensus tree with BioPython: {e}")
                logger.debug("Falling back to extraction method")
                # Fall back to extraction method
                posterior_probs = self._extract_mrbayes_posterior_probs(newick_str)
                
            return posterior_probs
            
        except Exception as e:
            logger.error(f"Error parsing MrBayes posterior probabilities: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return {}

    def _manual_parse_mrbayes_tree(self, newick_str):
        """
        Manually parse MrBayes consensus tree to extract clades and posterior probabilities.
        
        MrBayes format: (taxon1,taxon2)[0.95]
        
        Args:
            newick_str: Newick string with posterior probabilities in square brackets
            
        Returns:
            Dictionary mapping frozensets of taxa to posterior probabilities
        """
        import re
        
        posterior_probs = {}
        
        # First, let's use a more robust approach to parse the tree
        # We'll recursively parse the tree structure
        
        def parse_clade(s, pos=0):
            """Recursively parse a clade from the Newick string."""
            taxa_in_clade = set()
            
            # Skip whitespace
            while pos < len(s) and s[pos].isspace():
                pos += 1
            
            if pos >= len(s):
                return taxa_in_clade, pos
            
            if s[pos] == '(':
                # This is an internal node
                pos += 1  # Skip '('
                
                # Parse all children
                while pos < len(s) and s[pos] != ')':
                    child_taxa, pos = parse_clade(s, pos)
                    taxa_in_clade.update(child_taxa)
                    
                    # Skip whitespace and commas
                    while pos < len(s) and (s[pos].isspace() or s[pos] == ','):
                        pos += 1
                
                if pos < len(s) and s[pos] == ')':
                    pos += 1  # Skip ')'
                    
                    # Check for posterior probability
                    while pos < len(s) and s[pos].isspace():
                        pos += 1
                    
                    # Skip any metadata in square brackets (but not posterior probability)
                    # We'll handle posterior probability parsing separately
                    while pos < len(s) and s[pos] == '[':
                        bracket_depth = 1
                        pos += 1
                        while pos < len(s) and bracket_depth > 0:
                            if s[pos] == '[':
                                bracket_depth += 1
                            elif s[pos] == ']':
                                bracket_depth -= 1
                            pos += 1
                    
                    # Skip branch length if present
                    while pos < len(s) and s[pos].isspace():
                        pos += 1
                    
                    if pos < len(s) and s[pos] == ':':
                        pos += 1  # Skip ':'
                        # Skip the branch length
                        while pos < len(s) and (s[pos].isdigit() or s[pos] in '.eE-+'):
                            pos += 1
                
            else:
                # This is a leaf node (taxon name)
                taxon_start = pos
                while pos < len(s) and s[pos] not in '(),:;[]' and not s[pos].isspace():
                    pos += 1
                
                taxon = s[taxon_start:pos].strip()
                if taxon:
                    taxa_in_clade.add(taxon)
                    logger.debug(f"Found taxon: {taxon}")
                
                # Skip branch length if present
                while pos < len(s) and s[pos].isspace():
                    pos += 1
                
                if pos < len(s) and s[pos] == ':':
                    pos += 1  # Skip ':'
                    # Skip the branch length
                    while pos < len(s) and (s[pos].isdigit() or s[pos] in '.eE-+'):
                        pos += 1
            
            return taxa_in_clade, pos
        
        # Parse the entire tree
        try:
            all_taxa, _ = parse_clade(newick_str)
            logger.info(f"Parsed tree with {len(all_taxa)} taxa")
            logger.debug(f"All taxa: {sorted(all_taxa)}")
            
            # Filter out the root clade (contains all taxa)
            posterior_probs = {k: v for k, v in posterior_probs.items() if len(k) < len(all_taxa)}
            
            logger.info(f"Manual parsing found {len(posterior_probs)} clades with posterior probabilities")
            for clade, prob in sorted(posterior_probs.items(), key=lambda x: x[1], reverse=True)[:5]:
                logger.debug(f"  Top clade {sorted(clade)}: PP={prob:.3f}")
            
        except Exception as e:
            logger.error(f"Error in manual tree parsing: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
        
        return posterior_probs
    
    def _should_test_clade_wrapper(self, clade_obj, user_constraints):
        """Helper to check if a clade should be tested."""
        clade_taxa = [leaf.name for leaf in clade_obj.get_terminals()]
        total_taxa_count = len(self.ml_tree.get_terminals())
        
        # Skip trivial branches
        if len(clade_taxa) <= 1 or len(clade_taxa) >= total_taxa_count - 1:
            return False
            
        # Check constraint mode
        return self.should_test_clade(clade_taxa, user_constraints)
    
    def _get_box_chars(self):
        """Get box drawing characters based on output style."""
        if self.output_style == "unicode":
            return {
                'h': '─', 'v': '│', 'tl': '╭', 'tr': '╮', 'bl': '╰', 'br': '╯',
                'cross': '┼', 'hdown': '┬', 'hup': '┴', 'vright': '├', 'vleft': '┤',
                'h_thick': '═', 'v_thick': '║', 'tl_thick': '╔', 'tr_thick': '╗',
                'bl_thick': '╚', 'br_thick': '╝'
            }
        elif self.output_style == "ascii":
            return {
                'h': '-', 'v': '|', 'tl': '+', 'tr': '+', 'bl': '+', 'br': '+',
                'cross': '+', 'hdown': '+', 'hup': '+', 'vright': '+', 'vleft': '+',
                'h_thick': '=', 'v_thick': '|', 'tl_thick': '+', 'tr_thick': '+',
                'bl_thick': '+', 'br_thick': '+'
            }
        else:  # minimal
            return None
    
    def _format_table_row(self, values, widths, alignments=None):
        """Format a table row with proper alignment and spacing."""
        if alignments is None:
            alignments = ['<'] * len(values)  # Default left align
        
        formatted = []
        for val, width, align in zip(values, widths, alignments):
            if align == '>':  # Right align
                formatted.append(str(val).rjust(width))
            elif align == '^':  # Center align
                formatted.append(str(val).center(width))
            else:  # Left align
                formatted.append(str(val).ljust(width))
        
        return " │ ".join(formatted) if self.output_style == "unicode" else " | ".join(formatted)
    
    def _format_support_symbol(self, pvalue):
        """Convert p-value to support symbol."""
        if pvalue == 'N/A' or pvalue is None:
            return 'N/A'
        try:
            p = float(pvalue)
            if p < 0.001:
                return '***'
            elif p < 0.01:
                return '**'
            elif p < 0.05:
                return '*'
            else:
                return 'ns'
        except:
            return 'N/A'
    
    def _format_tree_annotation(self, clade_id, annotation_dict, style="compact"):
        """Format tree node annotations based on style."""
        if self.output_style == "minimal":
            # Simple format
            parts = []
            if clade_id:
                parts.append(clade_id)
            for key, val in annotation_dict.items():
                if val is not None:
                    parts.append(f"{key}:{val}")
            return " ".join(parts) if parts else None
        
        if style == "compact":
            # Compact bracket notation: Clade_5[AU=0.023,ΔlnL=4.57,BD=34.02]
            if not annotation_dict:
                return clade_id if clade_id else None
            
            formatted_values = []
            # Order matters for readability
            order = ['AU', 'ΔlnL', 'BD', 'PS', 'ES', 'PD', 'PP', 'BS']
            for key in order:
                if key in annotation_dict and annotation_dict[key] is not None:
                    val = annotation_dict[key]
                    if isinstance(val, float):
                        if key in ['AU', 'PP']:
                            formatted_values.append(f"{key}={val:.3f}")
                        elif key in ['ΔlnL', 'BD', 'ES']:
                            formatted_values.append(f"{key}={val:.2f}")
                        elif key == 'PS':
                            formatted_values.append(f"{key}={val:.6f}")
                        else:
                            formatted_values.append(f"{key}={val}")
                    else:
                        formatted_values.append(f"{key}={val}")
            
            if clade_id and formatted_values:
                return f"{clade_id}[{','.join(formatted_values)}]"
            elif formatted_values:
                return f"[{','.join(formatted_values)}]"
            else:
                return clade_id
        
        elif style == "symbols":
            # Symbol format with separators
            symbols = {
                'AU': '✓', 'ΔlnL': '△', 'BD': '◆', 'PS': '§', 'ES': '※', 
                'PD': '#', 'PP': '●', 'BS': '◯'
            }
            parts = []
            if clade_id:
                parts.append(clade_id)
            
            for key, symbol in symbols.items():
                if key in annotation_dict and annotation_dict[key] is not None:
                    val = annotation_dict[key]
                    parts.append(f"{symbol}{key}:{val}")
            
            return "[" + "|".join(parts) + "]" if parts else None
        
        else:  # original
            # Original pipe-separated format
            parts = []
            if clade_id:
                parts = [clade_id, "-"]
            for key, val in annotation_dict.items():
                if val is not None:
                    parts.append(f"{key}:{val}")
            return " ".join(parts) if len(parts) > 2 else clade_id
    
    def _get_display_path(self, path):
        """Get a display-friendly path (relative if possible, otherwise absolute)."""
        try:
            return str(path.relative_to(Path.cwd()))
        except ValueError:
            return str(path)
    
    def _format_progress_bar(self, current, total, width=30, elapsed_time=None):
        """Format a progress bar with percentage and time estimate."""
        if self.output_style == "minimal":
            return f"{current}/{total}"
        
        percent = current / total if total > 0 else 0
        filled = int(width * percent)
        
        if self.output_style == "unicode":
            bar = "█" * filled + "░" * (width - filled)
        else:
            bar = "#" * filled + "-" * (width - filled)
        
        progress_str = f"[{bar}] {percent*100:.0f}% | {current}/{total}"
        
        if elapsed_time and current > 0:
            avg_time = elapsed_time / current
            remaining = avg_time * (total - current)
            if remaining > 60:
                progress_str += f" | Est. time remaining: {remaining/60:.0f}m"
            else:
                progress_str += f" | Est. time remaining: {remaining:.0f}s"
        
        return progress_str
    
    def _format_progress_box(self, title, content_lines, width=78):
        """Format a progress box with title and content using simple dashed style."""
        if self.output_style == "minimal":
            return f"\n{title}\n" + "\n".join(content_lines) + "\n"
        
        output = []
        
        # Title line with dashes
        output.append(f"--- {title} ---")
        
        # Content lines (no padding needed for simple style)
        for line in content_lines:
            if line == "---":  # Skip separator lines in content
                continue
            output.append(line)
        
        # Bottom dashes (match the longest line)
        max_len = max(len(line) for line in output)
        output.append("-" * max_len)
        
        return "\n".join(output)
    
    def _extract_mrbayes_posterior_probs(self, newick_str):
        """
        Extract posterior probabilities from MrBayes extended NEXUS format.
        This format has [&prob=X.XXXe+00,...] annotations.
        """
        import re
        
        posterior_probs = {}
        
        try:
            # Add timeout protection
            import time
            start_time = time.time()
            max_time = self.mrbayes_parse_timeout if self.mrbayes_parse_timeout > 0 else float('inf')
            
            # Check if this is a large tree and warn
            if len(newick_str) > 1_000_000:  # 1MB
                logger.warning(f"Large consensus tree ({len(newick_str)/1_000_000:.1f}MB), parsing may take time...")
            
            # Remove the trailing semicolon
            newick_str = newick_str.rstrip(';')
            
            # Pattern to extract taxon names - they appear before [ or : or , or )
            taxon_pattern = r'[\(,]([^,\(\)\[\]:]+?)(?=\[|:|,|\))'
            
            # Pattern to extract prob values from annotations
            prob_value_pattern = r'&prob=([0-9.eE+-]+)'
            
            # First, extract all taxa names to identify terminals vs internals
            all_taxa = set()
            for match in re.finditer(taxon_pattern, newick_str):
                taxon = match.group(1).strip()
                if taxon:
                    all_taxa.add(taxon)
            
            logger.debug(f"Found {len(all_taxa)} taxa in tree")
            
            # Now parse clades by tracking parentheses and their prob values
            # Strategy: find each closing ) followed by [&prob=...]
            clade_pattern = r'\)(\[&[^\]]+\])?'
            
            # Track position in string and parse clades
            pos = 0
            clade_stack = []  # Stack of sets of taxa
            nodes_processed = 0
            
            i = 0
            while i < len(newick_str):
                # Check timeout
                if time.time() - start_time > max_time:
                    logger.warning(f"Posterior probability extraction timed out after {max_time}s")
                    break
                
                # Progress logging
                if nodes_processed > 0 and nodes_processed % 1000 == 0:
                    elapsed = time.time() - start_time
                    logger.info(f"  Processed {nodes_processed} nodes in {elapsed:.1f}s...")
                    
                char = newick_str[i]
                
                if char == '(':
                    # Start of a new clade
                    clade_stack.append(set())
                    i += 1
                    
                elif char == ')':
                    # End of a clade - check for posterior probability
                    nodes_processed += 1
                    if clade_stack:
                        current_clade = clade_stack.pop()
                        
                        # Look ahead for [&prob=...]
                        j = i + 1
                        while j < len(newick_str) and newick_str[j].isspace():
                            j += 1
                            
                        if j < len(newick_str) and newick_str[j] == '[':
                            # Find the closing ]
                            k = j + 1
                            bracket_depth = 1
                            while k < len(newick_str) and bracket_depth > 0:
                                if newick_str[k] == '[':
                                    bracket_depth += 1
                                elif newick_str[k] == ']':
                                    bracket_depth -= 1
                                k += 1
                            
                            if k <= len(newick_str):
                                annotation = newick_str[j:k]
                                # Extract prob value
                                prob_match = re.search(prob_value_pattern, annotation)
                                if prob_match and len(current_clade) > 1:
                                    # Only store multi-taxa clades (not terminals)
                                    prob_value = float(prob_match.group(1))
                                    clade_key = frozenset(current_clade)
                                    posterior_probs[clade_key] = prob_value
                                    
                                i = k  # Skip past the annotation
                                continue
                        
                        # Add this clade's taxa to parent clade if any
                        if clade_stack:
                            clade_stack[-1].update(current_clade)
                    
                    i += 1
                    
                elif char not in '[]():,':
                    # Possible start of a taxon name
                    taxon_start = i
                    while i < len(newick_str) and newick_str[i] not in '[]():,':
                        i += 1
                    
                    taxon = newick_str[taxon_start:i].strip()
                    if taxon and not taxon[0].isdigit():  # Skip branch lengths
                        # Add to current clade
                        if clade_stack:
                            clade_stack[-1].add(taxon)
                        
                        # Skip any following annotations
                        while i < len(newick_str) and newick_str[i] == '[':
                            # Skip annotation
                            bracket_depth = 1
                            i += 1
                            while i < len(newick_str) and bracket_depth > 0:
                                if newick_str[i] == '[':
                                    bracket_depth += 1
                                elif newick_str[i] == ']':
                                    bracket_depth -= 1
                                i += 1
                else:
                    i += 1
            
            logger.debug(f"Extracted posterior probabilities for {len(posterior_probs)} clades")
            
            # Debug: show some extracted values
            if self.debug and posterior_probs:
                for clade, prob in list(posterior_probs.items())[:3]:
                    taxa_list = sorted(list(clade))[:3]
                    logger.debug(f"  Clade {','.join(taxa_list)}{'...' if len(clade) > 3 else ''}: PP={prob:.3f}")
            
        except Exception as e:
            logger.warning(f"Failed to extract posterior probabilities: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
        
        return posterior_probs

    def _parse_mrbayes_convergence_diagnostics(self, nexus_file_path, output_prefix):
        """
        Parse convergence diagnostics from MrBayes output files.
        
        Args:
            nexus_file_path: Path to the MrBayes nexus file
            output_prefix: Output file prefix (e.g., 'unc', 'c_3')
            
        Returns:
            Dictionary with convergence metrics:
            {
                'min_ess': minimum ESS across all parameters,
                'max_psrf': maximum PSRF across all parameters,
                'asdsf': average standard deviation of split frequencies,
                'converged': boolean indicating if convergence criteria met,
                'warnings': list of warning messages
            }
        """
        convergence_data = {
            'min_ess': None,
            'max_psrf': None,
            'asdsf': None,
            'converged': False,
            'warnings': []
        }
        
        try:
            # Parse .pstat file for ESS and PSRF
            pstat_file = self.temp_path / f"{nexus_file_path.name}.pstat"
            if pstat_file.exists():
                pstat_content = pstat_file.read_text()
                
                # Parse ESS and PSRF values
                ess_values = []
                psrf_values = []
                
                for line in pstat_content.splitlines():
                    if line.startswith("#") or not line.strip():
                        continue
                    
                    parts = line.split()
                    if len(parts) >= 4:  # Parameter, Mean, Variance, ESS, PSRF
                        try:
                            # ESS is typically in column 3 or 4 (0-indexed)
                            if len(parts) > 3 and parts[3].replace('.','').replace('-','').isdigit():
                                ess = float(parts[3])
                                ess_values.append(ess)
                            
                            # PSRF is typically in column 4 or 5
                            if len(parts) > 4 and parts[4].replace('.','').replace('-','').isdigit():
                                psrf = float(parts[4])
                                psrf_values.append(psrf)
                        except (ValueError, IndexError):
                            continue
                
                if ess_values:
                    convergence_data['min_ess'] = min(ess_values)
                    if convergence_data['min_ess'] < self.min_ess:
                        convergence_data['warnings'].append(
                            f"Low ESS detected: {convergence_data['min_ess']:.0f} < {self.min_ess}"
                        )
                
                if psrf_values:
                    convergence_data['max_psrf'] = max(psrf_values)
                    if convergence_data['max_psrf'] > self.max_psrf:
                        convergence_data['warnings'].append(
                            f"High PSRF detected: {convergence_data['max_psrf']:.3f} > {self.max_psrf}"
                        )
            
            # Parse .mcmc file or stdout for ASDSF
            mcmc_file = self.temp_path / f"{nexus_file_path.name}.mcmc"
            if mcmc_file.exists():
                mcmc_content = mcmc_file.read_text()
                
                # Look for ASDSF in the last few lines
                for line in reversed(mcmc_content.splitlines()[-20:]):
                    if "Average standard deviation of split frequencies:" in line:
                        try:
                            asdsf_str = line.split(":")[-1].strip()
                            convergence_data['asdsf'] = float(asdsf_str)
                            if convergence_data['asdsf'] > self.max_asdsf:
                                convergence_data['warnings'].append(
                                    f"High ASDSF: {convergence_data['asdsf']:.6f} > {self.max_asdsf}"
                                )
                            break
                        except ValueError:
                            pass
            
            # Determine overall convergence
            convergence_data['converged'] = (
                (convergence_data['min_ess'] is None or convergence_data['min_ess'] >= self.min_ess) and
                (convergence_data['max_psrf'] is None or convergence_data['max_psrf'] <= self.max_psrf) and
                (convergence_data['asdsf'] is None or convergence_data['asdsf'] <= self.max_asdsf)
            )
            
            return convergence_data
            
        except Exception as e:
            logger.error(f"Error parsing convergence diagnostics: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return convergence_data
    
    def _check_mrbayes_convergence(self, convergence_data, output_prefix):
        """
        Check convergence criteria and log appropriate warnings.
        
        Args:
            convergence_data: Dictionary from _parse_mrbayes_convergence_diagnostics
            output_prefix: Run identifier for logging
            
        Returns:
            Boolean indicating if run should be considered valid
        """
        if not self.check_convergence:
            return True
        
        # Log convergence metrics at debug level (main info shown in results box)
        logger.debug(f"Convergence diagnostics for {output_prefix}:")
        if convergence_data['min_ess'] is not None:
            logger.debug(f"  Minimum ESS: {convergence_data['min_ess']:.0f}")
        if convergence_data['max_psrf'] is not None:
            logger.debug(f"  Maximum PSRF: {convergence_data['max_psrf']:.3f}")
        if convergence_data['asdsf'] is not None:
            logger.debug(f"  Final ASDSF: {convergence_data['asdsf']:.6f}")
        
        # Log warnings at warning level
        for warning in convergence_data['warnings']:
            logger.warning(f"  ⚠️  {warning}")
        
        # If strict mode and not converged, treat as failure
        if self.convergence_strict and not convergence_data['converged']:
            logger.error(f"Convergence criteria not met for {output_prefix} (strict mode enabled)")
            return False
        
        # Otherwise just warn
        if not convergence_data['converged']:
            logger.warning(f"Convergence criteria not met for {output_prefix}, but continuing (strict mode disabled)")
            logger.warning("Consider increasing --bayes-ngen or adjusting MCMC parameters")
        
        return True

    def _identify_testable_branches(self):
        """
        Identify all internal branches in the tree that can be tested.
        Returns a list of clade objects.
        """
        if not self.ml_tree:
            logger.error("No tree available for branch identification")
            return []
        
        internal_clades = [cl for cl in self.ml_tree.get_nonterminals() if cl and cl.clades]
        return internal_clades
    
    def run_parsimony_decay_analysis(self):
        """
        Run parsimony analysis to calculate traditional Bremer support values.
        
        Returns:
            Dictionary mapping clade IDs to parsimony decay values
        """
        logger.info("Running parsimony decay analysis...")
        
        # Build parsimony tree if not already done
        PARS_TREE_FN = "pars_tree.tre"
        PARS_SCORE_FN = "pars_score.txt"
        PARS_NEX_FN = "pars_search.nex"
        PARS_LOG_FN = "paup_pars.log"
        
        # Build initial parsimony tree
        script_cmds = [
            f"execute {NEXUS_ALIGNMENT_FN};",
            f"set criterion=parsimony;",
            f"hsearch start=stepwise addseq=random nreps=10 swap=tbr multrees=yes;",
            f"savetrees file={PARS_TREE_FN} replace=yes;",
            f"pscores 1 / scorefile={PARS_SCORE_FN} replace=yes;"
        ]
        
        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        pars_cmd_path = self.temp_path / PARS_NEX_FN
        pars_cmd_path.write_text(paup_script_content)
        
        # Initialize original_tree at the start to avoid scope issues
        original_tree = self.ml_tree
        
        try:
            self._run_paup_command_file(PARS_NEX_FN, PARS_LOG_FN)
            
            # Parse parsimony score
            score_path = self.temp_path / PARS_SCORE_FN
            pars_score = None
            if score_path.exists():
                score_content = score_path.read_text()
                logger.debug(f"Parsimony score file content:\n{score_content}")
                # Parse parsimony score from PAUP output
                # Try different patterns that PAUP might use
                for line in score_content.splitlines():
                    # Pattern 1: "Length = 123"
                    if "Length" in line and "=" in line:
                        try:
                            score_str = line.split("=")[1].strip().split()[0]
                            pars_score = int(score_str)
                            break
                        except (ValueError, IndexError):
                            pass
                    # Pattern 2: "Tree    Length"
                    # Next line: "1       123"
                    elif line.strip() and line.split()[0] == "1":
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            try:
                                pars_score = int(parts[1])
                                break
                            except ValueError:
                                pass
                    # Pattern 3: Original pattern "Length 123"
                    elif "Length" in line:
                        parts = line.strip().split()
                        for i, part in enumerate(parts):
                            if part == "Length" and i+1 < len(parts):
                                try:
                                    pars_score = int(parts[i+1])
                                    break
                                except (ValueError, IndexError):
                                    continue
                        if pars_score:
                            break
            
            if pars_score is None:
                logger.error("Could not parse parsimony score")
                return {}
                
            logger.info(f"Initial parsimony score: {pars_score}")
            
            # Now calculate decay for each clade
            parsimony_decay = {}
            
            # Load the parsimony tree for clade identification
            pars_tree_path = self.temp_path / PARS_TREE_FN
            if not pars_tree_path.exists():
                logger.error("Parsimony tree file not found")
                return {}
            
            # Temporarily use parsimony tree for clade identification if no ML tree
            # (original_tree already initialized above)
            if not self.ml_tree:
                try:
                    self.ml_tree = Phylo.read(str(pars_tree_path), 'newick')
                    logger.info("Using parsimony tree for clade identification")
                except Exception as e:
                    logger.error(f"Failed to load parsimony tree: {e}")
                    return {}
            
            branches = self._identify_testable_branches()
            total_taxa_count = len(self.ml_tree.get_terminals())
            
            # Parse user constraints if constraint mode is not "all"
            user_constraints = []
            if self.constraint_mode != "all":
                user_constraints = self.parse_constraints()
                if not user_constraints and self.constraint_mode == "specific":
                    logger.warning("Constraint mode is 'specific' but no constraints were provided. No branches will be tested.")
                    return {}
                logger.info(f"Parsed {len(user_constraints)} user-defined constraints for parsimony analysis")
            
            # Count testable branches (same logic as ML analysis)
            testable_branches = []
            for i, clade_obj in enumerate(branches):
                clade_log_idx = i + 1
                clade_taxa = [leaf.name for leaf in clade_obj.get_terminals()]
                
                if len(clade_taxa) <= 1 or len(clade_taxa) >= total_taxa_count - 1:
                    continue
                if not self.should_test_clade(clade_taxa, user_constraints):
                    continue
                testable_branches.append((i, clade_obj, clade_log_idx, clade_taxa))
            
            logger.info(f"Testing {len(testable_branches)} branches for parsimony decay...")
            
            # Process testable branches
            for branch_num, (i, clade_obj, clade_log_idx, clade_taxa) in enumerate(testable_branches, 1):
                clade_id = f"Clade_{clade_log_idx}"
                
                # Display progress box
                taxa_sample = ", ".join(clade_taxa[:2])
                if len(clade_taxa) > 2:
                    taxa_sample += "..."
                
                box_content = [
                    f"Testing clade {clade_log_idx} ({branch_num} of {len(testable_branches)}) • {len(clade_taxa)} taxa",
                    f"Testing constraint on: {taxa_sample}",
                    "Running PAUP* with converse constraint"
                ]
                logger.info(self._format_progress_box("Parsimony Analysis", box_content))
                
                # Create constraint forcing non-monophyly
                # Format taxa names for PAUP* constraint syntax
                formatted_clade_taxa = [self._format_taxon_for_paup(t) for t in clade_taxa]
                # Use same constraint specification format as ML analysis
                clade_spec = "((" + ", ".join(formatted_clade_taxa) + "));"
                constraint_cmds = [
                    f"execute {NEXUS_ALIGNMENT_FN};",
                    f"set criterion=parsimony;",
                    f"constraint broken_clade (MONOPHYLY) = {clade_spec}",
                    f"hsearch start=stepwise addseq=random nreps=10 swap=tbr multrees=yes enforce=yes converse=yes constraints=broken_clade;",
                    f"savetrees file=pars_constraint_{clade_log_idx}.tre replace=yes;",
                    f"pscores 1 / scorefile=pars_constraint_score_{clade_log_idx}.txt replace=yes;"
                ]
                
                constraint_script = f"#NEXUS\nbegin paup;\n" + "\n".join(constraint_cmds) + "\nquit;\nend;\n"
                constraint_path = self.temp_path / f"pars_constraint_{clade_log_idx}.nex"
                constraint_path.write_text(constraint_script)
                
                try:
                    self._run_paup_command_file(f"pars_constraint_{clade_log_idx}.nex", f"paup_pars_constraint_{clade_log_idx}.log")
                    
                    # Parse constrained score
                    constrained_score_path = self.temp_path / f"pars_constraint_score_{clade_log_idx}.txt"
                    constrained_score = None
                    
                    if constrained_score_path.exists():
                        score_content = constrained_score_path.read_text()
                        logger.debug(f"Constraint {clade_log_idx} parsimony score file content:\n{score_content}")
                        # Use same parsing logic as initial parsimony score
                        for line in score_content.splitlines():
                            # Pattern 1: "Length = 123"
                            if "Length" in line and "=" in line:
                                try:
                                    score_str = line.split("=")[1].strip().split()[0]
                                    constrained_score = int(score_str)
                                    break
                                except (ValueError, IndexError):
                                    pass
                            # Pattern 2: "Tree    Length"
                            # Next line: "1       123"
                            elif line.strip() and line.split()[0] == "1":
                                parts = line.strip().split()
                                if len(parts) >= 2:
                                    try:
                                        constrained_score = int(parts[1])
                                        break
                                    except ValueError:
                                        pass
                            # Pattern 3: Original pattern "Length 123"
                            elif "Length" in line:
                                parts = line.strip().split()
                                for i, part in enumerate(parts):
                                    if part == "Length" and i+1 < len(parts):
                                        try:
                                            constrained_score = int(parts[i+1])
                                            break
                                        except (ValueError, IndexError):
                                            continue
                                if constrained_score:
                                    break
                    
                    if constrained_score is not None:
                        decay_value = constrained_score - pars_score
                        parsimony_decay[clade_id] = {
                            'pars_decay': decay_value,
                            'pars_score': pars_score,
                            'constrained_score': constrained_score,
                            'taxa': clade_taxa
                        }
                        logger.info(f"{clade_id}: Parsimony decay = {decay_value} (constrained: {constrained_score}, optimal: {pars_score})")
                    
                except Exception as e:
                    logger.error(f"Failed to calculate parsimony decay for {clade_id}: {e}")
                    continue
                    
        except Exception as e:
            logger.error(f"Parsimony decay analysis failed: {e}")
            # Restore original tree if we temporarily used parsimony tree
            if original_tree is not None:
                self.ml_tree = original_tree
            return {}
        
        finally:
            # Restore original tree if we temporarily used parsimony tree
            if original_tree is not None:
                self.ml_tree = original_tree
            
        return parsimony_decay

    def parse_constraints(self):
        """Parse constraints from various sources and return a list of taxon sets."""
        constraints = []
        
        # Parse constraints from config file [constraints] section
        if self.config_constraints:
            for key, value in self.config_constraints.items():
                taxa = [t.strip() for t in value.split(',')]
                constraints.append(taxa)
                logger.debug(f"Added constraint from config [{key}]: {taxa}")
        
        # Parse constraints from --test-branches argument
        if self.test_branches:
            if self.test_branches.startswith('@'):
                # Read from file
                constraint_file = self.test_branches[1:]
                try:
                    with open(constraint_file, 'r') as f:
                        for line in f:
                            line = line.strip()
                            if line and not line.startswith('#'):
                                taxa = [t.strip() for t in line.split(',')]
                                constraints.append(taxa)
                                logger.debug(f"Added constraint from file: {taxa}")
                except Exception as e:
                    logger.error(f"Failed to read constraints from file {constraint_file}: {e}")
            else:
                # Parse semicolon-separated clades
                clades = self.test_branches.split(';')
                for clade in clades:
                    taxa = [t.strip() for t in clade.split(',')]
                    constraints.append(taxa)
                    logger.debug(f"Added constraint from command line: {taxa}")
        
        # Parse constraints from --constraint-file
        if self.constraint_file:
            try:
                with open(self.constraint_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            taxa = [t.strip() for t in line.split(',')]
                            constraints.append(taxa)
                            logger.debug(f"Added constraint from constraint file: {taxa}")
            except Exception as e:
                logger.error(f"Failed to read constraint file {self.constraint_file}: {e}")
        
        return constraints
    
    def should_test_clade(self, clade_taxa_names, user_constraints):
        """Determine if a clade should be tested based on constraint mode."""
        if self.constraint_mode == "all":
            return True
        
        # For specific mode, test only if clade matches a constraint
        if self.constraint_mode == "specific":
            for constraint_taxa in user_constraints:
                # Check if the clade taxa match the constraint
                # Allow for subset matching (constraint is subset of clade)
                if set(constraint_taxa).issubset(set(clade_taxa_names)):
                    return True
            return False
        
        # For exclude mode, test only if clade doesn't match any constraint
        if self.constraint_mode == "exclude":
            for constraint_taxa in user_constraints:
                if set(constraint_taxa).issubset(set(clade_taxa_names)):
                    return False
            return True
        
        return True
    
    def calculate_decay_indices(self, perform_site_analysis=False):
        """
        Calculate decay indices based on the selected analysis mode.
        
        Note: ML and parsimony analyses are performed separately even though both use PAUP*.
        This is because:
        1. They use different optimality criteria (likelihood vs. parsimony score)
        2. The optimal tree under one criterion may not be optimal under the other
        3. Tree searches are guided by the active criterion (likelihood-based vs. parsimony-based swapping)
        4. This design allows flexible analysis combinations (ML-only, parsimony-only, etc.)
        
        Future optimization: When both analyses are requested, we could find the constrained
        tree under one criterion and then score it under both criteria in a single PAUP* run.
        
        Args:
            perform_site_analysis: Whether to perform site-specific analysis (ML only)
            
        Returns:
            Dictionary of decay indices
        """
        # For any analysis mode, we need a tree to identify clades
        # Build ML tree if needed for ML analysis or if no tree exists
        if self.do_ml or not self.ml_tree:
            if not self.ml_tree:
                logger.info("Building tree to identify clades...")
                try:
                    self.build_ml_tree()
                except Exception as e:
                    logger.error(f"Failed to build tree: {e}")
                    return {}

        if not self.ml_tree:
            logger.error("Tree is missing. Cannot identify clades for decay analysis.")
            return {}
            
        # Initialize results dictionary
        self.decay_indices = {}
        
        # Calculate ML decay indices if requested
        if self.do_ml:
            if self.ml_likelihood is None:
                logger.error("ML likelihood is missing. Cannot calculate ML decay indices.")
                if not self.do_bayesian and not self.do_parsimony:
                    return {}
                # Continue with other analyses
            else:
                logger.info("Calculating ML decay indices...")
                self.decay_indices = self._calculate_ml_decay_indices(perform_site_analysis)
        
        # Calculate Bayesian decay indices if requested
        if self.do_bayesian:
            logger.info("Calculating Bayesian decay indices...")
            bayesian_results = self._calculate_bayesian_decay_indices()
            
            # Merge Bayesian results with existing results
            if bayesian_results:
                for clade_id, bayes_data in bayesian_results.items():
                    if clade_id in self.decay_indices:
                        # Add Bayesian fields to existing results
                        self.decay_indices[clade_id].update({
                            'bayes_unconstrained_ml': bayes_data.get('unconstrained_ml'),
                            'bayes_constrained_ml': bayes_data.get('constrained_ml'),
                            'bayes_decay': bayes_data.get('bayes_decay'),
                        })
                    else:
                        # Create new entry for Bayesian-only results
                        self.decay_indices[clade_id] = bayes_data
        
        # Calculate parsimony decay indices if requested
        if self.do_parsimony:
            logger.info("Calculating parsimony decay indices...")
            parsimony_results = self.run_parsimony_decay_analysis()
            
            # Merge parsimony results with existing results
            if parsimony_results:
                for clade_id, pars_data in parsimony_results.items():
                    if clade_id in self.decay_indices:
                        # Add parsimony fields to existing results
                        self.decay_indices[clade_id].update({
                            'pars_decay': pars_data.get('pars_decay'),
                            'pars_score': pars_data.get('pars_score'),
                            'pars_constrained_score': pars_data.get('constrained_score')
                        })
                    else:
                        # Create new entry for parsimony-only results
                        self.decay_indices[clade_id] = pars_data
        
        if not self.decay_indices:
            logger.warning("No branch support values were calculated.")
        else:
            logger.info(f"Calculated support values for {len(self.decay_indices)} branches.")
            
        return self.decay_indices
    
    def _calculate_ml_decay_indices(self, perform_site_analysis=False):
        """Calculate ML decay indices for all internal branches of the ML tree."""
        if not self.ml_tree or self.ml_likelihood is None:
            logger.error("ML tree or its likelihood is missing. Cannot calculate ML decay indices.")
            return {}

        logger.info("Calculating branch support (decay indices)...")
        all_tree_files_rel = [ML_TREE_FN] # ML tree is first
        constraint_info_map = {} # Maps clade_id_str to its info

        internal_clades = [cl for cl in self.ml_tree.get_nonterminals() if cl and cl.clades] # Biphasic, non-empty
        logger.info(f"ML tree has {len(internal_clades)} internal branches to test.")
        if not internal_clades:
            logger.warning("ML tree has no testable internal branches. No decay indices calculated.")
            return {}
        
        # Parse user constraints if constraint mode is not "all"
        user_constraints = []
        if self.constraint_mode != "all":
            user_constraints = self.parse_constraints()
            if not user_constraints and self.constraint_mode == "specific":
                logger.warning("Constraint mode is 'specific' but no constraints were provided. No branches will be tested.")
                return {}
            logger.info(f"Parsed {len(user_constraints)} user-defined constraints")

        # Count testable branches
        testable_branches = []
        for i, clade_obj in enumerate(internal_clades):
            clade_log_idx = i + 1
            clade_taxa_names = [leaf.name for leaf in clade_obj.get_terminals()]
            total_taxa_count = len(self.ml_tree.get_terminals())
            
            if len(clade_taxa_names) <= 1 or len(clade_taxa_names) >= total_taxa_count - 1:
                continue
            if not self.should_test_clade(clade_taxa_names, user_constraints):
                continue
            testable_branches.append((i, clade_obj, clade_log_idx, clade_taxa_names))
        
        logger.info(f"Testing {len(testable_branches)} branches for ML decay...")
        
        # Track time for progress estimation
        import time
        start_time = time.time()
        
        for branch_num, (i, clade_obj, clade_log_idx, clade_taxa_names) in enumerate(testable_branches, 1):
            elapsed = time.time() - start_time
            
            # Show progress
            # Display progress box
            taxa_sample = ", ".join(clade_taxa_names[:2])
            if len(clade_taxa_names) > 2:
                taxa_sample += "..."
            
            box_content = [
                f"Testing clade {clade_log_idx} ({branch_num} of {len(testable_branches)}) • {len(clade_taxa_names)} taxa",
                f"Testing constraint on: {taxa_sample}",
                "Running PAUP* with reverse constraint"
            ]
            logger.info(self._format_progress_box("ML Analysis", box_content))
            rel_constr_tree_fn, constr_lnl = self._generate_and_score_constraint_tree(clade_taxa_names, clade_log_idx)

            if rel_constr_tree_fn: # Successfully generated and scored (even if LNL is None)
                all_tree_files_rel.append(rel_constr_tree_fn)
                clade_id_str = f"Clade_{clade_log_idx}"

                lnl_diff = (constr_lnl - self.ml_likelihood) if constr_lnl is not None and self.ml_likelihood is not None else None
                if constr_lnl is None: logger.warning(f"{clade_id_str}: Constrained LNL is None.")

                constraint_info_map[clade_id_str] = {
                    'taxa': clade_taxa_names,
                    'paup_tree_index': len(all_tree_files_rel), # 1-based index for PAUP*
                    'constrained_lnl': constr_lnl,
                    'lnl_diff': lnl_diff,
                    'tree_filename': rel_constr_tree_fn  # Store tree filename for site analysis
                }
            else:
                logger.warning(f"Failed to generate/score constraint tree for branch {clade_log_idx}. It will be excluded.")

        if not constraint_info_map:
            logger.warning("No valid constraint trees were generated. Skipping AU test.")
            self.decay_indices = {}
            return self.decay_indices

        # Perform site-specific likelihood analysis if requested
        if perform_site_analysis:
            logger.info("Performing site-specific likelihood analysis for each branch...")

            for clade_id, cdata in list(constraint_info_map.items()):
                rel_constr_tree_fn = cdata.get('tree_filename')

                if rel_constr_tree_fn:
                    tree_files = [ML_TREE_FN, rel_constr_tree_fn]
                    site_analysis_result = self._calculate_site_likelihoods(tree_files, clade_id)

                    if site_analysis_result:
                        # Store all site analysis data
                        constraint_info_map[clade_id].update(site_analysis_result)

                        # Log key results
                        supporting = site_analysis_result.get('supporting_sites', 0)
                        conflicting = site_analysis_result.get('conflicting_sites', 0)
                        ratio = site_analysis_result.get('support_ratio', 0.0)
                        weighted = site_analysis_result.get('weighted_support_ratio', 0.0)

                        logger.info(f"Branch {clade_id}: {supporting} supporting sites, {conflicting} conflicting sites, ratio: {ratio:.2f}, weighted ratio: {weighted:.2f}")

        logger.info(f"Running AU test on {len(all_tree_files_rel)} trees (1 ML + {len(constraint_info_map)} constrained).")
        au_test_results = self.run_au_test(all_tree_files_rel)

        self.decay_indices = {}
        # Populate with LNL diffs first, then add AU results
        for cid, cdata in constraint_info_map.items():
            self.decay_indices[cid] = {
                'taxa': cdata['taxa'],
                'lnl_diff': cdata['lnl_diff'],
                'constrained_lnl': cdata['constrained_lnl'],
                'AU_pvalue': None,
                'significant_AU': None
            }

            # Add site analysis data if available
            if 'site_data' in cdata:
                # Copy all the site analysis fields
                for key in ['site_data', 'supporting_sites', 'conflicting_sites', 'neutral_sites',
                           'support_ratio', 'sum_supporting_delta', 'sum_conflicting_delta',
                           'weighted_support_ratio']:
                    if key in cdata:
                        self.decay_indices[cid][key] = cdata[key]

        if au_test_results:
            # Update ML likelihood if AU test scored it differently (should be rare)
            if 1 in au_test_results and self.ml_likelihood is not None:
                if abs(au_test_results[1]['lnL'] - self.ml_likelihood) > 1e-3: # Tolerate small diffs
                    logger.info(f"ML likelihood updated from AU test: {self.ml_likelihood} -> {au_test_results[1]['lnL']}")
                    self.ml_likelihood = au_test_results[1]['lnL']
                    # Need to recalculate all lnl_diffs if ML_LNL changed
                    for cid_recalc in self.decay_indices:
                        constr_lnl_recalc = self.decay_indices[cid_recalc]['constrained_lnl']
                        if constr_lnl_recalc is not None:
                            self.decay_indices[cid_recalc]['lnl_diff'] = constr_lnl_recalc - self.ml_likelihood


            for cid, cdata in constraint_info_map.items():
                paup_idx = cdata['paup_tree_index']
                if paup_idx in au_test_results:
                    au_res_for_tree = au_test_results[paup_idx]
                    self.decay_indices[cid]['AU_pvalue'] = au_res_for_tree['AU_pvalue']
                    if au_res_for_tree['AU_pvalue'] is not None:
                        self.decay_indices[cid]['significant_AU'] = au_res_for_tree['AU_pvalue'] < 0.05

                    # Update constrained LNL from AU test if different
                    current_constr_lnl = self.decay_indices[cid]['constrained_lnl']
                    au_constr_lnl = au_res_for_tree['lnL']
                    if current_constr_lnl is None or abs(current_constr_lnl - au_constr_lnl) > 1e-3:
                        if current_constr_lnl is not None: # Log if it changed significantly
                            logger.info(f"Constrained LNL for {cid} updated by AU test: {current_constr_lnl} -> {au_constr_lnl}")
                        self.decay_indices[cid]['constrained_lnl'] = au_constr_lnl
                        if self.ml_likelihood is not None: # Recalculate diff
                            self.decay_indices[cid]['lnl_diff'] = au_constr_lnl - self.ml_likelihood
                else:
                    logger.warning(f"No AU test result for PAUP tree index {paup_idx} (Clade: {cid}).")
        else:
            logger.warning("AU test failed or returned no results. Decay indices will lack AU p-values.")


        if not self.decay_indices:
            logger.warning("No branch support values were calculated.")
        else:
            logger.info(f"Calculated support values for {len(self.decay_indices)} branches.")

        return self.decay_indices
    
    def _calculate_normalized_bd_metrics(self, raw_bd, unconstrained_ml, site_data=None):
        """
        Calculate normalized Bayesian Decay metrics for better cross-study comparisons.
        
        Args:
            raw_bd: Raw Bayesian Decay value (unconstrained_ml - constrained_ml)
            unconstrained_ml: Unconstrained marginal likelihood
            site_data: Optional site-specific data for signal-to-noise calculation
            
        Returns:
            Dictionary with normalized metrics
        """
        normalized_metrics = {}
        
        # Per-site BD (BD / alignment_length)
        if "per_site" in self.bd_normalization_methods:
            bd_per_site = raw_bd / self.alignment_length if self.alignment_length > 0 else 0
            normalized_metrics["bd_per_site"] = bd_per_site
            
        # Relative BD (BD / |unconstrained_ML|)
        if "relative" in self.bd_normalization_methods:
            if unconstrained_ml != 0:
                bd_relative = raw_bd / abs(unconstrained_ml)
                normalized_metrics["bd_relative"] = bd_relative
            else:
                normalized_metrics["bd_relative"] = 0
                
        # Signal-to-noise ratio (if site data available)
        if "signal_to_noise" in self.bd_normalization_methods and site_data:
            # Calculate signal-to-noise as ratio of supporting to conflicting sites
            supporting_sites = site_data.get('supporting_sites', 0)
            conflicting_sites = site_data.get('conflicting_sites', 0)
            
            if conflicting_sites > 0:
                signal_to_noise = supporting_sites / conflicting_sites
            elif supporting_sites > 0:
                signal_to_noise = float('inf')  # All supporting, no conflicting
            else:
                signal_to_noise = 0  # No supporting sites
                
            normalized_metrics["signal_to_noise"] = signal_to_noise
        
        # Effect Size calculations (like Cohen's d) - requires site data
        if site_data:
            signal_std = site_data.get('signal_std', 0)
            signal_mad = site_data.get('signal_mad', 0)
            
            # Standard Effect Size: BD / SD(site_signals) - like Cohen's d
            if "effect_size" in self.bd_normalization_methods:
                if signal_std > 0:
                    effect_size = raw_bd / signal_std
                    normalized_metrics["effect_size"] = effect_size
                else:
                    normalized_metrics["effect_size"] = 0  # No variance, effect size undefined
            
            # Robust Effect Size: BD / MAD(site_signals) - resistant to outliers
            if "effect_size_robust" in self.bd_normalization_methods:
                if signal_mad > 0:
                    effect_size_robust = raw_bd / signal_mad
                    normalized_metrics["effect_size_robust"] = effect_size_robust
                else:
                    normalized_metrics["effect_size_robust"] = 0
            
            # Weighted Effect Size: weighted by site information content
            if "effect_size_weighted" in self.bd_normalization_methods:
                # For now, use MAD-based but could be enhanced with site weights
                if signal_mad > 0:
                    # Placeholder: Could weight by site conservation, information content, etc.
                    effect_size_weighted = raw_bd / signal_mad
                    normalized_metrics["effect_size_weighted"] = effect_size_weighted
                else:
                    normalized_metrics["effect_size_weighted"] = 0
            
        return normalized_metrics
    
    def _get_bd_interpretation_scale(self, bd_per_site):
        """
        Get interpretation scale for per-site BD values.
        
        Args:
            bd_per_site: Per-site Bayesian Decay value
            
        Returns:
            Tuple of (interpretation_level, symbol)
        """
        if bd_per_site is None:
            return ("unknown", "?")
        
        # Empirically reasonable thresholds for per-site BD
        if bd_per_site >= 0.01:
            return ("very_strong", "***")
        elif bd_per_site >= 0.005:
            return ("strong", "**")
        elif bd_per_site >= 0.001:
            return ("moderate", "*")
        else:
            return ("weak", "ns")
    
    def _get_effect_size_interpretation_scale(self, effect_size):
        """
        Get interpretation scale for effect sizes based on Cohen's d conventions.
        
        Args:
            effect_size: Effect size value (BD / SD or BD / MAD)
            
        Returns:
            Tuple of (interpretation_level, symbol, description)
        """
        if effect_size is None:
            return ("unknown", "?", "Unknown")
        
        # Cohen's d-based thresholds adapted for phylogenetics
        abs_effect_size = abs(effect_size)
        if abs_effect_size >= 1.2:
            return ("very_large", "***", "Very large effect")
        elif abs_effect_size >= 0.8:
            return ("large", "**", "Large effect") 
        elif abs_effect_size >= 0.5:
            return ("medium", "*", "Medium effect")
        elif abs_effect_size >= 0.2:
            return ("small", "~", "Small effect")
        else:
            return ("negligible", "ns", "Negligible effect")
    
    def _calculate_bayesian_decay_indices(self):
        """
        Calculate Bayesian decay indices for all internal branches.
        
        Returns:
            Dictionary of decay indices with Bayesian metrics
        """
        # Run Bayesian decay analysis
        bayesian_results = self.run_bayesian_decay_analysis()
        
        if not bayesian_results:
            logger.warning("No Bayesian decay indices were calculated.")
            return {}
            
        # Convert to standard decay_indices format
        converted_results = {}
        
        for clade_id, bayes_data in bayesian_results.items():
            converted_results[clade_id] = {
                'taxa': bayes_data['taxa'],
                'bayes_unconstrained_ml': bayes_data['unconstrained_ml'],
                'bayes_constrained_ml': bayes_data['constrained_ml'],
                'bayes_decay': bayes_data['bayes_decay'],
                # ML fields are None for Bayesian-only analysis
                'lnl_diff': None,
                'constrained_lnl': None,
                'AU_pvalue': None,
                'significant_AU': None
            }
            
            # Add normalized BD metrics if present
            for key, value in bayes_data.items():
                if key in ['bd_per_site', 'bd_relative', 'signal_to_noise', 
                          'effect_size', 'effect_size_robust', 'effect_size_weighted']:
                    converted_results[clade_id][key] = value
            
        logger.info(f"Calculated Bayesian decay indices for {len(converted_results)} branches.")
        return converted_results
    
    def _calculate_combined_decay_indices(self, perform_site_analysis=False):
        """
        Calculate both ML and Bayesian decay indices.
        
        Args:
            perform_site_analysis: Whether to perform site-specific analysis for ML
            
        Returns:
            Dictionary of decay indices with both ML and Bayesian metrics
        """
        logger.info("Calculating combined ML and Bayesian decay indices...")
        
        # First calculate ML decay indices
        ml_results = self._calculate_ml_decay_indices(perform_site_analysis)
        
        # Then run Bayesian analysis
        logger.info("Starting Bayesian analysis phase...")
        bayesian_results = self.run_bayesian_decay_analysis()
        
        # Merge results
        logger.info(f"Bayesian analysis returned {len(bayesian_results) if bayesian_results else 0} results")
        if bayesian_results:
            for clade_id in ml_results:
                if clade_id in bayesian_results:
                    # Add Bayesian fields to existing ML results
                    bayes_data = bayesian_results[clade_id]
                    ml_results[clade_id].update({
                        'bayes_unconstrained_ml': bayes_data['unconstrained_ml'],
                        'bayes_constrained_ml': bayes_data['constrained_ml'],
                        'bayes_decay': bayes_data['bayes_decay'],
                    })
                    
                    # Add normalized BD metrics if present
                    for key, value in bayes_data.items():
                        if key in ['bd_per_site', 'bd_relative', 'signal_to_noise', 
                                  'effect_size', 'effect_size_robust', 'effect_size_weighted']:
                            ml_results[clade_id][key] = value
                else:
                    # No Bayesian results for this clade
                    ml_results[clade_id].update({
                        'bayes_unconstrained_ml': None,
                        'bayes_constrained_ml': None,
                        'bayes_decay': None
                    })
        else:
            logger.warning("Bayesian analysis failed; results will contain ML metrics only.")
            
        # Add posterior probabilities if available
        if hasattr(self, 'posterior_probs') and self.posterior_probs:
            logger.info(f"Adding posterior probabilities for {len(self.posterior_probs)} clades")
            for clade_id in ml_results:
                # Get taxa set for this clade
                clade_taxa = frozenset(ml_results[clade_id]['taxa'])
                if clade_taxa in self.posterior_probs:
                    ml_results[clade_id]['posterior_prob'] = self.posterior_probs[clade_taxa]
                else:
                    ml_results[clade_id]['posterior_prob'] = None
            
        self.decay_indices = ml_results
        return self.decay_indices

    def _calculate_site_likelihoods(self, tree_files_list, branch_id):
        """
        Calculate site-specific likelihoods for ML tree vs constrained tree.

        Args:
            tree_files_list: List with [ml_tree_file, constrained_tree_file]
            branch_id: Identifier for the branch being analyzed

        Returns:
            Dictionary with site-specific likelihood differences or None if failed
        """
        if len(tree_files_list) != 2:
            logger.warning(f"Site analysis for branch {branch_id} requires exactly 2 trees (ML and constrained).")
            return None

        site_lnl_file = f"site_lnl_{branch_id}.txt"
        site_script_file = f"site_analysis_{branch_id}.nex"
        site_log_file = f"site_analysis_{branch_id}.log"

        # Create PAUP* script for site likelihood calculation
        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", self._get_paup_model_setup_cmds()]

        # Get both trees (ML and constrained)
        script_cmds.append(f"gettrees file={tree_files_list[0]} mode=3 storebrlens=yes;")
        script_cmds.append(f"gettrees file={tree_files_list[1]} mode=7 storebrlens=yes;")

        # Calculate site likelihoods for both trees
        script_cmds.append(f"lscores 1-2 / sitelikes=yes scorefile={site_lnl_file} replace=yes;")

        # Write PAUP* script
        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        script_path = self.temp_path / site_script_file
        script_path.write_text(paup_script_content)
        if self.debug:
            logger.debug(f"Site analysis script for {branch_id}:\n{paup_script_content}")

        try:
            # Run PAUP* to calculate site likelihoods
            self._run_paup_command_file(site_script_file, site_log_file, timeout_sec=600)

            # Parse the site likelihood file
            site_lnl_path = self.temp_path / site_lnl_file
            if not site_lnl_path.exists():
                logger.warning(f"Site likelihood file not found for branch {branch_id}.")
                return None

            # Read the site likelihoods file
            site_lnl_content = site_lnl_path.read_text()

            # Initialize dictionaries for tree likelihoods
            tree1_lnl = {}
            tree2_lnl = {}

            # Define patterns to extract data from the file
            # First pattern: Match the header line for each tree section
            tree_header_pattern = r'(\d+)\t([-\d\.]+)\t-\t-'

            # Second pattern: Match site and likelihood lines (indented with tabs)
            site_lnl_pattern = r'\t\t(\d+)\t([-\d\.]+)'

            # Find all tree headers
            tree_headers = list(re.finditer(tree_header_pattern, site_lnl_content))

            # Make sure we found at least 2 tree headers (Tree 1 and Tree 2)
            if len(tree_headers) < 2:
                logger.warning(f"Could not find enough tree headers in site likelihood file for branch {branch_id}")
                if self.debug:
                    logger.debug(f"Site likelihood file content (first 500 chars):\n{site_lnl_content[:500]}...")
                return None

            # Process each tree section
            for i, header_match in enumerate(tree_headers[:2]):  # Only process the first two trees
                tree_num = int(header_match.group(1))

                # If there's a next header, read up to it; otherwise, read to the end
                if i < len(tree_headers) - 1:
                    section_text = site_lnl_content[header_match.end():tree_headers[i+1].start()]
                else:
                    section_text = site_lnl_content[header_match.end():]

                # Find all site and likelihood entries
                site_matches = re.finditer(site_lnl_pattern, section_text)

                # Store data in appropriate dictionary
                for site_match in site_matches:
                    site_num = int(site_match.group(1))
                    lnl_val = float(site_match.group(2))

                    if tree_num == 1:
                        tree1_lnl[site_num] = lnl_val
                    else:
                        tree2_lnl[site_num] = lnl_val

            # Check if we have data for both trees
            if not tree1_lnl:
                logger.warning(f"No data found for Tree 1 in site likelihood file for branch {branch_id}")
                return None

            if not tree2_lnl:
                logger.warning(f"No data found for Tree 2 in site likelihood file for branch {branch_id}")
                return None

            # Create the site_data dictionary with differences
            site_data = {}
            all_sites = sorted(set(tree1_lnl.keys()) & set(tree2_lnl.keys()))

            for site_num in all_sites:
                ml_lnl = tree1_lnl[site_num]
                constrained_lnl = tree2_lnl[site_num]
                delta_lnl = ml_lnl - constrained_lnl

                site_data[site_num] = {
                    'lnL_ML': ml_lnl,
                    'lnL_constrained': constrained_lnl,
                    'delta_lnL': delta_lnl,
                    'supports_branch': delta_lnl < 0  # Negative delta means site supports ML branch
                }

            # Calculate summary statistics
            if site_data:
                deltas = [site_info['delta_lnL'] for site_info in site_data.values()]

                # Calculate signal statistics for effect size calculations
                signal_mean = np.mean(deltas)
                signal_std = np.std(deltas, ddof=1) if len(deltas) > 1 else 0  # Sample standard deviation
                signal_median = np.median(deltas)
                # Median Absolute Deviation (MAD) - more robust to outliers
                signal_mad = np.median(np.abs(np.array(deltas) - signal_median)) if len(deltas) > 0 else 0

                supporting_sites = sum(1 for d in deltas if d < 0)
                conflicting_sites = sum(1 for d in deltas if d > 0)
                neutral_sites = sum(1 for d in deltas if abs(d) < 1e-6)

                # Calculate sum of likelihood differences
                sum_supporting_delta = sum(d for d in deltas if d < 0)  # Sum of negative deltas (supporting)
                sum_conflicting_delta = sum(d for d in deltas if d > 0)  # Sum of positive deltas (conflicting)

                # Calculate weighted support ratio
                weighted_support_ratio = abs(sum_supporting_delta) / sum_conflicting_delta if sum_conflicting_delta > 0 else float('inf')

                # Calculate standard support ratio
                support_ratio = supporting_sites / conflicting_sites if conflicting_sites > 0 else float('inf')

                logger.info(f"Extracted site likelihoods for {len(site_data)} sites for branch {branch_id}")
                logger.info(f"Branch {branch_id}: {supporting_sites} supporting sites, {conflicting_sites} conflicting sites")
                logger.info(f"Branch {branch_id}: {supporting_sites} supporting sites, {conflicting_sites} conflicting sites, ratio: {support_ratio:.2f}")
                logger.info(f"Branch {branch_id}: Sum supporting delta: {sum_supporting_delta:.4f}, sum conflicting: {sum_conflicting_delta:.4f}, weighted ratio: {weighted_support_ratio:.2f}")

                # Return a comprehensive dictionary with all info
                return {
                    'site_data': site_data,
                    'supporting_sites': supporting_sites,
                    'conflicting_sites': conflicting_sites,
                    'neutral_sites': neutral_sites,
                    'support_ratio': support_ratio,
                    'sum_supporting_delta': sum_supporting_delta,
                    'sum_conflicting_delta': sum_conflicting_delta,
                    'weighted_support_ratio': weighted_support_ratio,
                    # Signal statistics for effect size calculations
                    'signal_mean': signal_mean,
                    'signal_std': signal_std,
                    'signal_median': signal_median,
                    'signal_mad': signal_mad
                }
            else:
                logger.warning(f"No comparable site likelihoods found for branch {branch_id}")
                return None

        except Exception as e:
            logger.error(f"Failed to calculate site likelihoods for branch {branch_id}: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Traceback for site likelihood calculation error:\n{traceback.format_exc()}")
            return None

    def run_au_test(self, tree_filenames_relative: list):
        """
        Run the AU (Approximately Unbiased) test on multiple trees using PAUP*.

        Args:
            tree_filenames_relative: List of tree filenames (relative to self.temp_path)

        Returns:
            Dictionary mapping tree indices to results or None if test fails
        """
        if not tree_filenames_relative:
            logger.error("No tree files for AU test.")
            return None
        num_trees = len(tree_filenames_relative)
        if num_trees < 2: # AU test is meaningful for comparing multiple trees
            logger.warning(f"AU test needs >= 2 trees; {num_trees} provided. Skipping AU test.")
            # If it's just the ML tree, we can return its own info conventionally
            if num_trees == 1 and tree_filenames_relative[0] == ML_TREE_FN and self.ml_likelihood is not None:
                return {1: {'lnL': self.ml_likelihood, 'AU_pvalue': 1.0}} # Best tree p-val = 1
            return None

        script_cmds = [f"execute {NEXUS_ALIGNMENT_FN};", self._get_paup_model_setup_cmds()]

        # Load all trees: first in mode=3 (reset tree buffer), rest in mode=7 (add to buffer)
        first_tree = tree_filenames_relative[0]
        script_cmds.append(f"gettrees file={first_tree} mode=3 storebrlens=yes;")

        for rel_fn in tree_filenames_relative[1:]:
            script_cmds.append(f"gettrees file={rel_fn} mode=7 storebrlens=yes;")

        # Add debugging commands to see tree status
        if self.debug:
            script_cmds.append("treeinfo;")     # Show tree information

        # Make sure tree indices match our expectations (1-based indexing)
        if num_trees > 1:
            script_cmds.append(f"lscores 1-{num_trees} / scorefile=all_tree_scores.txt replace=yes;")

        # AU test command with additional options for improved reliability
        script_cmds.append(f"lscores 1-{num_trees} / autest=yes scorefile={AU_TEST_SCORE_FN} replace=yes;")

        # Save log with results and trees for reference
        script_cmds.append("log file=au_test_detailed.log replace=yes;")

        # FIX: Use explicit tree range instead of 'all'
        script_cmds.append(f"describe 1-{num_trees} / plot=none;")  # Show tree descriptions
        script_cmds.append(f"lscores 1-{num_trees};")              # Show scores again
        script_cmds.append("log stop;")

        paup_script_content = f"#NEXUS\nbegin paup;\n" + "\n".join(script_cmds) + "\nquit;\nend;\n"
        au_cmd_path = self.temp_path / AU_TEST_NEX_FN
        au_cmd_path.write_text(paup_script_content)
        if self.debug: logger.debug(f"AU test PAUP* script ({au_cmd_path}):\n{paup_script_content}")

        try:
            self._run_paup_command_file(AU_TEST_NEX_FN, AU_LOG_FN, timeout_sec=max(1800, 600 * num_trees / 10)) # Dynamic timeout

            # Parse results from the AU test results file
            return self._parse_au_results(self.temp_path / AU_LOG_FN)
        except Exception as e:
            logger.error(f"AU test execution failed: {e}")
            return None

    def _parse_au_results(self, au_log_path: Path):
        """
        Parse the results of an Approximately Unbiased (AU) test from PAUP* log file.

        Args:
            au_log_path: Path to the PAUP* log file containing AU test results

        Returns:
            Dictionary mapping tree index to dict with 'lnL' and 'AU_pvalue' keys, or None if parsing failed
        """
        if not au_log_path.exists():
            logger.warning(f"AU test log file not found: {au_log_path}")
            return None

        try:
            log_content = au_log_path.read_text()
            if self.debug: logger.debug(f"AU test log file content (excerpt):\n{log_content[:1000]}...")

            # Look for the AU test results section
            # First try to find the formatted table with header and rows that looks like:
            #    Tree         -ln L    Diff -ln L      AU
            # --------------------------------------------
            #       1    6303.66091        (best)
            #       7    6304.45629       0.79537  0.6069
            # ...etc

            au_results = {}

            # Use a multi-line pattern to find the table with AU test results
            au_table_pattern = r'Tree\s+-ln L\s+Diff[^-]*\n-+\n(.*?)(?=\n\s*\n|\n[^\d\s]|$)'
            au_match = re.search(au_table_pattern, log_content, re.DOTALL)

            if au_match:
                table_text = au_match.group(1)
                logger.debug(f"Found AU test table:\n{table_text}")

                # Parse each line of the table
                for line in table_text.strip().split('\n'):
                    # Skip empty lines
                    if not line.strip():
                        continue

                    # Format is typically: tree_num, -ln L, Diff -ln L, AU p-value
                    # Example:       1    6303.66091        (best)
                    #                7    6304.45629       0.79537  0.6069
                    parts = line.strip().split()

                    # Make sure we have at least the tree number and likelihood
                    if len(parts) >= 2 and parts[0].isdigit():
                        tree_idx = int(parts[0])
                        ln_l = float(parts[1])

                        # Get p-value - it might be "(best)" for the best tree or a number for others
                        # The p-value is typically the last element in the line
                        p_val = None
                        if len(parts) >= 3:
                            p_val_str = parts[-1]  # Take the last element as the p-value

                            # Handle special cases
                            if p_val_str == "(best)":
                                p_val = 1.0  # Best tree has p-value of 1.0
                            elif "~0" in p_val_str:
                                # Values like "~0*" mean extremely small p-values
                                p_val = 0.0001  # Use a small non-zero value
                            else:
                                # Normal p-values, remove any trailing asterisks
                                p_val = float(p_val_str.rstrip("*"))

                        au_results[tree_idx] = {
                            'lnL': ln_l,
                            'AU_pvalue': p_val
                        }

                if au_results:
                    logger.info(f"Successfully parsed AU test results for {len(au_results)} trees.")
                    for tree_idx, data in sorted(au_results.items()):
                        logger.debug(f"Tree {tree_idx}: lnL={data['lnL']:.4f}, AU p-value={data['AU_pvalue']}")
                    return au_results

            # If we couldn't find the AU test table, try an alternative approach
            # Look for the detailed AU test results in the log
            detailed_pattern = r'P values for.*?Tree\s+(-ln L)\s+Diff.*?\n.*?\n(.*?)(?=\n\s*\n|\n[^\d\s]|$)'
            detailed_match = re.search(detailed_pattern, log_content, re.DOTALL)

            if detailed_match:
                results_text = detailed_match.group(2)
                logger.debug(f"Found detailed AU test results table:\n{results_text}")

                for line in results_text.strip().split('\n'):
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[0].isdigit():
                        tree_idx = int(parts[0])
                        ln_l = float(parts[1])

                        # AU p-value is typically at position 3 (index 2)
                        p_val = None
                        if len(parts) > 3:
                            p_val_str = parts[3] if len(parts) > 3 else parts[-1]

                            if p_val_str == "(best)":
                                p_val = 1.0
                            elif p_val_str.startswith("~"):
                                p_val = 0.0001
                            else:
                                p_val = float(p_val_str.rstrip("*"))

                        au_results[tree_idx] = {
                            'lnL': ln_l,
                            'AU_pvalue': p_val
                        }

                if au_results:
                    logger.info(f"Successfully parsed detailed AU test results for {len(au_results)} trees.")
                    return au_results

            # Third approach: try to parse AU scores from the direct output in the log
            # This pattern is more specific to the format observed in your logs
            au_direct_pattern = r'Tree\s+.*?-ln L\s+Diff -ln L\s+AU\n-+\n(.*?)(?=\n\s*\n|\Z)'
            direct_match = re.search(au_direct_pattern, log_content, re.DOTALL)

            if direct_match:
                direct_table = direct_match.group(1)
                logger.debug(f"Found direct AU test table:\n{direct_table}")

                for line in direct_table.strip().split('\n'):
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[0].isdigit():
                        tree_idx = int(parts[0])
                        ln_l = float(parts[1])

                        # Handle the p-value which is in the last column
                        p_val = None
                        if len(parts) >= 4:  # We need at least 4 parts for tree, lnL, diff, and p-value
                            p_val_str = parts[-1]

                            if p_val_str == "(best)":
                                p_val = 1.0
                            elif "~0" in p_val_str:
                                p_val = 0.0001
                            else:
                                try:
                                    p_val = float(p_val_str.rstrip("*"))
                                except ValueError:
                                    # If conversion fails, check if it's just due to an asterisk
                                    if p_val_str.endswith("*"):
                                        try:
                                            p_val = float(p_val_str[:-1])
                                        except ValueError:
                                            p_val = None

                        au_results[tree_idx] = {
                            'lnL': ln_l,
                            'AU_pvalue': p_val
                        }

                if au_results:
                    logger.info(f"Successfully parsed direct AU test results for {len(au_results)} trees.")
                    return au_results

            # Try to parse from scorefile if results not found in log
            au_test_score_fn = "au_test_results.txt"
            score_file_path = self.temp_path / au_test_score_fn
            if score_file_path.exists():
                return self._parse_au_results_from_scorefile(score_file_path)

            logger.warning("Failed to find AU test results in log file formats. Checking for other sources.")
            return None

        except Exception as e:
            logger.error(f"Error parsing AU test results: {e}")
            if self.debug:
                import traceback
                logger.debug(f"AU test parsing traceback: {traceback.format_exc()}")
            return None

    def _parse_au_results_from_scorefile(self, score_file_path: Path):
        """
        Parse AU test results directly from the score file produced by PAUP*.
        This serves as a backup to parsing the log file.

        Args:
            score_file_path: Path to the AU test score file

        Returns:
            Dictionary mapping tree index to dict with 'lnL' and 'AU_pvalue' keys, or None if parsing failed
        """
        if not score_file_path.exists():
            logger.warning(f"AU test score file not found: {score_file_path}")
            return None

        try:
            file_content = score_file_path.read_text()
            if self.debug: logger.debug(f"AU test score file content (excerpt):\n{file_content[:1000]}...")

            # AU test score files typically have headers like:
            # Tree     -lnL     Diff    P-value
            header_pattern = r'^\s*Tree\s+\-lnL\s+(?:Diff\s+)?P\-?value'

            # Find where results start
            lines = file_content.splitlines()
            results_start = -1
            for i, line in enumerate(lines):
                if re.match(header_pattern, line, re.IGNORECASE):
                    results_start = i + 1
                    break

            if results_start == -1:
                logger.warning("Could not find AU test results header in score file.")
                return None

            # Parse results
            au_results = {}
            for i in range(results_start, len(lines)):
                line = lines[i].strip()
                if not line or 'tree' in line.lower():  # Skip empty lines or new headers
                    continue

                # Try to parse - expecting tree_num, lnL, [diff], p_value
                parts = line.split()
                if len(parts) < 3:  # Need at least tree, lnL, p-value
                    continue

                try:
                    tree_idx = int(parts[0])
                    lnl = float(parts[1])

                    # The p-value is either the last or third column depending on format
                    p_val = float(parts[-1]) if len(parts) >= 3 else None

                    au_results[tree_idx] = {
                        'lnL': lnl,
                        'AU_pvalue': p_val
                    }
                except (ValueError, IndexError) as e:
                    logger.debug(f"Failed to parse AU result line '{line}': {e}")

            if au_results:
                logger.info(f"Successfully parsed {len(au_results)} AU test results from score file.")
                return au_results
            else:
                logger.warning("No AU test results could be parsed from score file.")
                return None

        except Exception as e:
            logger.error(f"Error parsing AU test score file: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Score file parsing error: {traceback.format_exc()}")
            return None

    def annotate_trees(self, output_dir: Path, base_filename: str = "annotated_tree"):
        """
        Create annotated trees with different support values:
        1. A tree with AU p-values as branch labels
        2. A tree with log-likelihood differences as branch labels
        3. A combined tree with both values as FigTree-compatible branch labels
        4. A tree with bootstrap values if bootstrap analysis was performed
        5. A comprehensive tree with bootstrap, AU, and LnL values if bootstrap was performed

        Args:
            output_dir: Directory to save the tree files
            base_filename: Base name for the tree files (without extension)

        Returns:
            Dict with paths to the created tree files
        """
        if not self.ml_tree or not self.decay_indices:
            logger.warning("ML tree or decay indices missing. Cannot annotate trees.")
            return {}

        output_dir.mkdir(parents=True, exist_ok=True)
        tree_files = {}

        try:
            # Create AU p-value annotated tree
            au_tree_path = output_dir / f"{base_filename}_au.nwk"
            try:
                # Work on a copy to avoid modifying self.ml_tree
                temp_tree_for_au = self.temp_path / f"ml_tree_for_au_annotation.nwk"
                Phylo.write(self.ml_tree, str(temp_tree_for_au), "newick")
                cleaned_tree_path = self._clean_newick_tree(temp_tree_for_au)
                au_tree = Phylo.read(str(cleaned_tree_path), "newick")

                annotated_nodes_count = 0
                for node in au_tree.get_nonterminals():
                    if not node or not node.clades: continue
                    node_taxa_set = set(leaf.name for leaf in node.get_terminals())

                    # Find matching entry in decay_indices by taxa set
                    matched_data = None
                    matched_clade_id = None
                    for decay_id_str, decay_info in self.decay_indices.items():
                        if 'taxa' in decay_info and set(decay_info['taxa']) == node_taxa_set:
                            matched_data = decay_info
                            matched_clade_id = decay_id_str
                            break

                    node.confidence = None  # Default
                    if matched_data and 'AU_pvalue' in matched_data and matched_data['AU_pvalue'] is not None:
                        au_pvalue = matched_data['AU_pvalue']
                        # Create clear separation between clade name and AU p-value
                        node.name = f"{matched_clade_id} - AU:{au_pvalue:.4f}"
                        annotated_nodes_count += 1

                Phylo.write(au_tree, str(au_tree_path), "newick")
                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {self._get_display_path(au_tree_path)} (type: au).")
                tree_files['au'] = au_tree_path
            except Exception as e:
                logger.error(f"Failed to create AU tree: {e}")

            # Create log-likelihood difference annotated tree
            lnl_tree_path = output_dir / f"{base_filename}_delta_lnl.nwk"
            try:
                temp_tree_for_lnl = self.temp_path / f"ml_tree_for_lnl_annotation.nwk"
                Phylo.write(self.ml_tree, str(temp_tree_for_lnl), "newick")
                cleaned_tree_path = self._clean_newick_tree(temp_tree_for_lnl)
                lnl_tree = Phylo.read(str(cleaned_tree_path), "newick")

                annotated_nodes_count = 0
                for node in lnl_tree.get_nonterminals():
                    if not node or not node.clades: continue
                    node_taxa_set = set(leaf.name for leaf in node.get_terminals())

                    matched_data = None
                    matched_clade_id = None
                    for decay_id_str, decay_info in self.decay_indices.items():
                        if 'taxa' in decay_info and set(decay_info['taxa']) == node_taxa_set:
                            matched_data = decay_info
                            matched_clade_id = decay_id_str
                            break

                    node.confidence = None  # Default
                    if matched_data and 'lnl_diff' in matched_data and matched_data['lnl_diff'] is not None:
                        lnl_diff = abs(matched_data['lnl_diff'])
                        # Create clear separation between clade name and LnL difference
                        node.name = f"{matched_clade_id} - ΔlnL:{lnl_diff:.4f}"
                        annotated_nodes_count += 1

                Phylo.write(lnl_tree, str(lnl_tree_path), "newick")
                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {self._get_display_path(lnl_tree_path)} (type: delta_lnl).")
                tree_files['lnl'] = lnl_tree_path
            except Exception as e:
                logger.error(f"Failed to create LNL tree: {e}")

            # Create combined annotation tree for FigTree
            combined_tree_path = output_dir / f"{base_filename}_combined.nwk"
            try:
                # For the combined approach, we'll directly modify the Newick string
                # First, get both trees as strings
                temp_tree_for_combined = self.temp_path / f"ml_tree_for_combined_annotation.nwk"
                Phylo.write(self.ml_tree, str(temp_tree_for_combined), "newick")

                # Create a mapping from node taxa sets to combined annotation strings
                node_annotations = {}

                # If bootstrap analysis was performed, get bootstrap values first
                bootstrap_values = {}
                if hasattr(self, 'bootstrap_tree') and self.bootstrap_tree:
                    for node in self.bootstrap_tree.get_nonterminals():
                        if node.confidence is not None:
                            taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                            bootstrap_values[taxa_set] = node.confidence

                for node in self.ml_tree.get_nonterminals():
                    if not node or not node.clades: continue
                    node_taxa_set = frozenset(leaf.name for leaf in node.get_terminals())

                    # Initialize annotation parts
                    annotation_parts = []
                    clade_id = None  # Store the matched clade_id

                    # Add bootstrap value if available
                    if bootstrap_values and node_taxa_set in bootstrap_values:
                        bs_val = bootstrap_values[node_taxa_set]
                        annotation_parts.append(f"BS:{int(bs_val)}")

                    # Add AU and LnL values if available
                    for decay_id_str, decay_info in self.decay_indices.items():
                        if 'taxa' in decay_info and frozenset(decay_info['taxa']) == node_taxa_set:
                            clade_id = decay_id_str  # Save clade ID for later use
                            au_val = decay_info.get('AU_pvalue')
                            lnl_val = decay_info.get('lnl_diff')
                            bayes_decay = decay_info.get('bayes_decay')

                            if au_val is not None:
                                annotation_parts.append(f"AU:{au_val:.4f}")

                            if lnl_val is not None:
                                annotation_parts.append(f"ΔlnL:{abs(lnl_val):.4f}")
                                
                            if bayes_decay is not None:
                                annotation_parts.append(f"BD:{bayes_decay:.4f}")
                                
                            # Add BD/site normalization if available
                            bd_per_site = decay_info.get('bd_per_site')
                            if bd_per_site is not None:
                                annotation_parts.append(f"PS:{bd_per_site:.6f}")
                                
                            # Add parsimony decay if available
                            pars_decay = decay_info.get('pars_decay')
                            if pars_decay is not None:
                                annotation_parts.append(f"PD:{pars_decay}")
                                
                            # Add posterior probability if available
                            post_prob = decay_info.get('posterior_prob')
                            if post_prob is not None:
                                annotation_parts.append(f"PP:{post_prob:.2f}")
                                
                            # Add effect size if available
                            effect_size = decay_info.get('effect_size') or decay_info.get('effect_size_robust') or decay_info.get('effect_size_weighted')
                            if effect_size is not None:
                                annotation_parts.append(f"ES:{effect_size:.2f}")
                            break

                    # Format annotations using the new method
                    if clade_id or annotation_parts:
                        # Convert annotation_parts to dict format
                        ann_dict = {}
                        for part in annotation_parts:
                            if ':' in part:
                                key, val = part.split(':', 1)
                                ann_dict[key] = val
                        
                        # Use compact format for combined tree
                        formatted = self._format_tree_annotation(clade_id, ann_dict, style="compact")
                        if formatted:
                            node_annotations[node_taxa_set] = formatted

                # Now, we'll manually construct a combined tree by using string replacement on the base tree
                # First, make a working copy of the ML tree
                cleaned_tree_path = self._clean_newick_tree(temp_tree_for_combined)
                combined_tree = Phylo.read(str(cleaned_tree_path), "newick")

                # Add our custom annotations
                annotated_nodes_count = 0
                for node in combined_tree.get_nonterminals():
                    if not node or not node.clades: continue
                    node_taxa_set = frozenset(leaf.name for leaf in node.get_terminals())

                    if node_taxa_set in node_annotations:
                        # We need to use string values for combined annotation
                        # Save our combined annotation as a string in .name instead of .confidence
                        # This is a hack that works with some tree viewers including FigTree
                        node.name = node_annotations[node_taxa_set]
                        annotated_nodes_count += 1

                # Write the modified tree
                Phylo.write(combined_tree, str(combined_tree_path), "newick")

                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {self._get_display_path(combined_tree_path)} (type: combined).")
                tree_files['combined'] = combined_tree_path
            except Exception as e:
                logger.error(f"Failed to create combined tree: {e}")
                import traceback
                logger.debug(f"Traceback: {traceback.format_exc()}")

            # Handle bootstrap tree if bootstrap analysis was performed
            if hasattr(self, 'bootstrap_tree') and self.bootstrap_tree:
                # 1. Create a bootstrap tree using ML tree topology with bootstrap values
                bootstrap_tree_path = output_dir / f"{base_filename}_bootstrap.nwk"
                try:
                    # Create a copy of the ML tree for bootstrap annotation
                    temp_tree_for_bootstrap = self.temp_path / f"ml_tree_for_bootstrap_annotation.nwk"
                    Phylo.write(self.ml_tree, str(temp_tree_for_bootstrap), "newick")
                    cleaned_tree_path = self._clean_newick_tree(temp_tree_for_bootstrap)
                    ml_tree_for_bootstrap = Phylo.read(str(cleaned_tree_path), "newick")
                    
                    # Extract bootstrap values from the consensus tree
                    bootstrap_values = {}
                    for node in self.bootstrap_tree.get_nonterminals():
                        if node.confidence is not None:
                            taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                            bootstrap_values[taxa_set] = node.confidence
                    
                    # Apply bootstrap values to the ML tree
                    annotated_nodes_count = 0
                    for node in ml_tree_for_bootstrap.get_nonterminals():
                        if not node or not node.clades: continue
                        node_taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                        
                        if node_taxa_set in bootstrap_values:
                            bs_val = bootstrap_values[node_taxa_set]
                            node.name = f"{int(bs_val)}"
                            annotated_nodes_count += 1
                    
                    # Write the ML tree with bootstrap values
                    Phylo.write(ml_tree_for_bootstrap, str(bootstrap_tree_path), "newick")
                    logger.info(f"Bootstrap tree (ML topology with bootstrap values) written to {self._get_display_path(bootstrap_tree_path)}")
                    tree_files['bootstrap'] = bootstrap_tree_path
                except Exception as e:
                    logger.error(f"Failed to write bootstrap tree: {e}")

                # 2. Create a comprehensive tree with bootstrap, AU and LnL values
                comprehensive_tree_path = output_dir / f"{base_filename}_comprehensive.nwk"
                try:
                    temp_tree_for_comprehensive = self.temp_path / f"ml_tree_for_comprehensive_annotation.nwk"
                    Phylo.write(self.ml_tree, str(temp_tree_for_comprehensive), "newick")
                    cleaned_tree_path = self._clean_newick_tree(temp_tree_for_comprehensive)
                    comprehensive_tree = Phylo.read(str(cleaned_tree_path), "newick")

                    # Get bootstrap values for each clade
                    bootstrap_values = {}
                    for node in self.bootstrap_tree.get_nonterminals():
                        if node.confidence is not None:
                            taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                            bootstrap_values[taxa_set] = node.confidence

                    # Create comprehensive annotations
                    node_annotations = {}
                    for node in self.ml_tree.get_nonterminals():
                        if not node or not node.clades: continue
                        node_taxa_set = frozenset(leaf.name for leaf in node.get_terminals())

                        # Find matching decay info
                        matched_data = None
                        matched_clade_id = None
                        for decay_id_str, decay_info in self.decay_indices.items():
                            if 'taxa' in decay_info and frozenset(decay_info['taxa']) == node_taxa_set:
                                matched_data = decay_info
                                matched_clade_id = decay_id_str
                                break

                        # Combine all values
                        annotation_parts = []

                        # Add bootstrap value if available
                        if node_taxa_set in bootstrap_values:
                            bs_val = bootstrap_values[node_taxa_set]
                            annotation_parts.append(f"BS:{int(bs_val)}")

                        # Add AU and LnL values if available
                        if matched_data:
                            au_val = matched_data.get('AU_pvalue')
                            lnl_val = matched_data.get('lnl_diff')
                            bayes_decay = matched_data.get('bayes_decay')

                            if au_val is not None:
                                annotation_parts.append(f"AU:{au_val:.4f}")

                            if lnl_val is not None:
                                annotation_parts.append(f"ΔlnL:{abs(lnl_val):.4f}")
                                
                            if bayes_decay is not None:
                                annotation_parts.append(f"BD:{bayes_decay:.4f}")
                                
                            # Add BD/site normalization if available
                            bd_per_site = matched_data.get('bd_per_site')
                            if bd_per_site is not None:
                                annotation_parts.append(f"PS:{bd_per_site:.6f}")
                                
                            # Add parsimony decay if available
                            pars_decay = matched_data.get('pars_decay')
                            if pars_decay is not None:
                                annotation_parts.append(f"PD:{pars_decay}")
                                
                            # Add posterior probability if available
                            post_prob = matched_data.get('posterior_prob')
                            if post_prob is not None:
                                annotation_parts.append(f"PP:{post_prob:.2f}")
                                
                            # Add effect size if available
                            effect_size = matched_data.get('effect_size') or matched_data.get('effect_size_robust') or matched_data.get('effect_size_weighted')
                            if effect_size is not None:
                                annotation_parts.append(f"ES:{effect_size:.2f}")

                        if annotation_parts:
                            # For comprehensive trees, add clear separation between clade and metrics if we have a clade ID
                            if matched_clade_id:
                                metrics_part = "|".join(annotation_parts)
                                node_annotations[node_taxa_set] = f"{matched_clade_id} - {metrics_part}"
                            else:
                                node_annotations[node_taxa_set] = "|".join(annotation_parts)

                    # Apply annotations to tree
                    annotated_nodes_count = 0
                    for node in comprehensive_tree.get_nonterminals():
                        if not node or not node.clades: continue
                        node_taxa_set = frozenset(leaf.name for leaf in node.get_terminals())

                        if node_taxa_set in node_annotations:
                            node.name = node_annotations[node_taxa_set]
                            annotated_nodes_count += 1

                    # Write the tree
                    Phylo.write(comprehensive_tree, str(comprehensive_tree_path), "newick")
                    logger.info(f"Comprehensive tree with {annotated_nodes_count} branch values written to {self._get_display_path(comprehensive_tree_path)}")
                    tree_files['comprehensive'] = comprehensive_tree_path
                except Exception as e:
                    logger.error(f"Failed to create comprehensive tree: {e}")
                    if self.debug:
                        import traceback
                        logger.debug(f"Traceback: {traceback.format_exc()}")

            # Create Bayesian-specific trees if Bayesian results are available
            has_bayesian = any(d.get('bayes_decay') is not None for d in self.decay_indices.values())
            if has_bayesian:
                # Create Bayes decay annotated tree
                bayes_decay_tree_path = output_dir / f"{base_filename}_bayes_decay.nwk"
                try:
                    temp_tree_for_bd = self.temp_path / f"ml_tree_for_bd_annotation.nwk"
                    Phylo.write(self.ml_tree, str(temp_tree_for_bd), "newick")
                    cleaned_tree_path = self._clean_newick_tree(temp_tree_for_bd)
                    bd_tree = Phylo.read(str(cleaned_tree_path), "newick")

                    annotated_nodes_count = 0
                    for node in bd_tree.get_nonterminals():
                        if not node or not node.clades: continue
                        node_taxa_set = set(leaf.name for leaf in node.get_terminals())

                        matched_data = None
                        matched_clade_id = None
                        for decay_id_str, decay_info in self.decay_indices.items():
                            if 'taxa' in decay_info and set(decay_info['taxa']) == node_taxa_set:
                                matched_data = decay_info
                                matched_clade_id = decay_id_str
                                break

                        node.confidence = None  # Default
                        if matched_data and 'bayes_decay' in matched_data and matched_data['bayes_decay'] is not None:
                            bayes_decay_val = matched_data['bayes_decay']
                            
                            # Create annotation with normalized values if available
                            if self.normalize_bd and 'bd_per_site' in matched_data:
                                bd_per_site = matched_data['bd_per_site']
                                annotation_parts = [f"PS:{bd_per_site:.6f}"]
                                
                                # Add effect size if available
                                effect_size = matched_data.get('effect_size') or matched_data.get('effect_size_robust') or matched_data.get('effect_size_weighted')
                                if effect_size is not None:
                                    annotation_parts.append(f"ES:{effect_size:.2f}")
                                
                                annotation_str = ", ".join(annotation_parts)
                                node.name = f"{matched_clade_id} - BD:{bayes_decay_val:.2f} ({annotation_str})"
                            else:
                                node.name = f"{matched_clade_id} - BD:{bayes_decay_val:.4f}"
                            annotated_nodes_count += 1

                    Phylo.write(bd_tree, str(bayes_decay_tree_path), "newick")
                    logger.info(f"Annotated tree with {annotated_nodes_count} Bayes decay values written to {self._get_display_path(bayes_decay_tree_path)}")
                    tree_files['bayes_decay'] = bayes_decay_tree_path
                except Exception as e:
                    logger.error(f"Failed to create Bayes decay tree: {e}")


            return tree_files

        except Exception as e:
            logger.error(f"Failed to annotate trees: {e}")
            if hasattr(self, 'debug') and self.debug:
                import traceback
                logger.debug(f"Traceback: {traceback.format_exc()}")
            return tree_files  # Return any successfully created files

    def _write_support_table(self, f, has_ml, has_bayesian, has_parsimony, has_posterior, has_bootstrap):
        """Write the formatted support values table."""
        box = self._get_box_chars()
        
        # Check if effect size data is available
        has_effect_size = any(
            any(key.startswith('effect_size') for key in data.keys()) 
            for data in self.decay_indices.values()
        )
        
        # Build header structure
        if self.output_style == "unicode":
            # Top border
            f.write("┌──────────┬──────┬")
            if has_ml:
                f.write("────────────────────────────────┬")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    # Extended width for BD + BD/site + BD% + Effect Size
                    f.write("──────────────────────────────────────────────┬")
                elif self.normalize_bd:
                    f.write("───────────────────────────────┬")
                else:
                    f.write("────────┬")
            if has_parsimony:
                f.write("───────────────────────┬")
            f.write("\n")
            
            # Main headers
            f.write("│ Clade ID │ Taxa │")
            if has_ml:
                f.write("    Maximum Likelihood          │")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    f.write("                 Bayesian (w/ Effect Size)                │")
                elif self.normalize_bd:
                    f.write("         Bayesian (Normalized)        │")
                else:
                    f.write("Bayesian │")
            if has_parsimony:
                f.write("      Parsimony        │")
            f.write("\n")
            
            # Sub-headers
            f.write("│          │      ├")
            if has_ml:
                f.write("──────────┬──────────┬──────────┼")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    f.write("────────┬────────┬────────┬────────┼")
                elif self.normalize_bd:
                    f.write("────────┬────────┬────────┼")
                else:
                    f.write("────────┼")
            if has_parsimony:
                f.write("───────────────────────┤")
            f.write("\n")
            
            # Column names
            f.write("│          │      │")
            if has_ml:
                f.write(" ΔlnL     │ AU p-val │ Support  │")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    f.write(" BD     │ BD/site │ BD%    │ ES     │")
                elif self.normalize_bd:
                    f.write(" BD     │ BD/site │ BD%    │")
                else:
                    f.write(" BD     │")
            if has_parsimony:
                f.write(" Decay │ Post.Prob     │" if has_posterior else " Decay                 │")
            f.write("\n")
            
            # Bottom border
            f.write("├──────────┼──────┼")
            if has_ml:
                f.write("──────────┼──────────┼──────────┼")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    f.write("────────┼────────┼────────┼────────┼")
                elif self.normalize_bd:
                    f.write("────────┼────────┼────────┼")
                else:
                    f.write("────────┼")
            if has_parsimony:
                f.write("───────┼───────────────┤" if has_posterior else "───────────────────────┤")
            f.write("\n")
        else:
            # ASCII version
            header_parts = ["Clade ID", "Taxa"]
            if has_ml:
                header_parts.extend(["ΔlnL", "AU p-val", "Support"])
            if has_bayesian:
                header_parts.append("BD")
            if has_parsimony:
                header_parts.append("P.Decay")
                if has_posterior:
                    header_parts.append("Post.Prob")
            
            f.write(self._format_table_row(header_parts, [10, 6] + [10] * (len(header_parts) - 2)) + "\n")
            f.write("-" * (sum([10, 6] + [10] * (len(header_parts) - 2)) + 3 * (len(header_parts) - 1)) + "\n")
        
        # Data rows
        for clade_id, data in sorted(self.decay_indices.items()):
            taxa_list = sorted(data.get('taxa', []))
            num_taxa = len(taxa_list)
            
            row_values = [clade_id, str(num_taxa)]
            
            if has_ml:
                lnl_diff = data.get('lnl_diff')
                au_pval = data.get('AU_pvalue')
                
                if lnl_diff is not None:
                    row_values.append(f"{lnl_diff:.3f}")
                else:
                    row_values.append("N/A")
                    
                if au_pval is not None:
                    row_values.append(f"{au_pval:.4f}")
                    row_values.append(self._format_support_symbol(au_pval))
                else:
                    row_values.extend(["N/A", "N/A"])
            
            if has_bayesian:
                bd = data.get('bayes_decay')
                
                if bd is not None:
                    row_values.append(f"{bd:.2f}")
                else:
                    row_values.append("N/A")
                
                # Add normalized metrics if enabled
                if self.normalize_bd:
                    bd_per_site = data.get('bd_per_site')
                    bd_relative = data.get('bd_relative')
                    
                    if bd_per_site is not None:
                        row_values.append(f"{bd_per_site:.6f}")
                    else:
                        row_values.append("N/A")
                        
                    if bd_relative is not None:
                        row_values.append(f"{bd_relative*100:.3f}")
                    else:
                        row_values.append("N/A")
                    
                    # Add effect size if available
                    if has_effect_size:
                        # Prefer standard effect size, but fall back to robust if available
                        effect_size = data.get('effect_size') or data.get('effect_size_robust') or data.get('effect_size_weighted')
                        if effect_size is not None:
                            row_values.append(f"{effect_size:.2f}")
                        else:
                            row_values.append("N/A")
                    
            
            if has_parsimony:
                pd = data.get('pars_decay')
                row_values.append(str(pd) if pd is not None else "N/A")
                
                if has_posterior:
                    pp = data.get('posterior_prob')
                    row_values.append(f"{pp:.2f}" if pp is not None else "N/A")
            
            if self.output_style == "unicode":
                f.write("│ " + " │ ".join(row_values) + " │\n")
            else:
                f.write(self._format_table_row(row_values, [10, 6] + [10] * (len(row_values) - 2)) + "\n")
        
        # Bottom border
        if self.output_style == "unicode":
            f.write("└──────────┴──────┴")
            if has_ml:
                f.write("──────────┴──────────┴──────────┴")
            if has_bayesian:
                if self.normalize_bd and has_effect_size:
                    f.write("────────┴────────┴────────┴────────┴")
                elif self.normalize_bd:
                    f.write("────────┴────────┴────────┴")
                else:
                    f.write("────────┴")
            if has_parsimony:
                f.write("───────┴───────────────┘" if has_posterior else "───────────────────────┘")
            f.write("\n")
    
    def write_formatted_results(self, output_path: Path):
        """Write results in formatted table style based on output_style setting."""
        if self.output_style == "minimal" or not self.decay_indices:
            # Fall back to original method for minimal style
            return self.write_results(output_path)
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Check which types of results we have
        has_ml = any(d.get('lnl_diff') is not None for d in self.decay_indices.values())
        has_bayesian = any(d.get('bayes_decay') is not None for d in self.decay_indices.values())
        has_parsimony = any(d.get('pars_decay') is not None for d in self.decay_indices.values())
        has_posterior = any(d.get('posterior_prob') is not None for d in self.decay_indices.values())
        has_bootstrap = hasattr(self, 'bootstrap_tree') and self.bootstrap_tree
        
        box = self._get_box_chars()
        
        with output_path.open('w') as f:
            # Header section
            if self.output_style == "unicode":
                f.write("═" * 100 + "\n")
                f.write(" " * 30 + "panDecay Branch Support Analysis Results" + " " * 30 + "\n")
                f.write("═" * 100 + "\n\n")
            else:
                f.write("=" * 100 + "\n")
                f.write(" " * 30 + "panDecay Branch Support Analysis Results" + " " * 30 + "\n")
                f.write("=" * 100 + "\n\n")
            
            # Analysis summary
            f.write("Analysis Summary\n")
            f.write("─" * 16 + "\n" if self.output_style == "unicode" else "-" * 16 + "\n")
            
            if self.ml_likelihood is not None:
                f.write(f"• ML tree log-likelihood: {self.ml_likelihood:.3f}\n")
            
            analysis_types = []
            if self.do_ml:
                analysis_types.append("Maximum Likelihood")
            if self.do_bayesian:
                analysis_types.append(f"Bayesian ({self.bayesian_software})")
            if self.do_parsimony:
                analysis_types.append("Parsimony")
            
            f.write(f"• Analysis types: {' + '.join(analysis_types)}\n")
            f.write(f"• Total clades analyzed: {len(self.decay_indices)}\n\n")
            
            # Branch support table
            f.write("Branch Support Values\n")
            f.write("─" * 21 + "\n" if self.output_style == "unicode" else "-" * 21 + "\n")
            
            # Write the formatted table
            self._write_support_table(f, has_ml, has_bayesian, has_parsimony, has_posterior, has_bootstrap)
            
            # Support legend
            f.write("\nSupport Legend: *** = p < 0.001, ** = p < 0.01, * = p < 0.05, ns = not significant\n")
            f.write("BD = Bayes Decay (log scale), Post.Prob = Posterior Probability\n")
            
            # Clade details section
            f.write("\nClade Details\n")
            f.write("─" * 13 + "\n" if self.output_style == "unicode" else "-" * 13 + "\n")
            
            for clade_id, data in sorted(self.decay_indices.items()):
                taxa_list = sorted(data.get('taxa', []))
                support_levels = []
                
                if has_ml and data.get('AU_pvalue') is not None:
                    p = data.get('AU_pvalue', 1.0)
                    if p < 0.05:
                        support_levels.append("ML")
                
                if has_bayesian and data.get('bayes_decay') is not None:
                    bd = data.get('bayes_decay', 0)
                    if bd > 5:  # Strong support threshold for BD
                        support_levels.append("Bayes")
                        
                if has_parsimony and data.get('pars_decay') is not None:
                    pd = data.get('pars_decay', 0)
                    if pd > 3:
                        support_levels.append("Parsimony")
                
                support_str = f"Strong support: {'/'.join(support_levels)}" if support_levels else "Weak support"
                
                f.write(f"→ {clade_id} ({support_str})\n")
                
                # Format taxa list with wrapping
                taxa_str = "  Taxa: "
                line_len = len(taxa_str)
                for i, taxon in enumerate(taxa_list):
                    if i > 0:
                        if line_len + len(taxon) + 2 > 80:  # Wrap at 80 chars
                            f.write(",\n        ")
                            line_len = 8
                        else:
                            f.write(", ")
                            line_len += 2
                    f.write(taxon)
                    line_len += len(taxon)
                f.write("\n\n")
        
        logger.info(f"Results written to {self._get_display_path(output_path)}")
    
    def write_results(self, output_path: Path):
        if not self.decay_indices:
            logger.warning("No branch support results to write.")
            # Create an empty or minimal file? For now, just return.
            try:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with output_path.open('w') as f:
                    f.write("No branch support results calculated.\n")
                    if self.ml_likelihood is not None:
                        f.write(f"ML tree log-likelihood: {self.ml_likelihood:.6f}\n")
                return
            except Exception as e_write:
                logger.error(f"Failed to write empty results file {output_path}: {e_write}")
                return

        # Check if bootstrap analysis was performed
        has_bootstrap = hasattr(self, 'bootstrap_tree') and self.bootstrap_tree
        
        # Check which types of results we have
        has_ml = any(d.get('lnl_diff') is not None for d in self.decay_indices.values())
        has_bayesian = any(d.get('bayes_decay') is not None for d in self.decay_indices.values())
        has_parsimony = any(d.get('pars_decay') is not None for d in self.decay_indices.values())
        has_posterior = any(d.get('posterior_prob') is not None for d in self.decay_indices.values())

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('w') as f:
            f.write("panDecay Branch Support Analysis\n")
            f.write("=" * 30 + "\n\n")
            
            # Write appropriate header based on analysis type
            analysis_types = []
            if self.do_ml:
                analysis_types.append("Maximum Likelihood")
            if self.do_bayesian:
                analysis_types.append(f"Bayesian ({self.bayesian_software})")
            if self.do_parsimony:
                analysis_types.append("Parsimony")
            
            f.write(f"Analysis mode: {' + '.join(analysis_types)}\n")
            f.write("\n")
            
            # ML tree likelihood (if available)
            if self.ml_likelihood is not None:
                f.write(f"ML tree log-likelihood: {self.ml_likelihood:.6f}\n\n")
            
            f.write("Branch Support Values:\n")
            f.write("-" * 120 + "\n")

            # Build dynamic header based on available data
            header_parts = ["Clade_ID", "Num_Taxa"]
            
            if has_ml:
                header_parts.extend(["Constrained_lnL", "Delta_LnL", "AU_p-value", "Significant_AU (p<0.05)"])
            if has_parsimony:
                header_parts.append("Pars_Decay")
            if has_bayesian:
                header_parts.append("Bayes_ML_Diff")
                # Add normalized BD metrics if enabled
                if self.normalize_bd:
                    if "per_site" in self.bd_normalization_methods:
                        header_parts.append("BD_Per_Site")
                    if "relative" in self.bd_normalization_methods:
                        header_parts.append("BD_Relative")
                    if "signal_to_noise" in self.bd_normalization_methods:
                        header_parts.append("Signal_To_Noise")
                    if "effect_size" in self.bd_normalization_methods:
                        header_parts.append("Effect_Size")
                    if "effect_size_robust" in self.bd_normalization_methods:
                        header_parts.append("Effect_Size_Robust")
                    if "effect_size_weighted" in self.bd_normalization_methods:
                        header_parts.append("Effect_Size_Weighted")
                if has_posterior:
                    header_parts.append("Posterior_Prob")
            if has_bootstrap:
                header_parts.append("Bootstrap")
            header_parts.append("Taxa_List")
            
            f.write("\t".join(header_parts) + "\n")

            # Create mapping of taxa sets to bootstrap values if bootstrap analysis was performed
            bootstrap_values = {}
            if has_bootstrap:
                for node in self.bootstrap_tree.get_nonterminals():
                    if node.confidence is not None:
                        taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                        bootstrap_values[taxa_set] = node.confidence

            for clade_id, data in sorted(self.decay_indices.items()): # Sort for consistent output
                taxa_list = sorted(data.get('taxa', []))
                taxa_str = ",".join(taxa_list)
                num_taxa = len(taxa_list)

                row_parts = [clade_id, str(num_taxa)]
                
                # Add ML fields if present
                if has_ml:
                    c_lnl = data.get('constrained_lnl', 'N/A')
                    if isinstance(c_lnl, float): c_lnl = f"{c_lnl:.4f}"
                    elif c_lnl is None: c_lnl = 'N/A'
                    
                    lnl_d = data.get('lnl_diff', 'N/A')
                    if isinstance(lnl_d, float): lnl_d = f"{lnl_d:.4f}"
                    elif lnl_d is None: lnl_d = 'N/A'
                    
                    au_p = data.get('AU_pvalue', 'N/A')
                    if isinstance(au_p, float): au_p = f"{au_p:.4f}"
                    elif au_p is None: au_p = 'N/A'
                    
                    sig_au = data.get('significant_AU', 'N/A')
                    if isinstance(sig_au, bool): sig_au = "Yes" if sig_au else "No"
                    elif sig_au is None: sig_au = 'N/A'
                    
                    row_parts.extend([c_lnl, lnl_d, au_p, sig_au])
                
                # Add parsimony decay if present
                if has_parsimony:
                    pars_decay = data.get('pars_decay', 'N/A')
                    if isinstance(pars_decay, (int, float)): pars_decay = str(pars_decay)
                    elif pars_decay is None: pars_decay = 'N/A'
                    row_parts.append(pars_decay)
                
                # Add Bayesian fields if present
                if has_bayesian:
                    bayes_decay = data.get('bayes_decay', 'N/A')
                    if isinstance(bayes_decay, float): bayes_decay = f"{bayes_decay:.4f}"
                    elif bayes_decay is None: bayes_decay = 'N/A'
                    
                    row_parts.append(bayes_decay)
                    
                    # Add normalized BD metrics if enabled
                    if self.normalize_bd:
                        if "per_site" in self.bd_normalization_methods:
                            bd_per_site = data.get('bd_per_site', 'N/A')
                            if isinstance(bd_per_site, float): bd_per_site = f"{bd_per_site:.6f}"
                            elif bd_per_site is None: bd_per_site = 'N/A'
                            row_parts.append(bd_per_site)
                        
                        if "relative" in self.bd_normalization_methods:
                            bd_relative = data.get('bd_relative', 'N/A')
                            if isinstance(bd_relative, float): bd_relative = f"{bd_relative:.6f}"
                            elif bd_relative is None: bd_relative = 'N/A'
                            row_parts.append(bd_relative)
                        
                        if "signal_to_noise" in self.bd_normalization_methods:
                            signal_to_noise = data.get('signal_to_noise', 'N/A')
                            if isinstance(signal_to_noise, float): signal_to_noise = f"{signal_to_noise:.4f}"
                            elif signal_to_noise is None: signal_to_noise = 'N/A'
                            row_parts.append(signal_to_noise)
                        
                        if "effect_size" in self.bd_normalization_methods:
                            effect_size = data.get('effect_size', 'N/A')
                            if isinstance(effect_size, float): effect_size = f"{effect_size:.4f}"
                            elif effect_size is None: effect_size = 'N/A'
                            row_parts.append(effect_size)
                        
                        if "effect_size_robust" in self.bd_normalization_methods:
                            effect_size_robust = data.get('effect_size_robust', 'N/A')
                            if isinstance(effect_size_robust, float): effect_size_robust = f"{effect_size_robust:.4f}"
                            elif effect_size_robust is None: effect_size_robust = 'N/A'
                            row_parts.append(effect_size_robust)
                        
                        if "effect_size_weighted" in self.bd_normalization_methods:
                            effect_size_weighted = data.get('effect_size_weighted', 'N/A')
                            if isinstance(effect_size_weighted, float): effect_size_weighted = f"{effect_size_weighted:.4f}"
                            elif effect_size_weighted is None: effect_size_weighted = 'N/A'
                            row_parts.append(effect_size_weighted)
                    
                    # Add posterior probability if present
                    if has_posterior:
                        post_prob = data.get('posterior_prob', 'N/A')
                        if isinstance(post_prob, float): post_prob = f"{post_prob:.2f}"
                        elif post_prob is None: post_prob = 'N/A'
                        row_parts.append(post_prob)

                # Add bootstrap value if available
                if has_bootstrap:
                    taxa_set = frozenset(taxa_list)
                    bs_val = bootstrap_values.get(taxa_set, "N/A")
                    # Convert any numeric type to string
                    if bs_val != "N/A" and bs_val is not None:
                        try:
                            bs_val = f"{int(float(bs_val))}"
                        except (ValueError, TypeError):
                            bs_val = str(bs_val)
                    elif bs_val is None:
                        bs_val = "N/A"
                    row_parts.append(bs_val)

                row_parts.append(taxa_str)
                # Ensure all items are strings before joining
                row_parts = [str(item) for item in row_parts]
                f.write("\t".join(row_parts) + "\n")

        logger.info(f"Results written to {self._get_display_path(output_path)}")

    def generate_detailed_report(self, output_path: Path):
        # Basic check
        if not self.decay_indices and self.ml_likelihood is None and not hasattr(self, 'bayes_marginal_likelihood'):
            logger.warning("No results available to generate detailed report.")
            try:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with output_path.open('w') as f: f.write("# ML-Decay Report\n\nNo analysis results to report.\n")
                return
            except Exception as e_write:
                 logger.error(f"Failed to write empty detailed report {output_path}: {e_write}")
                 return

        # Check which types of results we have
        has_ml = any(d.get('lnl_diff') is not None for d in self.decay_indices.values()) if self.decay_indices else False
        has_bayesian = any(d.get('bayes_decay') is not None for d in self.decay_indices.values()) if self.decay_indices else False
        has_bootstrap = hasattr(self, 'bootstrap_tree') and self.bootstrap_tree

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('w') as f:
            # Title based on analysis type
            if has_ml and has_bayesian:
                f.write(f"# panDecay Branch Support Analysis Report (v{VERSION})\n\n")
            elif has_bayesian:
                f.write(f"# panDecay Bayesian Branch Support Analysis Report (v{VERSION})\n\n")
            else:
                f.write(f"# panDecay ML Branch Support Analysis Report (v{VERSION})\n\n")
                
            f.write(f"Date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            f.write("## Analysis Parameters\n\n")
            f.write(f"- Alignment file: `{self.alignment_file.name}`\n")
            f.write(f"- Data type: `{self.data_type}`\n")
            f.write(f"- Analysis mode: `{self.analysis_mode}`\n")
            
            # ML parameters
            if has_ml or self.do_ml:
                if self.user_paup_block:
                    f.write("- ML Model: User-defined PAUP* block\n")
                else:
                    f.write(f"- ML Model string: `{self.model_str}`\n")
                    f.write(f"- PAUP* `lset` command: `{self.paup_model_cmds}`\n")
            
            # Bayesian parameters
            if has_bayesian or self.do_bayesian:
                f.write(f"- Bayesian software: `{self.bayesian_software}`\n")
                f.write(f"- Bayesian model: `{self.bayes_model}`\n")
                f.write(f"- MCMC generations: `{self.bayes_ngen}`\n")
                f.write(f"- Burnin fraction: `{self.bayes_burnin}`\n")
                # Report the actual method being used
                ml_method = "stepping-stone" if self.marginal_likelihood == "ss" else "harmonic mean"
                f.write(f"- Marginal likelihood method: `{ml_method}`\n")
                
            if has_bootstrap:
                f.write("- Bootstrap analysis: Performed\n")

            f.write("\n## Summary Statistics\n\n")
            
            # ML statistics
            if has_ml:
                ml_l = self.ml_likelihood if self.ml_likelihood is not None else "N/A"
                if isinstance(ml_l, float): ml_l = f"{ml_l:.6f}"
                f.write(f"- ML tree log-likelihood: **{ml_l}**\n")
                
            # Bayesian statistics
            if has_bayesian and hasattr(self, 'bayes_marginal_likelihood'):
                bayes_ml = self.bayes_marginal_likelihood if self.bayes_marginal_likelihood is not None else "N/A"
                if isinstance(bayes_ml, float): bayes_ml = f"{bayes_ml:.6f}"
                f.write(f"- Bayesian marginal likelihood: **{bayes_ml}**\n")
                
            f.write(f"- Number of internal branches tested: {len(self.decay_indices)}\n")

            if self.decay_indices:
                # ML-specific statistics
                if has_ml:
                    lnl_diffs = [d['lnl_diff'] for d in self.decay_indices.values() if d.get('lnl_diff') is not None]
                    if lnl_diffs:
                        f.write(f"- Avg ML log-likelihood difference (constrained vs ML): {np.mean(lnl_diffs):.4f}\n")
                        f.write(f"- Min ML log-likelihood difference: {min(lnl_diffs):.4f}\n")
                        f.write(f"- Max ML log-likelihood difference: {max(lnl_diffs):.4f}\n")

                    au_pvals = [d['AU_pvalue'] for d in self.decay_indices.values() if d.get('AU_pvalue') is not None]
                    if au_pvals:
                        sig_au_count = sum(1 for p in au_pvals if p < 0.05)
                        f.write(f"- Branches with significant AU support (p < 0.05): {sig_au_count} / {len(au_pvals)} evaluated\n")
                
                # Bayesian-specific statistics
                if has_bayesian:
                    bayes_decays = [d['bayes_decay'] for d in self.decay_indices.values() if d.get('bayes_decay') is not None]
                    if bayes_decays:
                        f.write(f"- Avg Bayesian decay (marginal lnL difference): {np.mean(bayes_decays):.4f}\n")
                        f.write(f"- Min Bayesian decay: {min(bayes_decays):.4f}\n")
                        f.write(f"- Max Bayesian decay: {max(bayes_decays):.4f}\n")
                        
                        # Check for negative values and add warning
                        negative_count = sum(1 for bd in bayes_decays if bd < 0)
                        if negative_count > 0:
                            f.write(f"\n**⚠️ WARNING**: {negative_count}/{len(bayes_decays)} branches have negative Bayes Decay values.\n")
                            f.write("This suggests potential issues with:\n")
                            f.write("- MCMC convergence (consider increasing --bayes-ngen)\n")
                            f.write("- Marginal likelihood estimation reliability\n")
                            f.write("- Model specification\n\n")
                        
                        # Report on BD distribution (basic statistics only)
                        f.write(f"\n**Bayesian Decay Distribution**:\n")
                        f.write(f"- Mean BD: {np.mean(bayes_decays):.3f}\n")
                        f.write(f"- Median BD: {np.median(bayes_decays):.3f}\n")
                        f.write(f"- Standard deviation: {np.std(bayes_decays):.3f}\n")
                        f.write(f"- Min BD: {min(bayes_decays):.3f}\n")
                        f.write(f"- Max BD: {max(bayes_decays):.3f}\n")
                        f.write(f"- Range: {max(bayes_decays) - min(bayes_decays):.3f}\n")
                        if negative_count > 0:
                            f.write(f"- Negative BD values: {negative_count} branches (see warning above)\n")
                        
                        # Add percentile information for relative comparison
                        sorted_bds = sorted([bd for bd in bayes_decays if bd >= 0])  # Exclude negative values from percentiles
                        if sorted_bds:
                            p75 = np.percentile(sorted_bds, 75)
                            p25 = np.percentile(sorted_bds, 25)
                            f.write(f"- 75th percentile: {p75:.3f} (branches above this have relatively high support)\n")
                            f.write(f"- 25th percentile: {p25:.3f} (branches below this have relatively low support)\n")
                        
                        f.write(f"\n**Note**: BD values closely approximating ML log-likelihood differences is expected behavior in phylogenetic topology testing.\n")
                        f.write(f"**Important**: BD values scale with dataset characteristics. Compare BD values only within this analysis, not across studies.\n\n")
                        
                        # Add normalized BD statistics if enabled
                        if self.normalize_bd:
                            f.write(f"**Normalized BD Statistics** (for cross-study comparisons):\n")
                            
                            if "per_site" in self.bd_normalization_methods:
                                bd_per_site_values = [d['bd_per_site'] for d in self.decay_indices.values() if d.get('bd_per_site') is not None]
                                if bd_per_site_values:
                                    f.write(f"- BD/site statistics:\n")
                                    f.write(f"  - Mean: {np.mean(bd_per_site_values):.6f}\n")
                                    f.write(f"  - Median: {np.median(bd_per_site_values):.6f}\n")
                                    f.write(f"  - Min: {min(bd_per_site_values):.6f}\n")
                                    f.write(f"  - Max: {max(bd_per_site_values):.6f}\n")
                                    
                            if "relative" in self.bd_normalization_methods:
                                bd_relative_values = [d['bd_relative'] for d in self.decay_indices.values() if d.get('bd_relative') is not None]
                                if bd_relative_values:
                                    f.write(f"- Avg BD%: {np.mean(bd_relative_values)*100:.3f}%\n")
                                    f.write(f"- Min BD%: {min(bd_relative_values)*100:.3f}%\n")
                                    f.write(f"- Max BD%: {max(bd_relative_values)*100:.3f}%\n")
                            
                            # Effect Size Statistics
                            effect_size_methods = [method for method in self.bd_normalization_methods if method.startswith("effect_size")]
                            if effect_size_methods:
                                for method in effect_size_methods:
                                    effect_size_values = [d[method] for d in self.decay_indices.values() if d.get(method) is not None]
                                    if effect_size_values:
                                        method_name = method.replace("_", " ").title()
                                        f.write(f"- {method_name} distribution:\n")
                                        f.write(f"  - Mean: {np.mean(effect_size_values):.3f}\n")
                                        f.write(f"  - Min: {min(effect_size_values):.3f}\n")
                                        f.write(f"  - Max: {max(effect_size_values):.3f}\n")
                                        
                                        # Add Cohen's d interpretation distribution
                                        interpretations = {}
                                        for value in effect_size_values:
                                            _, _, interpretation = self._get_effect_size_interpretation_scale(abs(value))
                                            if interpretation not in interpretations:
                                                interpretations[interpretation] = 0
                                            interpretations[interpretation] += 1
                                        
                                        f.write(f"  - Effect size interpretation distribution:\n")
                                        for interp, count in sorted(interpretations.items()):
                                            f.write(f"    - {interp}: {count} branches\n")
                                
                                # Add effect size explanation
                                f.write(f"- Effect Size Interpretation (Cohen's d framework):\n")
                                f.write(f"  - Effect size = BD / SD(site signals) - measures signal-to-noise ratio\n")
                                f.write(f"  - Enables cross-study comparison by normalizing for dataset variability\n")
                                f.write(f"  - Small effect (0.2-0.5): Weak phylogenetic signal\n")
                                f.write(f"  - Medium effect (0.5-0.8): Moderate phylogenetic signal\n")
                                f.write(f"  - Large effect (0.8-1.2): Strong phylogenetic signal\n")
                                f.write(f"  - Very large effect (≥1.2): Very strong phylogenetic signal\n")
                            
                            f.write(f"\n")
                        
                    # Note: We provide relative BD statistics rather than misleading absolute thresholds
                
                # Add convergence diagnostics section if available
                if has_bayesian and hasattr(self, 'convergence_diagnostics'):
                    f.write("\n## Bayesian Convergence Diagnostics\n\n")
                    
                    # Summary across all runs
                    all_ess = []
                    all_psrf = []
                    all_asdsf = []
                    convergence_issues = []
                    
                    for run_id, conv_data in self.convergence_diagnostics.items():
                        if conv_data['min_ess'] is not None:
                            all_ess.append(conv_data['min_ess'])
                        if conv_data['max_psrf'] is not None:
                            all_psrf.append(conv_data['max_psrf'])
                        if conv_data['asdsf'] is not None:
                            all_asdsf.append(conv_data['asdsf'])
                        if not conv_data['converged']:
                            convergence_issues.append(run_id)
                    
                    if all_ess:
                        f.write(f"- Minimum ESS across all runs: {min(all_ess):.0f} (threshold: {self.min_ess})\n")
                    if all_psrf:
                        f.write(f"- Maximum PSRF across all runs: {max(all_psrf):.3f} (threshold: {self.max_psrf})\n")
                    if all_asdsf:
                        f.write(f"- Final ASDSF range: {min(all_asdsf):.6f} - {max(all_asdsf):.6f} (threshold: {self.max_asdsf})\n")
                    
                    if convergence_issues:
                        f.write(f"\n**⚠️ WARNING**: {len(convergence_issues)} runs did not meet convergence criteria:\n")
                        for run_id in convergence_issues[:5]:  # Show first 5
                            f.write(f"  - {run_id}\n")
                        if len(convergence_issues) > 5:
                            f.write(f"  - ... and {len(convergence_issues) - 5} more\n")
                        f.write("\nConsider:\n")
                        f.write("- Increasing MCMC generations (--bayes-ngen)\n")
                        f.write("- Running longer burnin (--bayes-burnin)\n")
                        f.write("- Using more chains (--bayes-chains)\n")
                        f.write("- Checking model specification\n")

            f.write("\n## Detailed Branch Support Results\n\n")

            # Build dynamic table header based on available data
            header_parts = ["| Clade ID | Taxa Count "]
            separator_parts = ["|----------|------------ "]
            
            if has_ml:
                header_parts.extend(["| Constrained lnL | ΔlnL (from ML) | AU p-value | Significant (AU) "])
                separator_parts.extend(["|-----------------|------------------|------------|-------------------- "])
            
            if has_bayesian:
                header_parts.append("| Bayes Decay ")
                separator_parts.append("|------------- ")
                
            if has_bootstrap:
                header_parts.append("| Bootstrap ")
                separator_parts.append("|----------- ")
                
            header_parts.append("| Included Taxa (sample) |\n")
            separator_parts.append("|--------------------------|\n")
            
            f.write("".join(header_parts))
            f.write("".join(separator_parts))

            # Get bootstrap values if bootstrap analysis was performed
            bootstrap_values = {}
            if has_bootstrap:
                for node in self.bootstrap_tree.get_nonterminals():
                    if node.confidence is not None:
                        taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                        bootstrap_values[taxa_set] = node.confidence

            for clade_id, data in sorted(self.decay_indices.items()):
                taxa_list = sorted(data.get('taxa', []))
                taxa_count = len(taxa_list)
                taxa_sample = ", ".join(taxa_list[:3]) + ('...' if taxa_count > 3 else '')

                # Build the table row
                row_parts = [f"| {clade_id} | {taxa_count} "]
                
                # ML fields
                if has_ml:
                    c_lnl = data.get('constrained_lnl', 'N/A')
                    if isinstance(c_lnl, float): c_lnl = f"{c_lnl:.4f}"
                    
                    lnl_d = data.get('lnl_diff', 'N/A')
                    if isinstance(lnl_d, float): lnl_d = f"{lnl_d:.4f}"
                    
                    au_p = data.get('AU_pvalue', 'N/A')
                    if isinstance(au_p, float): au_p = f"{au_p:.4f}"
                    
                    sig_au = data.get('significant_AU', 'N/A')
                    if isinstance(sig_au, bool): sig_au = "**Yes**" if sig_au else "No"
                    
                    row_parts.append(f"| {c_lnl} | {lnl_d} | {au_p} | {sig_au} ")
                
                # Bayesian fields
                if has_bayesian:
                    bayes_d = data.get('bayes_decay', 'N/A')
                    if isinstance(bayes_d, float): bayes_d = f"{bayes_d:.4f}"
                    
                    row_parts.append(f"| {bayes_d} ")

                # Bootstrap column if available
                if has_bootstrap:
                    taxa_set = frozenset(taxa_list)
                    bs_val = bootstrap_values.get(taxa_set, "N/A")
                    # Convert any numeric type to string
                    if bs_val != "N/A" and bs_val is not None:
                        try:
                            bs_val = f"{int(float(bs_val))}"
                        except (ValueError, TypeError):
                            bs_val = str(bs_val)
                    elif bs_val is None:
                        bs_val = "N/A"
                    row_parts.append(f"| {bs_val} ")

                row_parts.append(f"| {taxa_sample} |\n")
                f.write("".join(row_parts))

            f.write("\n## Interpretation Guide\n\n")
            
            if has_ml:
                f.write("### ML Analysis\n")
                f.write("- **ΔlnL (from ML)**: Log-likelihood difference between the constrained tree (without the clade) and the ML tree. Calculated as: constrained_lnL - ML_lnL. Larger positive values indicate stronger support for the clade's presence in the ML tree.\n")
                f.write("- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.\n\n")
                
            if has_bayesian:
                f.write("### Bayesian Analysis\n")
                f.write("- **Bayes Decay (BD)**: Marginal log-likelihood difference (unconstrained - constrained). This is the primary metric for Bayesian support.\n")
                f.write("  - **Key insight**: In phylogenetic topology testing, BD values typically closely approximate ML log-likelihood differences\n")
                f.write("  - **Important**: BD values scale with dataset characteristics (alignment length, sequence diversity, substitution rates)\n")
                f.write("  - **Interpretation**: Compare BD values only within this analysis. Higher BD values indicate stronger relative support\n")
                f.write("  - **Cross-study comparisons**: Use normalized metrics (effect sizes) rather than raw BD values\n")
                f.write("  - **Why BD ≈ ΔlnL**: When comparing models that differ only by a topological constraint, the marginal likelihood is dominated by the likelihood component\n")
                
                # Add normalized BD documentation if enabled
                if self.normalize_bd:
                    f.write("  - **Normalized BD Metrics** (for cross-study comparisons):\n")
                    if "per_site" in self.bd_normalization_methods:
                        f.write("    - **BD/site**: Per-site Bayesian Decay (BD ÷ alignment length) - enables comparison across different alignment lengths\n")
                    if "relative" in self.bd_normalization_methods:
                        f.write("    - **BD%**: Relative BD as percentage of unconstrained marginal likelihood\n")
                        f.write("      - Provides context relative to total likelihood magnitude\n")
                        f.write("      - Useful for comparing datasets with different likelihood scales\n")
                    f.write("    - **Why normalize?**: Raw BD values scale with dataset size, making cross-study comparisons misleading\n")
                    f.write("    - **Recommendation**: Use normalized values for publications and meta-analyses\n")
                f.write("  - **Negative values** suggest the constrained analysis had higher marginal likelihood, which may indicate:\n")
                if self.marginal_likelihood == "ss":
                    f.write("    - Poor chain convergence or insufficient MCMC sampling\n")
                    f.write("    - Complex posterior distribution requiring more steps\n")
                    f.write("    - Genuine lack of support for the clade\n")
                else:
                    f.write("    - Harmonic mean estimator limitations (notoriously unreliable)\n")
                    f.write("    - Poor MCMC convergence (try increasing generations)\n")
                    f.write("    - Genuine lack of support for the clade\n")
                    f.write("    - **Consider using stepping-stone sampling (--marginal-likelihood ss) for more reliable estimates**\n")
                
            if has_bootstrap:
                f.write("- **Bootstrap**: Bootstrap support value (percentage of bootstrap replicates in which the clade appears). Higher values (e.g., > 70) suggest stronger support for the clade.\n")
            
            # Add detailed explanation about BD vs ML differences when both analyses are present
            if has_ml and has_bayesian:
                f.write("\n## Understanding BD ≈ ΔlnL in Phylogenetics\n\n")
                f.write("You may notice that Bayesian Decay (BD) values closely approximate the ML log-likelihood differences (ΔlnL). ")
                f.write("This is **expected behavior** in phylogenetic topology testing, not an anomaly. Here's why:\n\n")
                f.write("1. **Identical Models**: The constrained and unconstrained analyses use the same substitution model, ")
                f.write("differing only in whether a specific clade is allowed to exist.\n\n")
                f.write("2. **Likelihood Dominance**: When data strongly support a topology, the marginal likelihood ")
                f.write("(which integrates over all parameters) becomes dominated by the likelihood at the optimal parameter values.\n\n")
                f.write("3. **Minimal Prior Effects**: Since both analyses explore nearly identical parameter spaces ")
                f.write("(same model parameters, only different tree topologies), the prior's influence is minimal.\n\n")
                f.write("**Practical Implications**:\n")
                f.write("- Similar BD and ΔlnL values confirm your analyses are working correctly\n")
                f.write("- BD values scale with dataset characteristics - do not compare absolute values across studies\n")
                f.write("- For cross-study comparisons, use normalized metrics (effect sizes) instead of raw BD values\n")
                f.write("- Compare relative BD values across branches within this analysis to identify well-supported clades\n")
                
        logger.info(f"Detailed report written to {self._get_display_path(output_path)}")

    def write_site_analysis_results(self, output_dir: Path):
        """
        Write site-specific likelihood analysis results to files.

        Args:
            output_dir: Directory to save the site analysis files
        """
        if not self.decay_indices:
            logger.warning("No decay indices available for site analysis output.")
            return

        # Check if any clade has site data
        has_site_data = any('site_data' in data for data in self.decay_indices.values())
        if not has_site_data:
            logger.warning("No site-specific analysis data available to write.")
            return

        output_dir.mkdir(parents=True, exist_ok=True)

        # Create a summary file for all branches
        summary_path = output_dir / "site_analysis_summary.txt"
        with summary_path.open('w') as f:
            f.write("Branch Site Analysis Summary\n")
            f.write("=========================\n\n")
            f.write("Clade_ID\tSupporting_Sites\tConflicting_Sites\tNeutral_Sites\tSupport_Ratio\tSum_Supporting_Delta\tSum_Conflicting_Delta\tWeighted_Support_Ratio\n")

            for clade_id, data in sorted(self.decay_indices.items()):
                if 'site_data' not in data:
                    continue

                supporting = data.get('supporting_sites', 0)
                conflicting = data.get('conflicting_sites', 0)
                neutral = data.get('neutral_sites', 0)
                ratio = data.get('support_ratio', 0.0)
                sum_supporting = data.get('sum_supporting_delta', 0.0)
                sum_conflicting = data.get('sum_conflicting_delta', 0.0)
                weighted_ratio = data.get('weighted_support_ratio', 0.0)

                if ratio == float('inf'):
                    ratio_str = "Inf"
                else:
                    ratio_str = f"{ratio:.4f}"

                if weighted_ratio == float('inf'):
                    weighted_ratio_str = "Inf"
                else:
                    weighted_ratio_str = f"{weighted_ratio:.4f}"

                f.write(f"{clade_id}\t{supporting}\t{conflicting}\t{neutral}\t{ratio_str}\t{sum_supporting:.4f}\t{sum_conflicting:.4f}\t{weighted_ratio_str}\n")

        logger.info(f"Site analysis summary written to {self._get_display_path(summary_path)}")

        # For each branch, write detailed site data
        for clade_id, data in self.decay_indices.items():
            if 'site_data' not in data:
                continue

            site_data_path = output_dir / f"site_data_{clade_id}.txt"
            with site_data_path.open('w') as f:
                f.write(f"Site-Specific Likelihood Analysis for {clade_id}\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Supporting sites: {data.get('supporting_sites', 0)}\n")
                f.write(f"Conflicting sites: {data.get('conflicting_sites', 0)}\n")
                f.write(f"Neutral sites: {data.get('neutral_sites', 0)}\n")
                f.write(f"Support ratio: {data.get('support_ratio', 0.0):.4f}\n")
                f.write(f"Sum of supporting deltas: {data.get('sum_supporting_delta', 0.0):.4f}\n")
                f.write(f"Sum of conflicting deltas: {data.get('sum_conflicting_delta', 0.0):.4f}\n")
                f.write(f"Weighted support ratio: {data.get('weighted_support_ratio', 0.0):.4f}\n\n")
                f.write("Site\tML_Tree_lnL\tConstrained_lnL\tDelta_lnL\tSupports_Branch\n")

                # Make sure site_data is a dictionary with entries for each site
                site_data = data.get('site_data', {})
                if isinstance(site_data, dict) and site_data:
                    for site_num, site_info in sorted(site_data.items()):
                        # Safely access each field with a default
                        ml_lnl = site_info.get('lnL_ML', 0.0)
                        constrained_lnl = site_info.get('lnL_constrained', 0.0)
                        delta_lnl = site_info.get('delta_lnL', 0.0)
                        supports = site_info.get('supports_branch', False)

                        f.write(f"{site_num}\t{ml_lnl:.6f}\t{constrained_lnl:.6f}\t{delta_lnl:.6f}\t{supports}\n")

            logger.info(f"Detailed site data for {clade_id} written to {self._get_display_path(site_data_path)}")

        # Generate site analysis visualizations
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np

            # Get visualization options
            viz_format = getattr(self, 'viz_format', 'png')

            for clade_id, data in self.decay_indices.items():
                if 'site_data' not in data:
                    continue

                # Extract data for plotting
                site_data = data.get('site_data', {})
                if not site_data:
                    continue

                site_nums = sorted(site_data.keys())
                deltas = [site_data[site]['delta_lnL'] for site in site_nums if 'delta_lnL' in site_data[site]]

                if not deltas:
                    continue

                # Get taxa in this clade for visualization
                clade_taxa = data.get('taxa', [])

                # Prepare taxa list for title display
                if len(clade_taxa) <= 3:
                    taxa_display = ", ".join(clade_taxa)
                else:
                    taxa_display = f"{', '.join(sorted(clade_taxa)[:3])}... (+{len(clade_taxa)-3} more)"

                # Create standard site analysis plot
                fig = plt.figure(figsize=(12, 6))
                ax_main = fig.add_subplot(111)

                # Create the main bar plot
                bar_colors = ['green' if d < 0 else 'red' for d in deltas]
                ax_main.bar(range(len(deltas)), deltas, color=bar_colors, alpha=0.7)

                # Add x-axis ticks at reasonable intervals
                if len(site_nums) > 50:
                    tick_interval = max(1, len(site_nums) // 20)
                    tick_positions = range(0, len(site_nums), tick_interval)
                    tick_labels = [site_nums[i] for i in tick_positions if i < len(site_nums)]
                    ax_main.set_xticks(tick_positions)
                    ax_main.set_xticklabels(tick_labels, rotation=45)
                else:
                    ax_main.set_xticks(range(len(site_nums)))
                    ax_main.set_xticklabels(site_nums, rotation=45)

                # Add reference line at y=0
                ax_main.axhline(y=0, color='black', linestyle='-', alpha=0.3)

                # Add title that includes some taxa information
                ax_main.set_title(f"Site-Specific Likelihood Differences for {clade_id} ({taxa_display})")
                ax_main.set_xlabel("Site Position")
                ax_main.set_ylabel("Delta lnL (ML - Constrained)")

                # Add summary info text box
                support_ratio = data.get('support_ratio', 0.0)
                weighted_ratio = data.get('weighted_support_ratio', 0.0)

                ratio_text = "Inf" if support_ratio == float('inf') else f"{support_ratio:.2f}"
                weighted_text = "Inf" if weighted_ratio == float('inf') else f"{weighted_ratio:.2f}"

                info_text = (
                    f"Supporting sites: {data.get('supporting_sites', 0)}\n"
                    f"Conflicting sites: {data.get('conflicting_sites', 0)}\n"
                    f"Support ratio: {ratio_text}\n"
                    f"Weighted ratio: {weighted_text}"
                )

                # Add text box with summary info
                ax_main.text(
                    0.02, 0.95, info_text,
                    transform=ax_main.transAxes,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
                )

                plt.tight_layout()

                # Save plot in the requested format
                plot_path = output_dir / f"site_plot_{clade_id}.{viz_format}"
                plt.savefig(str(plot_path), dpi=150, format=viz_format)
                plt.close(fig)

                logger.info(f"Site-specific likelihood plot for {clade_id} saved to {plot_path}")

                # Optional: Create a histogram of delta values
                plt.figure(figsize=(10, 5))
                sns.histplot(deltas, kde=True, bins=30)
                plt.axvline(x=0, color='black', linestyle='--')
                plt.title(f"Distribution of Site Likelihood Differences for {clade_id}")
                plt.xlabel("Delta lnL (ML - Constrained)")
                plt.tight_layout()

                hist_path = output_dir / f"site_hist_{clade_id}.{viz_format}"
                plt.savefig(str(hist_path), dpi=150, format=viz_format)
                plt.close()

                logger.info(f"Site likelihood histogram for {clade_id} saved to {hist_path}")

                if not self.debug and not self.keep_files:
                    # Clean up tree files
                    for file_path in output_dir.glob("tree_*.nwk*"):
                        try:
                            file_path.unlink()
                            logger.debug(f"Deleted tree file: {file_path}")
                        except Exception as e:
                            logger.warning(f"Failed to delete tree file {file_path}: {e}")

        except ImportError:
            logger.warning("Matplotlib/seaborn not available for site analysis visualization.")
        except Exception as e:
            logger.error(f"Error creating site analysis visualizations: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Visualization error traceback: {traceback.format_exc()}")

    def visualize_support_distribution(self, output_path: Path, value_type="au", **kwargs):
        if not self.decay_indices: logger.warning("No data for support distribution plot."); return
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns # numpy is usually a dependency of seaborn/matplotlib

            vals = []
            for data in self.decay_indices.values():
                if value_type == "au" and data.get('AU_pvalue') is not None: vals.append(data['AU_pvalue'])
                elif value_type == "lnl" and data.get('lnl_diff') is not None: vals.append(abs(data['lnl_diff']))
            if not vals: logger.warning(f"No '{value_type}' values for distribution plot."); return

            plt.figure(figsize=(kwargs.get('width',10), kwargs.get('height',6)))
            sns.histplot(vals, kde=True)
            title, xlabel = "", ""
            if value_type == "au":
                plt.axvline(0.05, color='r', linestyle='--', label='p=0.05 threshold')
                title, xlabel = 'Distribution of AU Test p-values', 'AU p-value'
            else: # lnl
                mean_val = np.mean(vals)
                plt.axvline(mean_val, color='g', linestyle='--', label=f'Mean diff ({mean_val:.2f})')
                title, xlabel = 'Distribution of abs(Log-Likelihood Differences)', 'abs(LNL Difference)'
            plt.title(title); plt.xlabel(xlabel); plt.ylabel('Frequency'); plt.legend(); plt.tight_layout()

            output_path.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_path), format=kwargs.get('format',"png"), dpi=300); plt.close()
            logger.info(f"Support distribution plot saved to {output_path}")
        except ImportError: logger.error("Matplotlib/Seaborn not found for visualization.")
        except Exception as e: logger.error(f"Failed support distribution plot: {e}")


    def cleanup_intermediate_files(self):
        """
        Clean up intermediate files that are not needed for final output.
        This includes temporary .cleaned tree files and other intermediate files.
        """
        if self.debug or self.keep_files:
            logger.info("Skipping intermediate file cleanup due to debug or keep_files flag")
            return

        logger.info("Cleaning up intermediate files...")

        # Delete files explicitly marked for cleanup
        for file_path in self._files_to_cleanup:
            try:
                if file_path.exists():
                    file_path.unlink()
                    logger.debug(f"Deleted intermediate file: {file_path}")
            except Exception as e:
                logger.warning(f"Failed to delete intermediate file {file_path}: {e}")

        # Clean up other known intermediate files
        intermediate_patterns = [
            "*.cleaned",  # Cleaned tree files
            "constraint_tree_*.tre",  # Constraint trees
            "site_lnl_*.txt",  # Site likelihood files
            "ml_tree_for_*_annotation.nwk",  # Temporary annotation tree files
        ]

        for pattern in intermediate_patterns:
            for file_path in self.temp_path.glob(pattern):
                try:
                    file_path.unlink()
                    logger.debug(f"Deleted intermediate file: {file_path}")
                except Exception as e:
                    logger.warning(f"Failed to delete intermediate file {file_path}: {e}")


# --- Main Execution Logic ---
def get_display_path(path):
    """Get a display-friendly path (relative if possible, otherwise absolute)."""
    try:
        return str(Path(path).relative_to(Path.cwd()))
    except ValueError:
        return str(path)

def print_runtime_parameters(args_ns, model_str_for_print):
    """Prints a summary of runtime parameters."""
    # (args_ns is the namespace from argparse.ArgumentParser)
    print("\n" + "=" * 80)
    print(f"panDecay: Phylogenetic Analysis using Decay Indices v{VERSION}")
    print("=" * 80)
    print("\nRUNTIME PARAMETERS:")
    print(f"  Alignment file:     {args_ns.alignment}") # Original string path is fine for print
    print(f"  Alignment format:   {args_ns.format}")
    print(f"  Data type:          {args_ns.data_type}")
    if args_ns.paup_block:
        print(f"  PAUP* settings:     User-provided block from '{args_ns.paup_block}'")
    else:
        print(f"  Model string:       {model_str_for_print}")
        # Further model details can be extracted from args_ns if needed
    print(f"\n  PAUP* executable:   {args_ns.paup}")
    print(f"  Threads for PAUP*:  {args_ns.threads}")
    if args_ns.starting_tree:
        print(f"  Starting tree:      {args_ns.starting_tree}")
    output_p = Path(args_ns.output) # Use Path for consistent name generation
    print("\nOUTPUT SETTINGS:")
    print(f"  Results file:       {output_p}")
    print(f"  Annotated trees:    {args_ns.tree}_au.nwk, {args_ns.tree}_delta_lnl.nwk, {args_ns.tree}_combined.nwk")
    print(f"  Detailed report:    {output_p.with_suffix('.md')}")
    if args_ns.temp: print(f"  Temp directory:     {args_ns.temp}")
    if args_ns.debug: print(f"  Debug mode:         Enabled (log: mldecay_debug.log, if configured)")
    if args_ns.keep_files: print(f"  Keep temp files:    Enabled")
    if args_ns.visualize:
        print("\nVISUALIZATIONS:")
        print(f"  Enabled, format:    {args_ns.viz_format}")
        print(f"  Tree plot:          {output_p.parent / (output_p.stem + '_tree.' + args_ns.viz_format)}")
    print("\n" + "=" * 80 + "\n")

    @staticmethod
    def read_paup_block(paup_block_file_path: Path):
        if not paup_block_file_path.is_file():
            logger.error(f"PAUP block file not found: {paup_block_file_path}")
            return None
        try:
            content = paup_block_file_path.read_text()
            # Regex captures content *between* "BEGIN PAUP;" and "END;" (case-insensitive)
            match = re.search(r'BEGIN\s+PAUP\s*;(.*?)\s*END\s*;', content, re.DOTALL | re.IGNORECASE)
            if match:
                paup_cmds = match.group(1).strip()
                if not paup_cmds: logger.warning(f"PAUP block in {paup_block_file_path} is empty.")
                return paup_cmds
            else:
                logger.error(f"No valid PAUP block (BEGIN PAUP; ... END;) in {paup_block_file_path}")
                return None
        except Exception as e:
            logger.error(f"Error reading PAUP block file {paup_block_file_path}: {e}")
            return None


def generate_config_template(filepath):
    """Generate a template configuration file with all options and comments."""
    template = '''# panDecay Configuration File
# Lines starting with # are comments and will be ignored
# This file allows you to specify all parameters for a panDecay analysis
# Format: parameter = value (spaces around = are optional)

# ==============================================================================
# BASIC INPUT/OUTPUT SETTINGS
# ==============================================================================

# Input alignment file (required)
alignment = my_sequences.fasta

# Alignment format (default: fasta)
# Options: fasta, phylip, nexus, clustal, stockholm
format = fasta

# Data type (default: dna)
# Options: dna, protein, discrete
data_type = dna

# Output file prefix (default: pan_decay_indices)
# This will generate files like: myanalysis.txt, myanalysis_tree.nwk, etc.
output = pan_decay_indices.txt

# Base name for annotated tree files (default: annotated_tree)
tree = annotated_tree

# Keep temporary files after analysis (default: false)
# Options: true, false
keep_files = false

# Enable debug mode with detailed logging (default: false)
debug = false

# Custom temporary directory (optional)
# If not specified, uses system temp directory
# temp = /path/to/temp

# ==============================================================================
# ANALYSIS MODE SETTINGS
# ==============================================================================

# Analysis type (default: ml)
# Options: ml, bayesian, parsimony, ml+parsimony, bayesian+parsimony, ml+bayesian, all
# Note: 'all' runs all three analysis types (ML + Bayesian + Parsimony)
# Use '+' to combine multiple analyses (e.g., ml+parsimony)
analysis = ml

# Perform bootstrap analysis (default: false)
bootstrap = false

# Number of bootstrap replicates (default: 100)
bootstrap_reps = 100

# Perform site-specific analysis (default: false)
site_analysis = false

# Generate visualizations (default: false)
# Requires matplotlib and seaborn
visualize = false

# Visualization format (default: png)
# Options: png, pdf, svg
viz_format = png

# Annotation type for visualization (default: lnl)
# Options: au, lnl
annotation = lnl


# ==============================================================================
# MODEL SETTINGS (for ML and Bayesian analyses)
# ==============================================================================

# Base substitution model
# DNA options: JC, K80, HKY, GTR
# Protein options: JTT, WAG, LG, Dayhoff
# Discrete options: Mk
model = GTR

# Add gamma rate heterogeneity (default: false)
gamma = false

# Fixed gamma shape parameter (optional)
# If not specified, will be estimated
# gamma_shape = 0.5

# Add proportion of invariable sites (default: false)
invariable = false

# Fixed proportion of invariable sites (optional)
# If not specified, will be estimated
# prop_invar = 0.2

# Base frequencies (optional)
# Options: equal, estimate, empirical
# base_freq = empirical

# Number of substitution types for DNA (optional)
# Options: 1, 2, 6
# nst = 6

# Protein-specific model (overrides base model for protein data)
# protein_model = WAG

# ==============================================================================
# COMPUTATIONAL SETTINGS
# ==============================================================================

# Number of threads for PAUP* (default: auto)
# Options: auto (uses cores-2), all, or specific number (e.g., 8)
threads = auto

# Path to PAUP* executable (default: paup)
# Specify full path if not in system PATH
paup = paup

# Starting tree file in Newick format (optional)
# starting_tree = my_starting_tree.nwk

# Custom PAUP* commands file (optional)
# This overrides most model settings above
# paup_block = custom_paup_commands.txt

# ==============================================================================
# BAYESIAN-SPECIFIC SETTINGS
# ==============================================================================

# Bayesian software (required if analysis includes bayesian)
# Options: mrbayes
bayesian_software = mrbayes

# Path to MrBayes executable (default: mb)
mrbayes_path = mb


# Model for Bayesian analysis (optional)
# If not specified, uses same as ML model
# bayes_model = GTR+G+I

# Number of MCMC generations (default: 1000000)
bayes_ngen = 1000000

# Burnin fraction 0-1 (default: 0.25)
bayes_burnin = 0.25

# Number of MCMC chains (default: 4)
bayes_chains = 4

# Sample frequency (default: 1000)
bayes_sample_freq = 1000

# Marginal likelihood estimation method (default: ss)
# Options: ss (stepping-stone), ps (path sampling), hm (harmonic mean)
marginal_likelihood = ss

# Stepping-stone sampling parameters
ss_alpha = 0.4
ss_nsteps = 50

# ==============================================================================
# CONVERGENCE CHECKING (MrBayes only)
# ==============================================================================

# Check MCMC convergence diagnostics (default: true)
check_convergence = true

# Minimum ESS (Effective Sample Size) threshold (default: 200)
# Values below this indicate poor mixing
min_ess = 200

# Maximum PSRF (Potential Scale Reduction Factor) threshold (default: 1.01)
# Values above this indicate lack of convergence between chains
max_psrf = 1.01

# Maximum ASDSF (Average Standard Deviation of Split Frequencies) (default: 0.01)
# Values above this indicate lack of convergence between independent runs
max_asdsf = 0.01

# Fail analysis if convergence criteria not met (default: false)
# If false, warnings are issued but analysis continues
convergence_strict = false

# ==============================================================================
# PARALLEL PROCESSING (MrBayes only)
# ==============================================================================

# Use MPI version of MrBayes (default: false)
use_mpi = false

# Number of MPI processors (default: number of chains)
# mpi_processors = 8

# Path to mpirun executable (default: mpirun)
mpirun_path = mpirun

# Use BEAGLE library for acceleration (default: false)
use_beagle = false

# BEAGLE device (default: auto)
# Options: cpu, gpu, auto
beagle_device = auto

# BEAGLE precision (default: double)
# Options: single, double
beagle_precision = double

# BEAGLE scaling (default: dynamic)
# Options: none, dynamic, always
beagle_scaling = dynamic

# ==============================================================================
# CONSTRAINT SETTINGS
# ==============================================================================

# Constraint mode (default: all)
# Options: all (test all branches), specific (test listed branches), exclude
constraint_mode = all

# Specify branches to test (optional)
# Format: semicolon-separated clades, each clade is comma-separated taxa
# Example: Homo_sapiens,Pan_troglodytes;Mus_musculus,Rattus_norvegicus
# test_branches = 

# Constraint file (optional)
# File containing constraint definitions (one per line)
# constraint_file = constraints.txt

# ==============================================================================
# CONSTRAINT DEFINITIONS
# ==============================================================================
# Define specific clades to test when constraint_mode = specific
# Format: clade_name = comma-separated list of taxa
# Note: These are only used when constraint_mode = specific

[constraints]
# Example constraints (uncomment and modify as needed):
# clade1 = Homo_sapiens,Pan_troglodytes,Gorilla_gorilla
# clade2 = Rattus_norvegicus,Mus_musculus
# clade3 = Gallus_gallus,Taeniopygia_guttata

# ==============================================================================
# ADVANCED SETTINGS
# ==============================================================================

# Parsimony model for discrete data (default: true for discrete)
# parsmodel = true

# Site rate variation (overrides gamma setting)
# Options: equal, gamma
# rates = gamma

# ==============================================================================
# NOTES
# ==============================================================================
# 1. Command-line arguments override settings in this file
# 2. Paths can be absolute or relative to the working directory
# 3. Boolean values can be: true/false, yes/no, on/off, 1/0
# 4. Comments can appear on their own line or after values
# 5. Section headers like [constraints] are optional for organization
# 6. Empty lines are ignored
'''
    
    try:
        with open(filepath, 'w') as f:
            f.write(template)
        logger.info(f"Template configuration file generated: {filepath}")
        logger.info("Edit this file with your parameters and run:")
        logger.info(f"  python panDecay.py --config {filepath}")
    except Exception as e:
        logger.error(f"Failed to generate config template: {e}")
        sys.exit(1)


def parse_config(config_file, args):
    """Parse configuration file and update args namespace with values."""
    config = configparser.ConfigParser()
    
    # Define how to convert string values to appropriate types
    def str_to_bool(value):
        """Convert string to boolean."""
        if value.lower() in ('true', 'yes', 'on', '1'):
            return True
        elif value.lower() in ('false', 'no', 'off', '0'):
            return False
        else:
            raise ValueError(f"Cannot convert '{value}' to boolean")
    
    try:
        config.read(config_file)
        
        # Process main section (unnamed section at top of file)
        if config.has_section('DEFAULT') or len(config.sections()) > 0 or len(config.defaults()) > 0:
            # Get all items from DEFAULT section and any named sections
            items = dict(config.defaults())
            
            # Map config file parameters to argparse arguments
            param_map = {
                'alignment': 'alignment',
                'format': 'format',
                'model': 'model',
                'gamma': 'gamma',
                'invariable': 'invariable',
                'paup': 'paup',
                'output': 'output',
                'tree': 'tree',
                'site_analysis': 'site_analysis',
                'data_type': 'data_type',
                'gamma_shape': 'gamma_shape',
                'prop_invar': 'prop_invar',
                'base_freq': 'base_freq',
                'rates': 'rates',
                'protein_model': 'protein_model',
                'nst': 'nst',
                'parsmodel': 'parsmodel',
                'threads': 'threads',
                'starting_tree': 'starting_tree',
                'paup_block': 'paup_block',
                'temp': 'temp',
                'keep_files': 'keep_files',
                'debug': 'debug',
                'analysis': 'analysis',
                'bootstrap': 'bootstrap',
                'bootstrap_reps': 'bootstrap_reps',
                'bayesian_software': 'bayesian_software',
                'mrbayes_path': 'mrbayes_path',
                'bayes_model': 'bayes_model',
                'bayes_ngen': 'bayes_ngen',
                'bayes_burnin': 'bayes_burnin',
                'bayes_chains': 'bayes_chains',
                'bayes_sample_freq': 'bayes_sample_freq',
                'marginal_likelihood': 'marginal_likelihood',
                'ss_alpha': 'ss_alpha',
                'ss_nsteps': 'ss_nsteps',
                'use_mpi': 'use_mpi',
                'mpi_processors': 'mpi_processors',
                'mpirun_path': 'mpirun_path',
                'use_beagle': 'use_beagle',
                'beagle_device': 'beagle_device',
                'beagle_precision': 'beagle_precision',
                'beagle_scaling': 'beagle_scaling',
                'visualize': 'visualize',
                'viz_format': 'viz_format',
                'annotation': 'annotation',
                'constraint_mode': 'constraint_mode',
                'test_branches': 'test_branches',
                'constraint_file': 'constraint_file',
                'check_convergence': 'check_convergence',
                'min_ess': 'min_ess',
                'max_psrf': 'max_psrf',
                'max_asdsf': 'max_asdsf',
                'convergence_strict': 'convergence_strict',
                'normalize_bd': 'normalize_bd',
                'bd_normalization_methods': 'bd_normalization_methods',
            }
            
            # Process each parameter
            for config_param, arg_param in param_map.items():
                if config_param in items:
                    value = items[config_param]
                    
                    # Skip if command line already provided this argument
                    if hasattr(args, arg_param):
                        current_value = getattr(args, arg_param)
                        # Check if it's a default value or explicitly set
                        if arg_param == 'alignment':
                            continue  # Special case - required positional arg
                        
                    # Convert types as needed
                    if arg_param in ['gamma', 'invariable', 'keep_files', 'debug', 
                                     'site_analysis', 'bootstrap', 'use_mpi', 'use_beagle',
                                     'visualize', 'check_convergence', 'convergence_strict',
                                     'normalize_bd']:
                        value = str_to_bool(value)
                    elif arg_param in ['gamma_shape', 'prop_invar', 'bayes_burnin', 'ss_alpha',
                                       'max_psrf', 'max_asdsf']:
                        value = float(value)
                    elif arg_param in ['nst', 'bootstrap_reps', 'bayes_ngen', 'bayes_chains',
                                       'bayes_sample_freq', 'ss_nsteps', 'mpi_processors', 'min_ess']:
                        value = int(value)
                    elif arg_param == 'bd_normalization_methods':
                        # Handle list parameter - split by comma/space
                        if isinstance(value, str):
                            value = value.replace(',', ' ').split()
                    
                    setattr(args, arg_param, value)
        
        # Handle constraints section if present
        if config.has_section('constraints'):
            constraints = {}
            for key, value in config.items('constraints'):
                if key.startswith('clade'):
                    constraints[key] = value
            
            # Store constraints for later use
            if constraints:
                args.config_constraints = constraints
        
        logger.info(f"Loaded configuration from: {config_file}")
        
    except Exception as e:
        logger.error(f"Error parsing configuration file: {e}")
        sys.exit(1)


def find_alignment_files(input_dir, file_pattern):
    """Find alignment files matching the pattern in the input directory."""
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    # Use glob to find matching files
    pattern_path = input_path / file_pattern
    files = glob.glob(str(pattern_path))
    
    # Filter to only include files (not directories)
    alignment_files = [Path(f) for f in files if Path(f).is_file()]
    
    if not alignment_files:
        raise FileNotFoundError(f"No files found matching pattern '{file_pattern}' in {input_dir}")
    
    logger.info(f"Found {len(alignment_files)} alignment files for batch processing")
    return sorted(alignment_files)


def process_single_alignment(args_dict):
    """Process a single alignment file. Used by multiprocessing pool."""
    alignment_file, args, output_dir = args_dict
    
    # Create a copy of args with the specific alignment file
    args_copy = argparse.Namespace(**vars(args))
    args_copy.alignment = str(alignment_file)
    
    # Set output file based on alignment filename
    alignment_name = alignment_file.stem
    output_file = output_dir / f"{alignment_name}_decay_indices.txt"
    args_copy.output = str(output_file)
    
    # Set tree output name
    args_copy.tree = alignment_name + "_tree"
    
    try:
        logger.info(f"Processing {alignment_file.name}...")
        
        # Create panDecayIndices instance and run analysis
        decay_calc = create_decay_calc_from_args(args_copy)
        
        # Run the analysis
        decay_calc.build_ml_tree()
        
        if decay_calc.ml_tree and decay_calc.ml_likelihood is not None:
            # Run bootstrap if requested
            if args_copy.bootstrap:
                decay_calc.run_bootstrap_analysis(num_replicates=args_copy.bootstrap_reps)
            
            decay_calc.calculate_decay_indices(perform_site_analysis=args_copy.site_analysis)
            
            # Write results
            decay_calc.write_formatted_results(Path(args_copy.output))
            
            # Generate report
            report_path = Path(args_copy.output).with_suffix(".md")
            decay_calc.generate_detailed_report(report_path)
            
            # Create annotated trees
            tree_files = decay_calc.annotate_trees(output_dir, args_copy.tree)
            
            # Handle visualizations and site analysis
            if args_copy.visualize:
                try:
                    import matplotlib, seaborn
                    viz_out_dir = output_dir
                    viz_base_name = alignment_name
                    viz_kwargs = {'width': 10, 'height': 6, 'format': args_copy.viz_format}
                    
                    decay_calc.visualize_support_distribution(
                        viz_out_dir / f"{viz_base_name}_dist_{args_copy.annotation}.{args_copy.viz_format}",
                        value_type=args_copy.annotation, **viz_kwargs)
                except ImportError:
                    logger.warning(f"Matplotlib/Seaborn not available for {alignment_file.name}")
            
            if args_copy.site_analysis and hasattr(decay_calc, 'decay_indices') and decay_calc.decay_indices:
                site_output_dir = output_dir / f"{alignment_name}_site_analysis"
                decay_calc.write_site_analysis_results(site_output_dir)
            
            decay_calc.cleanup_intermediate_files()
            
            # Return summary statistics for batch report
            return {
                'file': alignment_file.name,
                'status': 'success',
                'num_branches': len(decay_calc.decay_indices) if hasattr(decay_calc, 'decay_indices') else 0,
                'ml_likelihood': decay_calc.ml_likelihood,
                'output_file': output_file.name
            }
        else:
            return {
                'file': alignment_file.name,
                'status': 'failed',
                'error': 'ML tree construction failed'
            }
            
    except Exception as e:
        logger.error(f"Error processing {alignment_file.name}: {e}")
        return {
            'file': alignment_file.name,
            'status': 'failed',
            'error': str(e)
        }


def create_decay_calc_from_args(args):
    """Create panDecayIndices instance from argument namespace."""
    # Convert effective model string
    effective_model_str = args.model
    if args.gamma: effective_model_str += "+G"
    if args.invariable: effective_model_str += "+I"
    
    # Handle PAUP block if specified
    paup_block_content = None
    if args.paup_block:
        paup_block_content = panDecayIndices.read_paup_block(Path(args.paup_block))
    
    # Convert paths
    temp_dir_path = Path(args.temp) if args.temp else None
    starting_tree_path = Path(args.starting_tree) if args.starting_tree else None
    
    return panDecayIndices(
        alignment_file=args.alignment,
        alignment_format=args.format,
        model=effective_model_str,
        temp_dir=temp_dir_path,
        paup_path=args.paup,
        threads=args.threads,
        starting_tree=starting_tree_path,
        data_type=args.data_type,
        debug=args.debug,
        keep_files=args.keep_files,
        gamma_shape=args.gamma_shape, prop_invar=args.prop_invar,
        base_freq=args.base_freq, rates=args.rates,
        protein_model=args.protein_model, nst=args.nst,
        parsmodel=args.parsmodel,
        paup_block=paup_block_content,
        analysis_mode=args.analysis,
        bayesian_software=args.bayesian_software,
        mrbayes_path=args.mrbayes_path,
        bayes_model=args.bayes_model,
        bayes_ngen=args.bayes_ngen,
        bayes_burnin=args.bayes_burnin,
        bayes_chains=args.bayes_chains,
        bayes_sample_freq=args.bayes_sample_freq,
        marginal_likelihood=args.marginal_likelihood,
        ss_alpha=args.ss_alpha,
        ss_nsteps=args.ss_nsteps,
        use_mpi=args.use_mpi,
        mpi_processors=args.mpi_processors,
        mpirun_path=args.mpirun_path,
        use_beagle=args.use_beagle,
        beagle_device=args.beagle_device,
        beagle_precision=args.beagle_precision,
        beagle_scaling=args.beagle_scaling,
        constraint_mode=args.constraint_mode,
        test_branches=args.test_branches,
        constraint_file=args.constraint_file,
        config_constraints=getattr(args, 'config_constraints', None),
        check_convergence=args.check_convergence,
        min_ess=args.min_ess,
        max_psrf=args.max_psrf,
        max_asdsf=args.max_asdsf,
        convergence_strict=args.convergence_strict,
        mrbayes_parse_timeout=args.mrbayes_parse_timeout,
        output_style=args.output_style,
        normalize_bd=args.normalize_bd,
        bd_normalization_methods=args.bd_normalization_methods
    )


def run_batch_analysis(args):
    """Run batch analysis on multiple alignment files."""
    logger.info("Starting batch processing mode...")
    
    # Validate batch arguments
    if not args.input_dir:
        logger.error("--input-dir is required for batch processing")
        sys.exit(1)
    
    # Find alignment files
    try:
        alignment_files = find_alignment_files(args.input_dir, args.file_pattern)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
    
    # Set up output directory
    if args.batch_output_dir:
        output_dir = Path(args.batch_output_dir)
    else:
        output_dir = Path.cwd()
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine number of parallel jobs
    max_jobs = min(args.batch_jobs, multiprocessing.cpu_count(), len(alignment_files))
    logger.info(f"Running batch analysis with {max_jobs} parallel jobs")
    
    # Prepare arguments for each file
    job_args = [(f, args, output_dir) for f in alignment_files]
    
    # Run parallel processing
    start_time = time.time()
    
    if max_jobs == 1:
        # Sequential processing
        results = []
        for job_arg in job_args:
            result = process_single_alignment(job_arg)
            results.append(result)
    else:
        # Parallel processing
        with multiprocessing.Pool(processes=max_jobs) as pool:
            results = pool.map(process_single_alignment, job_args)
    
    end_time = time.time()
    
    # Generate batch summary
    write_batch_summary(results, output_dir / args.batch_summary, end_time - start_time)
    
    # Report results
    successful = sum(1 for r in results if r['status'] == 'success')
    failed = len(results) - successful
    
    logger.info(f"Batch processing completed:")
    logger.info(f"  Successful: {successful}/{len(results)}")
    logger.info(f"  Failed: {failed}/{len(results)}")
    logger.info(f"  Total time: {end_time - start_time:.1f} seconds")
    logger.info(f"  Summary: {output_dir / args.batch_summary}")


def write_batch_summary(results, summary_file, total_time):
    """Write batch processing summary to file."""
    with open(summary_file, 'w') as f:
        f.write("# panDecay Batch Processing Summary\n")
        f.write(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Total processing time: {total_time:.1f} seconds\n\n")
        
        # Summary statistics
        successful = [r for r in results if r['status'] == 'success']
        failed = [r for r in results if r['status'] == 'failed']
        
        f.write(f"Files processed: {len(results)}\n")
        f.write(f"Successful: {len(successful)}\n")
        f.write(f"Failed: {len(failed)}\n\n")
        
        # Detailed results table
        f.write("File\tStatus\tBranches\tML_Likelihood\tOutput_File\tError\n")
        
        for result in results:
            file_name = result['file']
            status = result['status']
            branches = result.get('num_branches', 'N/A')
            likelihood = result.get('ml_likelihood', 'N/A')
            output_file = result.get('output_file', 'N/A')
            error = result.get('error', '')
            
            f.write(f"{file_name}\t{status}\t{branches}\t{likelihood}\t{output_file}\t{error}\n")


def main():
    parser = argparse.ArgumentParser(
        description=f"panDecay v{VERSION}: Calculate phylogenetic decay indices (ML, Bayesian, and parsimony).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter # Shows defaults in help
    )
    # Arguments (similar to original, ensure help messages are clear)
    parser.add_argument("alignment", nargs='?', help="Input alignment file path (can be specified in config file).")
    parser.add_argument("--format", default="fasta", help="Alignment format.")
    parser.add_argument("--model", default="GTR", help="Base substitution model (e.g., GTR, HKY, JC). Combine with --gamma and --invariable.")
    parser.add_argument("--gamma", action="store_true", help="Add Gamma rate heterogeneity (+G) to model.")
    parser.add_argument("--invariable", action="store_true", help="Add invariable sites (+I) to model.")

    parser.add_argument("--paup", default="paup", help="Path to PAUP* executable.")
    parser.add_argument("--output", default="pan_decay_indices.txt", help="Output file for summary results.")
    parser.add_argument("--tree", default="annotated_tree", help="Base name for annotated tree files. Three trees will be generated with suffixes: _au.nwk (AU p-values), _delta_lnl.nwk (likelihood differences), and _combined.nwk (both values).")
    parser.add_argument("--site-analysis", action="store_true", help="Perform site-specific likelihood analysis to identify supporting/conflicting sites for each branch.")
    parser.add_argument("--data-type", default="dna", choices=["dna", "protein", "discrete"], help="Type of sequence data.")
    
    # Batch processing options
    batch_opts = parser.add_argument_group('Batch Processing Options')
    batch_opts.add_argument("--batch", action="store_true", help="Enable batch processing mode for multiple alignment files.")
    batch_opts.add_argument("--input-dir", help="Directory containing alignment files for batch processing (use with --batch).")
    batch_opts.add_argument("--file-pattern", default="*", help="File pattern for batch processing (e.g., '*.fasta', '*.nex'). Default: '*' (all files)")
    batch_opts.add_argument("--batch-output-dir", help="Output directory for batch results (default: current directory).")
    batch_opts.add_argument("--batch-jobs", type=int, default=1, help="Number of parallel batch jobs (default: 1, max: number of CPU cores).")
    batch_opts.add_argument("--batch-summary", default="batch_summary.txt", help="Summary file for batch analysis results.")
    
    # Model parameter overrides
    mparams = parser.add_argument_group('Model Parameter Overrides (optional)')
    mparams.add_argument("--gamma-shape", type=float, help="Fixed Gamma shape value (default: estimate if +G).")
    mparams.add_argument("--prop-invar", type=float, help="Fixed proportion of invariable sites (default: estimate if +I).")
    mparams.add_argument("--base-freq", choices=["equal", "estimate", "empirical"], help="Base/state frequencies (default: model-dependent). 'empirical' uses observed frequencies.")
    mparams.add_argument("--rates", choices=["equal", "gamma"], help="Site rate variation model (overrides --gamma flag if specified).")
    mparams.add_argument("--protein-model", help="Specific protein model (e.g., JTT, WAG; overrides base --model for protein data).")
    mparams.add_argument("--nst", type=int, choices=[1, 2, 6], help="Number of substitution types (DNA; overrides model-based nst).")
    mparams.add_argument("--parsmodel", action=argparse.BooleanOptionalAction, default=None, help="Use parsimony-based branch lengths (discrete data; default: yes for discrete). Use --no-parsmodel to disable.")

    run_ctrl = parser.add_argument_group('Runtime Control')
    run_ctrl.add_argument("--threads", default="auto", help="Number of threads for PAUP* (e.g., 4 or 'auto').")
    run_ctrl.add_argument("--starting-tree", help="Path to a user-provided starting tree file (Newick).")
    run_ctrl.add_argument("--paup-block", help="Path to file with custom PAUP* commands for model/search setup (overrides most model args).")
    run_ctrl.add_argument("--temp", help="Custom directory for temporary files (default: system temp).")
    run_ctrl.add_argument("--keep-files", action="store_true", help="Keep temporary files after analysis.")
    run_ctrl.add_argument("--debug", action="store_true", help="Enable detailed debug logging (implies --keep-files).")

    # Analysis mode selection
    analysis_mode = parser.add_argument_group('Analysis Mode')
    analysis_mode.add_argument("--analysis", 
                              choices=["ml", "bayesian", "parsimony", "ml+parsimony", "bayesian+parsimony", "ml+bayesian", "all"], 
                              default="ml",
                              help="Type of decay analysis to perform (default: ml). "
                                   "Options: ml, bayesian, parsimony, ml+parsimony, bayesian+parsimony, ml+bayesian, all")
    
    # Add bootstrap options
    bootstrap_opts = parser.add_argument_group('Bootstrap Analysis (optional)')
    bootstrap_opts.add_argument("--bootstrap", action="store_true", help="Perform bootstrap analysis to calculate support values.")
    bootstrap_opts.add_argument("--bootstrap-reps", type=int, default=100, help="Number of bootstrap replicates (default: 100)")
    
    # Bayesian-specific options
    bayesian_opts = parser.add_argument_group('Bayesian Analysis Options')
    bayesian_opts.add_argument("--bayesian-software", choices=["mrbayes"], 
                              default="mrbayes", help="Bayesian software to use (default: mrbayes)")
    bayesian_opts.add_argument("--mrbayes-path", default="mb", help="Path to MrBayes executable")
    bayesian_opts.add_argument("--bayes-model", help="Model for Bayesian analysis (if different from ML model)")
    bayesian_opts.add_argument("--bayes-ngen", type=int, default=1000000, help="Number of MCMC generations")
    bayesian_opts.add_argument("--bayes-burnin", type=float, default=0.25, help="Burnin fraction (0-1)")
    bayesian_opts.add_argument("--bayes-chains", type=int, default=4, help="Number of MCMC chains")
    bayesian_opts.add_argument("--bayes-sample-freq", type=int, default=1000, help="Sample frequency for MCMC")
    bayesian_opts.add_argument("--marginal-likelihood", choices=["ss", "ps", "hm"], default="ss",
                              help="Marginal likelihood estimation method: ss=stepping-stone, ps=path sampling, hm=harmonic mean")
    bayesian_opts.add_argument("--ss-alpha", type=float, default=0.4, help="Alpha parameter for stepping-stone sampling")
    bayesian_opts.add_argument("--ss-nsteps", type=int, default=50, help="Number of steps for stepping-stone sampling")
    
    # Normalization options
    normalization_opts = parser.add_argument_group('Bayesian Decay Normalization Options')
    normalization_opts.add_argument("--normalize-bd", action=argparse.BooleanOptionalAction, default=True,
                                   help="Calculate normalized BD metrics for cross-study comparisons (default: enabled)")
    normalization_opts.add_argument("--bd-normalization-methods", nargs="+", 
                                   choices=["per_site", "relative", "signal_to_noise", 
                                           "effect_size", "effect_size_robust", "effect_size_weighted"], 
                                   default=["per_site", "relative"],
                                   help="Which normalization methods to calculate: per_site (BD/sites), relative (BD%%), signal_to_noise (supporting/conflicting ratio), effect_size (BD/SD like Cohen's d), effect_size_robust (BD/MAD), effect_size_weighted (weighted by site information) (default: per_site relative)")
    
    # MPI and BEAGLE options
    parallel_opts = parser.add_argument_group('Parallel Processing Options (MrBayes)')
    parallel_opts.add_argument("--use-mpi", action="store_true", help="Use MPI version of MrBayes")
    parallel_opts.add_argument("--mpi-processors", type=int, help="Number of MPI processors (default: number of chains)")
    parallel_opts.add_argument("--mpirun-path", default="mpirun", help="Path to mpirun executable")
    parallel_opts.add_argument("--use-beagle", action="store_true", help="Enable BEAGLE library for GPU/CPU acceleration")
    parallel_opts.add_argument("--beagle-device", choices=["cpu", "gpu", "auto"], default="auto", 
                             help="BEAGLE device preference")
    parallel_opts.add_argument("--beagle-precision", choices=["single", "double"], default="double",
                             help="BEAGLE precision mode")
    parallel_opts.add_argument("--beagle-scaling", choices=["none", "dynamic", "always"], default="dynamic",
                             help="BEAGLE scaling frequency")
    
    # Convergence checking options
    convergence_opts = parser.add_argument_group('Convergence Checking Options (MrBayes)')
    convergence_opts.add_argument("--check-convergence", action=argparse.BooleanOptionalAction, default=True,
                                 help="Check MCMC convergence diagnostics")
    convergence_opts.add_argument("--min-ess", type=int, default=200,
                                 help="Minimum ESS (Effective Sample Size) threshold (default: 200)")
    convergence_opts.add_argument("--max-psrf", type=float, default=1.01,
                                 help="Maximum PSRF (Potential Scale Reduction Factor) threshold (default: 1.01)")
    convergence_opts.add_argument("--max-asdsf", type=float, default=0.01,
                                 help="Maximum ASDSF (Average Standard Deviation of Split Frequencies) threshold (default: 0.01)")
    convergence_opts.add_argument("--convergence-strict", action="store_true",
                                 help="Fail analysis if convergence criteria not met (default: warn only)")
    convergence_opts.add_argument("--mrbayes-parse-timeout", type=float, default=30.0,
                                 help="Timeout for parsing MrBayes consensus trees in seconds (0 to disable, default: 30)")

    viz_opts = parser.add_argument_group('Visualization Output (optional)')
    viz_opts.add_argument("--visualize", action="store_true", help="Generate static visualization plots (requires matplotlib, seaborn).")
    viz_opts.add_argument("--viz-format", default="png", choices=["png", "pdf", "svg"], help="Format for static visualizations.")
    viz_opts.add_argument("--annotation", default="lnl", choices=["au", "lnl"], help="Type of support values to visualize in distribution plots (au=AU p-values, lnl=likelihood differences).")
    viz_opts.add_argument("--output-style", choices=["unicode", "ascii", "minimal"], default="unicode",
                         help="Output formatting style: unicode (modern), ascii (compatible), minimal (basic)")
    
    # Configuration file and constraint selection options
    config_opts = parser.add_argument_group('Configuration and Constraint Options')
    config_opts.add_argument("--config", help="Read parameters from configuration file (INI format)")
    config_opts.add_argument("--generate-config", help="Generate a template configuration file at the specified path and exit")
    config_opts.add_argument("--constraint-mode", choices=["all", "specific", "exclude"], default="all",
                           help="Branch selection mode: all (test all branches), specific (test only specified), exclude (test all except specified)")
    config_opts.add_argument("--test-branches", 
                           help="Specify branches to test. Format: 'taxon1,taxon2,taxon3;taxon4,taxon5' for clades, '1,3,5' for branch IDs, or '@file.txt' to read from file")
    config_opts.add_argument("--constraint-file", help="File containing constraint definitions (one per line)")
    
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    
    # Handle --generate-config first
    if args.generate_config:
        generate_config_template(args.generate_config)
        sys.exit(0)
    
    # Handle --config if provided
    if args.config:
        parse_config(args.config, args)
    
    # Handle batch processing mode
    if args.batch:
        # For batch mode, alignment file is not required as individual arguments
        # but input_dir is required
        if not args.input_dir:
            parser.error("--input-dir is required when using --batch mode")
        
        # Validate batch job count
        if args.batch_jobs > multiprocessing.cpu_count():
            logger.warning(f"Requested {args.batch_jobs} jobs, but only {multiprocessing.cpu_count()} CPU cores available. Using {multiprocessing.cpu_count()} jobs.")
            args.batch_jobs = multiprocessing.cpu_count()
        
        # Run batch processing and exit
        run_batch_analysis(args)
        return
    
    # Validate that alignment is provided (either via command line or config) for single file mode
    if not args.alignment:
        parser.error("Alignment file must be specified either as a command-line argument or in the config file")
    
    # Convert data_type to lowercase to handle case-insensitive input
    args.data_type = args.data_type.lower()
    
    # Validate Bayesian analysis arguments - bayesian_software now has a default
    # No validation needed since default is "mrbayes"
    
    if args.analysis == "ml" and args.bayesian_software:
        logger.warning("Bayesian software specified but analysis mode is ML-only. Bayesian options will be ignored.")

    if args.debug:
        logger.setLevel(logging.DEBUG)
        # Setup a dedicated debug file handler
        debug_log_path = Path.cwd() / "mldecay_debug.log" # Or in temp_path once it's known
        fh = logging.FileHandler(debug_log_path, mode='w') # Overwrite for each run
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.info(f"Debug logging enabled. Detailed log: {debug_log_path}")
        args.keep_files = True # Debug implies keeping files

    # Construct full model string for display and internal use if not using paup_block
    effective_model_str = args.model
    if args.gamma: effective_model_str += "+G"
    if args.invariable: effective_model_str += "+I"
    # Adjust for protein/discrete if base model is not specific enough
    if args.data_type == "protein" and not args.protein_model and not any(pm in args.model.upper() for pm in ["JTT", "WAG", "LG", "DAYHOFF"]):
        logger.info(f"Protein data with generic model '{args.model}'. Effective model might default to JTT within PAUP* settings.")
        # effective_model_str might be JTT+G+I if G/I are on
    elif args.data_type == "discrete" and "MK" not in args.model.upper():
        logger.info(f"Discrete data with non-Mk model '{args.model}'. Effective model might default to Mk within PAUP* settings.")

    paup_block_content = None
    if args.paup_block:
        pbf_path = Path(args.paup_block)
        logger.info(f"Reading PAUP block from: {pbf_path}")
        paup_block_content = panDecayIndices.read_paup_block(pbf_path)
        if paup_block_content is None: # Handles not found or invalid block
            logger.error("Failed to read or validate PAUP block file. Exiting.")
            sys.exit(1)

    print_runtime_parameters(args, effective_model_str)

    try:
        # Use shared function to create decay calculator
        decay_calc = create_decay_calc_from_args(args)

        decay_calc.build_ml_tree() # Can raise exceptions

        if decay_calc.ml_tree and decay_calc.ml_likelihood is not None:
            # Run bootstrap analysis if requested
            if args.bootstrap:
                logger.info(f"Running bootstrap analysis with {args.bootstrap_reps} replicates...")
                decay_calc.run_bootstrap_analysis(num_replicates=args.bootstrap_reps)

            decay_calc.calculate_decay_indices(perform_site_analysis=args.site_analysis)

            # Add this new code snippet here
            if hasattr(decay_calc, 'decay_indices') and decay_calc.decay_indices:
                for clade_id, data in decay_calc.decay_indices.items():
                    if 'site_data' in data:
                        site_output_dir = Path(args.output).parent / f"{Path(args.output).stem}_site_analysis"
                        decay_calc.write_site_analysis_results(site_output_dir)
                        logger.info(f"Site-specific analysis results written to {get_display_path(site_output_dir)}")
                        break  # Only need to do this once if any site_data exists

            output_main_path = Path(args.output)
            decay_calc.write_formatted_results(output_main_path)

            report_path = output_main_path.with_suffix(".md")
            decay_calc.generate_detailed_report(report_path)

            output_dir = output_main_path.resolve().parent
            tree_base_name = args.tree  # Use the tree argument directly as the base name
            tree_files = decay_calc.annotate_trees(output_dir, tree_base_name)

            if tree_files:
                logger.info(f"Successfully created {len(tree_files)} annotated trees.")
                for tree_type, path in tree_files.items():
                    logger.info(f"  - {tree_type} tree: {get_display_path(path)}")
            else:
                logger.warning("Failed to create annotated trees.")

            if args.visualize:
                viz_out_dir = output_main_path.resolve().parent # Ensure absolute path for parent
                viz_base_name = output_main_path.stem
                viz_kwargs = {'width': 10, 'height': 6, 'format': args.viz_format}

                # Check for viz library availability early
                try: import matplotlib, seaborn
                except ImportError:
                    logger.warning("Matplotlib/Seaborn not installed. Skipping static visualizations.")
                    args.visualize = False # Disable further attempts

                if args.visualize:
                    decay_calc.visualize_support_distribution(
                        viz_out_dir / f"{viz_base_name}_dist_{args.annotation}.{args.viz_format}",
                        value_type=args.annotation, **viz_kwargs)

                if args.site_analysis:
                    # Pass visualization preferences to the panDecayIndices instance
                    if args.visualize:
                        decay_calc.viz_format = args.viz_format

                    site_output_dir = output_main_path.parent / f"{output_main_path.stem}_site_analysis"
                    decay_calc.write_site_analysis_results(site_output_dir)
                    logger.info(f"Site-specific analysis results written to {get_display_path(site_output_dir)}")

            decay_calc.cleanup_intermediate_files()
            logger.info("panDecay analysis completed successfully.")
        else:
            logger.error("ML tree construction failed or likelihood missing. Halting.")
            sys.exit(1) # Ensure exit if ML tree is critical and failed

    except Exception as e:
        logger.error(f"panDecay analysis terminated with an error: {e}")
        if args.debug: # Print traceback in debug mode
            import traceback
            logger.debug("Full traceback:\n%s", traceback.format_exc())
        sys.exit(1)

    finally:
        # This block executes whether try succeeds or fails.
        # If decay_calc was initialized and keep_files is false, __del__ will handle cleanup.
        # If __init__ failed before self.temp_path was set, no specific cleanup here yet.
        if 'decay_calc' in locals() and (args.debug or args.keep_files):
            logger.info(f"Temporary files are preserved in: {decay_calc.temp_path}")


if __name__ == "__main__":
    main()
