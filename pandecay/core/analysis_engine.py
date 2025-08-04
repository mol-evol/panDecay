#!/usr/bin/env python3
"""
Core Analysis Engine for panDecay.

Contains the main panDecayIndices class extracted from the original
monolithic system to preserve exact behavior and functionality.
"""

import os
import sys
import json
import csv
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
from pathlib import Path

from pandecay.core.utils import ProgressIndicator, OverwritingProgress


class AnalysisEngineError(Exception):
    """Custom exception for analysis engine errors."""
    pass


class ExternalToolError(AnalysisEngineError):
    """Exception for external tool execution errors."""
    pass
from pandecay.core.constants import (
    VERSION,
    NEXUS_ALIGNMENT_FN,
    ML_TREE_FN,
    ML_SCORE_FN,
    ML_SEARCH_NEX_FN,
    ML_LOG_FN,
    AU_TEST_NEX_FN,
    AU_TEST_SCORE_FN,
    AU_LOG_FN,
    BOOTSTRAP_TREES_FN,
    BOOTSTRAP_SEARCH_NEX_FN,
    BOOTSTRAP_LOG_FN,
    PARSIMONY_TREES_FN,
    PARSIMONY_SCORE_FN,
    PARSIMONY_SEARCH_NEX_FN,
    PARSIMONY_LOG_FN,
    CONSTRAINT_TREE_PREFIX,
    CONSTRAINT_SCORE_PREFIX,
    CONSTRAINT_SEARCH_PREFIX,
    CONSTRAINT_LOG_PREFIX,
    SITE_ANALYSIS_PREFIX,
    SITE_LNL_PREFIX,
    MRBAYES_LOG_FN,
    ANNOTATED_TREE_AU,
    ANNOTATED_TREE_DELTA_LNL,
    ANNOTATED_TREE_COMBINED,
    ANNOTATED_TREE_BOOTSTRAP,
    ANNOTATED_TREE_COMPREHENSIVE
)

logger = logging.getLogger(__name__)


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
                 mrbayes_parse_timeout=30.0, output_style="unicode"):

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
            logger.info(f"Loaded alignment: {len(self.alignment)} sequences, {self.alignment.get_alignment_length()} sites.")
        except FileNotFoundError:
            msg = f"Alignment file not found: '{self.alignment_file}'"
            logger.error(msg)
            if self._temp_dir_obj: self._temp_dir_obj.cleanup() # Manual cleanup if init fails early
            raise AnalysisEngineError(msg)
        except Exception as e:
            msg = f"Failed to load alignment '{self.alignment_file}': {e}"
            logger.error(msg)
            if self._temp_dir_obj: self._temp_dir_obj.cleanup() # Manual cleanup if init fails early
            raise AnalysisEngineError(msg)

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
        progress = ProgressIndicator()
        progress.start("Building maximum likelihood tree...")
        
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
                patterns = [r'-ln\s*L\s*=\s*([0-9.]+)', r'likelihood\s*=\s*([0-9.]+)', r'score\s*=\s*([0-9.]+)']
                for p in patterns:
                    m = re.findall(p, paup_result.stdout, re.IGNORECASE)
                    if m: self.ml_likelihood = float(m[-1]); break
                if not self.ml_likelihood: logger.warning("Could not extract ML likelihood from PAUP* log.")

            ml_tree_path = self.temp_path / ML_TREE_FN
            if ml_tree_path.exists() and ml_tree_path.stat().st_size > 0:
                # Clean the tree file if it has metadata after semicolon
                cleaned_tree_path = self._clean_newick_tree(ml_tree_path)
                self.ml_tree = Phylo.read(str(cleaned_tree_path), "newick")
                
                # Consolidated final message
                lk_display = f"{self.ml_likelihood:.3f}" if self.ml_likelihood is not None else "N/A"
                progress.stop(f"ML tree built → Log-likelihood: {lk_display}")
                
                if self.ml_likelihood is None:
                    logger.error("ML tree built, but likelihood could not be determined. Analysis may be compromised.")
            else:
                progress.stop()
                logger.error(f"ML tree file {ml_tree_path} not found or is empty after PAUP* run.")
                raise FileNotFoundError(f"ML tree file missing or empty: {ml_tree_path}")
        except Exception as e:
            progress.stop()
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
            
        progress = ProgressIndicator()
        progress.start(f"Running Bayesian analysis ({self.bayes_ngen:,} gen, {self.bayes_chains} chains)...")
        
        # Log parallel processing settings
        if self.use_mpi:
            n_procs = self.mpi_processors if self.mpi_processors else self.bayes_chains
            logger.info(f"MPI enabled: Using {n_procs} processors with {self.mpirun_path}")
        if self.use_beagle:
            logger.info(f"BEAGLE enabled: device={self.beagle_device}, precision={self.beagle_precision}, scaling={self.beagle_scaling}")
        
        # First, run unconstrained Bayesian analysis
        progress.update("Running unconstrained analysis...")
        
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
            progress.stop()
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
        
        try:
            for branch_num, (i, clade_obj, clade_log_idx, clade_taxa) in enumerate(testable_branches, 1):
                clade_id = f"Clade_{clade_log_idx}"
                
                # Update progress with current branch
                progress.update(f"Testing branch {branch_num}/{len(testable_branches)} ({clade_id})...")
                
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
                    
                    bayesian_results[clade_id] = {
                        'unconstrained_ml': unconstrained_ml,
                        'constrained_ml': constrained_ml,
                        'bayes_decay': bayes_decay,
                        'taxa': clade_taxa
                    }
                    
                    # Results will be displayed in final summary
                    if bayes_decay < 0:
                        logger.warning(f"WARNING: {clade_id} has negative Bayes Decay ({bayes_decay:.4f}), suggesting potential convergence or estimation issues")
                else:
                    logger.warning(f"Constrained analysis failed for {clade_id}")
        
            progress.stop(f"Bayesian analysis completed → {len(bayesian_results)} branches analyzed")
        except Exception as e:
            progress.stop()
            raise
        
        # Display consolidated results for each branch with overwriting progress
        if bayesian_results:
            overwrite_progress = OverwritingProgress()
            try:
                branch_list = list(bayesian_results.items())
                for i, (clade_id, data) in enumerate(branch_list, 1):
                    bayes_decay = data['bayes_decay']
                    overwrite_progress.update(f"Processing results: {clade_id} ({i}/{len(branch_list)})")
                overwrite_progress.finish(f"Bayesian results: {len(branch_list)} branches completed")
            except Exception:
                overwrite_progress.finish(f"Bayesian results: {len(bayesian_results)} branches completed")
            
            # Check for negative Bayes Decay values and issue summary warning
            negative_clades = [cid for cid, data in bayesian_results.items() if data['bayes_decay'] < 0]
            if negative_clades:
                logger.warning(f"\nWARNING: {len(negative_clades)}/{len(bayesian_results)} clades have negative Bayes Decay values!")
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
            order = ['AU', 'ΔlnL', 'BD', 'PD', 'PP', 'BS']
            for key in order:
                if key in annotation_dict and annotation_dict[key] is not None:
                    val = annotation_dict[key]
                    if isinstance(val, float):
                        if key in ['AU', 'PP']:
                            formatted_values.append(f"{key}={val:.3f}")
                        elif key in ['ΔlnL', 'BD']:
                            formatted_values.append(f"{key}={val:.2f}")
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
                'AU': '✓', 'ΔlnL': '△', 'BD': '◆', 
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
            logger.warning(f"  WARNING: {warning}")
        
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
        progress = ProgressIndicator()
        progress.start("Building parsimony tree...")
        
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
            
            progress.update(f"Testing branches for parsimony decay ({len(testable_branches)} branches)...")
            
            # Process testable branches
            try:
                results_to_display = []
                for branch_num, (i, clade_obj, clade_log_idx, clade_taxa) in enumerate(testable_branches, 1):
                    clade_id = f"Clade_{clade_log_idx}"
                    
                    # Update progress with current branch
                    progress.update(f"Testing branch {branch_num}/{len(testable_branches)} ({clade_id})...")
                    
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
                            # Store result for final display
                            results_to_display.append((clade_id, decay_value, pars_score, constrained_score))
                        
                    except Exception as e:
                        logger.error(f"Failed to calculate parsimony decay for {clade_id}: {e}")
                        continue
                    
                progress.stop(f"Parsimony analysis completed → {len(parsimony_decay)} branches analyzed")
                
                # Display consolidated results with overwriting progress
                if results_to_display:
                    overwrite_progress = OverwritingProgress()
                    try:
                        for i, (clade_id, decay_value, pars_score, constrained_score) in enumerate(results_to_display, 1):
                            overwrite_progress.update(f"Processing results: {clade_id} ({i}/{len(results_to_display)})")
                        overwrite_progress.finish(f"Parsimony results: {len(results_to_display)} branches completed")
                    except Exception:
                        overwrite_progress.finish(f"Parsimony results: {len(results_to_display)} branches completed")
                    
            except Exception as inner_e:
                progress.stop()
                raise inner_e
                
        except Exception as e:
            progress.stop()
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
                            'bayes_decay': bayes_data.get('bayes_decay')
                        })
                    else:
                        # Create new entry for Bayesian-only results
                        self.decay_indices[clade_id] = bayes_data
        
        # Calculate parsimony decay indices if requested
        if self.do_parsimony:
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
        
        progress_ml = ProgressIndicator()
        progress_ml.start(f"Testing {len(testable_branches)} branches for ML decay...")
        
        try:
            for branch_num, (i, clade_obj, clade_log_idx, clade_taxa_names) in enumerate(testable_branches, 1):
                # Update progress with current branch
                progress_ml.update(f"Testing branch {branch_num}/{len(testable_branches)} (Clade_{clade_log_idx})...")
                
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

            progress_ml.stop(f"ML constraint analysis completed → {len(constraint_info_map)} branches analyzed")
        except Exception as e:
            progress_ml.stop()
            raise
            
        if not constraint_info_map:
            logger.warning("No valid constraint trees were generated. Skipping AU test.")
            self.decay_indices = {}
            return self.decay_indices

        # Perform site-specific likelihood analysis if requested
        if perform_site_analysis:
            progress = ProgressIndicator()
            progress.start("Analyzing site-specific likelihoods...")
            
            try:
                for clade_id, cdata in list(constraint_info_map.items()):
                    rel_constr_tree_fn = cdata.get('tree_filename')

                    if rel_constr_tree_fn:
                        tree_files = [ML_TREE_FN, rel_constr_tree_fn]
                        site_analysis_result = self._calculate_site_likelihoods(tree_files, clade_id)

                        if site_analysis_result:
                            # Store all site analysis data
                            constraint_info_map[clade_id].update(site_analysis_result)
                            
                progress.stop("Site analysis completed")
                
                # Display site analysis results
                # Use actual alignment length as total sites (constant for all branches)
                total_sites = self.alignment.get_alignment_length()
                
                for clade_id, cdata in constraint_info_map.items():
                    if 'supporting_sites' in cdata:  # Has site analysis data
                        supporting = cdata.get('supporting_sites', 0)
                        conflicting = cdata.get('conflicting_sites', 0)
                        neutral = cdata.get('neutral_sites', 0)
                        
                        # Verify site counting adds up (debugging info)
                        calculated_total = supporting + conflicting + neutral
                        if calculated_total != total_sites:
                            logger.debug(f"Site count mismatch for {clade_id}: {calculated_total} calculated vs {total_sites} alignment length")
                        
                        ratio = cdata.get('support_ratio', 0.0)
                        weighted = cdata.get('weighted_support_ratio', 0.0)
                        sum_supporting = cdata.get('sum_supporting_delta', 0.0)
                        sum_conflicting = cdata.get('sum_conflicting_delta', 0.0)
                        
                        supporting_pct = (supporting / total_sites * 100) if total_sites > 0 else 0
                        conflicting_pct = (conflicting / total_sites * 100) if total_sites > 0 else 0
                        
                        logger.info(f"Branch {clade_id}: {total_sites} sites → {supporting} supporting ({supporting_pct:.1f}%), {conflicting} conflicting ({conflicting_pct:.1f}%) | Support ratio: {ratio:.2f} | Weighted: {weighted:.2f} (Δ support: {sum_supporting:.2f}, Δ conflict: {sum_conflicting:.2f})")
            except Exception as e:
                progress.stop()
                raise

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
                        'bayes_decay': bayes_data['bayes_decay']
                    })
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

                supporting_sites = sum(1 for d in deltas if d < 0)
                conflicting_sites = sum(1 for d in deltas if d > 0)
                # Calculate neutral sites as remainder to ensure all sites are accounted for
                total_alignment_sites = self.alignment.get_alignment_length()
                sites_analyzed = len(deltas)
                neutral_sites = sites_analyzed - supporting_sites - conflicting_sites
                
                # If we have fewer analyzed sites than alignment length, note the difference
                if sites_analyzed != total_alignment_sites:
                    logger.debug(f"Branch {branch_id}: Analyzed {sites_analyzed} sites, alignment has {total_alignment_sites} sites")

                # Calculate sum of likelihood differences
                sum_supporting_delta = sum(d for d in deltas if d < 0)  # Sum of negative deltas (supporting)
                sum_conflicting_delta = sum(d for d in deltas if d > 0)  # Sum of positive deltas (conflicting)

                # Calculate weighted support ratio
                weighted_support_ratio = abs(sum_supporting_delta) / sum_conflicting_delta if sum_conflicting_delta > 0 else float('inf')

                # Calculate standard support ratio
                support_ratio = supporting_sites / conflicting_sites if conflicting_sites > 0 else float('inf')

                logger.info(f"Extracted site likelihoods for {len(site_data)} sites for branch {branch_id}")
                # Verbose logging removed - consolidated output handled by caller

                # Return a comprehensive dictionary with all info
                return {
                    'site_data': site_data,
                    'supporting_sites': supporting_sites,
                    'conflicting_sites': conflicting_sites,
                    'neutral_sites': neutral_sites,
                    'support_ratio': support_ratio,
                    'sum_supporting_delta': sum_supporting_delta,
                    'sum_conflicting_delta': sum_conflicting_delta,
                    'weighted_support_ratio': weighted_support_ratio
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
                logger.info(f"AU tree saved: {au_tree_path.name} ({annotated_nodes_count} branches)")
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
                logger.info(f"Delta-LnL tree saved: {lnl_tree_path.name} ({annotated_nodes_count} branches)")
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
                                
                            # Add parsimony decay if available
                            pars_decay = decay_info.get('pars_decay')
                            if pars_decay is not None:
                                annotation_parts.append(f"PD:{pars_decay}")
                                
                            # Add posterior probability if available
                            post_prob = decay_info.get('posterior_prob')
                            if post_prob is not None:
                                annotation_parts.append(f"PP:{post_prob:.2f}")
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

                logger.info(f"Combined tree saved: {combined_tree_path.name} ({annotated_nodes_count} branches)")
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
                    logger.info(f"Bootstrap tree saved: {bootstrap_tree_path.name}")
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
                                
                            # Add parsimony decay if available
                            pars_decay = decay_info.get('pars_decay')
                            if pars_decay is not None:
                                annotation_parts.append(f"PD:{pars_decay}")
                                
                            # Add posterior probability if available
                            post_prob = decay_info.get('posterior_prob')
                            if post_prob is not None:
                                annotation_parts.append(f"PP:{post_prob:.2f}")

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
                    logger.info(f"Comprehensive tree saved: {comprehensive_tree_path.name} ({annotated_nodes_count} branches)")
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
                            node.name = f"{matched_clade_id} - BD:{bayes_decay_val:.4f}"
                            annotated_nodes_count += 1

                    Phylo.write(bd_tree, str(bayes_decay_tree_path), "newick")
                    logger.info(f"Bayes decay tree saved: {bayes_decay_tree_path.name} ({annotated_nodes_count} branches)")
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
        
        # Build header structure
        if self.output_style == "unicode":
            # Top border
            f.write("┌──────────┬──────┬")
            if has_ml:
                f.write("────────────────────────────────┬")
            if has_bayesian:
                f.write("───────────────────────┬")
            if has_parsimony:
                f.write("───────────────────────┬")
            f.write("\n")
            
            # Main headers
            f.write("│ Clade ID │ Taxa │")
            if has_ml:
                f.write("    Maximum Likelihood          │")
            if has_bayesian:
                f.write("       Bayesian        │")
            if has_parsimony:
                f.write("      Parsimony        │")
            f.write("\n")
            
            # Sub-headers
            f.write("│          │      ├")
            if has_ml:
                f.write("──────────┬──────────┬──────────┼")
            if has_bayesian:
                f.write("───────────────────────┼")
            if has_parsimony:
                f.write("───────────────────────┤")
            f.write("\n")
            
            # Column names
            f.write("│          │      │")
            if has_ml:
                f.write(" ΔlnL     │ AU p-val │ Support  │")
            if has_bayesian:
                f.write(" BD                     │")
            if has_parsimony:
                f.write(" Decay │ Post.Prob     │" if has_posterior else " Decay                 │")
            f.write("\n")
            
            # Bottom border
            f.write("├──────────┼──────┼")
            if has_ml:
                f.write("──────────┼──────────┼──────────┼")
            if has_bayesian:
                f.write("───────────────────────┼")
            if has_parsimony:
                f.write("───────┼───────────────┤" if has_posterior else "───────────────────────┤")
            f.write("\n")
        else:
            # ASCII version
            header_parts = ["Clade ID", "Taxa"]
            if has_ml:
                header_parts.extend(["ΔlnL", "AU p-val", "Support"])
            if has_bayesian:
                header_parts.extend(["BD"])
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
                f.write("───────────────────────┴")
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
                
                f.write(f"→ {clade_id}\n")
                
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
        
        logger.info(f"Results saved: {output_path.name}")
    
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

        logger.info(f"Results saved: {output_path.name}")

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
        has_parsimony = any(d.get('pars_decay') is not None for d in self.decay_indices.values()) if self.decay_indices else False
        has_bootstrap = hasattr(self, 'bootstrap_tree') and self.bootstrap_tree
        has_site_data = any('site_data' in data for data in self.decay_indices.values()) if self.decay_indices else False
        
        # Debugging information
        logger.debug(f"Markdown report generation - decay_indices entries: {len(self.decay_indices) if self.decay_indices else 0}")
        logger.debug(f"Markdown report flags: has_ml={has_ml}, has_bayesian={has_bayesian}, has_parsimony={has_parsimony}, has_bootstrap={has_bootstrap}, has_site_data={has_site_data}")
        if self.decay_indices:
            for cid, data in list(self.decay_indices.items())[:2]:  # Show first 2 entries for debugging
                logger.debug(f"Sample data for {cid}: keys={list(data.keys())}")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
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
                                f.write(f"\n**WARNING**: {negative_count}/{len(bayes_decays)} branches have negative Bayes Decay values.\n")
                                f.write("This suggests potential issues with:\n")
                                f.write("- MCMC convergence (consider increasing --bayes-ngen)\n")
                                f.write("- Marginal likelihood estimation reliability\n")
                                f.write("- Model specification\n\n")
                            
                            f.write(f"\n**Note**: BD values closely approximating ML log-likelihood differences is expected behavior in phylogenetic topology testing.\n\n")
                            
                        # Note: We no longer report traditional BF thresholds as they don't apply well to phylogenetics
                    
                    # Parsimony-specific statistics
                    if has_parsimony:
                        pars_decays = [d['pars_decay'] for d in self.decay_indices.values() if d.get('pars_decay') is not None]
                        if pars_decays:
                            f.write(f"- Avg Parsimony decay (step difference): {np.mean(pars_decays):.1f}\n")
                            f.write(f"- Min Parsimony decay: {min(pars_decays)}\n")
                            f.write(f"- Max Parsimony decay: {max(pars_decays)}\n")
                            
                            # Count branches with no parsimony support
                            zero_support = sum(1 for pd in pars_decays if pd == 0)
                            if zero_support > 0:
                                f.write(f"- Branches with no parsimony support (decay = 0): {zero_support} / {len(pars_decays)}\n")
                            
                            f.write("\n")
                
                    # Site analysis statistics
                    if has_site_data:
                        site_ratios = [d.get('support_ratio', 0.0) for d in self.decay_indices.values() if 'site_data' in d]
                        weighted_ratios = [d.get('weighted_support_ratio', 0.0) for d in self.decay_indices.values() if 'site_data' in d]
                        total_sites = self.alignment.get_alignment_length() if hasattr(self.alignment, 'get_alignment_length') else 'N/A'
                        
                        if site_ratios:
                            f.write(f"- Total alignment sites analyzed: {total_sites}\n")
                            f.write(f"- Avg site support ratio: {np.mean(site_ratios):.2f}\n")
                            f.write(f"- Min site support ratio: {min(site_ratios):.2f}\n")
                            f.write(f"- Max site support ratio: {max(site_ratios):.2f}\n")
                            
                            if weighted_ratios:
                                f.write(f"- Avg weighted support ratio: {np.mean(weighted_ratios):.2f}\n")
                            
                            # Count clades with strong site support (ratio > 2.0)
                            strong_support = sum(1 for ratio in site_ratios if ratio > 2.0)
                            f.write(f"- Clades with strong site support (ratio > 2.0): {strong_support} / {len(site_ratios)}\n")
                            f.write("\n")
                
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
                            f.write(f"\n**WARNING**: {len(convergence_issues)} runs did not meet convergence criteria:\n")
                            for run_id in convergence_issues[:5]:  # Show first 5
                                f.write(f"  - {run_id}\n")
                            if len(convergence_issues) > 5:
                                f.write(f"  - ... and {len(convergence_issues) - 5} more\n")
                            f.write("\nConsider:\n")
                            f.write("- Increasing MCMC generations (--bayes-ngen)\n")
                            f.write("- Running longer burnin (--bayes-burnin)\n")
                            f.write("- Using more chains (--bayes-chains)\n")
                            f.write("- Checking model specification\n")
                    else:
                        f.write("Convergence diagnostics not available. This may occur when:\n")
                        f.write("- Using harmonic mean marginal likelihood estimation\n")
                        f.write("- Convergence monitoring was disabled\n")
                        f.write("- Analysis completed before diagnostics were calculated\n\n")
                        f.write("For detailed convergence assessment, consider using stepping-stone sampling with `--marginal-likelihood ss`.\n\n")

                f.write("\n## Detailed Branch Support Results\n\n")

                # Ensure we have data to display
                if not self.decay_indices:
                    f.write("No branch support data available.\n\n")
                    logger.warning("No decay_indices data available for detailed report table")
                else:
                    logger.debug(f"Generating table for {len(self.decay_indices)} clades")

                # Build dynamic table header based on available data
                header_parts = ["| Clade ID | Taxa Count "]
                separator_parts = ["|----------|------------ "]
                
                # Always include basic columns, add analysis-specific ones as available  
                if has_ml:
                    header_parts.extend(["| Constrained lnL | ΔlnL (from ML) | AU p-value | Significant (AU) "])
                    separator_parts.extend(["|-----------------|------------------|------------|-------------------- "])
                    logger.debug("Added ML columns to report table")
                
                if has_bayesian:
                    header_parts.extend(["| Bayes Decay "])
                    separator_parts.extend(["|------------- "])
                    logger.debug("Added Bayesian columns to report table")
                    
                if has_parsimony:
                    header_parts.extend(["| Parsimony Decay "])
                    separator_parts.extend(["|----------------- "])
                    logger.debug("Added Parsimony columns to report table")
                    
                if has_bootstrap:
                    header_parts.append("| Bootstrap ")
                    separator_parts.append("|----------- ")
                    logger.debug("Added Bootstrap columns to report table")
                    
                header_parts.append("| Included Taxa (sample) |\n")
                separator_parts.append("|--------------------------|\n")
                
                # Log final table structure
                logger.debug(f"Table header: {''.join(header_parts).strip()}")
                
                f.write("".join(header_parts))
                f.write("".join(separator_parts))

                # Get bootstrap values if bootstrap analysis was performed
                bootstrap_values = {}
                if has_bootstrap:
                    for node in self.bootstrap_tree.get_nonterminals():
                        if node.confidence is not None:
                            taxa_set = frozenset(leaf.name for leaf in node.get_terminals())
                            bootstrap_values[taxa_set] = node.confidence

                logger.debug(f"Starting table generation loop for {len(self.decay_indices)} clades")
                row_count = 0
                
                try:
                    for clade_id, data in sorted(self.decay_indices.items()):
                        row_count += 1
                        logger.debug(f"Processing row {row_count}: {clade_id}")
                        
                        taxa_list = sorted(data.get('taxa', []))
                        taxa_count = len(taxa_list)
                        taxa_sample = ", ".join(taxa_list[:3]) + ('...' if taxa_count > 3 else '')

                        # Build the table row
                        row_parts = [f"| {clade_id} | {taxa_count} "]
                        
                        # ML fields
                        if has_ml:
                            try:
                                c_lnl = data.get('constrained_lnl', 'N/A')
                                if isinstance(c_lnl, float): c_lnl = f"{c_lnl:.4f}"
                                
                                lnl_d = data.get('lnl_diff', 'N/A')
                                if isinstance(lnl_d, float): lnl_d = f"{lnl_d:.4f}"
                                
                                au_p = data.get('AU_pvalue', 'N/A')
                                if isinstance(au_p, float): au_p = f"{au_p:.4f}"
                                
                                sig_au = data.get('significant_AU', 'N/A')
                                if isinstance(sig_au, bool): sig_au = "**Yes**" if sig_au else "No"
                                
                                row_parts.append(f"| {c_lnl} | {lnl_d} | {au_p} | {sig_au} ")
                            except Exception as e:
                                logger.error(f"Error processing ML data for {clade_id}: {e}")
                                row_parts.append("| N/A | N/A | N/A | N/A ")
                        
                        # Bayesian fields
                        if has_bayesian:
                            try:
                                bayes_d = data.get('bayes_decay', 'N/A')
                                if isinstance(bayes_d, float): bayes_d = f"{bayes_d:.4f}"
                                
                                row_parts.append(f"| {bayes_d} ")
                            except Exception as e:
                                logger.error(f"Error processing Bayesian data for {clade_id}: {e}")
                                row_parts.append("| N/A ")
                        
                        # Parsimony fields
                        if has_parsimony:
                            try:
                                pars_d = data.get('pars_decay', 'N/A')
                                if isinstance(pars_d, (int, float)): pars_d = str(int(pars_d))
                                
                                row_parts.append(f"| {pars_d} ")
                            except Exception as e:
                                logger.error(f"Error processing Parsimony data for {clade_id}: {e}")
                                row_parts.append("| N/A ")

                        # Bootstrap column if available
                        if has_bootstrap:
                            try:
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
                            except Exception as e:
                                logger.error(f"Error processing Bootstrap data for {clade_id}: {e}")
                                row_parts.append("| N/A ")

                        row_parts.append(f"| {taxa_sample} |\n")
                        
                        try:
                            row_content = "".join(row_parts)
                            f.write(row_content)
                            logger.debug(f"Successfully wrote row for {clade_id}")
                        except Exception as e:
                            logger.error(f"Error writing row for {clade_id}: {e}")
                            logger.debug(f"Row parts: {row_parts}")
                            # Write a basic fallback row
                            fallback_row = f"| {clade_id} | {taxa_count} | ERROR | {taxa_sample} |\n"
                            f.write(fallback_row)
                
                except Exception as e:
                    logger.error(f"Major error in table generation loop: {e}")
                    f.write(f"\nError generating detailed results table: {e}\n\n")

                logger.debug(f"Completed table generation for {row_count} rows")
                
                # Add site analysis summary table
                f.write("\n## Site Analysis Summary\n\n")
                
                # Use the has_site_data variable defined at the top of the function
                if has_site_data and self.decay_indices:
                    f.write("Site-specific likelihood analysis showing supporting vs. conflicting sites for each clade:\n\n")
                    
                    # Create site analysis table header
                    f.write("| Clade ID | Supporting Sites | Conflicting Sites | Neutral Sites | Support Ratio | Sum Supporting Δ | Sum Conflicting Δ | Weighted Support Ratio |\n")
                    f.write("|----------|------------------|-------------------|---------------|---------------|------------------|-------------------|------------------------|\n")
                    
                    # Generate table rows for each clade with site data
                    for clade_id in sorted(self.decay_indices.keys()):
                        data = self.decay_indices[clade_id]
                        if 'site_data' not in data:
                            continue
                        
                        supporting = data.get('supporting_sites', 0)
                        conflicting = data.get('conflicting_sites', 0)
                        neutral = data.get('neutral_sites', 0)
                        ratio = data.get('support_ratio', 0.0)
                        sum_supporting = data.get('sum_supporting_delta', 0.0)
                        sum_conflicting = data.get('sum_conflicting_delta', 0.0)
                        weighted_ratio = data.get('weighted_support_ratio', 0.0)
                        
                        f.write(f"| {clade_id} | {supporting} | {conflicting} | {neutral} | {ratio:.4f} | {sum_supporting:.4f} | {sum_conflicting:.4f} | {weighted_ratio:.4f} |\n")
                    
                    f.write("\n**Note**: Site analysis shows which alignment positions support or conflict with each clade. ")
                    f.write("Supporting sites have positive likelihood differences when the clade is present, ")
                    f.write("while conflicting sites favor alternative topologies.\n\n")
                else:
                    f.write("Site-specific analysis data not available.\n\n")
                
                # Add visualizations section
                f.write("\n## Visualizations\n\n")
                
                # Check for main visualization files
                visualization_dir = output_path.parent / "visualizations"
                site_analysis_dir = output_path.parent / "site_analysis"
                
                if visualization_dir.exists():
                    # Main distribution and correlation plots
                    support_dist_plot = visualization_dir / f"{self.alignment_file.stem}_support_distribution.png"
                    support_corr_plot = visualization_dir / f"{self.alignment_file.stem}_support_correlation.png"
                    
                    if support_dist_plot.exists():
                        f.write(f"### Support Distribution Plot\n")
                        f.write(f"![Support Distribution](./visualizations/{support_dist_plot.name})\n\n")
                        f.write("Distribution of support values across all analytical methods.\n\n")
                    
                    if support_corr_plot.exists():
                        f.write(f"### Support Correlation Plot\n")
                        f.write(f"![Support Correlation](./visualizations/{support_corr_plot.name})\n\n")
                        f.write("Correlation between different support measures (ML, Bayesian, Parsimony).\n\n")
                
                # Site-specific plots
                if site_analysis_dir.exists() and has_site_data:
                    f.write("### Site-Specific Analysis Plots\n\n")
                    f.write("Individual clade site analysis visualizations:\n\n")
                    
                    for clade_id in sorted(self.decay_indices.keys()):
                        if 'site_data' in self.decay_indices[clade_id]:
                            site_hist = site_analysis_dir / f"site_hist_{clade_id}.png"
                            site_plot = site_analysis_dir / f"site_plot_{clade_id}.png"
                            
                            if site_hist.exists() or site_plot.exists():
                                f.write(f"#### {clade_id}\n")
                                
                                if site_hist.exists():
                                    f.write(f"![{clade_id} Site Histogram](./site_analysis/{site_hist.name})\n")
                                
                                if site_plot.exists():
                                    f.write(f"![{clade_id} Site Plot](./site_analysis/{site_plot.name})\n")
                                
                                f.write("\n")
                    
                    f.write("Site plots show likelihood differences across alignment positions, with positive values indicating support for the clade.\n\n")
                
                f.write("\n## Interpretation Guide\n\n")
                
                if has_ml:
                    f.write("### ML Analysis\n")
                    f.write("- **ΔlnL (from ML)**: Log-likelihood difference between the constrained tree (without the clade) and the ML tree. Calculated as: constrained_lnL - ML_lnL. Larger positive values indicate stronger support for the clade's presence in the ML tree.\n")
                    f.write("- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.\n\n")
                    
                if has_bayesian:
                    f.write("### Bayesian Analysis\n")
                    f.write("- **Bayes Decay (BD)**: Marginal log-likelihood difference (unconstrained - constrained). This is the primary metric for Bayesian support.\n")
                    f.write("  - **Key insight**: In phylogenetic topology testing, BD values typically closely approximate ML log-likelihood differences\n")
                    f.write("  - **Note**: BD values of 30-50 or higher are common when data support a clade and should not be considered anomalous\n")
                    f.write("  - **Why BD ≈ ΔlnL**: When comparing models that differ only by a topological constraint, the marginal likelihood is dominated by the likelihood component\n")
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
                    
                if has_parsimony:
                    f.write("### Parsimony Analysis\n")
                    f.write("- **Parsimony Decay**: The difference in parsimony steps between the unconstrained most parsimonious tree and the constrained tree (without the clade). Calculated as: constrained_steps - unconstrained_steps.\n")
                    f.write("  - **Interpretation**: Higher positive values indicate that removing the clade requires more evolutionary steps, suggesting stronger parsimony support for the clade.\n")
                    f.write("  - **Zero values**: Indicate that the clade can be removed without increasing the number of parsimony steps, suggesting weak parsimony support.\n")
                    f.write("  - **Negative values**: Should not occur in proper decay analysis (would indicate an error).\n\n")
                    
                if has_bootstrap:
                    f.write("### Bootstrap Analysis\n")
                    f.write("- **Bootstrap**: Bootstrap support value (percentage of bootstrap replicates in which the clade appears). Higher values (e.g., > 70) suggest stronger support for the clade.\n\n")
                
                # Add detailed explanation about BD vs ML differences when both analyses are present
                if has_ml and has_bayesian:
                    f.write("\n## Understanding BD ≈ ΔlnL in Phylogenetics\n\n")
                    f.write("You may notice that Bayesian Decay (BD) values closely approximate the ML log-likelihood differences (ΔlnL). ")
                    f.write("This is **expected behavior** in phylogenetic topology testing, not an anomaly. Here's why:\n\n")
                    f.write("1. **Identical Models**: The constrained and unconstrained analyses use the same substitution model, ")
                    f.write("differing only in whether a specific clade is allowed to exist.\n\n")
                    f.write("2. **Likelihood Dominance**: When data support a topology, the marginal likelihood ")
                    f.write("(which integrates over all parameters) becomes dominated by the likelihood at the optimal parameter values.\n\n")
                    f.write("3. **Minimal Prior Effects**: Since both analyses explore nearly identical parameter spaces ")
                    f.write("(same model parameters, only different tree topologies), the prior's influence is minimal.\n\n")
                    f.write("**Practical Implications**:\n")
                    f.write("- Similar BD and ΔlnL values confirm your analyses are working correctly\n")
                    f.write("- Large BD values (30-50+) simply reflect strong data signal, not \"astronomical\" support\n")
                    f.write("- Use the phylogenetic-specific BD interpretation scale provided above\n")
                    f.write("- Compare relative BD values across branches rather than focusing on absolute values\n")
                
                # Add file references section
                f.write("\n## Additional Files and References\n\n")
                f.write("This analysis generates several output files for detailed examination:\n\n")
                
                # Main data files
                f.write("### Data Files\n")
                csv_file = output_path.parent / f"{self.alignment_file.stem}_data.csv"
                summary_file = output_path.parent / f"{self.alignment_file.stem}_summary.txt"
                
                if csv_file.exists():
                    f.write(f"- **CSV Export**: [`{csv_file.name}`](./{csv_file.name}) - Complete numerical results in spreadsheet format\n")
                if summary_file.exists():
                    f.write(f"- **Text Summary**: [`{summary_file.name}`](./{summary_file.name}) - Formatted summary with symbols and tables\n")
                
                # Tree files
                trees_dir = output_path.parent / "trees"
                if trees_dir.exists():
                    f.write(f"\n### Tree Files\n")
                    f.write(f"- **Tree Directory**: [`./trees/`](./trees/) - Contains phylogenetic trees with support values:\n")
                    for tree_type in ['au', 'delta_lnl', 'bayes_decay', 'combined']:
                        tree_file = trees_dir / f"{self.alignment_file.stem}_{tree_type}.nwk"
                        if tree_file.exists():
                            f.write(f"  - `{tree_file.name}` - {tree_type.replace('_', ' ').title()} annotated tree\n")
                
                # Site analysis files
                if site_analysis_dir.exists():
                    f.write(f"\n### Site Analysis Files\n")
                    f.write(f"- **Site Analysis Directory**: [`./site_analysis/`](./site_analysis/) - Detailed site-by-site likelihood data\n")
                    f.write(f"  - `site_analysis_summary.txt` - Summary table of site statistics\n")
                    f.write(f"  - `site_data_Clade_X.txt` - Individual site likelihood differences for each clade\n")
                    f.write(f"  - `site_hist_Clade_X.png` - Histogram visualizations of site support\n")
                    f.write(f"  - `site_plot_Clade_X.png` - Detailed site-by-site support plots\n")
                
                # Supplementary files
                supplementary_dir = output_path.parent / "supplementary"
                if supplementary_dir.exists():
                    f.write(f"\n### Configuration and Logs\n")
                    f.write(f"- **Supplementary Directory**: [`./supplementary/`](./supplementary/) - Analysis configuration and metadata\n")
                    config_file = supplementary_dir / f"{self.alignment_file.stem}_analysis_config.txt"
                    if config_file.exists():
                        f.write(f"  - `{config_file.name}` - Complete analysis parameters and settings\n")
                
                    f.write(f"\n---\n")
                    f.write(f"*Generated by panDecay v1.1.0 on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")
                
        except Exception as e:
            logger.error(f"Error generating markdown report: {e}")
            # Write a minimal report on error
            try:
                with output_path.open('w') as f:
                    f.write(f"# panDecay Analysis Report - Error\n\n")
                    f.write(f"Report generation failed with error: {e}\n\n")
                    f.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            except Exception as fallback_error:
                logger.error(f"Failed to write fallback report: {fallback_error}")
                
        logger.info(f"Report saved: {output_path.name}")

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

        logger.info(f"Site summary saved: {summary_path.name}")

        # For each branch, write detailed site data with overwriting progress
        site_data_clades = [(clade_id, data) for clade_id, data in self.decay_indices.items() if 'site_data' in data]
        if site_data_clades:
            overwrite_progress = OverwritingProgress()
            try:
                for i, (clade_id, data) in enumerate(site_data_clades, 1):
                    overwrite_progress.update(f"Saving site data: {clade_id} ({i}/{len(site_data_clades)})")
                    
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
                
                overwrite_progress.finish(f"Site data saved: {len(site_data_clades)} files")
            except Exception:
                overwrite_progress.finish(f"Site data saved: {len(site_data_clades)} files")

        # Generate site analysis visualizations
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np

            # Get visualization options
            viz_format = getattr(self, 'viz_format', 'png')

            # Collect clades that have site data for plotting
            plot_clades = [(clade_id, data) for clade_id, data in self.decay_indices.items() if 'site_data' in data and data.get('site_data', {})]
            
            if plot_clades:
                overwrite_progress = OverwritingProgress()
                plot_count = 0
                try:
                    for i, (clade_id, data) in enumerate(plot_clades, 1):
                        overwrite_progress.update(f"Generating plots: {clade_id} ({i}/{len(plot_clades)})")
                        
                        # Extract data for plotting
                        site_data = data.get('site_data', {})

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

                        plot_count += 1

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
                        
                        plot_count += 1
                    
                    overwrite_progress.finish(f"Plots saved: {plot_count} files ({len(plot_clades)} plots + {len(plot_clades)} histograms)")
                except Exception:
                    overwrite_progress.finish(f"Plots saved: {plot_count} files")

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
            logger.info(f"Distribution plot saved: {output_path.name}")
        except ImportError: logger.error("Matplotlib/Seaborn not found for visualization.")
        except Exception as e: logger.error(f"Failed support distribution plot: {e}")

    def visualize_support_correlation(self, output_file: Path, **kwargs) -> None:
        """
        Create support value correlation plot.
        
        Args:
            output_file: Path to output plot file
            **kwargs: Additional plotting parameters
        """
        import logging
        local_logger = logging.getLogger(__name__)
        
        if not self.decay_indices:
            local_logger.warning("No data for support correlation plot.")
            return
        
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import pandas as pd
            
            # Collect support measures for correlation analysis
            correlation_data = []
            
            for clade_id, data in self.decay_indices.items():
                row = {'clade_id': clade_id}
                
                # Add different support measures (using correct key names)
                if 'AU_pvalue' in data:
                    row['AU_pvalue'] = data['AU_pvalue']
                if 'lnl_diff' in data:
                    row['Delta_LnL'] = abs(data['lnl_diff'])
                if 'bayes_decay' in data:
                    row['Bayesian_Decay'] = data['bayes_decay']
                if 'pars_decay' in data:
                    row['Parsimony_Decay'] = data['pars_decay']
                if 'bootstrap_support' in data:
                    row['Bootstrap_Support'] = data['bootstrap_support']
                
                if len(row) > 1:  # Has at least one support measure
                    correlation_data.append(row)
            
            if len(correlation_data) < 2:
                local_logger.warning("Insufficient data for correlation plot (need at least 2 clades with support measures).")
                return
            
            # Create DataFrame
            df = pd.DataFrame(correlation_data)
            
            # Select only numeric columns for correlation
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            if 'clade_id' in numeric_cols:
                numeric_cols.remove('clade_id')
            
            if len(numeric_cols) < 2:
                local_logger.warning("Insufficient numeric support measures for correlation plot.")
                return
            
            # Create correlation matrix
            corr_data = df[numeric_cols]
            correlation_matrix = corr_data.corr()
            
            # Create plot
            plt.figure(figsize=(kwargs.get('width', 10), kwargs.get('height', 8)))
            
            # Use seaborn heatmap for correlation matrix
            mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
            sns.heatmap(correlation_matrix, 
                       mask=mask,
                       annot=True, 
                       cmap='coolwarm', 
                       center=0,
                       square=True,
                       fmt='.2f',
                       cbar_kws={"shrink": .8})
            
            plt.title('Support Measure Correlation Matrix')
            plt.tight_layout()
            
            # Save plot
            output_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(output_file), format=kwargs.get('format', "png"), dpi=300)
            plt.close()
            
            local_logger.info(f"Support correlation plot saved: {output_file.name}")
            
        except ImportError:
            local_logger.error("Required packages (matplotlib/seaborn/pandas) not found for correlation visualization.")
        except Exception as e:
            local_logger.error(f"Failed to create support correlation plot: {e}")

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

    def export_results_csv(self, output_file: Path) -> None:
        """
        Export results to CSV format.
        
        Args:
            output_file: Path to CSV output file
        """
        try:
            with open(output_file, 'w', newline='') as csvfile:
                if not self.decay_indices:
                    csvfile.write("No results available\n")
                    return
                
                # Define organized fieldnames with proper support measure names
                fieldnames = ['clade_id', 'taxa_count', 'taxa']
                
                # Add support measures in logical order: PD, LD, BD, BS, AU
                # Using correct field names from actual data structure
                support_measures = [
                    ('PD', 'pars_decay', None),
                    ('LD', 'lnl_diff', None), 
                    ('BD', 'bayes_decay', None),
                    ('BS', 'bootstrap_support', None),
                    ('AU', 'AU_pvalue', None)
                ]
                
                # Add site analysis measures
                site_analysis_measures = [
                    ('Supporting_Sites', 'supporting_sites'),
                    ('Conflicting_Sites', 'conflicting_sites'), 
                    ('Neutral_Sites', 'neutral_sites'),
                    ('Site_Support_Ratio', 'support_ratio'),
                    ('Weighted_Support_Ratio', 'weighted_support_ratio'),
                    ('Sum_Supporting_Delta', 'sum_supporting_delta'),
                    ('Sum_Conflicting_Delta', 'sum_conflicting_delta')
                ]
                
                # Check which measures are available and add columns
                available_measures = []
                for short_name, primary_key, fallback_key in support_measures:
                    has_measure = any(
                        primary_key in d or (fallback_key and fallback_key in d) 
                        for d in self.decay_indices.values()
                    )
                    if has_measure:
                        fieldnames.append(short_name)
                        available_measures.append((short_name, primary_key, fallback_key))
                
                # Check for site analysis measures
                available_site_measures = []
                for short_name, primary_key in site_analysis_measures:
                    has_measure = any(
                        primary_key in d for d in self.decay_indices.values()
                    )
                    if has_measure:
                        fieldnames.append(short_name)
                        available_site_measures.append((short_name, primary_key))
                
                # Add any remaining fields (legacy support)
                all_keys = set()
                for decay_data in self.decay_indices.values():
                    all_keys.update(decay_data.keys())
                
                # Exclude internal/technical fields from CSV export
                exclude_keys = {'clade_id', 'taxa_count', 'taxa', 'analysis_types', 
                               'pars_decay', 'lnl_diff', 'bayes_decay', 'bootstrap_support', 'AU_pvalue',
                               'site_data', 'supporting_sites', 'conflicting_sites', 'neutral_sites',
                               'support_ratio', 'weighted_support_ratio', 'sum_supporting_delta', 
                               'sum_conflicting_delta', 'significant_AU', 'constrained_lnl', 
                               'bayes_constrained_ml', 'bayes_unconstrained_ml', 'pars_constrained_score',
                               'pars_score', 'posterior_prob'}
                
                for key in sorted(all_keys):
                    if key not in fieldnames and key not in exclude_keys:
                        fieldnames.append(key)
                
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for clade_id, decay_data in self.decay_indices.items():
                    row_data = {'clade_id': clade_id}
                    taxa_list = decay_data.get('taxa', [])
                    row_data['taxa_count'] = len(taxa_list)
                    row_data['taxa'] = ', '.join(taxa_list)
                    
                    # Add support measures with proper names
                    for short_name, primary_key, fallback_key in available_measures:
                        value = decay_data.get(primary_key)
                        if value is None and fallback_key:
                            value = decay_data.get(fallback_key)
                        row_data[short_name] = value if value is not None else 0
                    
                    # Add site analysis measures
                    for short_name, primary_key in available_site_measures:
                        value = decay_data.get(primary_key, 0)
                        row_data[short_name] = value
                    
                    # Add any remaining fields
                    for key in fieldnames[3 + len(available_measures) + len(available_site_measures):]:  # Skip all processed fields
                        if key in decay_data:
                            value = decay_data[key]
                            if isinstance(value, list):
                                value = ', '.join(map(str, value))
                            row_data[key] = value
                    writer.writerow(row_data)
                
            logger.info(f"Results exported to CSV: {output_file}")
            
        except Exception as e:
            logger.error(f"Failed to export CSV: {e}")
            raise Exception(f"CSV export failed: {e}")


# --- Main Execution Logic ---