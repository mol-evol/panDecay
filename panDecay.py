#!/usr/bin/env python3

import os
import sys
import argparse
import configparser
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
    Implements ML-based phylogenetic decay indices (Bremer support) using PAUP*.
    Calculates support by comparing optimal tree likelihood with constrained trees,
    using PAUP*'s backbone constraint followed by reverse constraints and AU test.
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
                 max_psrf=1.01, max_asdsf=0.01, convergence_strict=False):

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
        self._convert_to_nexus() # Writes to self.nexus_file_path

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

        logger.info(f"Running PAUP* command file: {paup_cmd_filename_str} (Log: {log_filename_str})")

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
        if self.debug: logger.debug(f"ML search PAUP* script ({ml_search_cmd_path}):\n{paup_script_content}")

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

                logger.info(f"Cleaned tree file {tree_path} - removed metadata after semicolon")
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
            
            logger.info(f"Running MrBayes: {' '.join(cmd)} in directory {self.temp_path}")
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
            
            # Log successful completion
            logger.info(f"MrBayes completed successfully for {output_prefix}")
            
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
                    logger.info(f"Parsing posterior probabilities from {con_tree_path}")
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
        logger.info("Running unconstrained Bayesian analysis...")
        
        # Create NEXUS file with MrBayes block
        nexus_content = self.nexus_file_path.read_text()
        mrbayes_block = self._generate_mrbayes_nexus()
        combined_nexus = nexus_content + "\n" + mrbayes_block
        
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
        
        for i, clade_obj in enumerate(internal_clades):
            clade_log_idx = i + 1  # For filenames and logging (1-based)
            clade_taxa = [leaf.name for leaf in clade_obj.get_terminals()]
            total_taxa_count = len(self.ml_tree.get_terminals())
            
            # Skip trivial branches (consistent with ML analysis)
            if len(clade_taxa) <= 1 or len(clade_taxa) >= total_taxa_count - 1:
                logger.info(f"Skipping trivial branch {clade_log_idx} (taxa: {len(clade_taxa)}/{total_taxa_count}).")
                continue
            
            # Check if this clade should be tested based on constraint mode
            if not self.should_test_clade(clade_taxa, user_constraints):
                logger.info(f"Skipping branch {clade_log_idx} based on constraint mode '{self.constraint_mode}'")
                continue
            
            clade_id = f"Clade_{clade_log_idx}"
            logger.info(f"Running Bayesian constraint analysis for {clade_id} (taxa: {len(clade_taxa)})")
            
            # Create constrained NEXUS file
            mrbayes_block = self._generate_mrbayes_nexus(
                clade_taxa=clade_taxa, 
                constraint_id=clade_id
            )
            combined_nexus = nexus_content + "\n" + mrbayes_block
            
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
                
                # Calculate Bayes Factor with reasonable bounds for display
                # Note: log(BF) = Bayes Decay
                # Cap display at BF = 10^6 (decay ~13.8) for readability
                MAX_DISPLAY_DECAY = 13.8  # log(10^6)
                
                if bayes_decay > MAX_DISPLAY_DECAY:
                    bayes_factor = 10**6  # Display cap
                    bayes_factor_display = ">10^6"
                    if bayes_decay > 10:
                        logger.warning(f"{clade_id}: Large Bayes Decay ({bayes_decay:.2f}) may reflect model dimension differences")
                elif bayes_decay < -MAX_DISPLAY_DECAY:
                    bayes_factor = 10**-6  # Display cap
                    bayes_factor_display = "<10^-6"
                else:
                    bayes_factor = np.exp(bayes_decay)
                    if bayes_factor > 1000:
                        bayes_factor_display = f"{bayes_factor:.2e}"
                    else:
                        bayes_factor_display = f"{bayes_factor:.2f}"
                
                bayesian_results[clade_id] = {
                    'unconstrained_ml': unconstrained_ml,
                    'constrained_ml': constrained_ml,
                    'bayes_decay': bayes_decay,
                    'bayes_factor': bayes_factor,
                    'bayes_factor_display': bayes_factor_display,
                    'taxa': clade_taxa
                }
                
                logger.info(f"{clade_id}: Bayes decay = {bayes_decay:.4f}, BF = {bayes_factor_display}")
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
                logger.info(f"Calculated stepping-stone marginal likelihood for {output_prefix}: {ml_value}")
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
            posterior_probs = {}
            
            # Read the consensus tree file
            tree_content = con_tree_path.read_text()
            
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
            
            # Remove [&U] or other tree attributes
            if newick_str.startswith("["):
                end_bracket = newick_str.find("]")
                if end_bracket != -1:
                    newick_str = newick_str[end_bracket+1:].strip()
            
            # Parse the tree to extract clades and their posterior probabilities
            # This is a simplified parser - for production use, consider using Bio.Phylo
            # or another proper Newick parser that handles posterior probabilities
            
            # For now, log that we found the tree
            logger.info(f"Found consensus tree with posterior probabilities")
            
            # Parse using BioPython
            from io import StringIO
            tree_io = StringIO(newick_str)
            
            try:
                tree = Phylo.read(tree_io, "newick")
                
                # Extract posterior probabilities from internal nodes
                for node in tree.get_nonterminals():
                    if node.confidence is not None and node != tree.root:
                        # Get taxa in this clade
                        clade_taxa = frozenset(term.name for term in node.get_terminals())
                        posterior_probs[clade_taxa] = float(node.confidence)
                        
                logger.info(f"Extracted posterior probabilities for {len(posterior_probs)} clades")
                
            except Exception as e:
                logger.warning(f"Could not parse consensus tree with BioPython: {e}")
                # Could implement manual parsing here if needed
                
            return posterior_probs
            
        except Exception as e:
            logger.error(f"Error parsing MrBayes posterior probabilities: {e}")
            if self.debug:
                import traceback
                logger.debug(traceback.format_exc())
            return {}

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
        
        # Log convergence metrics
        logger.info(f"Convergence diagnostics for {output_prefix}:")
        if convergence_data['min_ess'] is not None:
            logger.info(f"  Minimum ESS: {convergence_data['min_ess']:.0f}")
        if convergence_data['max_psrf'] is not None:
            logger.info(f"  Maximum PSRF: {convergence_data['max_psrf']:.3f}")
        if convergence_data['asdsf'] is not None:
            logger.info(f"  Final ASDSF: {convergence_data['asdsf']:.6f}")
        
        # Log warnings
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
            
            # Process branches with same logic as ML analysis
            for i, clade_obj in enumerate(branches):
                clade_log_idx = i + 1  # For filenames and logging (1-based)
                clade_taxa = [leaf.name for leaf in clade_obj.get_terminals()]
                clade_id = f"Clade_{clade_log_idx}"
                
                # Skip trivial branches (consistent with ML analysis)
                if len(clade_taxa) <= 1 or len(clade_taxa) >= total_taxa_count - 1:
                    logger.info(f"Skipping trivial branch {clade_log_idx} (taxa: {len(clade_taxa)}/{total_taxa_count}).")
                    continue
                
                # Check if this clade should be tested based on constraint mode
                if not self.should_test_clade(clade_taxa, user_constraints):
                    logger.info(f"Skipping branch {clade_log_idx} based on constraint mode '{self.constraint_mode}'")
                    continue
                
                logger.info(f"Calculating parsimony decay for {clade_id} ({len(clade_taxa)} taxa)")
                
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
                        logger.debug(f"Constraint {idx} parsimony score file content:\n{score_content}")
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
                            'bayes_factor': bayes_data.get('bayes_factor')
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

        for i, clade_obj in enumerate(internal_clades):
            clade_log_idx = i + 1 # For filenames and logging (1-based)
            clade_taxa_names = [leaf.name for leaf in clade_obj.get_terminals()]
            total_taxa_count = len(self.ml_tree.get_terminals())

            if len(clade_taxa_names) <= 1 or len(clade_taxa_names) >= total_taxa_count -1:
                logger.info(f"Skipping trivial branch {clade_log_idx} (taxa: {len(clade_taxa_names)}/{total_taxa_count}).")
                continue
            
            # Check if this clade should be tested based on constraint mode
            if not self.should_test_clade(clade_taxa_names, user_constraints):
                logger.info(f"Skipping branch {clade_log_idx} based on constraint mode '{self.constraint_mode}'")
                continue

            logger.info(f"Processing branch {clade_log_idx}/{len(internal_clades)} (taxa: {len(clade_taxa_names)})")
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
    
    def _calculate_bayesian_decay_indices(self):
        """
        Calculate Bayesian decay indices for all internal branches.
        
        Returns:
            Dictionary of decay indices with Bayesian metrics
        """
        logger.info("Calculating Bayesian decay indices...")
        
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
                'bayes_factor': bayes_data['bayes_factor'],
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
                        'bayes_decay': bayes_data['bayes_decay'],
                        'bayes_factor': bayes_data['bayes_factor']
                    })
                else:
                    # No Bayesian results for this clade
                    ml_results[clade_id].update({
                        'bayes_unconstrained_ml': None,
                        'bayes_constrained_ml': None,
                        'bayes_decay': None,
                        'bayes_factor': None
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
                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {au_tree_path} (type: au).")
                tree_files['au'] = au_tree_path
            except Exception as e:
                logger.error(f"Failed to create AU tree: {e}")

            # Create log-likelihood difference annotated tree
            lnl_tree_path = output_dir / f"{base_filename}_lnl.nwk"
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
                        node.name = f"{matched_clade_id} - LnL:{lnl_diff:.4f}"
                        annotated_nodes_count += 1

                Phylo.write(lnl_tree, str(lnl_tree_path), "newick")
                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {lnl_tree_path} (type: lnl).")
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
                            bayes_factor = decay_info.get('bayes_factor')

                            if au_val is not None:
                                annotation_parts.append(f"AU:{au_val:.4f}")

                            if lnl_val is not None:
                                annotation_parts.append(f"LnL:{abs(lnl_val):.4f}")
                                
                            if bayes_decay is not None:
                                annotation_parts.append(f"BD:{bayes_decay:.4f}")
                                
                            # Use display string for Bayes Factor if available
                            bayes_factor_display = decay_info.get('bayes_factor_display')
                            if bayes_factor_display:
                                annotation_parts.append(f"BF:{bayes_factor_display}")
                            elif bayes_factor is not None:
                                # Fallback for older results
                                if bayes_factor > 10**6:
                                    annotation_parts.append(f"BF:>10^6")
                                elif bayes_factor > 1000:
                                    annotation_parts.append(f"BF:{bayes_factor:.2e}")
                                else:
                                    annotation_parts.append(f"BF:{bayes_factor:.2f}")
                                
                            # Add parsimony decay if available
                            pars_decay = decay_info.get('pars_decay')
                            if pars_decay is not None:
                                annotation_parts.append(f"PD:{pars_decay}")
                                
                            # Add posterior probability if available
                            post_prob = decay_info.get('posterior_prob')
                            if post_prob is not None:
                                annotation_parts.append(f"PP:{post_prob:.2f}")
                            break

                    # Add the clade ID to the beginning of annotation if we have one
                    if clade_id and annotation_parts:
                        # Create clear separation between clade name and metrics
                        clade_part = clade_id
                        metrics_part = "|".join(annotation_parts)
                        node_annotations[node_taxa_set] = f"{clade_part} - {metrics_part}"
                    elif annotation_parts:
                        node_annotations[node_taxa_set] = "|".join(annotation_parts)

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

                logger.info(f"Annotated tree with {annotated_nodes_count} branch values written to {combined_tree_path} (type: combined).")
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
                    logger.info(f"Bootstrap tree (ML topology with bootstrap values) written to {bootstrap_tree_path}")
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
                            bayes_factor = matched_data.get('bayes_factor')

                            if au_val is not None:
                                annotation_parts.append(f"AU:{au_val:.4f}")

                            if lnl_val is not None:
                                annotation_parts.append(f"LnL:{abs(lnl_val):.4f}")
                                
                            if bayes_decay is not None:
                                annotation_parts.append(f"BD:{bayes_decay:.4f}")
                                
                            # Use display string for Bayes Factor if available
                            bayes_factor_display = decay_info.get('bayes_factor_display')
                            if bayes_factor_display:
                                annotation_parts.append(f"BF:{bayes_factor_display}")
                            elif bayes_factor is not None:
                                # Fallback for older results
                                if bayes_factor > 10**6:
                                    annotation_parts.append(f"BF:>10^6")
                                elif bayes_factor > 1000:
                                    annotation_parts.append(f"BF:{bayes_factor:.2e}")
                                else:
                                    annotation_parts.append(f"BF:{bayes_factor:.2f}")
                                
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
                    logger.info(f"Comprehensive tree with {annotated_nodes_count} branch values written to {comprehensive_tree_path}")
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
                    logger.info(f"Annotated tree with {annotated_nodes_count} Bayes decay values written to {bayes_decay_tree_path}")
                    tree_files['bayes_decay'] = bayes_decay_tree_path
                except Exception as e:
                    logger.error(f"Failed to create Bayes decay tree: {e}")

                # Create Bayes factor annotated tree
                bayes_factor_tree_path = output_dir / f"{base_filename}_bayes_factor.nwk"
                try:
                    temp_tree_for_bf = self.temp_path / f"ml_tree_for_bf_annotation.nwk"
                    Phylo.write(self.ml_tree, str(temp_tree_for_bf), "newick")
                    cleaned_tree_path = self._clean_newick_tree(temp_tree_for_bf)
                    bf_tree = Phylo.read(str(cleaned_tree_path), "newick")

                    annotated_nodes_count = 0
                    for node in bf_tree.get_nonterminals():
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
                        if matched_data and 'bayes_factor' in matched_data and matched_data['bayes_factor'] is not None:
                            # Use display string if available
                            bayes_factor_display = matched_data.get('bayes_factor_display')
                            if bayes_factor_display:
                                node.name = f"{matched_clade_id} - BF:{bayes_factor_display}"
                            else:
                                bayes_factor_val = matched_data['bayes_factor']
                                if bayes_factor_val > 10**6:
                                    node.name = f"{matched_clade_id} - BF:>10^6"
                                elif bayes_factor_val > 1000:
                                    node.name = f"{matched_clade_id} - BF:{bayes_factor_val:.2e}"
                                else:
                                    node.name = f"{matched_clade_id} - BF:{bayes_factor_val:.2f}"
                            annotated_nodes_count += 1

                    Phylo.write(bf_tree, str(bayes_factor_tree_path), "newick")
                    logger.info(f"Annotated tree with {annotated_nodes_count} Bayes factor values written to {bayes_factor_tree_path}")
                    tree_files['bayes_factor'] = bayes_factor_tree_path
                except Exception as e:
                    logger.error(f"Failed to create Bayes factor tree: {e}")

            return tree_files

        except Exception as e:
            logger.error(f"Failed to annotate trees: {e}")
            if hasattr(self, 'debug') and self.debug:
                import traceback
                logger.debug(f"Traceback: {traceback.format_exc()}")
            return tree_files  # Return any successfully created files

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
                header_parts.extend(["Constrained_lnL", "LnL_Diff_from_ML", "AU_p-value", "Significant_AU (p<0.05)"])
            if has_parsimony:
                header_parts.append("Pars_Decay")
            if has_bayesian:
                header_parts.extend(["Bayes_ML_Diff", "Bayes_Factor"])
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
                    
                    bayes_factor = data.get('bayes_factor', 'N/A')
                    if isinstance(bayes_factor, float): 
                        # Use scientific notation for large Bayes factors
                        if bayes_factor > 1000:
                            bayes_factor = f"{bayes_factor:.2e}"
                        else:
                            bayes_factor = f"{bayes_factor:.2f}"
                    elif bayes_factor is None: bayes_factor = 'N/A'
                    
                    row_parts.extend([bayes_decay, bayes_factor])
                    
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

        logger.info(f"Results written to {output_path}")

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
                        
                        # Check for extreme positive values
                        extreme_count = sum(1 for bd in bayes_decays if bd > 10)
                        if extreme_count > 0:
                            f.write(f"\n**⚠️ NOTE**: {extreme_count}/{len(bayes_decays)} branches have Bayes Decay values > 10.\n")
                            f.write("Very large values may partially reflect:\n")
                            f.write("- Model dimension differences between constrained/unconstrained trees\n")
                            f.write("- Prior distribution effects\n")
                            f.write("- Genuine very strong support for the clade\n\n")
                        
                    bayes_factors = [d['bayes_factor'] for d in self.decay_indices.values() if d.get('bayes_factor') is not None]
                    if bayes_factors:
                        strong_bf_count = sum(1 for bf in bayes_factors if bf > 10)
                        f.write(f"- Branches with strong Bayes Factor support (BF > 10): {strong_bf_count} / {len(bayes_factors)} evaluated\n")
                
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
                header_parts.extend(["| Constrained lnL | LnL Diff from ML | AU p-value | Significant (AU) "])
                separator_parts.extend(["|-----------------|------------------|------------|-------------------- "])
            
            if has_bayesian:
                header_parts.extend(["| Bayes Decay | Bayes Factor | BF Support "])
                separator_parts.extend(["|-------------|--------------|------------ "])
                
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
                    
                    bayes_f_raw = data.get('bayes_factor', 'N/A')
                    
                    # Interpret support based on Bayes Decay
                    bf_support = 'N/A'
                    bayes_d_val = data.get('bayes_decay')
                    if isinstance(bayes_d_val, float):
                        if bayes_d_val > 5:
                            bf_support = "**Very Strong**"
                        elif bayes_d_val > 3:
                            bf_support = "**Strong**"
                        elif bayes_d_val > 1:
                            bf_support = "Positive"
                        elif bayes_d_val > 0:
                            bf_support = "Weak"
                        else:
                            bf_support = "None"
                    
                    # Format for display
                    bayes_f_display = data.get('bayes_factor_display')
                    if bayes_f_display:
                        bayes_f = bayes_f_display
                    elif isinstance(bayes_f_raw, float):
                        if bayes_f_raw > 10**6:
                            bayes_f = ">10^6"
                        elif bayes_f_raw > 1000:
                            bayes_f = f"{bayes_f_raw:.2e}"
                        else:
                            bayes_f = f"{bayes_f_raw:.2f}"
                    else:
                        bayes_f = bayes_f_raw
                    
                    row_parts.append(f"| {bayes_d} | {bayes_f} | {bf_support} ")

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
                f.write("- **LnL Diff from ML**: Log-likelihood of the best tree *without* the clade minus ML tree's log-likelihood. More negative (larger absolute difference) implies stronger support for the clade's presence in the ML tree.\n")
                f.write("- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.\n\n")
                
            if has_bayesian:
                f.write("### Bayesian Analysis\n")
                f.write("- **Bayes Decay (BD)**: Marginal log-likelihood difference (unconstrained - constrained). This is the primary metric for Bayesian support.\n")
                f.write("  - **Note**: Bayes Decay = log(Bayes Factor), so BD is more interpretable than BF for large values\n")
                f.write("  - Interpretation scale:\n")
                f.write("    - BD 0-1: Weak evidence for clade\n")
                f.write("    - BD 1-3: Positive evidence\n")
                f.write("    - BD 3-5: Strong evidence\n")
                f.write("    - BD > 5: Very strong evidence\n")
                f.write("  - **⚠️ Important**: BD values > 10 may partially reflect model dimension differences between constrained/unconstrained trees\n")
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
                f.write("\n- **Bayes Factor (BF)**: Exponential of the Bayes Decay value (BF = e^BD). Represents the ratio of marginal likelihoods.\n")
                f.write("  - Traditional interpretations:\n")
                f.write("    - BF > 100: Very strong support\n")
                f.write("    - BF > 10: Strong support\n")
                f.write("    - BF > 3: Moderate support\n")
                f.write("    - BF > 1: Weak support\n")
                f.write("    - BF ≤ 1: No support\n")
                f.write("  - **Note**: BF values are capped at 10^6 for display to avoid astronomical numbers\n\n")
                
            if has_bootstrap:
                f.write("- **Bootstrap**: Bootstrap support value (percentage of bootstrap replicates in which the clade appears). Higher values (e.g., > 70) suggest stronger support for the clade.\n")
                
        logger.info(f"Detailed report written to {output_path}")

    def write_site_analysis_results(self, output_dir: Path, keep_tree_files=False):
        """
        Write site-specific likelihood analysis results to files.

        Args:
            output_dir: Directory to save the site analysis files
            keep_tree_files: Whether to keep the Newick files used for HTML visualization
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

        logger.info(f"Site analysis summary written to {summary_path}")

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

            logger.info(f"Detailed site data for {clade_id} written to {site_data_path}")

        # Generate site analysis visualizations
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np

            # Get visualization options
            viz_format = getattr(self, 'viz_format', 'png')
            generate_html = getattr(self, 'generate_html', True)
            js_cdn = getattr(self, 'js_cdn', True)

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

                # Create interactive HTML tree visualization if enabled
                if generate_html and clade_taxa:
                    # Create HTML tree visualization
                    html_path = self.create_interactive_tree_html(output_dir, clade_id, clade_taxa)
                    if html_path:
                        logger.info(f"Interactive tree visualization for {clade_id} created at {html_path}")

                if not keep_tree_files and not self.debug and not self.keep_files:
                    for file_path in output_dir.glob("tree_*.nwk"):
                        try:
                            file_path.unlink()
                            logger.debug(f"Deleted tree file for HTML: {file_path}")
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

    def create_interactive_tree_html(self, output_dir, clade_id, highlight_taxa):
        """
        Create an interactive HTML tree visualization for a specific clade.

        Args:
            output_dir: Directory to save the HTML file
            clade_id: Identifier for the clade
            highlight_taxa: List of taxa names to highlight in the tree

        Returns:
            Path to the created HTML file or None if creation failed
        """
        try:
            import json
            from Bio import Phylo

            # Ensure output directory exists
            output_dir.mkdir(parents=True, exist_ok=True)

            # Output filenames
            html_path = output_dir / f"tree_{clade_id}.html"
            tree_path = output_dir / f"tree_{clade_id}.nwk"

            # Write the ML tree to a Newick file (needed for the HTML to load)
            Phylo.write(self.ml_tree, str(tree_path), "newick")

            # Clean up the tree file if needed
            cleaned_tree_path = self._clean_newick_tree(tree_path)

            # Get tree statistics
            total_taxa = len(self.ml_tree.get_terminals())
            highlight_ratio = len(highlight_taxa) / total_taxa if total_taxa > 0 else 0

            # Get site analysis data if available
            site_data = None
            if hasattr(self, 'decay_indices') and clade_id in self.decay_indices:
                clade_data = self.decay_indices[clade_id]
                if 'site_data' in clade_data:
                    supporting = clade_data.get('supporting_sites', 0)
                    conflicting = clade_data.get('conflicting_sites', 0)
                    support_ratio = clade_data.get('support_ratio', None)
                    if support_ratio == float('inf'):
                        support_ratio_str = "Infinity"
                    elif support_ratio is not None:
                        support_ratio_str = f"{support_ratio:.2f}"
                    else:
                        support_ratio_str = "N/A"

                    weighted_ratio = clade_data.get('weighted_support_ratio', None)
                    if weighted_ratio == float('inf'):
                        weighted_ratio_str = "Infinity"
                    elif weighted_ratio is not None:
                        weighted_ratio_str = f"{weighted_ratio:.2f}"
                    else:
                        weighted_ratio_str = "N/A"

                    site_data = {
                        'supporting': supporting,
                        'conflicting': conflicting,
                        'support_ratio': support_ratio_str,
                        'weighted_ratio': weighted_ratio_str
                    }

            # Get AU test and likelihood data
            au_data = None
            if hasattr(self, 'decay_indices') and clade_id in self.decay_indices:
                clade_data = self.decay_indices[clade_id]
                au_pvalue = clade_data.get('AU_pvalue', None)
                lnl_diff = clade_data.get('lnl_diff', None)

                if au_pvalue is not None or lnl_diff is not None:
                    au_data = {
                        'au_pvalue': f"{au_pvalue:.4f}" if au_pvalue is not None else "N/A",
                        'lnl_diff': f"{lnl_diff:.4f}" if lnl_diff is not None else "N/A",
                        'significant': clade_data.get('significant_AU', False)
                    }

            # Get bootstrap data if available
            bootstrap_value = None
            if hasattr(self, 'bootstrap_tree') and self.bootstrap_tree:
                # Find the corresponding node in the bootstrap tree
                for node in self.bootstrap_tree.get_nonterminals():
                    node_taxa = set(leaf.name for leaf in node.get_terminals())
                    if node_taxa == set(highlight_taxa):
                        bootstrap_value = int(node.confidence) if node.confidence is not None else None
                        break

            # Format taxa for title display
            if len(highlight_taxa) <= 5:
                taxa_display = ", ".join(highlight_taxa)
            else:
                taxa_display = f"{', '.join(sorted(highlight_taxa)[:5])}... (+{len(highlight_taxa)-5} more)"

            # Start building HTML content in parts to avoid f-string backslash issues
            # Basic structure
            html_parts = []

            # Header
            html_parts.append("<!DOCTYPE html>")
            html_parts.append("<html lang=\"en\">")
            html_parts.append("<head>")
            html_parts.append("    <meta charset=\"UTF-8\">")
            html_parts.append("    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">")
            html_parts.append(f"    <title>panDecay - Interactive Tree for {clade_id}</title>")

            # CSS
            html_parts.append("""    <style>
            body {
                font-family: Arial, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
            }
            h1, h2, h3 {
                color: #2c3e50;
            }
            .container {
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
            }
            .tree-section {
                flex: 1;
                min-width: 500px;
            }
            .info-section {
                flex: 1;
                min-width: 300px;
                background-color: #f8f9fa;
                padding: 15px;
                border-radius: 5px;
            }
            #tree_container {
                width: 100%;
                height: 600px;
                border: 1px solid #ddd;
                border-radius: 4px;
                overflow: hidden;
            }
            .table {
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 15px;
            }
            .table th, .table td {
                padding: 8px;
                text-align: left;
                border-bottom: 1px solid #ddd;
            }
            .table th {
                background-color: #f2f2f2;
            }
            .highlight {
                color: #e74c3c;
                font-weight: bold;
            }
            .significant {
                background-color: #d4edda;
            }
            .not-significant {
                background-color: #f8d7da;
            }
            .buttons {
                margin: 10px 0;
            }
            button {
                background-color: #4CAF50;
                color: white;
                padding: 8px 12px;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                margin-right: 5px;
            }
            button:hover {
                background-color: #45a049;
            }
            .phylotree-node circle {
                fill: #999;
            }
            .highlighted-node circle {
                fill: #e74c3c !important;
                r: 5 !important;
            }
            .highlighted-node text {
                fill: #e74c3c !important;
                font-weight: bold !important;
            }
            .highlighted-branch {
                stroke: #e74c3c !important;
                stroke-width: 3px !important;
            }
            .legend {
                margin-top: 10px;
                font-size: 0.9em;
            }
            .legend-item {
                display: inline-block;
                margin-right: 15px;
            }
            .legend-color {
                display: inline-block;
                width: 12px;
                height: 12px;
                margin-right: 5px;
                vertical-align: middle;
            }
            .download-links {
                margin-top: 20px;
            }
            .download-links a {
                display: inline-block;
                margin-right: 10px;
                padding: 5px 10px;
                background-color: #f8f9fa;
                border: 1px solid #ddd;
                border-radius: 3px;
                text-decoration: none;
                color: #333;
            }
            .download-links a:hover {
                background-color: #e9ecef;
            }
        </style>""")

            # JavaScript imports based on CDN preference
            use_cdn = getattr(self, 'js_cdn', True)
            if use_cdn:
                html_parts.append("""    <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.8.5/d3.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/phylotree@1.0.0-alpha.3/dist/phylotree.js"></script>""")
            else:
                html_parts.append("""    <!-- Embedded D3 and Phylotree libraries would go here -->
        <!-- This would make the file much larger -->""")

            html_parts.append("</head>")
            html_parts.append("<body>")

            # Body content
            html_parts.append(f"    <h1>panDecay - Interactive Tree Visualization</h1>")
            html_parts.append(f"    <h2>Clade: {clade_id} - {taxa_display}</h2>")

            html_parts.append("    <div class=\"container\">")
            html_parts.append("        <div class=\"tree-section\">")
            html_parts.append("            <div class=\"buttons\">")
            html_parts.append("                <button onclick=\"tree.spacing_x(100).spacing_y(20).update()\">Expand Tree</button>")
            html_parts.append("                <button onclick=\"tree.spacing_x(30).spacing_y(10).update()\">Compact Tree</button>")
            html_parts.append("                <button onclick=\"resetTree()\">Reset View</button>")
            html_parts.append("                <button onclick=\"toggleLabels()\">Toggle Labels</button>")
            html_parts.append("            </div>")
            html_parts.append("            <div id=\"tree_container\"></div>")
            html_parts.append("            <div class=\"legend\">")
            html_parts.append("                <div class=\"legend-item\">")
            html_parts.append("                    <span class=\"legend-color\" style=\"background-color: #e74c3c;\"></span>")
            html_parts.append("                    <span>Highlighted Clade</span>")
            html_parts.append("                </div>")
            html_parts.append("                <div class=\"legend-item\">")
            html_parts.append("                    <span class=\"legend-color\" style=\"background-color: #999;\"></span>")
            html_parts.append("                    <span>Other Nodes</span>")
            html_parts.append("                </div>")
            html_parts.append("            </div>")
            html_parts.append("        </div>")

            # Info section
            html_parts.append("        <div class=\"info-section\">")
            html_parts.append("            <h3>Clade Information</h3>")
            html_parts.append("            <table class=\"table\">")
            html_parts.append(f"                <tr><th>Number of Taxa:</th><td>{len(highlight_taxa)} of {total_taxa} ({highlight_ratio:.1%})</td></tr>")

            taxa_sample = ", ".join(sorted(highlight_taxa)[:10])
            if len(highlight_taxa) > 10:
                taxa_sample += " ..."
            html_parts.append(f"                <tr><th>Taxa:</th><td>{taxa_sample}</td></tr>")

            # Bootstrap value if available
            if bootstrap_value is not None:
                html_parts.append(f"                <tr><th>Bootstrap Support:</th><td>{bootstrap_value}%</td></tr>")

            html_parts.append("            </table>")

            # Branch support section if available
            if au_data:
                html_parts.append("            <h3>Branch Support</h3>")
                html_parts.append("            <table class=\"table\">")

                significance_class = 'significant' if au_data.get('significant') else 'not-significant'
                significance_text = '(significant)' if au_data.get('significant') else '(not significant)'
                html_parts.append(f"                <tr><th>AU Test p-value:</th><td class=\"{significance_class}\">{au_data['au_pvalue']} {significance_text}</td></tr>")
                html_parts.append(f"                <tr><th>Log-Likelihood Difference:</th><td>{au_data['lnl_diff']}</td></tr>")

                html_parts.append("            </table>")

            # Site analysis section if available
            if site_data:
                html_parts.append("            <h3>Site Analysis</h3>")
                html_parts.append("            <table class=\"table\">")
                html_parts.append(f"                <tr><th>Supporting Sites:</th><td>{site_data['supporting']}</td></tr>")
                html_parts.append(f"                <tr><th>Conflicting Sites:</th><td>{site_data['conflicting']}</td></tr>")
                html_parts.append(f"                <tr><th>Support Ratio:</th><td>{site_data['support_ratio']}</td></tr>")
                html_parts.append(f"                <tr><th>Weighted Support Ratio:</th><td>{site_data['weighted_ratio']}</td></tr>")
                html_parts.append("            </table>")

            # Download section
            html_parts.append("            <h3>Downloads</h3>")
            html_parts.append("            <div class=\"download-links\">")
            html_parts.append(f"                <a href=\"{tree_path.name}\" download>Download Newick Tree</a>")
            html_parts.append("                <a href=\"#\" onclick=\"saveSvg()\">Download SVG</a>")
            html_parts.append("            </div>")
            html_parts.append("        </div>")
            html_parts.append("    </div>")

            # JavaScript section
            html_parts.append("    <script>")
            html_parts.append(f"        // Taxa to highlight")
            html_parts.append(f"        const highlightTaxa = {json.dumps(list(highlight_taxa))};")
            html_parts.append("")
            html_parts.append("        // Create tree")
            tree_path_js = str(tree_path.name).replace('\\', '/')
            html_parts.append(f"        let tree = new phylotree.phylotree(\"{tree_path_js}\");")
            html_parts.append("        let showLabels = true;")
            html_parts.append("")
            html_parts.append("        // Initialize visualization on page load")
            html_parts.append("        document.addEventListener(\"DOMContentLoaded\", function() {")
            html_parts.append("            loadAndDisplayTree();")
            html_parts.append("        });")
            html_parts.append("")
            html_parts.append("        function loadAndDisplayTree() {")
            html_parts.append(f"            fetch(\"{tree_path_js}\").then(response => {{")
            html_parts.append("                if (response.ok) {")
            html_parts.append("                    return response.text();")
            html_parts.append("                }")
            html_parts.append("                throw new Error('Tree file not found');")
            html_parts.append("            })")
            html_parts.append("            .then(treeData => {")
            html_parts.append("                // Set up tree visualization")
            html_parts.append("                tree = new phylotree.phylotree(treeData);")
            html_parts.append("")
            html_parts.append("                // Configure tree display settings")
            html_parts.append("                tree.branch_length(null)  // Use branch lengths from the tree")
            html_parts.append("                    .branch_name(function(node) {")
            html_parts.append("                        return node.data.name;")
            html_parts.append("                    })")
            html_parts.append("                    .node_span(function(node) {")
            html_parts.append("                        return showLabels ? 5 : 2;")
            html_parts.append("                    })")
            html_parts.append("                    .node_circle_size(function(node) {")
            html_parts.append("                        return isHighlighted(node) ? 5 : 3;")
            html_parts.append("                    })")
            html_parts.append("                    .font_size(14)")
            html_parts.append("                    .scale_bar_font_size(12)")
            html_parts.append("                    .node_styler(nodeStyler)")
            html_parts.append("                    .branch_styler(branchStyler)")
            html_parts.append("                    .layout_handler(phylotree.layout_handlers.radial)")
            html_parts.append("                    .spacing_x(40) // Controls horizontal spacing")
            html_parts.append("                    .spacing_y(15) // Controls vertical spacing")
            html_parts.append("                    .size([550, 550])")
            html_parts.append("                    .radial(false); // Start with rectangular layout")
            html_parts.append("")
            html_parts.append("                // Get the container")
            html_parts.append("                let container = document.getElementById('tree_container');")
            html_parts.append("")
            html_parts.append("                // Render the tree")
            html_parts.append("                tree.render(\"#tree_container\");")
            html_parts.append("            })")
            html_parts.append("            .catch(error => {")
            html_parts.append("                console.error(\"Error loading tree data:\", error);")
            html_parts.append("                document.getElementById('tree_container').innerHTML =")
            html_parts.append("                    \"<p style='color:red;padding:20px;'>Error loading tree data. Check console for details.</p>\";")
            html_parts.append("            });")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Style branches that belong to the highlighted clade")
            html_parts.append("        function branchStyler(dom_element, link_data) {")
            html_parts.append("            if (isHighlighted(link_data.target)) {")
            html_parts.append("                dom_element.style.stroke = \"#e74c3c\";")
            html_parts.append("                dom_element.style.strokeWidth = \"3px\";")
            html_parts.append("                dom_element.classList.add(\"highlighted-branch\");")
            html_parts.append("            }")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Style nodes that belong to the highlighted clade")
            html_parts.append("        function nodeStyler(dom_element, node_data) {")
            html_parts.append("            if (isHighlighted(node_data)) {")
            html_parts.append("                dom_element.classList.add(\"highlighted-node\");")
            html_parts.append("            }")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Check if a node belongs to the highlighted clade")
            html_parts.append("        function isHighlighted(node) {")
            html_parts.append("            if (!node.data) return false;")
            html_parts.append("")
            html_parts.append("            // Directly check leaf nodes")
            html_parts.append("            if (node.children && node.children.length === 0) {")
            html_parts.append("                return highlightTaxa.includes(node.data.name);")
            html_parts.append("            }")
            html_parts.append("")
            html_parts.append("            // For internal nodes, check if all descendants are in the highlighted taxa")
            html_parts.append("            let leaves = getAllLeaves(node);")
            html_parts.append("            let leafNames = leaves.map(leaf => leaf.data.name);")
            html_parts.append("")
            html_parts.append("            if (leafNames.length === 0) return false;")
            html_parts.append("")
            html_parts.append("            // Check if this node represents exactly our clade")
            html_parts.append("            // or is contained within our clade")
            html_parts.append("            return leafNames.every(name => highlightTaxa.includes(name));")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Get all leaves (terminal nodes) descending from a node")
            html_parts.append("        function getAllLeaves(node) {")
            html_parts.append("            if (!node.children || node.children.length === 0) {")
            html_parts.append("                return [node];")
            html_parts.append("            }")
            html_parts.append("")
            html_parts.append("            let leaves = [];")
            html_parts.append("            for (let child of node.children) {")
            html_parts.append("                leaves = leaves.concat(getAllLeaves(child));")
            html_parts.append("            }")
            html_parts.append("            return leaves;")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Reset tree to default view")
            html_parts.append("        function resetTree() {")
            html_parts.append("            tree.spacing_x(40)")
            html_parts.append("                .spacing_y(15)")
            html_parts.append("                .radial(false)")
            html_parts.append("                .update();")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Toggle display of leaf labels")
            html_parts.append("        function toggleLabels() {")
            html_parts.append("            showLabels = !showLabels;")
            html_parts.append("            tree.node_span(function(node) {")
            html_parts.append("                return showLabels ? 5 : 2;")
            html_parts.append("            }).update();")
            html_parts.append("        }")
            html_parts.append("")
            html_parts.append("        // Save tree as SVG")
            html_parts.append("        function saveSvg() {")
            html_parts.append("            let svg = document.querySelector(\"#tree_container svg\");")
            html_parts.append("            let serializer = new XMLSerializer();")
            html_parts.append("            let source = serializer.serializeToString(svg);")
            html_parts.append("")
            html_parts.append("            // Add name spaces")
            html_parts.append("            if (!source.match(/^<svg[^>]+xmlns=\"http:\\/\\/www\\.w3\\.org\\/2000\\/svg\"/)) {")
            html_parts.append("                source = source.replace(/^<svg/, '<svg xmlns=\"http://www.w3.org/2000/svg\"');")
            html_parts.append("            }")
            html_parts.append("            if (!source.match(/^<svg[^>]+\"http:\\/\\/www\\.w3\\.org\\/1999\\/xlink\"/)) {")
            html_parts.append("                source = source.replace(/^<svg/, '<svg xmlns:xlink=\"http://www.w3.org/1999/xlink\"');")
            html_parts.append("            }")
            html_parts.append("")
            html_parts.append("            // Add XML declaration")
            html_parts.append("            source = '<?xml version=\"1.0\" standalone=\"no\"?>\\r\\n' + source;")
            html_parts.append("")
            html_parts.append("            // Create download link")
            html_parts.append("            let downloadLink = document.createElement(\"a\");")
            html_parts.append("            downloadLink.href = \"data:image/svg+xml;charset=utf-8,\" + encodeURIComponent(source);")
            html_parts.append(f"            downloadLink.download = \"tree_{clade_id}.svg\";")
            html_parts.append("            document.body.appendChild(downloadLink);")
            html_parts.append("            downloadLink.click();")
            html_parts.append("            document.body.removeChild(downloadLink);")
            html_parts.append("        }")
            html_parts.append("    </script>")
            html_parts.append("</body>")
            html_parts.append("</html>")

            # Join all parts into the complete HTML
            html_content = "\n".join(html_parts)

            # Write the HTML file
            with open(html_path, 'w') as f:
                f.write(html_content)

            logger.info(f"Created interactive tree visualization for {clade_id}: {html_path}")
            return html_path

        except Exception as e:
            logger.error(f"Failed to create interactive tree visualization for {clade_id}: {e}")
            if self.debug:
                import traceback
                logger.debug(f"Traceback: {traceback.format_exc()}")
            return None

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
    print(f"  Annotated trees:    {args_ns.tree}_au.nwk, {args_ns.tree}_lnl.nwk, {args_ns.tree}_combined.nwk")
    print(f"  Detailed report:    {output_p.with_suffix('.md')}")
    if args_ns.temp: print(f"  Temp directory:     {args_ns.temp}")
    if args_ns.debug: print(f"  Debug mode:         Enabled (log: mldecay_debug.log, if configured)")
    if args_ns.keep_files: print(f"  Keep temp files:    Enabled")
    if args_ns.visualize:
        print("\nVISUALIZATIONS:")
        print(f"  Enabled, format:    {args_ns.viz_format}")
        print(f"  Tree plot:          {output_p.parent / (output_p.stem + '_tree.' + args_ns.viz_format)}")
        if args_ns.html_trees:
            print(f"  HTML trees:         Enabled (using {'CDN' if args_ns.js_cdn else 'embedded'} JavaScript)")
        else:
            print(f"  HTML trees:         Disabled")
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

# Output file prefix (default: ml_decay_indices)
# This will generate files like: myanalysis.txt, myanalysis_tree.nwk, etc.
output = ml_decay_indices.txt

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

# Generate interactive HTML trees (default: true)
html_trees = true

# Use CDN for JavaScript libraries (default: true)
js_cdn = true

# Keep tree files for HTML visualization (default: false)
keep_tree_files = false

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
                'html_trees': 'html_trees',
                'js_cdn': 'js_cdn',
                'keep_tree_files': 'keep_tree_files',
                'constraint_mode': 'constraint_mode',
                'test_branches': 'test_branches',
                'constraint_file': 'constraint_file',
                'check_convergence': 'check_convergence',
                'min_ess': 'min_ess',
                'max_psrf': 'max_psrf',
                'max_asdsf': 'max_asdsf',
                'convergence_strict': 'convergence_strict',
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
                                     'visualize', 'html_trees', 'js_cdn', 'keep_tree_files',
                                     'check_convergence', 'convergence_strict']:
                        value = str_to_bool(value)
                    elif arg_param in ['gamma_shape', 'prop_invar', 'bayes_burnin', 'ss_alpha',
                                       'max_psrf', 'max_asdsf']:
                        value = float(value)
                    elif arg_param in ['nst', 'bootstrap_reps', 'bayes_ngen', 'bayes_chains',
                                       'bayes_sample_freq', 'ss_nsteps', 'mpi_processors', 'min_ess']:
                        value = int(value)
                    
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
    parser.add_argument("--output", default="ml_decay_indices.txt", help="Output file for summary results.")
    parser.add_argument("--tree", default="annotated_tree", help="Base name for annotated tree files. Three trees will be generated with suffixes: _au.nwk (AU p-values), _lnl.nwk (likelihood differences), and _combined.nwk (both values).")
    parser.add_argument("--site-analysis", action="store_true", help="Perform site-specific likelihood analysis to identify supporting/conflicting sites for each branch.")
    parser.add_argument("--data-type", default="dna", choices=["dna", "protein", "discrete"], help="Type of sequence data.")
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
                                 help="Check MCMC convergence diagnostics (default: True)")
    convergence_opts.add_argument("--min-ess", type=int, default=200,
                                 help="Minimum ESS (Effective Sample Size) threshold (default: 200)")
    convergence_opts.add_argument("--max-psrf", type=float, default=1.01,
                                 help="Maximum PSRF (Potential Scale Reduction Factor) threshold (default: 1.01)")
    convergence_opts.add_argument("--max-asdsf", type=float, default=0.01,
                                 help="Maximum ASDSF (Average Standard Deviation of Split Frequencies) threshold (default: 0.01)")
    convergence_opts.add_argument("--convergence-strict", action="store_true",
                                 help="Fail analysis if convergence criteria not met (default: warn only)")

    viz_opts = parser.add_argument_group('Visualization Output (optional)')
    viz_opts.add_argument("--visualize", action="store_true", help="Generate static visualization plots (requires matplotlib, seaborn).")
    viz_opts.add_argument("--viz-format", default="png", choices=["png", "pdf", "svg"], help="Format for static visualizations.")
    viz_opts.add_argument("--annotation", default="lnl", choices=["au", "lnl"], help="Type of support values to visualize in distribution plots (au=AU p-values, lnl=likelihood differences).")
    viz_opts.add_argument("--html-trees", action=argparse.BooleanOptionalAction, default=True, help="Generate interactive HTML tree visualizations (default: True)")
    viz_opts.add_argument("--js-cdn", action="store_true", default=True, help="Use CDN for JavaScript libraries (faster but requires internet connection)")
    viz_opts.add_argument("--keep-tree-files", action="store_true", default=False, help="Keep Newick tree files used for HTML visualization (default: False)")
    
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
    
    # Validate that alignment is provided (either via command line or config)
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
        # Convert string paths from args to Path objects for panDecayIndices
        temp_dir_path = Path(args.temp) if args.temp else None
        starting_tree_path = Path(args.starting_tree) if args.starting_tree else None

        decay_calc = panDecayIndices(
            alignment_file=args.alignment, # Converted to Path in __init__
            alignment_format=args.format,
            model=effective_model_str, # Pass the constructed string
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
            parsmodel=args.parsmodel, # Pass the BooleanOptionalAction value
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
            convergence_strict=args.convergence_strict
        )

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
                        logger.info(f"Site-specific analysis results written to {site_output_dir}")
                        break  # Only need to do this once if any site_data exists

            output_main_path = Path(args.output)
            decay_calc.write_results(output_main_path)

            report_path = output_main_path.with_suffix(".md")
            decay_calc.generate_detailed_report(report_path)

            output_dir = output_main_path.resolve().parent
            tree_base_name = args.tree  # Use the tree argument directly as the base name
            tree_files = decay_calc.annotate_trees(output_dir, tree_base_name)

            if tree_files:
                logger.info(f"Successfully created {len(tree_files)} annotated trees.")
                for tree_type, path in tree_files.items():
                    logger.info(f"  - {tree_type} tree: {path}")
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
                        decay_calc.generate_html = args.html_trees
                        decay_calc.js_cdn = args.js_cdn
                        decay_calc.viz_format = args.viz_format

                    site_output_dir = output_main_path.parent / f"{output_main_path.stem}_site_analysis"
                    decay_calc.write_site_analysis_results(site_output_dir, keep_tree_files=args.keep_tree_files)
                    logger.info(f"Site-specific analysis results written to {site_output_dir}")

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
