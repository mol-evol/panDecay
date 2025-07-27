#!/usr/bin/env python3
"""
Main entry point for panDecay modular system.

This module provides the same command-line interface as the original
monolithic version while using the modular components.
"""

import sys
import argparse
import logging
from pathlib import Path

from .core.constants import (
    VERSION,
    DEFAULT_ALIGNMENT_FORMAT,
    DEFAULT_MODEL,
    DEFAULT_PAUP_PATH,
    DEFAULT_OUTPUT_FILE,
    DEFAULT_TREE_BASE,
    DEFAULT_DATA_TYPE,
    DEFAULT_THREADS,
    DEFAULT_ANALYSIS_MODE,
    DEFAULT_BOOTSTRAP_REPS,
    DEFAULT_MRBAYES_PATH,
    DEFAULT_BAYES_NGEN,
    DEFAULT_BAYES_BURNIN,
    DEFAULT_BAYES_CHAINS,
    DEFAULT_BAYES_SAMPLE_FREQ,
    DEFAULT_MARGINAL_LIKELIHOOD,
    DEFAULT_SS_ALPHA,
    DEFAULT_SS_NSTEPS,
    DEFAULT_MPIRUN_PATH,
    DEFAULT_BEAGLE_DEVICE,
    DEFAULT_BEAGLE_PRECISION,
    DEFAULT_BEAGLE_SCALING,
    DEFAULT_MIN_ESS,
    DEFAULT_MAX_PSRF,
    DEFAULT_MAX_ASDSF,
    DEFAULT_MRBAYES_PARSE_TIMEOUT,
    DEFAULT_VIZ_FORMAT,
    DEFAULT_ANNOTATION,
    DEFAULT_OUTPUT_STYLE,
    DEFAULT_CONSTRAINT_MODE
)
from .core.configuration import generate_config_template, parse_config
from .core.analysis_engine import panDecayIndices


def get_display_path(path):
    """Get a display-friendly path representation."""
    if isinstance(path, Path):
        path_obj = path
    else:
        path_obj = Path(path)
    
    if path_obj.is_absolute():
        try:
            # Try to get relative path from current directory
            rel_path = path_obj.relative_to(Path.cwd())
            # Use relative path if it's shorter and doesn't go up too many levels
            if len(str(rel_path)) < len(str(path_obj)) and not str(rel_path).startswith('../../../'):
                return str(rel_path)
        except ValueError:
            pass  # Path is not relative to current directory
    
    return str(path_obj)


def print_runtime_parameters(args_ns, model_str_for_print):
    """Print runtime parameters in a formatted box."""
    print("┌─" + "─" * 76 + "─┐")
    print("│" + f" panDecay v{VERSION} - Runtime Parameters".center(76) + " │")
    print("├─" + "─" * 76 + "─┤")
    print(f"│ {'Alignment:':<20} {get_display_path(args_ns.alignment):<54} │")
    print(f"│ {'Format:':<20} {args_ns.format:<54} │")
    print(f"│ {'Data type:':<20} {args_ns.data_type:<54} │")
    print(f"│ {'Model:':<20} {model_str_for_print:<54} │")
    print(f"│ {'Analysis mode:':<20} {args_ns.analysis:<54} │")
    print(f"│ {'Threads:':<20} {str(args_ns.threads):<54} │")
    print(f"│ {'Output file:':<20} {get_display_path(args_ns.output):<54} │")
    
    if args_ns.starting_tree:
        print(f"│ {'Starting tree:':<20} {get_display_path(args_ns.starting_tree):<54} │")
    
    if args_ns.site_analysis:
        print(f"│ {'Site analysis:':<20} {'Enabled':<54} │")
    
    if args_ns.bootstrap:
        print(f"│ {'Bootstrap:':<20} {f'{args_ns.bootstrap_reps} replicates':<54} │")
    
    if args_ns.visualize:
        print(f"│ {'Visualization:':<20} {f'{args_ns.viz_format.upper()} format':<54} │")
    
    print("└─" + "─" * 76 + "─┘")


def main():
    """Main entry point for panDecay."""
    parser = argparse.ArgumentParser(
        description=f"panDecay v{VERSION}: Calculate phylogenetic decay indices (ML, Bayesian, and parsimony).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Basic arguments
    parser.add_argument("alignment", nargs='?', help="Input alignment file path (can be specified in config file).")
    parser.add_argument("--format", default=DEFAULT_ALIGNMENT_FORMAT, help="Alignment format.")
    parser.add_argument("--model", default=DEFAULT_MODEL, help="Base substitution model (e.g., GTR, HKY, JC). Combine with --gamma and --invariable.")
    parser.add_argument("--gamma", action="store_true", help="Add Gamma rate heterogeneity (+G) to model.")
    parser.add_argument("--invariable", action="store_true", help="Add invariable sites (+I) to model.")
    
    parser.add_argument("--paup", default=DEFAULT_PAUP_PATH, help="Path to PAUP* executable.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_FILE, help="Output file for summary results.")
    parser.add_argument("--tree", default=DEFAULT_TREE_BASE, help="Base name for annotated tree files. Three trees will be generated with suffixes: _au.nwk (AU p-values), _delta_lnl.nwk (likelihood differences), and _combined.nwk (both values).")
    parser.add_argument("--site-analysis", action="store_true", help="Perform site-specific likelihood analysis to identify supporting/conflicting sites for each branch.")
    parser.add_argument("--data-type", default=DEFAULT_DATA_TYPE, choices=["dna", "protein", "discrete"], help="Type of sequence data.")
    
    # Model parameter overrides
    mparams = parser.add_argument_group('Model Parameter Overrides (optional)')
    mparams.add_argument("--gamma-shape", type=float, help="Fixed Gamma shape value (default: estimate if +G).")
    mparams.add_argument("--prop-invar", type=float, help="Fixed proportion of invariable sites (default: estimate if +I).")
    mparams.add_argument("--base-freq", choices=["equal", "estimate", "empirical"], help="Base/state frequencies (default: model-dependent). 'empirical' uses observed frequencies.")
    mparams.add_argument("--rates", choices=["equal", "gamma"], help="Site rate variation model (overrides --gamma flag if specified).")
    mparams.add_argument("--protein-model", help="Specific protein model (e.g., JTT, WAG; overrides base --model for protein data).")
    mparams.add_argument("--nst", type=int, choices=[1, 2, 6], help="Number of substitution types (DNA; overrides model-based nst).")
    mparams.add_argument("--parsmodel", action=argparse.BooleanOptionalAction, default=None, help="Use parsimony-based branch lengths (discrete data; default: yes for discrete). Use --no-parsmodel to disable.")
    
    # Runtime control
    run_ctrl = parser.add_argument_group('Runtime Control')
    run_ctrl.add_argument("--threads", default=DEFAULT_THREADS, help="Number of threads for PAUP* (e.g., 4 or 'auto').")
    run_ctrl.add_argument("--starting-tree", help="Path to a user-provided starting tree file (Newick).")
    run_ctrl.add_argument("--paup-block", help="Path to file with custom PAUP* commands for model/search setup (overrides most model args).")
    run_ctrl.add_argument("--temp", help="Custom directory for temporary files (default: system temp).")
    run_ctrl.add_argument("--keep-files", action="store_true", help="Keep temporary files after analysis.")
    run_ctrl.add_argument("--debug", action="store_true", help="Enable detailed debug logging (implies --keep-files).")
    
    # Analysis mode selection
    analysis_mode = parser.add_argument_group('Analysis Mode')
    analysis_mode.add_argument("--analysis", 
                              choices=["ml", "bayesian", "parsimony", "ml+parsimony", "bayesian+parsimony", "ml+bayesian", "all"], 
                              default=DEFAULT_ANALYSIS_MODE,
                              help="Type of decay analysis to perform (default: ml). "
                                   "Options: ml, bayesian, parsimony, ml+parsimony, bayesian+parsimony, ml+bayesian, all")
    
    # Bootstrap options
    bootstrap_opts = parser.add_argument_group('Bootstrap Analysis (optional)')
    bootstrap_opts.add_argument("--bootstrap", action="store_true", help="Perform bootstrap analysis to calculate support values.")
    bootstrap_opts.add_argument("--bootstrap-reps", type=int, default=DEFAULT_BOOTSTRAP_REPS, help=f"Number of bootstrap replicates (default: {DEFAULT_BOOTSTRAP_REPS})")
    
    # Bayesian-specific options
    bayesian_opts = parser.add_argument_group('Bayesian Analysis Options')
    bayesian_opts.add_argument("--bayesian-software", choices=["mrbayes"], 
                              default="mrbayes", help="Bayesian software to use (default: mrbayes)")
    bayesian_opts.add_argument("--mrbayes-path", default=DEFAULT_MRBAYES_PATH, help="Path to MrBayes executable")
    bayesian_opts.add_argument("--bayes-model", help="Model for Bayesian analysis (if different from ML model)")
    bayesian_opts.add_argument("--bayes-ngen", type=int, default=DEFAULT_BAYES_NGEN, help="Number of MCMC generations")
    bayesian_opts.add_argument("--bayes-burnin", type=float, default=DEFAULT_BAYES_BURNIN, help="Burnin fraction (0-1)")
    bayesian_opts.add_argument("--bayes-chains", type=int, default=DEFAULT_BAYES_CHAINS, help="Number of MCMC chains")
    bayesian_opts.add_argument("--bayes-sample-freq", type=int, default=DEFAULT_BAYES_SAMPLE_FREQ, help="Sample frequency for MCMC")
    bayesian_opts.add_argument("--marginal-likelihood", choices=["ss", "ps", "hm"], default=DEFAULT_MARGINAL_LIKELIHOOD,
                              help="Marginal likelihood estimation method: ss=stepping-stone, ps=path sampling, hm=harmonic mean")
    bayesian_opts.add_argument("--ss-alpha", type=float, default=DEFAULT_SS_ALPHA, help="Alpha parameter for stepping-stone sampling")
    bayesian_opts.add_argument("--ss-nsteps", type=int, default=DEFAULT_SS_NSTEPS, help="Number of steps for stepping-stone sampling")
    
    # MPI and BEAGLE options
    parallel_opts = parser.add_argument_group('Parallel Processing Options (MrBayes)')
    parallel_opts.add_argument("--use-mpi", action="store_true", help="Use MPI version of MrBayes")
    parallel_opts.add_argument("--mpi-processors", type=int, help="Number of MPI processors (default: number of chains)")
    parallel_opts.add_argument("--mpirun-path", default=DEFAULT_MPIRUN_PATH, help="Path to mpirun executable")
    parallel_opts.add_argument("--use-beagle", action="store_true", help="Enable BEAGLE library for GPU/CPU acceleration")
    parallel_opts.add_argument("--beagle-device", choices=["cpu", "gpu", "auto"], default=DEFAULT_BEAGLE_DEVICE, 
                             help="BEAGLE device preference")
    parallel_opts.add_argument("--beagle-precision", choices=["single", "double"], default=DEFAULT_BEAGLE_PRECISION,
                             help="BEAGLE precision mode")
    parallel_opts.add_argument("--beagle-scaling", choices=["none", "dynamic", "always"], default=DEFAULT_BEAGLE_SCALING,
                             help="BEAGLE scaling frequency")
    
    # Convergence checking options
    convergence_opts = parser.add_argument_group('Convergence Checking Options (MrBayes)')
    convergence_opts.add_argument("--check-convergence", action=argparse.BooleanOptionalAction, default=True,
                                 help="Check MCMC convergence diagnostics")
    convergence_opts.add_argument("--min-ess", type=int, default=DEFAULT_MIN_ESS,
                                 help=f"Minimum ESS (Effective Sample Size) threshold (default: {DEFAULT_MIN_ESS})")
    convergence_opts.add_argument("--max-psrf", type=float, default=DEFAULT_MAX_PSRF,
                                 help=f"Maximum PSRF (Potential Scale Reduction Factor) threshold (default: {DEFAULT_MAX_PSRF})")
    convergence_opts.add_argument("--max-asdsf", type=float, default=DEFAULT_MAX_ASDSF,
                                 help=f"Maximum ASDSF (Average Standard Deviation of Split Frequencies) threshold (default: {DEFAULT_MAX_ASDSF})")
    convergence_opts.add_argument("--convergence-strict", action="store_true",
                                 help="Fail analysis if convergence criteria not met (default: warn only)")
    convergence_opts.add_argument("--mrbayes-parse-timeout", type=float, default=DEFAULT_MRBAYES_PARSE_TIMEOUT,
                                 help=f"Timeout for parsing MrBayes consensus trees in seconds (0 to disable, default: {DEFAULT_MRBAYES_PARSE_TIMEOUT})")

    # Visualization options
    viz_opts = parser.add_argument_group('Visualization Output (optional)')
    viz_opts.add_argument("--visualize", action="store_true", help="Generate static visualization plots (requires matplotlib, seaborn).")
    viz_opts.add_argument("--viz-format", choices=["png", "pdf", "svg"], default=DEFAULT_VIZ_FORMAT, help="Format for visualization output.")
    viz_opts.add_argument("--annotation", choices=["au", "lnl"], default=DEFAULT_ANNOTATION, help="Annotation type for visualization: AU p-values or log-likelihood differences.")
    viz_opts.add_argument("--output-style", choices=["unicode", "ascii"], default=DEFAULT_OUTPUT_STYLE, help="Output formatting style for tables and progress.")
    
    # Constraint options
    constraint_opts = parser.add_argument_group('Constraint Testing Options')
    constraint_opts.add_argument("--constraint-mode", choices=["all", "specific", "exclude"], default=DEFAULT_CONSTRAINT_MODE,
                                help="Constraint testing mode: 'all' tests all branches, 'specific' tests only listed branches, 'exclude' tests all except listed.")
    constraint_opts.add_argument("--test-branches", help="Semicolon-separated clades to test (format: 'taxon1,taxon2;taxon3,taxon4').")
    constraint_opts.add_argument("--constraint-file", help="File containing constraint definitions (one per line).")
    
    # Configuration file option
    parser.add_argument("--config", help="Configuration file with analysis parameters (overrides command-line options).")
    parser.add_argument("--generate-config", help="Generate a template configuration file at the specified path and exit.")
    
    args = parser.parse_args()
    
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    logger = logging.getLogger(__name__)
    
    # Handle config generation
    if args.generate_config:
        generate_config_template(args.generate_config)
        return
    
    # Load configuration file if provided
    if args.config:
        args = parse_config(args.config, args)
    
    # Validate required arguments
    if not args.alignment:
        logger.error("Error: Alignment file is required (either as positional argument or in config file)")
        parser.print_help()
        sys.exit(1)
    
    # Convert data_type to lowercase
    args.data_type = args.data_type.lower()
    
    # Validate Bayesian analysis arguments
    if args.analysis == "ml" and args.bayesian_software:
        logger.warning("Bayesian software specified but analysis mode is ML-only. Bayesian options will be ignored.")
    
    # Set up debug logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
        debug_log_path = Path.cwd() / "mldecay_debug.log"
        fh = logging.FileHandler(debug_log_path, mode='w')
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        logger.info(f"Debug logging enabled. Detailed log: {debug_log_path}")
        args.keep_files = True
    
    # Construct full model string
    effective_model_str = args.model
    if args.gamma: 
        effective_model_str += "+G"
    if args.invariable: 
        effective_model_str += "+I"
    
    # Handle protein/discrete model adjustments
    if args.data_type == "protein" and not args.protein_model and not any(pm in args.model.upper() for pm in ["JTT", "WAG", "LG", "DAYHOFF"]):
        logger.info(f"Protein data with generic model '{args.model}'. Effective model might default to JTT within PAUP* settings.")
    elif args.data_type == "discrete" and "MK" not in args.model.upper():
        logger.info(f"Discrete data with non-Mk model '{args.model}'. Effective model might default to Mk within PAUP* settings.")
    
    # Read PAUP block if specified
    paup_block_content = None
    if args.paup_block:
        from .core.configuration import read_paup_block
        pbf_path = Path(args.paup_block)
        logger.info(f"Reading PAUP block from: {pbf_path}")
        paup_block_content = read_paup_block(pbf_path)
        if paup_block_content is None:
            logger.error("Failed to read or validate PAUP block file. Exiting.")
            sys.exit(1)
    
    # Print runtime parameters
    print_runtime_parameters(args, effective_model_str)
    
    try:
        # Convert string paths to Path objects
        temp_dir_path = Path(args.temp) if args.temp else None
        starting_tree_path = Path(args.starting_tree) if args.starting_tree else None
        
        # Create panDecayIndices instance
        decay_calc = panDecayIndices(
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
            gamma_shape=args.gamma_shape, 
            prop_invar=args.prop_invar,
            base_freq=args.base_freq, 
            rates=args.rates,
            protein_model=args.protein_model, 
            nst=args.nst,
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
            output_style=args.output_style
        )
        
        # Build ML tree
        decay_calc.build_ml_tree()
        
        if decay_calc.ml_tree and decay_calc.ml_likelihood is not None:
            # Run bootstrap analysis if requested
            if args.bootstrap:
                logger.info(f"Running bootstrap analysis with {args.bootstrap_reps} replicates...")
                decay_calc.run_bootstrap_analysis(num_replicates=args.bootstrap_reps)
            
            # Calculate decay indices
            decay_calc.calculate_decay_indices(perform_site_analysis=args.site_analysis)
            
            # Write site analysis results if available
            if hasattr(decay_calc, 'decay_indices') and decay_calc.decay_indices:
                for clade_id, data in decay_calc.decay_indices.items():
                    if 'site_data' in data:
                        site_output_dir = Path(args.output).parent / f"{Path(args.output).stem}_site_analysis"
                        decay_calc.write_site_analysis_results(site_output_dir)
                        logger.info(f"Site-specific analysis results written to {get_display_path(site_output_dir)}")
                        break
            
            # Write main results
            output_main_path = Path(args.output)
            decay_calc.write_formatted_results(output_main_path)
            
            # Generate detailed report
            report_path = output_main_path.with_suffix(".md")
            decay_calc.generate_detailed_report(report_path)
            
            # Create annotated trees
            output_dir = output_main_path.resolve().parent
            tree_base_name = args.tree
            tree_files = decay_calc.annotate_trees(output_dir, tree_base_name)
            
            if tree_files:
                logger.info(f"Successfully created {len(tree_files)} annotated trees.")
                for tree_type, path in tree_files.items():
                    logger.info(f"  - {tree_type} tree: {get_display_path(path)}")
            else:
                logger.warning("Failed to create annotated trees.")
            
            # Handle visualization
            if args.visualize:
                viz_out_dir = output_main_path.resolve().parent
                viz_base_name = output_main_path.stem
                viz_kwargs = {'width': 10, 'height': 6, 'format': args.viz_format}
                
                # Check for visualization libraries
                try: 
                    import matplotlib, seaborn
                except ImportError:
                    logger.warning("Matplotlib/Seaborn not installed. Skipping static visualizations.")
                    args.visualize = False
                
                if args.visualize:
                    decay_calc.visualize_support_distribution(
                        viz_out_dir / f"{viz_base_name}_dist_{args.annotation}.{args.viz_format}",
                        value_type=args.annotation, **viz_kwargs)
                
                if args.site_analysis:
                    if args.visualize:
                        decay_calc.viz_format = args.viz_format
                    
                    site_output_dir = output_main_path.parent / f"{output_main_path.stem}_site_analysis"
                    decay_calc.write_site_analysis_results(site_output_dir)
                    logger.info(f"Site-specific analysis results written to {get_display_path(site_output_dir)}")
            
            # Cleanup
            decay_calc.cleanup_intermediate_files()
            logger.info("panDecay analysis completed successfully.")
        else:
            logger.error("ML tree construction failed or likelihood missing. Halting.")
            sys.exit(1)
    
    except Exception as e:
        logger.error(f"panDecay analysis terminated with an error: {e}")
        if args.debug:
            import traceback
            logger.debug("Full traceback:\n%s", traceback.format_exc())
        sys.exit(1)
    
    finally:
        if 'decay_calc' in locals() and (args.debug or args.keep_files):
            logger.info(f"Temporary files are preserved in: {decay_calc.temp_path}")


if __name__ == "__main__":
    main()