![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg) ![Version](https://img.shields.io/badge/version-1.1.0-orange.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# panDecay: Phylogenetic Analysis using Decay Indices

panDecay is a Python command-line tool for calculating phylogenetic decay indices across multiple analysis frameworks. It can compute parsimony-based decay indices (traditional Bremer support), Maximum Likelihood (ML)-based decay indices, and Bayesian decay indices using MrBayes   

## Table of Contents

1.  [Background](#background)
    *   [What are Decay Indices / Bremer Support?](#what-are-decay-indices--bremer-support)
    *   [Why ML-based Decay Indices?](#why-ml-based-decay-indices)
    *   [Bayesian Decay Indices](#bayesian-decay-indices)
2.  [Features](#features)
3.  [Installation](#installation)
    *   [Dependencies](#dependencies)
    *   [Installing panDecay](#installing-pandecay)
4.  [Usage](#usage)
    *   [Basic Command](#basic-command)
    *   [Command-Line Arguments](#command-line-arguments)
5.  [Input Files](#input-files)
    *   [Sequence Alignment](#sequence-alignment)
    *   [Optional Starting Tree](#optional-starting-tree)
    *   [Optional PAUP\* Block File](#optional-paup-block-file)
6.  [Output Files](#output-files)
    *   [Main Results File (`pan_decay_indices.txt`)](#main-results-file-pan_decay_indicestxt)
    *   [Annotated Trees](#annotated-trees)
    *   [Detailed Markdown Report (`_detailed.md`)](#detailed-markdown-report-_detailedmd)
    *   [Site-Specific Analysis (Optional)](#site-specific-analysis-optional)
    *   [Visualizations (Optional)](#visualizations-optional)
    *   [Temporary Files (Debug/Keep)](#temporary-files-debugkeep)
7.  [Examples & Recipes](#examples--recipes)
    *   [Example 1: Basic DNA Analysis](#example-1-basic-dna-analysis)
    *   [Example 2: Protein Data with Specific Model](#example-2-protein-data-with-specific-model)
    *   [Example 3: Discrete Morphological Data](#example-3-discrete-morphological-data)
    *   [Example 4: Using a Starting Tree](#example-4-using-a-starting-tree)
    *   [Example 5: Advanced Control with PAUP\* Block](#example-5-advanced-control-with-paup-block)
    *   [Example 6: Site-Specific Analysis](#example-6-site-specific-analysis)
    *   [Example 7: Bootstrap Analysis](#example-7-bootstrap-analysis)
8.  [Interpreting Results](#interpreting-results)
9.  [Troubleshooting](#troubleshooting)
10. [Citations](#citations)
11. [License](#license)
12. [Contributing](#contributing)
13. [Contact](#contact)

## Background

### What are Decay Indices / Bremer Support?

In phylogenetics, assessing the support for individual branches (clades) in a tree is crucial. Traditional bootstrap methods resample characters to estimate support. Decay indices, originally developed for parsimony (Bremer, 1988; Bremer, 1994), measure how much worse a tree must be (e.g., how many extra steps in parsimony) to lose a particular clade. A higher decay value indicates stronger support for that clade.

### Why ML-based Decay Indices?

While parsimony decay indices are well-established, maximum likelihood (ML) is a statistically robust framework for phylogenetic inference. ML-based decay indices extend this concept to the likelihood framework. Instead of "extra steps," we look at the difference in log-likelihood scores between the optimal ML tree and the best tree where a specific clade is constrained to be non-monophyletic (i.e., the branch defining that clade is collapsed).

panDecay automates this process by:
1.  Finding the optimal ML tree and its likelihood score.
2.  For each internal branch in the ML tree:
    a.  Defining a constraint that forces the taxa in that clade to *not* form a monophyletic group (using PAUP\*'s `converse=yes` constraint).
    b.  Searching for the best ML tree under this reverse-constraint and recording its likelihood.
3.  Calculating the difference in log-likelihood between the unconstrained ML tree and each constrained tree.
4.  Performing an Approximately Unbiased (AU) test (Shimodaira, 2002) to statistically compare the unconstrained ML tree against all the constrained alternative trees. The p-value from the AU test indicates the significance of the support for the original clade.

A significantly worse likelihood for the constrained tree (and a low AU test p-value for that constrained tree) provides strong evidence in favor of the monophyly of the original clade.

### Bayesian Decay Indices

panDecay now supports Bayesian phylogenetic decay indices, extending the decay index concept to Bayesian inference. Instead of comparing log-likelihoods, Bayesian decay indices compare marginal likelihoods between:
1. An unconstrained Bayesian analysis where all topologies are explored
2. Constrained analyses where specific clades are forced to be non-monophyletic

The Bayesian decay index for a clade is calculated as:
- **Bayesian Decay = ln(ML_unconstrained) - ln(ML_constrained)**
- **Bayes Factor = exp(Bayesian Decay)**

Where ML represents the marginal likelihood (not to be confused with maximum likelihood). A positive Bayesian decay value indicates support for the clade, with larger values indicating stronger support.

**Important Note on Interpretation**: In phylogenetic applications, Bayesian decay values tend to closely approximate ML log-likelihood differences. This occurs because:
- The compared models differ only in topological constraints, not in substitution models or parameters
- When data strongly support a topology, the marginal likelihood is dominated by the likelihood component
- Traditional Bayes Factor interpretation scales (e.g., BF >10 = strong, >100 = decisive) were developed for comparing fundamentally different models and do not apply well to phylogenetic topology testing

**Interpreting Bayesian decay values**:
- BD values should be interpreted comparatively across branches rather than using absolute thresholds
- Compare BD values with other support metrics (ΔlnL, AU test p-values, bootstrap, parsimony decay)
- BD values may scale with alignment size and sequence diversity
- Strong support is best identified when multiple metrics concordantly indicate robust clade support

Note that BD values of 30-50 or higher are common when data strongly support a clade and should not be considered anomalous.

panDecay can perform Bayesian analyses using:
- **MrBayes**: Currently supported with stepping-stone sampling (default) or harmonic mean marginal likelihood estimation

## Features

### Analysis Types
*   **Parsimony Analysis**: Calculates traditional Bremer support values (parsimony decay indices)
*   **ML Analysis**: Calculates ML-based decay values using log-likelihood differences
*   **Bayesian Analysis**: Calculates Bayesian decay indices using marginal likelihood comparisons
*   **Combined Analysis**: Performs multiple analysis types in a single run (e.g., ML+Bayesian, or all three)
*   **Bootstrap Analysis**: Optional bootstrap support values alongside decay indices

### Core Capabilities
*   Performs the Approximately Unbiased (AU) test for statistical assessment of ML branch support
*   Supports DNA, Protein, and binary discrete morphological data
*   Optional site-specific likelihood analysis to identify which alignment positions support or conflict with each branch
*   Flexible model specification (e.g., GTR, HKY, JTT, WAG, Mk) with options for gamma-distributed rate heterogeneity (+G) and proportion of invariable sites (+I)
*   Allows fine-grained control over model parameters (gamma shape, pinvar, base frequencies, etc.)
*   Option to provide a custom PAUP\* block for complex model or search strategy definitions
*   Option to provide a starting tree for the initial ML search

### Bayesian Features
*   Support for MrBayes with automatic constraint generation
*   Marginal likelihood estimation using:
    *   **Stepping-stone sampling** (recommended, more accurate)
    *   **Harmonic mean** (faster but less reliable)
*   **Improved Bayes Factor reporting**:
    *   Primary focus on Bayes Decay (log Bayes Factor) for interpretability
    *   Bayes Factor values capped at 10^6 for display to avoid numerical issues
    *   Clear warnings about model dimension effects for extreme values
*   **MCMC Convergence Diagnostics**:
    *   Automatic checking of ESS (Effective Sample Size)
    *   PSRF (Potential Scale Reduction Factor) monitoring
    *   ASDSF (Average Standard Deviation of Split Frequencies) tracking
    *   Configurable convergence thresholds with strict mode option
*   Flexible MCMC parameters (generations, chains, burnin, sampling frequency)
*   **MPI support**: Run chains in parallel with MPI-enabled MrBayes
*   **BEAGLE support**: GPU/CPU acceleration for likelihood calculations

### Output Files
*   Tab-delimited results file with ML and/or Bayesian support values
*   Multiple Newick trees annotated with different support values
*   Detailed Markdown report summarizing the analysis and results
*   Comprehensive trees combining all support metrics when multiple analyses are performed
*   Optional static visualization (requires `matplotlib` and `seaborn`):
    *   Distribution of support values
    *   Site-specific support visualizations

### Technical Features
*   Multi-threaded PAUP\* execution (configurable)
*   Parallel execution of constraint analyses
*   Debug mode and option to keep temporary files
*   Robust error handling and recovery

### User Interface Features
*   **Professional Runtime Banner**: Clean, professional parameter display at startup
    *   Fixed 80-character width for consistent appearance
    *   Organized sections (Analysis Configuration, Runtime Settings, Bayesian Parameters, etc.)
    *   Includes software citation: "McInerney, J.O. (2025) panDecay: Phylogenetic Analysis with Decay Indices"
    *   No jagged edges or complex table formatting
*   **Enhanced Progress Display**: Clean, non-verbose progress tracking
    *   Unicode spinner animations with dynamic updates
    *   Consolidated messaging without repetitive output
    *   Progress indicators that overwrite cleanly
*   **Comprehensive Output Reports**: Enhanced markdown reports with complete analysis summaries
    *   Embedded site analysis summary tables
    *   Integrated parsimony analysis statistics  
    *   PNG visualization embedding with relative links
    *   Complete file navigation and cross-references
*   **Clean Console Output**: Minimal, professional terminal output
*   **Informative Error Messages**: Clear error reporting with helpful troubleshooting suggestions

## Installation

### Dependencies

panDecay requires Python 3.8 or higher and has several dependencies that can be easily installed using pip.

1. **PAUP\***: Required for ML analysis. You must have a working PAUP\* command-line executable installed and accessible in your system's PATH, or provide the full path to it. PAUP\* can be obtained from [phylosolutions.com](http://phylosolutions.com/paup-test/).

2. **MrBayes** (optional): Required for Bayesian analysis. Install MrBayes and ensure it's accessible as `mb` in your PATH, or specify the path with `--mrbayes-path`. MrBayes can be obtained from [nbisweden.github.io/MrBayes/](https://nbisweden.github.io/MrBayes/).

3. **Python Dependencies**: All Python dependencies can be installed using the provided `requirements.txt` file:

   ```bash
   pip install -r requirements.txt
   ```

   This will install all required packages including BioPython, NumPy, and the optional visualization packages Matplotlib and Seaborn.

### Installing panDecay

#### Option 1: Install from PyPI (Recommended)
```bash
pip install pandecay
```

#### Option 2: Install from GitHub (Latest Development Version)
```bash
pip install git+https://github.com/mol-evol/panDecay.git
```

#### Option 3: Development Installation
For developers or if you want to modify the code:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/mol-evol/panDecay.git
   cd panDecay
   ```

2. **Install in Development Mode:**
   ```bash
   pip install -e .
   ```

This installs panDecay in "editable" mode, so changes to the source code are immediately available.

## Usage

### Basic Command

```bash
pandecay <alignment_file> --model <model_name> [options...]
```

### Command-Line Arguments

```
usage: pandecay [-h] [--format FORMAT] [--model MODEL] [--gamma] [--invariable] [--paup PAUP] [--output OUTPUT] [--tree TREE]
                  [--data-type {dna,protein,discrete}] [--gamma-shape GAMMA_SHAPE] [--prop-invar PROP_INVAR] 
                  [--base-freq {equal,estimate,empirical}] [--rates {equal,gamma}] [--protein-model PROTEIN_MODEL] 
                  [--nst {1,2,6}] [--parsmodel | --no-parsmodel] [--threads THREADS] [--starting-tree STARTING_TREE] 
                  [--paup-block PAUP_BLOCK] [--temp TEMP] [--keep-files] [--debug] [--site-analysis] 
                  [--analysis {ml,bayesian,parsimony,ml+parsimony,bayesian+parsimony,ml+bayesian,all}] 
                  [--bayesian-software {mrbayes}] [--mrbayes-path MRBAYES_PATH] [--bayes-model BAYES_MODEL]
                  [--bayes-ngen BAYES_NGEN] [--bayes-burnin BAYES_BURNIN] [--bayes-chains BAYES_CHAINS]
                  [--bayes-sample-freq BAYES_SAMPLE_FREQ] [--marginal-likelihood {ss,ps,hm}]
                  [--ss-alpha SS_ALPHA] [--ss-nsteps SS_NSTEPS] [--use-mpi] [--mpi-processors MPI_PROCESSORS]
                  [--mpirun-path MPIRUN_PATH] [--use-beagle] [--beagle-device {auto,cpu,gpu}]
                  [--beagle-precision {single,double}] [--beagle-scaling {none,dynamic,always}]
                  [--check-convergence | --no-check-convergence] [--min-ess MIN_ESS] [--max-psrf MAX_PSRF]
                  [--max-asdsf MAX_ASDSF] [--convergence-strict] [--mrbayes-parse-timeout TIMEOUT]
                  [--bootstrap] [--bootstrap-reps BOOTSTRAP_REPS] [--visualize] [--viz-format {png,pdf,svg}]
                  [--annotation {au,lnl}] [--output-style {unicode,ascii,minimal}]
                  [--config CONFIG] [--generate-config GENERATE_CONFIG]
                  [--constraint-mode {all,specific,exclude}] [--test-branches TEST_BRANCHES] 
                  [--constraint-file CONSTRAINT_FILE] [-v]
                  [alignment]

panDecay v1.1.0: Calculate phylogenetic decay indices (ML, Bayesian, and parsimony).

positional arguments:
  alignment             Input alignment file path (can also be provided via config file).

options:
  -h, --help            show this help message and exit
  --format FORMAT       Alignment format. (default: fasta)
  --model MODEL         Base substitution model (e.g., GTR, HKY, JC). Combine with --gamma and --invariable. (default: GTR)
  --gamma               Add Gamma rate heterogeneity (+G) to model. (default: False)
  --invariable          Add invariable sites (+I) to model. (default: False)
  --paup PAUP           Path to PAUP* executable. (default: paup)
  --output OUTPUT       Output file for summary results. (default: pan_decay_indices.txt)
  --tree TREE           Base name for annotated tree files. Three trees will be generated with suffixes: _au.nwk (AU p-values), 
                        _delta_lnl.nwk (likelihood differences), and _combined.nwk (both values). (default: annotated_tree)
  --data-type {dna,protein,discrete}
                        Type of sequence data. (default: dna)
  --site-analysis       Perform site-specific likelihood analysis to identify supporting/conflicting sites for each branch. (default: False)
  -v, --version         show program's version number and exit

Model Parameter Overrides (optional):
  --gamma-shape GAMMA_SHAPE
                        Fixed Gamma shape value (default: estimate if +G).
  --prop-invar PROP_INVAR
                        Fixed proportion of invariable sites (default: estimate if +I).
  --base-freq {equal,estimate,empirical}
                        Base/state frequencies (default: model-dependent). 'empirical' uses observed frequencies.
  --rates {equal,gamma}
                        Site rate variation model (overrides --gamma flag if specified).
  --protein-model PROTEIN_MODEL
                        Specific protein model (e.g., JTT, WAG; overrides base --model for protein data).
  --nst {1,2,6}         Number of substitution types (DNA; overrides model-based nst).
  --parsmodel, --no-parsmodel
                        Use parsimony-based branch lengths (discrete data; default: yes for discrete). Use --no-parsmodel to disable. (default: None)

Runtime Control:
  --threads THREADS     Number of threads for PAUP* (e.g., 4, 'auto', or 'all'). 'auto' uses: total_cores - 2 (if cores > 2), 
                        total_cores - 1 (if cores > 1), or 1 core. Leaving some cores free is recommended for system stability. (default: auto)
  --starting-tree STARTING_TREE
                        Path to a user-provided starting tree file (Newick).
  --paup-block PAUP_BLOCK
                        Path to file with custom PAUP* commands for model/search setup (overrides most model args).
  --temp TEMP           Custom directory for temporary files (default: system temp).
  --keep-files          Keep temporary files after analysis. (default: False)
  --debug               Enable detailed debug logging (implies --keep-files). (default: False)

Analysis Mode:
  --analysis {ml,bayesian,parsimony,ml+parsimony,bayesian+parsimony,ml+bayesian,all}
                        Type of decay analysis to perform (default: ml). 
                        Options: ml, bayesian, parsimony, ml+parsimony, bayesian+parsimony, ml+bayesian, all

Bayesian Analysis Options:
  --bayesian-software {mrbayes}
                        Bayesian software to use (default: mrbayes)
  --mrbayes-path MRBAYES_PATH
                        Path to MrBayes executable (default: mb)
  --bayes-model BAYES_MODEL
                        Model for Bayesian analysis (if different from ML model)
  --bayes-ngen BAYES_NGEN
                        Number of MCMC generations (default: 1000000)
  --bayes-burnin BAYES_BURNIN
                        Burnin fraction (0-1) (default: 0.25)
  --bayes-chains BAYES_CHAINS
                        Number of MCMC chains (default: 4)
  --bayes-sample-freq BAYES_SAMPLE_FREQ
                        Sample frequency for MCMC (default: 1000)
  --marginal-likelihood {ss,ps,hm}
                        Marginal likelihood estimation method: ss=stepping-stone, ps=path sampling, hm=harmonic mean (default: ss)
  --ss-alpha SS_ALPHA   Alpha parameter for stepping-stone sampling (default: 0.4)
  --ss-nsteps SS_NSTEPS Number of steps for stepping-stone sampling (default: 50)

Parallel Processing Options (MrBayes):
  --use-mpi             Use MPI version of MrBayes for parallel chains (default: False)
  --mpi-processors MPI_PROCESSORS
                        Number of processors for MPI (default: number of chains)
  --mpirun-path MPIRUN_PATH
                        Path to mpirun executable (default: mpirun)
  --use-beagle          Use BEAGLE library for likelihood calculations (default: False)
  --beagle-device {auto,cpu,gpu}
                        BEAGLE device type (default: auto)
  --beagle-precision {single,double}
                        BEAGLE precision level (default: double)
  --beagle-scaling {none,dynamic,always}
                        BEAGLE scaling frequency (default: dynamic)

Convergence Checking Options (MrBayes):
  --check-convergence/--no-check-convergence
                        Check MCMC convergence diagnostics (default: True)
  --min-ess MIN_ESS     Minimum ESS (Effective Sample Size) threshold (default: 200)
  --max-psrf MAX_PSRF   Maximum PSRF (Potential Scale Reduction Factor) threshold (default: 1.01)
  --max-asdsf MAX_ASDSF Maximum ASDSF (Average Standard Deviation of Split Frequencies) threshold (default: 0.01)
  --convergence-strict  Fail analysis if convergence criteria not met (default: warn only)
  --mrbayes-parse-timeout TIMEOUT
                        Timeout for parsing MrBayes output files in seconds (default: 30.0)

Bootstrap Analysis (optional):
  --bootstrap           Perform bootstrap analysis to calculate support values. (default: False)
  --bootstrap-reps BOOTSTRAP_REPS
                        Number of bootstrap replicates (default: 100)

Visualization Output (optional):
  --visualize           Generate static visualization plots (requires matplotlib, seaborn). (default: False)
  --viz-format {png,pdf,svg}
                        Format for static visualizations. (default: png)
  --annotation {au,lnl} Type of support values to visualize in distribution plots (au=AU p-values, lnl=likelihood differences). (default: lnl)
  --output-style {unicode,ascii,minimal}
                        Output formatting style: unicode (modern), ascii (compatible), minimal (basic). (default: unicode)

Configuration and Constraint Options:
  --config CONFIG       Read parameters from configuration file (INI format)
  --generate-config GENERATE_CONFIG
                        Generate a template configuration file at the specified path and exit
  --constraint-mode {all,specific,exclude}
                        Branch selection mode: all (test all branches), specific (test only specified), 
                        exclude (test all except specified) (default: all)
  --test-branches TEST_BRANCHES
                        Specify branches to test. Format: 'taxon1,taxon2,taxon3;taxon4,taxon5' for clades, 
                        '1,3,5' for branch IDs, or '@file.txt' to read from file
  --constraint-file CONSTRAINT_FILE
                        File containing constraint definitions (one per line)
```

## Input Files

### Sequence Alignment
A multiple sequence alignment file.
*   **Formats:** FASTA, NEXUS, PHYLIP, Clustal, etc. (any format BioPython's `AlignIO` can read). Use the `--format` option if not FASTA.
*   **Content:** DNA, protein, or binary (0/1) discrete morphological characters. Use `--data-type` to specify.
    *   For discrete data, characters should be '0' or '1'. Missing data as '?' and gaps as '-' are also recognized.

### Optional Starting Tree
A Newick tree file specified with `--starting-tree <path_to_tree.nwk>`.
*   If provided, PAUP\* will use this tree as the initial tree for the ML search, potentially speeding up the search or helping to find a better likelihood peak. Branch lengths from the starting tree are typically re-optimized.

### Optional PAUP\* Block File
A text file specified with `--paup-block <path_to_block.nex>`.
*   This file should contain valid PAUP\* commands that will be inserted into the PAUP\* script.
*   It typically starts after `execute <alignment_file>;` and should define the model (`lset`), search strategy (`hsearch`), and potentially how trees/scores are saved.
*   This allows advanced users to have full control over PAUP\*'s behavior for model setup and tree searching.
*   **Format:** The content should be the commands that would normally go *between* `BEGIN PAUP;` and `END;` in a PAUP block. For example:
    ```paup
    lset nst=2 basefreq=empirical rates=gamma shape=estimate pinv=estimate;
    hsearch nreps=100 swap=tbr multrees=yes;
    savetrees file=my_custom_ml.tre replace=yes;
    lscores 1 /scorefile=my_custom_scores.txt replace=yes;
    ```
    panDecay will try to defensively add `savetrees` and `lscores` commands if they appear to be missing from the user's block when needed for its internal workflow.

### Configuration File (INI format)
A configuration file specified with `--config <path_to_config.ini>`.
*   **Format:** Standard INI format with key-value pairs and optional sections
*   **Purpose:** Allows specifying all command-line parameters in a file for reproducibility and convenience
*   **Template:** Generate a fully-commented template with `--generate-config template.ini`
*   **Sections:**
    *   Main section (no header): Contains most parameters like alignment, model, output settings
    *   `[constraints]`: Define specific clades to test when using `constraint_mode = specific`
*   **Example:**
    ```ini
    alignment = my_data.fas
    model = GTR
    gamma = true
    constraint_mode = specific
    
    [constraints]
    clade1 = taxonA,taxonB,taxonC
    clade2 = taxonD,taxonE
    ```
*   **Note:** Command-line arguments override configuration file values

### Constraint File
A text file specified with `--constraint-file <path>` or `--test-branches @<path>`.
*   **Format:** One constraint per line, taxa separated by commas
*   **Comments:** Lines starting with # are ignored
*   **Example:**
    ```
    # Primates clade
    Homo_sapiens,Pan_troglodytes,Gorilla_gorilla
    # Rodents clade  
    Mus_musculus,Rattus_norvegicus
    ```

## Output Files

Unless specified with `--output`, `--tree`, etc., output files are created in the current working directory.

### Main Results File (`pan_decay_indices.txt` by default)
A tab-delimited text file containing:
*   The log-likelihood of the best ML tree found (for ML analyses).
*   For each internal branch tested:
    *   `Clade_ID`: An internal identifier for the branch.
    *   `Num_Taxa`: Number of taxa in the clade defined by this branch.
    *   **ML Metrics** (when ML analysis is performed):
        *   `Constrained_lnL`: Log-likelihood of the best tree found when this clade was constrained to be non-monophyletic.
        *   `Delta_LnL`: Log-likelihood difference (ΔlnL) between the constrained tree and the ML tree (constrained_lnL - ML_lnL).
        *   `AU_p-value`: The p-value from the Approximately Unbiased test.
        *   `Significant_AU (p<0.05)`: "Yes" if AU p-value < 0.05, "No" otherwise.
    *   **Bayesian Metrics** (when Bayesian analysis is performed):
        *   `Bayes_ML_Diff`: Marginal likelihood difference (unconstrained - constrained).
        *   `Bayes_Factor`: Exponential of the Bayes_ML_Diff, indicating support strength.
    *   `Bootstrap` (if bootstrap analysis performed): Bootstrap support value for the clade.
    *   `Taxa_List`: A comma-separated list of taxa in the clade.

### Annotated Trees
panDecay generates several different annotated tree files:
* `<tree_base>_au.nwk`: Tree with AU test p-values as branch labels
* `<tree_base>_delta_lnl.nwk`: Tree with log-likelihood differences (ΔlnL) as branch labels
* `<tree_base>_combined.nwk`: Tree with both values as branch labels in the format "AU:0.95|ΔlnL:2.34"

If bootstrap analysis is performed, additional tree files:
* `<tree_base>_bootstrap.nwk`: Tree with bootstrap support values
* `<tree_base>_comprehensive.nwk`: Tree with bootstrap values, AU test p-values, and log-likelihood differences combined in format "BS:80|AU:0.95|ΔlnL:2.34"

These trees can be visualized in standard tree viewers like [FigTree](https://github.com/rambaut/figtree/), [Dendroscope](https://github.com/husonlab/dendroscope3), [iTOL](https://itol.embl.de/), etc. The combined tree is particularly suited for FigTree which handles string labels well.

#### Example Annotated Tree Visualization
The following example shows a comprehensive annotated tree from a combined analysis (ML + Bayesian + Parsimony + Bootstrap):

![Annotated tree example](./images/annotated_tree.png)

In this visualization:
- Branch labels show multiple support metrics in the format: `Clade_X - BS:XX|AU:X.XXXX|ΔlnL:XX.XX|BD:XX.XX|BF:X.XXe+XX|PD:XX`
- **BS**: Bootstrap support percentage (when available)
- **AU**: Approximately Unbiased test p-value (lower values = stronger support)
- **ΔlnL**: Log-likelihood difference (higher values = stronger support)
- **BD**: Bayesian Decay (higher values = stronger support)
- **BF**: Bayes Factor (exponential of BD)
- **PD**: Parsimony Decay (traditional Bremer support)
- Branch colors and exact label formatting may vary by tree viewer; this example uses FigTree's visualization

### Detailed Markdown Report (`<output_stem>.md`)
A Markdown file providing a more human-readable summary of the analysis parameters, summary statistics, and detailed branch support results in a table format. It also includes a brief interpretation guide. A good markdown viewer is [Joplin](https://joplinapp.org/) or [MarkdownLivePreview](https://markdownlivepreview.com/).

### Site-Specific Analysis (Optional)
If `--site-analysis` is used, additional output files are generated in a directory named `<output_stem>_site_analysis/`:

1. **`site_analysis_summary.txt`**: A summary of supporting vs. conflicting sites for each branch.
2. **`site_data_Clade_X.txt`**: For each branch, detailed site-by-site likelihood differences.
3. **`site_plot_Clade_X.png`**: Visualization of site-specific support/conflict (if matplotlib is available).
4. **`site_hist_Clade_X.png`**: Histogram showing the distribution of site likelihood differences.

This feature allows you to identify which alignment positions support or conflict with each branch in the tree.

### Visualizations (Optional)
If `--visualize` is used, static plots are generated (requires `matplotlib` and `seaborn`):
*   **Support Distribution Plot** (`<output_stem>_dist_au.<viz_format>` and `<output_stem>_dist_delta_lnl.<viz_format>`): Histograms showing the distribution of AU p-values and ΔlnL values across all tested branches.

# Understanding the Site Analysis Plots in panDecay

## What the Bar Colours Mean

In the site-specific likelihood plots generated by panDecay (such as `site_plot_Clade_X.png`): ![site_plot_Clade_X.png](./images/site_plot_Clade_X.png)

- **Green bars** represent sites that **support** the branch/clade being tested. These are alignment positions where the ML tree (with the clade present) has a better likelihood than the constrained tree (where the clade is forced to be non-monophyletic).

- **Red bars** represent sites that **conflict with** the branch/clade being tested. These are alignment positions where the constrained tree (without the clade) actually has a better likelihood than the ML tree.

## What "Delta lnL" Means

![site_hist_Clade_X.png](./images/site_hist_Clade_X.png)

"Delta lnL" (ΔlnL) refers to the difference in site-specific log-likelihood between the ML tree and the constrained tree for each site in your alignment. Specifically:

```
Delta lnL = lnL_ML - lnL_constrained
```

Where:
- **lnL_ML** is the log-likelihood of that specific site in the maximum likelihood tree (with the clade present)
- **lnL_constrained** is the log-likelihood of that site in the constrained tree (where the clade is forced to be non-monophyletic)

## Interpreting the Values

1. **Negative Delta lnL (green bars)**: 
   - When Delta lnL is negative, it means the ML tree has a better (less negative) likelihood for that site than the constrained tree
   - These sites provide evidence **supporting** the clade's existence in the tree
   - The more negative the value, the stronger the support from that site

2. **Positive Delta lnL (red bars)**:
   - When Delta lnL is positive, it means the constrained tree has a better likelihood for that site than the ML tree
   - These sites provide evidence **against** the clade's existence
   - The more positive the value, the stronger the conflict

3. **Values near zero**:
   - Sites with Delta lnL values very close to zero are effectively neutral regarding this particular branch
   - They don't strongly support or conflict with the branch

## Additional Information in the Plots

The site-specific analysis plots also contain a text box with summary statistics:
- **Supporting sites**: Total number of sites with negative Delta lnL (green bars)
- **Conflicting sites**: Total number of sites with positive Delta lnL (red bars)
- **Support ratio**: The ratio of supporting sites to conflicting sites
- **Weighted ratio**: The ratio of the sum of absolute values of supporting deltas to the sum of conflicting deltas

## Practical Significance

These visualizations allow you to identify which specific positions in your alignment are driving the support or conflict for a particular branch. This can be useful for:

1. Detecting potential alignment errors or problematic regions
2. Identifying sites that might be under different selective pressures
3. Finding evidence of recombination or horizontal gene transfer
4. Understanding the strength of evidence for contentious branches in your phylogeny

A branch with many strong green bars and few red bars has robust evidence across many sites. A branch with a more balanced mix of green and red bars, or with only a few strong green bars, has more tenuous support and might be less reliable.


### Temporary Files (Debug/Keep)
If `--debug` or `--keep-files` is used, a temporary directory (usually in `debug_runs/mldecay_<timestamp>/` or a user-specified path) will be retained. This directory contains:
*   `alignment.nex`: The alignment converted to NEXUS format.
*   `ml_search.nex`, `paup_ml.log`: PAUP\* script and log for the initial ML tree search.
*   `ml_tree.tre`, `ml_score.txt`: The best ML tree and its likelihood score file.
*   `constraint_search_*.nex`, `paup_constraint_*.log`: PAUP\* scripts and logs for each constrained search.
*   `constraint_tree_*.tre`, `constraint_score_*.txt`: Constrained trees and their score files.
*   `au_test.nex`, `paup_au.log`: PAUP\* script and log for the AU test.
*   `au_test_results.txt`: Score file from the AU test (though the log is primarily parsed).
*   `bootstrap_search.nex`, `paup_bootstrap.log`, `bootstrap_trees.tre` (if `--bootstrap` used): Bootstrap analysis files.
*   `site_analysis_*.nex`, `site_lnl_*.txt` (if `--site-analysis` used): Site-specific likelihood files.
*   `mldecay_debug.log` (in the main execution directory if `--debug` is on): Detailed script execution log.

## Examples & Recipes

Let [alignment.fas](./alignment.fas) be a FASTA DNA alignment, [proteins.phy](./proteins.phy) be a PHYLIP protein alignment and [morpho.nex](./morpho.nex) be a morphological dataset.

### Example 1: Basic DNA Analysis
Analyze a DNA alignment with GTR+G+I model, automatically estimating parameters.

```bash
pandecay alignment.fas --model GTR --gamma --invariable --data-type dna \
    --output dna_decay.txt --tree dna_annotated
```

### Example 2: Parsimony Analysis (Traditional Bremer Support)
Calculate traditional Bremer support values using parsimony analysis.

```bash
pandecay alignment.fas --analysis parsimony \
    --output parsimony_bremer.txt --tree parsimony_annotated
```

### Example 3: Protein Data with Specific Model
Analyze a protein alignment using the WAG model, fixed gamma shape, and estimating proportion of invariable sites.

```bash
pandecay proteins.phy --format phylip --data-type protein \
    --protein-model WAG --gamma --gamma-shape 0.85 --invariable \
    --output protein_decay.txt --tree protein_annotated --threads 8
```

### Example 4: Discrete Morphological Data
Analyze a binary (0/1) discrete morphological dataset (e.g., in NEXUS format `morpho.nex`) using the Mk+G model.

```bash
pandecay morpho.nex --format nexus --data-type discrete \
    --model Mk --gamma \
    --output morpho_decay.txt --tree morpho_annotated
```
*Note: For discrete data, ensure characters are '0' and '1'. `--parsmodel` (default for discrete) will use parsimony-like branch lengths.*

### Example 5: Using a Starting Tree
Perform a GTR+G analysis, but provide PAUP* with a starting tree to potentially speed up or refine the initial ML search.

```bash
pandecay alignment.fas --model GTR --gamma \
    --starting-tree my_start_tree.nwk \
    --output results_with_start_tree.txt
```

### Example 6: Advanced Control with PAUP\* Block
Use a custom PAUP\* block for complex settings. Assume `my_paup_commands.txt` contains:
```paup
lset nst=6 basefreq=empirical rates=gamma(categories=8) shape=estimate pinv=0.1;
hsearch nreps=50 swap=tbr addseq=random hold=1 multrees=yes;
```
Then run:
```bash
pandecay alignment.fas --paup-block my_paup_commands.txt \
    --output results_custom_block.txt
```
*(panDecay will still handle the constraint generation and AU test logic around your block.)*

### Example 7: Site-Specific Analysis
Analyze which sites in the alignment support or conflict with each clade:

```bash
pandecay alignment.fas --model GTR --gamma --site-analysis --visualize \
    --output site_analysis_results.txt
```

This will generate site-specific likelihood analyses in addition to the standard branch support results.

### Example 8: Bootstrap Analysis
Perform bootstrap analysis (100 replicates by default) alongside ML decay indices:

```bash
pandecay alignment.fas --model GTR --gamma --bootstrap \
    --output with_bootstrap.txt
```

For more bootstrap replicates:

```bash
pandecay alignment.fas --model GTR --gamma --bootstrap --bootstrap-reps 500 \
    --output bootstrap500.txt
```

This will produce additional tree files with bootstrap values and a comprehensive tree that combines bootstrap values with ML decay indices.

### Example 9: Bayesian Analysis Only
Perform only Bayesian decay analysis using MrBayes:

```bash
pandecay alignment.fas --analysis bayesian --bayesian-software mrbayes \
    --bayes-model GTR --bayes-ngen 500000 --output bayesian_only.txt
```

### Example 10: Combined ML and Bayesian Analysis
Run both ML and Bayesian analyses:

```bash
pandecay alignment.fas --model GTR --gamma --analysis ml+bayesian --bayesian-software mrbayes \
    --bayes-ngen 1000000 --output combined_analysis.txt
```

### Example 11: Combined ML and Parsimony Analysis
Run both ML and parsimony analyses to compare modern and traditional support values:

```bash
pandecay alignment.fas --model GTR --gamma --analysis ml+parsimony \
    --output ml_parsimony_analysis.txt
```

### Example 12: All Three Analysis Types
Run ML, Bayesian, and parsimony analyses in a single run:

```bash
pandecay alignment.fas --model GTR --gamma --analysis all \
    --bayesian-software mrbayes --bayes-ngen 1000000 --output complete_analysis.txt
```

### Example 13: Using MPI for Parallel MrBayes
If you have MPI-enabled MrBayes installed:

```bash
pandecay alignment.fas --analysis bayesian --bayesian-software mrbayes --use-mpi \
    --mpi-processors 8 --bayes-chains 4 --bayes-ngen 2000000
```

This runs 4 chains across 8 processors (2 chains per processor for better mixing).

### Example 14: Using BEAGLE for GPU Acceleration
If MrBayes is compiled with BEAGLE support:

```bash
pandecay alignment.fas --analysis ml+bayesian --bayesian-software mrbayes --use-beagle \
    --beagle-device gpu --beagle-precision single --bayes-ngen 5000000
```

For CPU-based BEAGLE acceleration:

```bash
pandecay alignment.fas --analysis ml+bayesian --bayesian-software mrbayes --use-beagle \
    --beagle-device cpu --beagle-precision double
```

### Example 15: Combined MPI and BEAGLE
For maximum performance with both MPI and BEAGLE:

```bash
pandecay large_alignment.fas --analysis bayesian --bayesian-software mrbayes \
    --use-mpi --mpi-processors 16 --use-beagle --beagle-device gpu \
    --bayes-chains 4 --bayes-ngen 10000000 --bayes-sample-freq 5000
```

## Installation Requirements for Parallel Processing

### For MPI Support
To use `--use-mpi`, you need MrBayes compiled with MPI support. Follow the MrBayes manual instructions to compile with `--enable-mpi=yes`.

### For BEAGLE Support  
To use `--use-beagle`, you need:
1. BEAGLE library installed (GPU or CPU version)
2. MrBayes compiled with `--with-beagle` flag

Example installation on macOS:
```bash
# Install BEAGLE
brew install beagle-lib

# Compile MrBayes with BEAGLE and MPI
./configure --with-beagle --enable-mpi=yes
make && sudo make install
```

### Example 16: Quick Bayesian Test
For a quick test with minimal MCMC generations:

```bash
pandecay alignment.fas --analysis ml+bayesian --bayesian-software mrbayes \
    --bayes-ngen 10000 --bayes-sample-freq 100 \
    --output quick_test.txt
```

### Example 17: Using Configuration Files
Generate a template configuration file and use it for analysis:

```bash
# Generate template
pandecay --generate-config my_analysis.ini

# Edit my_analysis.ini with your parameters, then run:
pandecay --config my_analysis.ini

# Override config file settings with command-line arguments
pandecay --config my_analysis.ini --threads 16 --output different_output.txt
```

### Example 18: Testing Specific Branches
Test only specific clades of interest:

```bash
# Test only clades containing specific taxa (semicolon-separated)
pandecay alignment.fas --constraint-mode specific \
    --test-branches "Homo_sapiens,Pan_troglodytes;Mus_musculus,Rattus_norvegicus"

# Test specific branch IDs from a previous analysis
pandecay alignment.fas --constraint-mode specific \
    --test-branches "1,3,5,7"

# Read constraints from a file
pandecay alignment.fas --constraint-mode specific \
    --test-branches "@my_constraints.txt"

# Test all branches EXCEPT specified ones
pandecay alignment.fas --constraint-mode exclude \
    --test-branches "Drosophila_melanogaster,Anopheles_gambiae"
```

### Example 19: Combined Config File with Constraints
Create a configuration file with constraint definitions:

```ini
# my_analysis.ini
alignment = vertebrates.fas
model = GTR
gamma = true
analysis = all
constraint_mode = specific

[constraints]
primates = Homo_sapiens,Pan_troglodytes,Gorilla_gorilla
rodents = Mus_musculus,Rattus_norvegicus
birds = Gallus_gallus,Taeniopygia_guttata
```

Then run:
```bash
pandecay --config my_analysis.ini
```

## Interpreting Results

*   **ML Tree Log-Likelihood:** The baseline score for your optimal tree.
*   **Constrained Log-Likelihood (`Constrained_lnL`):** The score of the best tree found when a particular clade was forced to be non-monophyletic. This score will typically be worse (more positive, since they are -lnL) than the ML tree's score.
*   **Log-Likelihood Difference (`Delta_LnL` or `ΔlnL`):**
    *   Calculated as `Constrained_lnL - ML_lnL`.
    *   A larger positive value (i.e., the constrained tree has much worse likelihood) indicates stronger support for the original clade. This is the "decay" value in the likelihood sense.
*   **AU p-value:**
    *   Tests the null hypothesis that the ML tree is not significantly better than the constrained alternative tree (where the clade is broken).
    *   A **low p-value (e.g., < 0.05)** leads to rejecting the null hypothesis. This means the constrained tree is significantly worse, providing statistical support for the original clade's monophyly.
    *   A **high p-value (e.g., > 0.05)** means we cannot reject the null hypothesis; the data do not provide strong statistical evidence to prefer the ML tree (with the clade) over the alternative (clade broken). This implies weaker support for that specific clade.
*   **Bootstrap Value (if bootstrap analysis performed):**
    *   Percentage of bootstrap replicates in which the clade appears.
    *   Higher values (e.g., > 70%) indicate stronger support.
    *   Bootstrap is a widely-used and well-understood method, providing a complementary measure of support to the AU test and ΔlnL values.
*   **Bayesian Decay (BD):**
    *   The primary metric for Bayesian support: marginal log-likelihood difference (unconstrained - constrained)
    *   **Key insight**: In phylogenetic topology testing, BD values typically closely approximate ML log-likelihood differences
    *   **Interpretation**: BD values should be interpreted comparatively rather than using absolute thresholds:
        *   Compare BD values across branches to identify relatively well-supported vs poorly-supported clades
        *   Evaluate BD values alongside other metrics (ΔlnL, AU test, bootstrap, parsimony decay)
        *   Consider that BD values may scale with alignment properties (size, diversity)
        *   The strongest evidence for clade support comes from concordance across multiple metrics
    *   **Note**: BD values of 30-50 or higher are common when data strongly support a clade and should not be considered anomalous
    *   **Why traditional BF scales don't apply**: Traditional Bayes factor thresholds were developed for comparing fundamentally different models, not for topology testing where models differ only by a single constraint
    *   **Negative values** suggest potential issues:
        *   Poor MCMC convergence (check convergence diagnostics)
        *   Marginal likelihood estimation problems
        *   Genuine lack of support for the clade
*   **Bayes Factor (BF):**
    *   The exponential of the Bayesian decay value (BF = e^BD)
    *   Displayed with cap at 10^6 to avoid astronomical numbers
    *   **Not recommended for interpretation in phylogenetics**: Use BD values instead, as traditional BF interpretation scales are misleading for topology testing

**Site-Specific Analysis Interpretation:**
* **Negative delta lnL values** indicate sites that support the branch (they become less likely when the branch is constrained to be absent).
* **Positive delta lnL values** indicate sites that conflict with the branch (they become more likely when the branch is removed).
* **Values near zero** indicate sites that are neutral regarding this branch.

Generally, clades with large positive `ΔlnL` values, low `AU_p-value`s, high bootstrap values, and many supporting sites are considered well-supported.

## Troubleshooting

*   **"PAUP\* not found"**: Ensure PAUP\* is installed and either in your system PATH or you are providing the full path via `--paup /path/to/paup`.
*   **"I/O operation on closed file" / PAUP\* crashes for some constraints**:
    *   This can occur if PAUP\* or the system is under extreme load. Ensure you are not using all CPU cores for PAUP\*. Use `--threads auto` (default) or specify a number of threads less than your total core count (e.g., `total_cores - 2`).
    *   Check the PAUP\* logs in the temporary directory (if kept with `--keep-files` or `--debug`) for specific PAUP\* error messages.
*   **Low Support Values / Unexpected Results**:
    *   Ensure your chosen evolutionary model is appropriate for your data.
    *   The heuristic search in PAUP\* (`hsearch`) may not always find the global optimum. More intensive search settings (e.g., more `nreps`, different `swap` algorithms) might be needed, potentially by using a custom `--paup-block`.
    *   The data itself may not contain strong signal for certain relationships.
*   **Python Errors**: Ensure all Python dependencies (BioPython, NumPy) are correctly installed for the Python interpreter you are using.
*   **Site-specific analysis errors**: If the site analysis fails but the main analysis succeeds, try running the analysis again with `--keep-files` and check the site-specific likelihood files in the temporary directory.
*   **Bootstrap fails but ML analysis succeeds**: Bootstrap analysis can be computationally intensive. Try using fewer bootstrap replicates (`--bootstrap-reps 50`) or allocate more processing time by increasing the timeout value in the code.
*   **MrBayes errors**:
    *   "Command not found": Ensure MrBayes is installed and accessible. Use `--mrbayes-path /path/to/mb` if it's not in your PATH.
    *   "Error in command 'Ss'": Your version of MrBayes doesn't support stepping-stone sampling. The program will automatically use harmonic mean estimation instead.
    *   Path errors with spaces: panDecay handles paths with spaces automatically by using relative paths for MrBayes execution.
    *   No Bayesian output: Check the debug log for specific errors. Ensure your alignment is compatible with MrBayes (e.g., taxon names without special characters).
*   **Bayesian analysis takes too long**: Reduce the number of generations (`--bayes-ngen 50000`) or increase sampling frequency (`--bayes-sample-freq 500`) for testing. Production runs typically need at least 1 million generations.

## Understanding Bayes Factors in Phylogenetic Topology Testing

Traditional Bayes factor interpretation guidelines were developed for comparing fundamentally different models (e.g., linear vs. polynomial regression) and do not apply well to phylogenetic topology testing. Here's why:

### Why Phylogenetic Bayes Factors Behave Differently

1. **Minimal Model Complexity Differences**
   - In topology testing, we compare models that differ only by a single topological constraint
   - All other model parameters (substitution model, rates, frequencies) remain identical
   - Traditional BF penalties for model complexity are largely irrelevant here

2. **Likelihood Dominates Marginal Likelihood**
   - When data strongly support a topology, the posterior concentrates around the ML tree
   - The marginal likelihood becomes dominated by the likelihood component
   - This explains why Bayesian decay values closely approximate ML log-likelihood differences

3. **Different Prior Effects**
   - Traditional BF applications involve priors over fundamentally different parameter spaces
   - In topology testing, the prior difference is only in tree space
   - When data signal is strong, this prior difference has minimal impact

### Key Insights for Phylogenetics

1. **BD ≈ ΔlnL**: The Bayesian decay value typically approximates the ML log-likelihood difference
   - This is expected behavior, not an anomaly
   - Both metrics primarily reflect the strength of data support for the topology

2. **Large BD Values Are Normal**: 
   - BD values of 30-50 or higher are common for well-supported clades
   - These don't indicate "astronomical" support as traditional BF scales would suggest
   - They simply reflect strong data signal for the clade

3. **Phylogenetic-Specific Interpretation**:
   - Focus on BD values rather than Bayes factors
   - Use phylogenetic-specific thresholds (BD: 0-2 weak, 2-5 moderate, 5-10 strong, >10 very strong)
   - Compare BD values across clades in your tree for relative support assessment

### Practical Implications

When interpreting panDecay results:
- **Don't be alarmed** by large Bayes factors or BD values
- **Compare BD and ML differences** - similar values confirm proper analysis
- **Use BD values** for interpretation, not traditional BF thresholds
- **Consider relative support** across branches rather than absolute thresholds

This understanding helps avoid misinterpretation and provides more accurate assessment of phylogenetic support.

## Citations

If you use panDecay in your research, please cite this GitHub repository. Additionally, consider citing the relevant methodological papers:

*   **PAUP\***:
    *   Swofford, D. L. (2003). PAUP\*. Phylogenetic Analysis Using Parsimony (\*and Other Methods). Version 4. Sinauer Associates, Sunderland, Massachusetts.
*   **Bremer Support (Decay Index) - Original Concept (Parsimony)**:
    *   Bremer, K. (1988). The limits of amino acid sequence data in angiosperm phylogenetic reconstruction. *Evolution*, 42(4), 795-803.
    *   Bremer, K. (1994). Branch support and tree stability. *Cladistics*, 10(3), 295-304.
*  **Approximately Unbiased (AU) Test**:
    *   Shimodaira, H. (2002). An approximately unbiased test of phylogenetic tree selection. *Systematic Biology*, 51(3), 492-508.
*   **General ML Phylogenetics**:
    *   Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368-376.
*   **Site-specific likelihood methods**:
    *   Goldman, N., Anderson, J. P., & Rodrigo, A. G. (2000). Likelihood-based tests of topologies in phylogenetics. *Systematic Biology*, 49(4), 652-670.
*   **Bootstrap Methods**:
    *   Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. *Evolution*, 39(4), 783-791.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to:
*   Open an issue on GitHub.
*   Submit a pull request with your changes.

## Contact

For questions or support, please open an issue on the GitHub repository.

Project Maintainer: James McInerney

