![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg) ![Version](https://img.shields.io/badge/version-1.0.0-orange.svg) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# MLDecay: Maximum Likelihood-based Phylogenetic Decay Indices

MLDecay is a Python command-line tool for calculating Maximum Likelihood (ML)-based phylogenetic decay indices, also known as ML-Bremer support or ML branch support. It leverages the phylogenetic software PAUP* to assess the robustness of clades in a phylogenetic tree by comparing the likelihood of the optimal ML tree with trees where specific clades are constrained to be non-monophyletic. The primary statistical measure used is the Approximately Unbiased (AU) test.

## Table of Contents

1.  [Background](#background)
    *   [What are Decay Indices / Bremer Support?](#what-are-decay-indices--bremer-support)
    *   [Why ML-based Decay Indices?](#why-ml-based-decay-indices)
2.  [Features](#features)
3.  [Installation](#installation)
    *   [Dependencies](#dependencies)
    *   [Installing MLDecay](#installing-mldecay)
4.  [Usage](#usage)
    *   [Basic Command](#basic-command)
    *   [Command-Line Arguments](#command-line-arguments)
5.  [Input Files](#input-files)
    *   [Sequence Alignment](#sequence-alignment)
    *   [Optional Starting Tree](#optional-starting-tree)
    *   [Optional PAUP\* Block File](#optional-paup-block-file)
6.  [Output Files](#output-files)
    *   [Main Results File (`ml_decay_indices.txt`)](#main-results-file-ml_decay_indicestxt)
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

MLDecay automates this process by:
1.  Finding the optimal ML tree and its likelihood score.
2.  For each internal branch in the ML tree:
    a.  Defining a constraint that forces the taxa in that clade to *not* form a monophyletic group (using PAUP\*'s `converse=yes` constraint).
    b.  Searching for the best ML tree under this "anti-constraint" and recording its likelihood.
3.  Calculating the difference in log-likelihood between the unconstrained ML tree and each constrained tree.
4.  Performing an Approximately Unbiased (AU) test (Shimodaira, 2002) to statistically compare the unconstrained ML tree against all the constrained alternative trees. The p-value from the AU test indicates the significance of the support for the original clade.

A significantly worse likelihood for the constrained tree (and a low AU test p-value for that constrained tree) provides strong evidence in favor of the monophyly of the original clade.

## Features

*   Calculates ML-based decay values using log-likelihood differences.
*   Performs the Approximately Unbiased (AU) test for statistical assessment of branch support.
*   Supports DNA, Protein, and binary discrete morphological data.
*   Optional site-specific likelihood analysis to identify which alignment positions support or conflict with each branch.
*   Flexible model specification (e.g., GTR, HKY, JTT, WAG, Mk) with options for gamma-distributed rate heterogeneity (+G) and proportion of invariable sites (+I).
*   Allows fine-grained control over model parameters (gamma shape, pinvar, base frequencies, etc.).
*   Option to provide a custom PAUP\* block for complex model or search strategy definitions.
*   Option to provide a starting tree for the initial ML search.
*   Generates comprehensive output:
    *   Tab-delimited results file.
    *   Multiple Newick trees annotated with different support values.
    *   Detailed Markdown report summarizing the analysis and results.
*   Optional static visualization (requires `matplotlib` and `seaborn`):
    *   Distribution of support values.
    *   Site-specific support visualizations.
*   Multi-threaded PAUP\* execution (configurable).
*   Debug mode and option to keep temporary files.

## Installation

### Dependencies

MLDecay requires Python 3.8 or higher and has several dependencies that can be easily installed using pip.

1. **PAUP\***: You must have a working PAUP\* command-line executable installed and accessible in your system's PATH, or provide the full path to it. PAUP\* can be obtained from [phylosolutions.com](http://phylosolutions.com/paup-test/).

2. **Python Dependencies**: All Python dependencies can be installed using the provided `requirements.txt` file:

   ```bash
   pip install -r requirements.txt
   ```

   This will install all required packages including BioPython, NumPy, and the optional visualization packages Matplotlib and Seaborn.

### Installing MLDecay

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/mol-evol/MLDecay.git
   cd MLDecay
   ```

2. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Make the script executable (optional, for convenience):**
   ```bash
   chmod +x MLDecay.py
   ```

4. You can then run the script directly:
   ```bash
   ./MLDecay.py [arguments...]
   ```
   or using the python interpreter:
   ```bash
   python3 MLDecay.py [arguments...]
   ```

5. **Optional: Make MLDecay available system-wide**
   
   Consider adding the MLDecay directory to your system's PATH or creating a symbolic link to `MLDecay.py` in a directory that is already in your PATH (e.g., `~/.local/bin/` or `/usr/local/bin/`).

   For example, to create a symbolic link:
   ```bash
   ln -s $(pwd)/MLDecay.py ~/.local/bin/mldecay
   ```

## Usage

### Basic Command

```bash
python3 MLDecay.py <alignment_file> --model <model_name> [options...]
```

### Command-Line Arguments

```
usage: MLDecay.py [-h] [--format FORMAT] [--model MODEL] [--gamma] [--invariable] [--paup PAUP] [--output OUTPUT] [--tree TREE]
                  [--data-type {dna,protein,discrete}] [--gamma-shape GAMMA_SHAPE] [--prop-invar PROP_INVAR] 
                  [--base-freq {equal,estimate,empirical}] [--rates {equal,gamma}] [--protein-model PROTEIN_MODEL] 
                  [--nst {1,2,6}] [--parsmodel | --no-parsmodel] [--threads THREADS] [--starting-tree STARTING_TREE] 
                  [--paup-block PAUP_BLOCK] [--temp TEMP] [--keep-files] [--debug] [--site-analysis] [--visualize]
                  [--viz-format {png,pdf,svg}] [-v]
                  alignment

MLDecay v0.2.1: Calculate ML-based phylogenetic decay indices using PAUP*.

positional arguments:
  alignment             Input alignment file path.

options:
  -h, --help            show this help message and exit
  --format FORMAT       Alignment format. (default: fasta)
  --model MODEL         Base substitution model (e.g., GTR, HKY, JC). Combine with --gamma and --invariable. (default: GTR)
  --gamma               Add Gamma rate heterogeneity (+G) to model. (default: False)
  --invariable          Add invariable sites (+I) to model. (default: False)
  --paup PAUP           Path to PAUP* executable. (default: paup)
  --output OUTPUT       Output file for summary results. (default: ml_decay_indices.txt)
  --tree TREE           Base name for annotated tree files. Three trees will be generated with suffixes: _au.nwk (AU p-values), 
                        _lnl.nwk (likelihood differences), and _combined.nwk (both values). (default: annotated_tree)
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
  --threads THREADS     Number of threads for PAUP* (e.g., 4, 'auto' for (total_cores - 2), or 'all'). Using 'auto' or leaving some cores free is
                        recommended for system stability. (default: auto)
  --starting-tree STARTING_TREE
                        Path to a user-provided starting tree file (Newick).
  --paup-block PAUP_BLOCK
                        Path to file with custom PAUP* commands for model/search setup (overrides most model args).
  --temp TEMP           Custom directory for temporary files (default: system temp).
  --keep-files          Keep temporary files after analysis. (default: False)
  --debug               Enable detailed debug logging (implies --keep-files). (default: False)

Visualization Output (optional):
  --visualize           Generate static visualization plots (requires matplotlib, seaborn). (default: False)
  --viz-format {png,pdf,svg}
                        Format for static visualizations. (default: png)
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
    MLDecay will try to defensively add `savetrees` and `lscores` commands if they appear to be missing from the user's block when needed for its internal workflow.

## Output Files

Unless specified with `--output`, `--tree`, etc., output files are created in the current working directory.

### Main Results File (`ml_decay_indices.txt` by default)
A tab-delimited text file containing:
*   The log-likelihood of the best ML tree found.
*   For each internal branch tested:
    *   `Clade_ID`: An internal identifier for the branch.
    *   `Num_Taxa`: Number of taxa in the clade defined by this branch.
    *   `Constrained_lnL`: Log-likelihood of the best tree found when this clade was constrained to be non-monophyletic.
    *   `LnL_Diff_from_ML`: Difference between `Constrained_lnL` and the ML tree's likelihood.
    *   `AU_p-value`: The p-value from the Approximately Unbiased test.
    *   `Significant_AU (p<0.05)`: "Yes" if AU p-value < 0.05, "No" otherwise.
    *   `Taxa_List`: A comma-separated list of taxa in the clade.

### Annotated Trees
MLDecay now generates three different annotated tree files:
* `<tree_base>_au.nwk`: Tree with AU test p-values as branch labels
* `<tree_base>_lnl.nwk`: Tree with log-likelihood differences as branch labels
* `<tree_base>_combined.nwk`: Tree with both values as branch labels in the format "AU:0.95|LnL:2.34"

These trees can be visualized in standard tree viewers like FigTree, Dendroscope, iTOL, etc. The combined tree is particularly suited for FigTree which handles string labels well.

### Detailed Markdown Report (`<output_stem>.md`)
A Markdown file providing a more human-readable summary of the analysis parameters, summary statistics, and detailed branch support results in a table format. It also includes a brief interpretation guide.

### Site-Specific Analysis (Optional)
If `--site-analysis` is used, additional output files are generated in a directory named `<output_stem>_site_analysis/`:

1. **`site_analysis_summary.txt`**: A summary of supporting vs. conflicting sites for each branch.
2. **`site_data_Clade_X.txt`**: For each branch, detailed site-by-site likelihood differences.
3. **`site_plot_Clade_X.png`**: Visualization of site-specific support/conflict (if matplotlib is available).

This feature allows you to identify which alignment positions support or conflict with each branch in the tree.


# Understanding the Site Analysis Plots in MLDecay

## What the Bar Colours Mean

In the site-specific likelihood plots generated by MLDecay (such as `site_plot_Clade_X.png` ![site_plot_Clade_X.png](:/site_plot_Clade_X.png)):

- **Green bars** represent sites that **support** the branch/clade being tested. These are alignment positions where the ML tree (with the clade present) has a better likelihood than the constrained tree (where the clade is forced to be non-monophyletic).

- **Red bars** represent sites that **conflict with** the branch/clade being tested. These are alignment positions where the constrained tree (without the clade) actually has a better likelihood than the ML tree.

## What "Delta lnL" Means

![site_hist_Clade_X.png](:/site_hist_Clade_X.png)

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



### Visualizations (Optional)
If `--visualize` is used, static plots are generated (requires `matplotlib` and `seaborn`):
*   **Support Distribution Plot** (`<output_stem>_dist_au.<viz_format>` and `<output_stem>_dist_lnl.<viz_format>`): Histograms showing the distribution of AU p-values and LNL differences across all tested branches.

### Temporary Files (Debug/Keep)
If `--debug` or `--keep-files` is used, a temporary directory (usually in `debug_runs/mldecay_<timestamp>/` or a user-specified path) will be retained. This directory contains:
*   `alignment.nex`: The alignment converted to NEXUS format.
*   `ml_search.nex`, `paup_ml.log`: PAUP\* script and log for the initial ML tree search.
*   `ml_tree.tre`, `ml_score.txt`: The best ML tree and its likelihood score file.
*   `constraint_search_*.nex`, `paup_constraint_*.log`: PAUP\* scripts and logs for each constrained search.
*   `constraint_tree_*.tre`, `constraint_score_*.txt`: Constrained trees and their score files.
*   `au_test.nex`, `paup_au.log`: PAUP\* script and log for the AU test.
*   `au_test_results.txt`: Score file from the AU test (though the log is primarily parsed).
*   `site_analysis_*.nex`, `site_lnl_*.txt` (if `--site-analysis` used): Site-specific likelihood files.
*   `mldecay_debug.log` (in the main execution directory if `--debug` is on): Detailed script execution log.

## Examples & Recipes

Let `alignment.fas` be a FASTA DNA alignment and `proteins.phy` be a PHYLIP protein alignment.

### Example 1: Basic DNA Analysis
Analyze a DNA alignment with GTR+G+I model, automatically estimating parameters.

```bash
python3 MLDecay.py alignment.fas --model GTR --gamma --invariable --data-type dna \
    --output dna_decay.txt --tree dna_annotated
```

### Example 2: Protein Data with Specific Model
Analyze a protein alignment using the WAG model, fixed gamma shape, and estimating proportion of invariable sites.

```bash
python3 MLDecay.py proteins.phy --format phylip --data-type protein \
    --protein-model WAG --gamma --gamma-shape 0.85 --invariable \
    --output protein_decay.txt --tree protein_annotated --threads 8
```

### Example 3: Discrete Morphological Data
Analyze a binary (0/1) discrete morphological dataset (e.g., in NEXUS format `morpho.nex`) using the Mk+G model.

```bash
python3 MLDecay.py morpho.nex --format nexus --data-type discrete \
    --model Mk --gamma \
    --output morpho_decay.txt --tree morpho_annotated
```
*Note: For discrete data, ensure characters are '0' and '1'. `--parsmodel` (default for discrete) will use parsimony-like branch lengths.*

### Example 4: Using a Starting Tree
Perform a GTR+G analysis, but provide PAUP* with a starting tree to potentially speed up or refine the initial ML search.

```bash
python3 MLDecay.py alignment.fas --model GTR --gamma \
    --starting-tree my_start_tree.nwk \
    --output results_with_start_tree.txt
```

### Example 5: Advanced Control with PAUP\* Block
Use a custom PAUP\* block for complex settings. Assume `my_paup_commands.txt` contains:
```paup
lset nst=6 basefreq=empirical rates=gamma(categories=8) shape=estimate pinv=0.1;
hsearch nreps=50 swap=tbr addseq=random hold=1 multrees=yes;
```
Then run:
```bash
python3 MLDecay.py alignment.fas --paup-block my_paup_commands.txt \
    --output results_custom_block.txt
```
*(MLDecay will still handle the constraint generation and AU test logic around your block.)*

### Example 6: Site-Specific Analysis
Analyze which sites in the alignment support or conflict with each clade:

```bash
python3 MLDecay.py alignment.fas --model GTR --gamma --site-analysis --visualize \
    --output site_analysis_results.txt
```

This will generate site-specific likelihood analyses in addition to the standard branch support results.

## Interpreting Results

*   **ML Tree Log-Likelihood:** The baseline score for your optimal tree.
*   **Constrained Log-Likelihood (`Constrained_lnL`):** The score of the best tree found when a particular clade was forced to be non-monophyletic. This score will typically be worse (more positive, since they are -lnL) than the ML tree's score.
*   **Log-Likelihood Difference (`LnL_Diff_from_ML`):**
    *   Calculated as `Constrained_lnL - ML_lnL`.
    *   A larger positive value (i.e., the constrained tree is much worse) indicates stronger support for the original clade. This is the "decay" value in the likelihood sense.
*   **AU p-value:**
    *   Tests the null hypothesis that the ML tree is not significantly better than the constrained alternative tree (where the clade is broken).
    *   A **low p-value (e.g., < 0.05)** leads to rejecting the null hypothesis. This means the constrained tree is significantly worse, providing statistical support for the original clade's monophyly.
    *   A **high p-value (e.g., > 0.05)** means we cannot reject the null hypothesis; the data do not provide strong statistical evidence to prefer the ML tree (with the clade) over the alternative (clade broken). This implies weaker support for that specific clade.

**Site-Specific Analysis Interpretation:**
* **Negative delta lnL values** indicate sites that support the branch (they become less likely when the branch is constrained to be absent).
* **Positive delta lnL values** indicate sites that conflict with the branch (they become more likely when the branch is removed).
* **Values near zero** indicate sites that are neutral regarding this branch.

Generally, clades with large positive `LnL_Diff_from_ML` values, low `AU_p-value`s, and many supporting sites are considered well-supported.

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

## Citations

If you use MLDecay in your research, please cite this GitHub repository. Additionally, consider citing the relevant methodological papers:

*   **PAUP\***:
    *   Swofford, D. L. (2003). PAUP\*. Phylogenetic Analysis Using Parsimony (\*and Other Methods). Version 4. Sinauer Associates, Sunderland, Massachusetts.
*   **Bremer Support (Decay Index) - Original Concept (Parsimony)**:
    *   Bremer, K. (1988). The limits of amino acid sequence data in angiosperm phylogenetic reconstruction. *Evolution*, 42(4), 795-803.
    *   Bremer, K. (1994). Branch support and tree stability. *Cladistics*, 10(3), 295-304.
*   **Approximately Unbiased (AU) Test**:
    *   Shimodaira, H. (2002). An approximately unbiased test of phylogenetic tree selection. *Systematic Biology*, 51(3), 492-508.
*   **General ML Phylogenetics**:
    *   Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of Molecular Evolution*, 17(6), 368-376.
*   **Site-specific likelihood methods**:
    *   Goldman, N., Anderson, J. P., & Rodrigo, A. G. (2000). Likelihood-based tests of topologies in phylogenetics. *Systematic Biology*, 49(4), 652-670.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions, bug reports, and feature requests are welcome! Please feel free to:
*   Open an issue on GitHub.
*   Submit a pull request with your changes.

## Contact

For questions or support, please open an issue on the GitHub repository.
Project Maintainer: James McInerney
