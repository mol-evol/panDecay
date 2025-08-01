# panDecay Output Guide

## Overview

panDecay generates a comprehensive set of output files organized in timestamp-based directories. This guide explains each output file type and how to interpret the results.

## Directory Structure

```
<analysis_name>_pandecay_results/
├── <alignment_name>.csv              # Main results (tab-delimited)
├── <alignment_name>_summary.txt      # Human-readable summary
├── <alignment_name>_report.md        # Detailed markdown report
├── annotated_trees/                  # Tree files with annotations
│   ├── tree_with_support.tre         # Support values annotated
│   ├── tree_with_clade_ids.tre       # Clade identifiers
│   └── tree_original.tre             # Original optimal tree
├── site_analysis/                    # Per-site likelihood analysis
│   ├── site_analysis_summary.txt     # Site support summary
│   ├── clade_<N>_sites.csv          # Detailed site data per clade
│   └── likelihood_profiles/          # Site-by-site likelihood profiles
└── temp_files/                       # Intermediate files (if --keep-files)
```

## Main Results File (CSV Format)

**File**: `<alignment_name>.csv`

This tab-delimited file contains quantitative support metrics for each analyzed clade.

### Column Descriptions

| Column | Description | Interpretation |
|--------|-------------|----------------|
| `clade_id` | Unique identifier (e.g., Clade_3) | References internal tree branches |
| `taxa_count` | Number of taxa in clade | Larger clades = more inclusive groups |
| `taxa` | Comma-separated taxa list | Specific organisms in this clade |
| `PD` | Parsimony Decay Index | Steps lost when clade constrained non-monophyletic |
| `LD` | Likelihood Decay Index (ΔlnL) | Log-likelihood difference for clade constraint |
| `BD` | Bayes Decay Index | Marginal likelihood difference (log scale) |
| `AU` | Approximately Unbiased p-value | Statistical significance of constraint |
| `Supporting_Sites` | Sites favoring clade monophyly | Number of alignment positions supporting clade |
| `Conflicting_Sites` | Sites opposing clade monophyly | Number of alignment positions opposing clade |
| `Neutral_Sites` | Sites with no preference | Alignment positions with equal likelihood |
| `Site_Support_Ratio` | Supporting/Conflicting ratio | Higher values = stronger site-level support |
| `Weighted_Support_Ratio` | Likelihood-weighted ratio | Accounts for strength of site preferences |
| `Sum_Supporting_Delta` | Total supporting likelihood | Sum of likelihood differences for supporting sites |
| `Sum_Conflicting_Delta` | Total conflicting likelihood | Sum of likelihood differences for opposing sites |

### Example Record
```csv
Clade_7,3,"Homo_sapiens, Pan_troglodytes, Gorilla_gorilla",19,45.50,40.48,0.0001,686,211,1,3.25,1.95,-93.42,47.92
```

**Interpretation**: 
- Human-chimp-gorilla clade with very strong support
- 19 parsimony steps lost when broken
- ΔlnL = 45.50 (strong likelihood support) 
- AU p-value = 0.0001 (highly significant)
- 686 sites support vs 211 oppose (3.25:1 ratio)

## Summary Report (TXT Format)

**File**: `<alignment_name>_summary.txt`

Human-readable table with formatted results and statistical significance indicators.

### Support Significance Levels
- `***` = p < 0.001 (highly significant)
- `**` = p < 0.01 (significant) 
- `*` = p < 0.05 (marginally significant)
- `ns` = not significant (p ≥ 0.05)

### Analysis Types Included
- **Maximum Likelihood**: ΔlnL and AU p-values
- **Bayesian**: Bayes Decay (BD) indices
- **Parsimony**: Traditional Bremer support (Decay steps)

## Detailed Markdown Report

**File**: `<alignment_name>_report.md`

Comprehensive analysis report including:
- Analysis parameters and settings
- Statistical summaries and interpretations
- Site analysis integration
- Methodological details
- References and citations

## Annotated Tree Files

**Directory**: `annotated_trees/`

### tree_with_support.tre
Newick format tree with support values as branch labels. Values correspond to the primary analysis type specified.

### tree_with_clade_ids.tre  
Newick format tree with clade identifiers (Clade_3, Clade_7, etc.) as branch labels for cross-referencing with results tables.

### tree_original.tre
Original optimal tree without annotations, suitable for visualization software.

## Site Analysis Files

**Directory**: `site_analysis/`

### site_analysis_summary.txt
Tab-delimited summary of site-level support for each clade:
- `Supporting_Sites`: Sites favoring clade monophyly
- `Conflicting_Sites`: Sites opposing clade monophyly  
- `Support_Ratio`: Simple ratio of supporting/conflicting sites
- `Weighted_Support_Ratio`: Likelihood-weighted support ratio
- `Sum_Supporting_Delta`: Total likelihood support from all sites
- `Sum_Conflicting_Delta`: Total likelihood conflict from all sites

### Individual Clade Files
**Pattern**: `clade_<N>_sites.csv`

Detailed per-site likelihood analysis for each clade showing:
- Site position in alignment
- Unconstrained likelihood
- Constrained likelihood  
- Likelihood difference (ΔlnL)
- Support classification (Supporting/Conflicting/Neutral)

## Understanding Support Metrics

### Parsimony Decay (PD)
- **Range**: 0 to maximum possible steps
- **Interpretation**: Higher values = stronger support
- **Threshold**: PD ≥ 1 considered minimal support

### Likelihood Decay (LD) 
- **Range**: 0 to unlimited (ΔlnL units)
- **Interpretation**: Higher values = stronger support  
- **Threshold**: LD ≥ 2 often considered meaningful

### AU Test p-values
- **Range**: 0.0 to 1.0
- **Interpretation**: Lower values = stronger rejection of constraint
- **Threshold**: p < 0.05 for statistical significance

### Site Support Ratios
- **Range**: 0 to unlimited
- **Interpretation**: Higher ratios = more sites favor the clade
- **Balanced**: Ratios near 1.0 indicate conflicted signal

## Quality Assessment

### Strong Support Indicators
- AU p-value < 0.01
- LD (ΔlnL) > 10
- Site support ratio > 2.0
- Weighted support ratio > 1.5

### Weak/Conflicted Support Indicators  
- AU p-value > 0.05
- LD (ΔlnL) < 2
- Site support ratio < 1.5
- High numbers of conflicting sites

### Interpretation Guidelines
1. **Multiple metrics agreement**: Strong when PD, LD, and site ratios all high
2. **Conflicting metrics**: Consider data quality and alignment regions
3. **Site-level analysis**: Examine individual site files for detailed patterns
4. **Statistical significance**: Always consider AU p-values for objective assessment

## File Organization Features

### Timestamp-Based Directories
All outputs organized in directories with timestamps (YYYYMMDD_HHMMSS format) for:
- **Version control**: Multiple analyses preserved
- **Comparison**: Easy to compare different parameter sets
- **Organization**: No overwriting of previous results

### Automatic Cleanup
- Temporary files removed unless `--keep-files` specified
- Only final results and essential intermediates retained
- Clean directory structure for easy navigation

## Troubleshooting Output Issues

### Missing Files
- Check console output for error messages
- Verify input file formats and paths
- Ensure external software (PAUP*, MrBayes) installation

### Empty or Incomplete Results
- Review analysis parameters for consistency
- Check alignment quality and taxon sampling
- Verify computational resources (memory, time limits)

### Inconsistent Support Values
- Normal for conflicted phylogenetic signal
- Consider alignment trimming or outlier removal
- Examine site-level patterns for systematic biases

## Citation

When publishing results from panDecay, please cite:

McInerney, J.O. (2025) panDecay: Phylogenetic Analysis with Decay Indices. http://github.com/mol-evol/panDecay/