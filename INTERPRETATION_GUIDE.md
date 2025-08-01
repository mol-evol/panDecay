# panDecay Results Interpretation Guide

## Overview

This guide provides detailed instructions for interpreting panDecay phylogenetic decay index results, understanding statistical significance, and making biological conclusions from the analysis output.

## Understanding Decay Indices

### What are Decay Indices?

Decay indices measure the amount of phylogenetic signal supporting a particular clade by quantifying how much worse the data fit becomes when that clade is constrained to be non-monophyletic.

**Key Concept**: Higher decay values indicate stronger phylogenetic support for the clade.

### Types of Decay Indices

#### 1. Parsimony Decay (PD)
- **Definition**: Number of extra steps required when clade is broken
- **Units**: Parsimony steps (integer values)
- **Range**: 0 to theoretical maximum (depends on data)
- **Interpretation**: 
  - PD = 0: No parsimony support
  - PD = 1-2: Minimal support  
  - PD = 3-10: Moderate support
  - PD ≥ 10: Strong support

#### 2. Likelihood Decay (LD) 
- **Definition**: Log-likelihood difference (ΔlnL) between optimal and constrained trees
- **Units**: Natural log units
- **Range**: 0 to unlimited
- **Interpretation**:
  - LD < 2: Weak support (difference not meaningful)
  - LD = 2-10: Moderate support
  - LD = 10-50: Strong support  
  - LD ≥ 50: Very strong support

#### 3. Bayes Decay (BD)
- **Definition**: Marginal likelihood difference between models
- **Units**: Natural log units (log Bayes factors)
- **Range**: 0 to unlimited
- **Interpretation** (following Kass & Raftery 1995):
  - BD < 2: Not worth more than a bare mention
  - BD = 2-6: Positive evidence
  - BD = 6-10: Strong evidence
  - BD ≥ 10: Very strong evidence

## Statistical Testing

### AU Test (Approximately Unbiased)
- **Purpose**: Tests whether constraint significantly worsens tree likelihood
- **Null Hypothesis**: Constrained tree is as good as optimal tree
- **p-value interpretation**:
  - p ≥ 0.05: Constraint not significantly worse (weak support for clade)
  - p < 0.05: Constraint significantly worse (significant support)
  - p < 0.01: Strong evidence against constraint
  - p < 0.001: Very strong evidence against constraint

### Significance Categories in Output
- `***` = p < 0.001 (highly significant)
- `**` = p < 0.01 (significant)
- `*` = p < 0.05 (marginally significant)  
- `ns` = p ≥ 0.05 (not significant)

## Site-Level Analysis

### Site Support Metrics

#### Supporting vs. Conflicting Sites
- **Supporting Sites**: Alignment positions that favor clade monophyly
- **Conflicting Sites**: Alignment positions that oppose clade monophyly
- **Neutral Sites**: Positions with no preference (rare in practice)

#### Site Support Ratio
```
Site Support Ratio = Supporting Sites / Conflicting Sites
```

**Interpretation**:
- Ratio = 1.0: Equal support and conflict
- Ratio > 1.0: More sites support than oppose  
- Ratio < 1.0: More sites oppose than support
- Ratio ≥ 2.0: Strong site-level support
- Ratio ≥ 3.0: Very strong site-level support

#### Weighted Support Ratio
Accounts for the strength of site preferences, not just counts.

```
Weighted Ratio = |Sum Supporting ΔlnL| / |Sum Conflicting ΔlnL|
```

**Interpretation**:
- Similar to site support ratio but weighted by likelihood differences
- More reliable than simple ratios for strong phylogenetic inference
- Values > 1.5 generally indicate robust support

## Interpreting Results Across Analysis Types

### Concordant Support (Recommended Pattern)
When multiple analysis types agree:

```
Example: Clade_7 (Human-Chimp-Gorilla)
- PD = 19 (strong parsimony support)
- LD = 45.5 (very strong likelihood support)  
- BD = 40.5 (very strong Bayesian support)
- AU p-value = 0.0001 (highly significant)
- Site ratio = 3.25 (strong site support)
```

**Conclusion**: Robust support across all methods - high confidence in clade

### Discordant Support (Requires Caution)
When analysis types disagree:

```
Example: Clade_8 (Human-Chimp)
- PD = 0 (no parsimony support)
- LD = 0.8 (very weak likelihood support)
- BD = 6.9 (moderate Bayesian support)
- AU p-value = 0.61 (not significant)
- Site ratio = 2.23 (moderate site support)
```

**Conclusion**: Conflicting evidence - examine data quality and consider alternative relationships

### Analysis-Specific Considerations

#### When Parsimony Disagrees with Likelihood/Bayesian
- **Common Causes**: 
  - Rate heterogeneity not modeled in parsimony
  - Long-branch attraction in parsimony
  - Model misspecification in likelihood/Bayesian
- **Recommendation**: Consider rate variation, examine long branches

#### When Likelihood Disagrees with Bayesian
- **Common Causes**:
  - Prior sensitivity in Bayesian analysis
  - MCMC convergence issues
  - Model adequacy differences
- **Recommendation**: Check MCMC diagnostics, try different priors

## Making Biological Conclusions

### Strong Support Criteria
A clade has strong support when:
1. **AU p-value < 0.01** (statistical significance)
2. **LD ≥ 10** (substantial likelihood support)
3. **Site ratio ≥ 2.0** (favorable site-level evidence)
4. **Concordance** across multiple analysis types

### Interpreting Weak Support
Weak support may indicate:
- **Insufficient phylogenetic signal** in the data
- **Conflicting evolutionary histories** (e.g., incomplete lineage sorting)
- **Rapid evolutionary radiations** with short internal branches
- **Data quality issues** (alignment errors, model misspecification)

### Guidelines for Publication

#### High Confidence Clades
- Multiple analysis types with LD ≥ 10, AU p < 0.01
- Site ratios ≥ 2.0 across most clades
- **Statement**: "This clade received strong support across all analyses"

#### Moderate Confidence Clades  
- LD = 2-10, AU p < 0.05, but some discordance
- **Statement**: "This clade received moderate support with some conflicting signal"

#### Low Confidence Clades
- LD < 2, AU p ≥ 0.05, or strong discordance
- **Statement**: "Support for this clade was weak/conflicting - relationships remain uncertain"

## Common Interpretation Mistakes

### 1. Overinterpreting Single Metrics
❌ **Wrong**: "LD = 15, so this clade is strongly supported"
✅ **Right**: "LD = 15 and AU p < 0.001 provide strong likelihood support, consistent with site-level analysis"

### 2. Ignoring Statistical Significance
❌ **Wrong**: "High decay value means strong support"  
✅ **Right**: "High decay value with significant AU test indicates strong support"

### 3. Not Considering Conflicting Evidence
❌ **Wrong**: "Bayesian analysis supports this clade"
✅ **Right**: "Bayesian analysis supports this clade, but parsimony and likelihood show conflicting patterns"

### 4. Misunderstanding Site Analysis
❌ **Wrong**: "More supporting sites means stronger support"
✅ **Right**: "Site support ratio and weighted ratios together indicate the balance of phylogenetic signal"

## Data Quality Assessment

### Red Flags for Data Quality
- **Very low site ratios** (< 1.2) across most clades
- **High numbers of conflicting sites** relative to supporting sites
- **Large discordance** between analysis types systematically
- **Many non-significant AU tests** despite reasonable decay values

### Recommended Actions for Poor Quality
1. **Alignment inspection**: Check for misaligned regions
2. **Taxon sampling**: Consider adding more taxa or removing problematic sequences
3. **Model selection**: Test alternative evolutionary models
4. **Data partitioning**: Consider gene-by-gene or codon position analyses

## Advanced Interpretation

### Examining Site-Level Patterns
Use individual clade site files (`clade_N_sites.csv`) to:
- **Identify regions** of strong support/conflict
- **Assess whether support is localized** to specific alignment regions
- **Detect potential alignment artifacts** or problematic sites

### Comparative Analysis
When analyzing multiple datasets:
- **Compare site ratios** across datasets for same clades
- **Look for consistent patterns** of support/conflict
- **Identify dataset-specific** vs. universal phylogenetic signal

### Time-Calibrated Interpretations
Consider geological/fossil context:
- **Recent divergences** may have weaker support due to incomplete lineage sorting
- **Ancient divergences** with weak support may indicate rapid radiations
- **Strong support for recent splits** may indicate strong selection or barriers

## Summary Checklist for Interpretation

### Before Drawing Conclusions:
- [ ] Check AU test significance levels
- [ ] Compare results across multiple analysis types  
- [ ] Examine site support ratios and weighted ratios
- [ ] Consider biological plausibility of relationships
- [ ] Assess overall data quality indicators
- [ ] Review potential confounding factors (long branches, missing data)

### For Publication:
- [ ] Report multiple support metrics (not just one)
- [ ] Acknowledge uncertainties where support is weak
- [ ] Discuss biological implications of strongly supported clades
- [ ] Consider alternative explanations for conflicting signal
- [ ] Provide appropriate caveats about data limitations

## Citation

When interpreting and publishing panDecay results, please cite:

McInerney, J.O. (2025) panDecay: Phylogenetic Analysis with Decay Indices. http://github.com/mol-evol/panDecay/