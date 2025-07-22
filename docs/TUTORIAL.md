# panDecay 101: A Complete Beginner's Guide

Welcome to panDecay! This tutorial will take you from zero to confidently running phylogenetic decay analyses. No prior experience with decay indices is required - we'll build up your understanding step by step.

## Table of Contents

1. [What is Branch Support?](#what-is-branch-support)
2. [Understanding Decay Indices](#understanding-decay-indices)
3. [Your First Analysis](#your-first-analysis)
4. [Understanding Your Results](#understanding-your-results)
5. [Choosing Analysis Methods](#choosing-analysis-methods)
6. [Progressive Examples](#progressive-examples)
7. [Interpreting Statistics](#interpreting-statistics)
8. [Best Practices](#best-practices)
9. [Troubleshooting](#troubleshooting)
10. [Next Steps](#next-steps)

## What is Branch Support?

Imagine you've built a family tree showing how different species are related. But how confident can you be that each branch in your tree is correct? **Branch support** measures tell you how much evidence supports each relationship.

### Why Do We Need Branch Support?

Consider this simple tree:

```
    ┌─ Human
────┤
    └─┬─ Chimp  
      └─ Gorilla
```

This tree suggests humans are more distantly related to chimps and gorillas than those two are to each other. But what if the data actually supports this alternative:

```
    ┌─ Chimp
────┤
    └─┬─ Human
      └─ Gorilla
```

**Branch support measures help you determine which tree the data actually supports**, and how strongly.

### Traditional Approaches

**Bootstrap Support** (you may be familiar with this):

- Randomly resamples your data many times
- Builds a tree from each resampled dataset
- Counts how often each branch appears
- Gives percentages (e.g., "95% bootstrap support")

**Decay Indices** (what panDecay calculates):

- Tests how much worse a tree becomes when you force a specific branch to be absent
- Measures the "cost" of breaking each relationship
- Gives quantitative measures of support strength

### The Key Insight

Decay indices ask: *"If I force this branch to be absent, how much worse does my tree become?"* The larger the cost, the stronger the support for that branch.

## Understanding Decay Indices

### The Basic Concept

Think of decay indices like asking: "How much would I have to pay in 'likelihood points' to make this relationship disappear?"

1. **Start with your best tree** (highest likelihood score)
2. **Force a specific relationship to break** (add a constraint)
3. **Find the best tree with that constraint** (it will be worse)
4. **Calculate the difference** (this is the decay index)



### Types of Decay Indices

panDecay calculates three types:

1. **ML Decay** (Maximum Likelihood)
   - Uses likelihood scores
   - Uses complex models
   - Statistically rigorous
   - **Recommended for most studies**

2. **Bayesian Decay** 
   - Uses Bayesian marginal likelihoods
   - Good for complex evolutionary models
   - Statistically rigorous
   - Computationally intensive

3. **Parsimony Decay** (Traditional Bremer Support)
   - Uses parsimony scores (tree length)
   - Historical standard
   - Fast but less sophisticated than likelihood/Bayesian methods

## Your First Analysis

Let's run your first decay analysis step by step. We'll start simple and build complexity gradually.

### Prerequisites

Before starting, make sure you have:
- panDecay installed (see [Installation Guide](../README.md#installation))
- A sequence alignment file
- Basic familiarity with command-line interfaces

### Step 1: Get Example Data

panDecay comes with example datasets. Let's use the primate dataset:

```bash
# Navigate to your panDecay directory
cd panDecay

# Check what example data is available
ls examples/data/
```

You should see files like `Primate.nex`, `morpho.nex`, and `proteins.phy`.

### Step 2: Your First Command

Let's start with the simplest possible analysis:

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml
```

**What this does:**

- `examples/data/Primate.nex` - Use the example primate dataset
- `--analysis ml` - Perform Maximum Likelihood decay analysis
- Uses default settings (GTR+Gamma model, auto-detect threads)

### Step 3: Watch It Run

You'll see output like:

```
panDecay v1.1.0 - Phylogenetic Decay Analysis
==============================================

Input: examples/data/Primate.nex (nexus format)
Analysis: Maximum Likelihood (ML)
Model: GTR+Gamma
Threads: 4 (auto-detected)

Performing ML tree search...
COMPLETED: ML tree found (log-likelihood: -1234.567)

Generating constraints for 8 internal branches...
Testing Branch 1/8: Clade_1 (4 taxa)
Testing Branch 2/8: Clade_2 (3 taxa)
...
COMPLETED: All constraint analyses completed

Performing AU test for statistical significance...
COMPLETED: AU test completed

Creating organized output directory: 20241015_143022_panDecay_alignment/
Writing results to: 20241015_143022_panDecay_alignment/results/pan_decay_indices.txt
Writing annotated trees to: 20241015_143022_panDecay_alignment/trees/annotated_tree_*.nwk
Writing detailed report to: 20241015_143022_panDecay_alignment/reports/pan_decay_indices.md

Analysis completed successfully!
```

**Congratulations!** You've just completed your first decay analysis.

### Step 4: Find Your Results

Look for the new organized directory structure:
```
20241015_143022_panDecay_alignment/
├── results/
│   └── pan_decay_indices.txt         # Main results table
├── trees/ 
│   ├── annotated_tree_au.nwk         # Tree with AU p-values
│   ├── annotated_tree_delta_lnl.nwk  # Tree with likelihood differences
│   └── annotated_tree_combined.nwk   # Tree with combined annotations
├── reports/
│   └── pan_decay_indices.md          # Human-readable report
└── logs/
    └── panDecay_debug.log            # Analysis log
```

## Understanding Your Results

Let's decode what panDecay found. Navigate to the results directory and open the `pan_decay_indices.txt` file:

```bash
cd 20241015_143022_panDecay_alignment/results/
```

### Main Results Table

```
Clade_ID    Taxa_Count    ΔlnL        AU_pvalue    Significant    Included_Taxa
Clade_1     4            15.234      0.001        ***           Human,Chimp,Gorilla,Orangutan
Clade_2     3            8.456       0.023        *             Human,Chimp,Gorilla  
Clade_3     2            22.789      0.000        ***           Human,Chimp
...
```

### What Each Column Means

**Clade_ID**: panDecay's name for each branch it tested
**Taxa_Count**: How many species are in this group
**ΔlnL**: The decay index (likelihood difference)
**AU_pvalue**: Statistical significance (lower = more significant)
**Significant**: Star rating system
- `***` = Very strong support (p < 0.001)
- `**` = Strong support (p < 0.01) 
- `*` = Moderate support (p < 0.05)
- `ns` = Not significant (p > 0.05)

### Interpreting the Numbers

**ΔlnL (Decay Index)**:
- Higher numbers = stronger support
- `15.234` means the tree got 15.234 log-likelihood units worse when this relationship was forced to be wrong

**AU p-value**:
- Tests if the constrained tree is significantly worse
- `0.001` means there's only a 0.1% chance this strong support occurred by random chance
- Lower p-values = more confident in the support

### Real-World Interpretation

For `Clade_1` above:
> "The relationship grouping Human, Chimp, Gorilla, and Orangutan together has very strong support. Forcing this group to be non-monophyletic makes the tree 15.234 log-likelihood units worse, and this difference is highly statistically significant (p = 0.001)."

### The Detailed Report

Navigate to the reports directory and open the detailed report:

```bash
cd ../reports/
open pan_decay_indices.md
```

This file contains a more complete explanation:

```markdown
# panDecay ML Branch Support Analysis Report

## Summary Statistics
- ML tree log-likelihood: **-1234.567**
- Number of branches tested: 8
- Branches with significant AU support: 6 / 8

## Detailed Results
### Clade_1 (Strong support: ML)
**Taxa**: Human, Chimp, Gorilla, Orangutan
**Support Evidence**: 
- ΔlnL = 15.234 (strong likelihood preference)
- AU p-value = 0.001 (highly significant)
- **Interpretation**: This grouping has very strong support
```

This report explains your results in plain language with biological context.

## Choosing Analysis Methods

panDecay offers three analysis frameworks. Here's when to use each:

### Maximum Likelihood (ML) - **RECOMMENDED FOR MOST STUDIES**

```bash
python3 panDecay.py alignment.fas --analysis ml
```

**Use when:**
- You want statistically rigorous results
- You have moderate to large datasets (>20 taxa)
- You want to use sophisticated evolutionary models
- You need results comparable to modern phylogenetic studies

**Advantages:**
- Most statistically sound
- Fast for most datasets
- Well-established interpretation
- Compatible with complex evolutionary models

### Bayesian Analysis

```bash
python3 panDecay.py alignment.fas --analysis bayesian
```

**Use when:**
- You want to incorporate prior information
- You have very complex evolutionary models
- You need to account for uncertainty in model parameters
- You have computational resources for longer analyses

**Advantages:**
- Accounts for parameter uncertainty
- Can incorporate prior knowledge
- Provides comprehensive uncertainty estimates

**Disadvantages:**
- Much slower than ML
- More complex to interpret
- Requires more computational resources

### Parsimony (Traditional Bremer Support)

```bash
python3 panDecay.py alignment.fas --analysis parsimony
```

**Use when:**
- You want results comparable to older studies
- You have morphological or discrete character data
- You need very fast analyses
- You're working with small datasets

**Advantages:**
- Very fast
- Simple interpretation
- Historical standard (comparable to literature)
- Good for morphological data

**Disadvantages:**
- Less statistically sophisticated
- Makes strong assumptions about evolution
- May be misleading for molecular data

### Combined Analyses

You can run multiple methods together:

```bash
# ML + Bayesian
python3 panDecay.py alignment.fas --analysis ml+bayesian

# All three methods
python3 panDecay.py alignment.fas --analysis all
```

**Use combined analyses when:**
- You want to compare methods
- You need comprehensive support measures
- You're working on an important study requiring multiple lines of evidence

## Progressive Examples

Let's build your skills with increasingly complex analyses.

### Example 1: Basic ML Analysis (START HERE)

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml
```

**What you learn:**
- Basic command structure
- Reading results
- Understanding decay indices

### Example 2: Add Site-Specific Visualization

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml --site-analysis --create-alignment-viz
```

**New concepts:**
- `--site-analysis` performs site-specific likelihood analysis
- `--create-alignment-viz` creates alignment visualization with support overlays
- Visual representation of supporting/conflicting sites for each branch

**Check your results:**
- Navigate to your timestamped results directory
- Look in `visualizations/` for `.png` distribution plots
- Look in `site_analysis/` for site-specific plots and histograms
- All plots are publication-ready static images

### Example 3: Protein Data with Model Selection

```bash
python3 panDecay.py examples/data/proteins.phy --analysis ml --model AUTO --data-type protein
```

**New concepts:**
- `--data-type protein` tells panDecay this is protein data
- `--model AUTO` automatically selects the best protein evolutionary model
- Different data types require different approaches

### Example 4: Adding Statistical Rigor

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml --bootstrap --bootstrap-reps 100
```

**New concepts:**
- `--bootstrap` adds traditional bootstrap support
- `--bootstrap-reps 100` runs 100 bootstrap replicates
- Comparing decay indices with bootstrap support

**Interpretation tip:**
- Strong branches should have both high decay indices AND high bootstrap support
- Disagreement between methods suggests caution

### Example 5: Site-by-Site Analysis

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml --site-analysis --viz-format both
```

**New concepts:**
- `--site-analysis` examines which parts of your alignment support each branch
- Identifies conflicting signal in your data
- Creates detailed site-specific plots

**What to look for:**
- Sites that support vs. conflict with each branch
- Regions of alignment with strong phylogenetic signal
- Potential alignment problems or recombination

### Example 6: Configuration Files for Complex Analyses

First, generate a configuration template:

```bash
python3 panDecay.py --generate-yaml-config my_analysis.yaml
```

Edit `my_analysis.yaml` to customize your analysis, then run:

```bash
python3 panDecay.py examples/data/Primate.nex --config my_analysis.yaml
```

**New concepts:**
- Configuration files for reproducible analyses
- YAML format for complex parameter specification
- Saving and sharing analysis settings

### Example 7: Performance Optimization

```bash
python3 panDecay.py examples/data/Primate.nex --analysis ml --async-constraints --max-async-workers 4 --threads 8
```

**New concepts:**
- `--async-constraints` runs constraint analyses in parallel
- `--max-async-workers 4` uses 4 parallel processes
- `--threads 8` uses 8 threads for each analysis
- Performance optimization for large datasets

## Interpreting Statistics

Understanding the statistics panDecay provides is crucial for proper interpretation.

### Decay Index Values (ΔlnL)

**What they measure:**
- How much likelihood decreases when a branch is forced to be wrong
- Larger values = stronger support

**Typical ranges:**
- `0-5`: Weak support
- `5-15`: Moderate support  
- `15+`: Strong support
- `50+`: Very strong support

**Important caveat:**
- These ranges depend on your dataset size and evolutionary divergence
- Always compare values within your study, not across studies
- Use statistical tests (AU p-values) for significance

### AU Test P-values

**What they test:**
- Whether the constrained tree is significantly worse than the ML tree
- H₀: "The constrained tree is not significantly worse"
- H₁: "The constrained tree is significantly worse"

**Interpretation:**
- `p < 0.001`: Very strong evidence against the constraint (strong support for branch)
- `p < 0.01`: Strong evidence against the constraint
- `p < 0.05`: Moderate evidence against the constraint
- `p > 0.05`: No significant evidence (weak or no support)

**Important notes:**
- Lower p-values = stronger support
- This is a statistical test - consider both effect size (ΔlnL) and significance (p-value)
- Multiple testing: with many branches, some low p-values may occur by chance

### Support Level Interpretation

panDecay uses this classification system:

```
***  = p < 0.001  = Very strong support
**   = p < 0.01   = Strong support  
*    = p < 0.05   = Moderate support
ns   = p > 0.05   = Not significant (weak/no support)
```

### Combining Multiple Lines of Evidence

**Strong branch support** should show:
- High decay index (ΔlnL > 10)
- Low AU p-value (p < 0.01)
- High bootstrap support (if calculated, >70%)
- Consistent support across different methods

**Weak or questionable support:**
- Low decay index (ΔlnL < 5)
- High AU p-value (p > 0.05)
- Low bootstrap support (<50%)
- Disagreement between methods

### Dataset-Relative Normalization

For within-study comparisons, panDecay offers normalized metrics:

```bash
python3 panDecay.py alignment.fas --normalization full
```

This provides:
- **Percentile ranks**: "This branch has stronger support than 85% of branches in this tree"
- **Z-scores**: "This branch is 2.3 standard deviations above average support"
- **Dataset-relative values**: "This branch has 0.8 relative support (0=weakest, 1=strongest in this dataset)"

## Best Practices

### Study Design

**Before you start:**
1. **Check your alignment quality**
   - Remove poorly aligned regions
   - Verify no obvious errors or contamination
   - Consider alignment trimming tools

2. **Choose appropriate models**
   - Use model selection (AIC, BIC) for DNA data
   - Consider `--model AUTO` for proteins
   - Don't use overly simple models for complex data

3. **Plan your computational resources**
   - Larger datasets need more time and memory
   - Consider `--async-constraints` for datasets with >10 taxa
   - Plan for multiple hours for Bayesian analyses

### Running Analyses

**Start simple:**
```bash
# Begin with basic ML analysis
python3 panDecay.py alignment.fas --analysis ml

# Add complexity gradually
python3 panDecay.py alignment.fas --analysis ml --bootstrap --viz-format both
```

**For important studies:**
```bash
# Run multiple methods for comparison
python3 panDecay.py alignment.fas --analysis all --bootstrap --site-analysis --viz-format both
```

**Performance optimization:**
```bash
# For large datasets
python3 panDecay.py alignment.fas --analysis ml --async-constraints --max-async-workers 4 --threads auto
```

### Interpreting Results

**Focus on biological questions:**
- Which relationships are well-supported?
- Which relationships need more data?
- Are there conflicting signals in your data?

**Be conservative:**
- Don't over-interpret borderline results
- Consider multiple lines of evidence
- Report uncertainty honestly

**Compare with other studies:**
- How do your results compare to previous work?
- Are strongly supported relationships consistent across studies?
- Where do you find novel or conflicting results?

### Reporting Results

**In your methods:**
```
"Branch support was assessed using ML-based decay indices implemented in panDecay v1.1 (citation). We performed maximum likelihood searches under the GTR+Γ model, followed by constrained analyses for each internal branch. Support was evaluated using approximately unbiased (AU) tests, with p < 0.05 considered significant."
```

**In your results:**
```
"The relationship between X and Y received strong support (ΔlnL = 15.2, AU p < 0.001), while the position of Z was poorly supported (ΔlnL = 2.1, AU p = 0.23)."
```

**In figures:**
- Annotate trees with both decay indices and AU significance
- Use consistent color coding (e.g., red = non-significant, green = significant)
- Include scale bars and model information

## Troubleshooting

### Common Issues and Solutions

**Problem**: "PAUP* not found"
```
Solution: Install PAUP* and ensure it's in your PATH
# Download from: https://phylosolutions.com/paup-test/
# Add to PATH or use --paup-path option
```

**Problem**: Analysis takes forever
```
Solutions:
1. Use --async-constraints for parallel processing
2. Reduce dataset size for testing
3. Check if alignment has errors
4. Consider simpler models for initial tests
```

**Problem**: All branches show weak support
```
Possible causes:
1. Insufficient phylogenetic signal (need more data)
2. Conflicting signal (recombination, incomplete lineage sorting)
3. Alignment problems
4. Inappropriate evolutionary model

Solutions:
1. Try --site-analysis to understand signal distribution
2. Check alignment quality
3. Consider different molecular markers
4. Use model selection tools
```

**Problem**: Results don't match bootstrap analysis
```
This is not necessarily a problem:
1. Different methods can give different results
2. Decay indices are more sensitive to likelihood differences
3. Bootstrap tests different aspects of uncertainty
4. Report both and discuss differences
```

**Problem**: "Constraint analysis failed"
```
Solutions:
1. Check that constraint is biologically meaningful
2. Verify taxa names match alignment exactly
3. Try simpler constraints first
4. Check PAUP* log files for specific errors
```

### Performance Optimization

**For large datasets:**
```bash
# Use all available optimizations
python3 panDecay.py large_alignment.fas \
  --analysis ml \
  --async-constraints \
  --max-async-workers 8 \
  --threads auto \
  --constraint-timeout 3600
```

**For quick testing:**
```bash
# Minimal analysis for testing
python3 panDecay.py alignment.fas \
  --analysis ml \
  --constraint-mode specific \
  --constraints "taxon1,taxon2;taxon3,taxon4"
```

### Getting Help

**Check these in order:**
1. **Error messages** - panDecay provides detailed error descriptions
2. **Log files** - Check temporary directory for PAUP*/MrBayes logs
3. **Example data** - Test with provided examples to isolate issues
4. **Documentation** - Check README.md and other guides
5. **GitHub Issues** - Search existing issues or create new one

**When reporting bugs:**
- Include your exact command
- Include error messages
- Mention your operating system
- Provide a minimal example that reproduces the problem

## Next Steps

Congratulations! You now understand the basics of decay index analysis. Here's how to continue building your expertise:

### Immediate Next Steps

1. **Practice with your own data**
   - Start with a small, well-understood dataset
   - Compare results with published bootstrap values
   - Try different analysis methods

2. **Explore advanced features**
   - Site-specific analysis for understanding conflicting signal
   - Combined analyses for comprehensive support assessment
   - Configuration files for reproducible research

3. **Read the literature**
   - Bremer (1994) - Original parsimony decay indices
   - Shimodaira (2002) - AU test methodology
   - Recent papers using ML-based decay indices

### Advanced Topics

**For complex evolutionary scenarios:**
- Multiple gene analyses
- Accounting for incomplete lineage sorting
- Dealing with horizontal gene transfer
- Time-calibrated phylogenies

**For methodological development:**
- Comparing support measures across methods
- Understanding the relationship between different support metrics
- Developing new approaches for specific data types

**For large-scale studies:**
- Phylogenomic datasets
- Automated pipeline development
- High-performance computing considerations

### Additional Resources

**Documentation:**
- [README.md](../README.md) - Comprehensive reference
- [DEVELOPER_GUIDE.md](DEVELOPER_GUIDE.md) - For extending panDecay
- [API_REFERENCE.md](API_REFERENCE.md) - For programmatic use

**Example datasets:**
- `examples/data/` - Practice datasets with different characteristics
- `examples/configs/` - Example configuration files

**Community:**
- GitHub Issues - Questions and bug reports
- GitHub Discussions - General questions and sharing

### Final Thoughts

Phylogenetic analysis is both an art and a science. Decay indices provide powerful tools for understanding branch support, but they're just one piece of the puzzle. Always consider:

- **Multiple lines of evidence** - No single metric tells the whole story
- **Biological context** - Do your results make biological sense?
- **Data quality** - Good analyses require good data
- **Method limitations** - Every approach has strengths and weaknesses

Most importantly, **don't be afraid to experiment**! panDecay is designed to be flexible and forgiving. Try different settings, compare methods, and build your intuition about what works best for your research questions.

Good luck with your phylogenetic analyses!