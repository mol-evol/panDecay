# panDecay Active Context

## Current Project Status

panDecay is a mature phylogenetic analysis tool that calculates decay indices across multiple frameworks (ML, Bayesian, Parsimony). The project recently completed a major enhancement implementing **Effect Size Approach for BD Normalization**.

## Recent Major Enhancement: Effect Size Approach

### Problem Addressed
Traditional Bayesian Decay (BD) values scale with dataset characteristics (alignment length, sequence diversity, substitution rates), making cross-study comparisons impossible and reducing biological interpretability.

### Solution Implemented
**Effect Size Approach** treating BD like Cohen's d statistical effect size:
- **Standard Effect Size**: `BD / SD(site signals)`
- **Robust Effect Size**: `BD / MAD(site signals)` 
- **Weighted Effect Size**: Accounts for site-specific information content

### Key Benefits
- Enables meaningful cross-study comparisons
- Provides standardized interpretation thresholds (Cohen's d framework)
- Accounts for dataset-specific variability
- Maintains biological relevance
- Creates signal-to-noise ratio measurements

## Implementation Status: ✅ COMPLETE

### Code Implementation ✅
- Extended site analysis to compute signal statistics (mean, SD, MAD)
- Added effect size calculation methods with Cohen's d framework
- Updated CLI options with effect size normalization methods
- Enhanced table formatting with dynamic Effect Size columns
- Updated tree annotations to show effect sizes
- Added comprehensive effect size statistics to reports
- Implemented Cohen's d interpretation thresholds

### Documentation ✅
- Added Effect Size Approach section to README.md
- Updated CLI documentation with all effect size options
- Added comprehensive examples (basic, multi-study, meta-analysis)
- Updated Output Files section with effect size documentation
- Enhanced Interpretation Guide with Cohen's d thresholds
- Added Technical Background section explaining statistical foundation

## Current Capabilities

### Analysis Types
- Maximum Likelihood (ML) decay indices
- Bayesian decay indices (MrBayes integration)
- Parsimony decay indices (traditional Bremer support)
- Combined analyses with comprehensive output

### Support Metrics
- Raw Bayesian Decay (BD) values
- **Effect Size normalization (NEW)**
- Per-site normalization
- Relative normalization
- Signal-to-noise ratios
- AU test p-values
- Bootstrap support
- Site-specific analysis

### Advanced Features
- Multi-threaded execution
- MPI and BEAGLE support for Bayesian analyses
- Site-specific likelihood analysis
- Visualization capabilities
- Configuration file support
- Flexible constraint testing

## File Structure
- `panDecay.py`: Main implementation (comprehensive tool)
- `README.md`: Complete documentation with effect size coverage
- `CHANGELOG.md`: Version history
- Example files: `Primate.nex`, `alignment.fas`, etc.
- Output examples: `pan_decay_indices.*` files

## Next Steps / Future Enhancements

### Immediate Opportunities
1. **Testing and Validation**
   - Test effect size calculations on diverse datasets
   - Validate Cohen's d thresholds for phylogenetic applications
   - Compare effect sizes across known phylogenetic relationships

2. **Performance Optimization**
   - Profile effect size calculations for large alignments
   - Optimize site-specific analysis memory usage
   - Parallel processing for normalization calculations

### Medium-term Enhancements
1. **Additional Effect Size Variants**
   - Bias-corrected effect sizes (Hedges' g)
   - Bootstrap confidence intervals for effect sizes
   - Bayesian effect size estimation

2. **Meta-Analysis Support**
   - Effect size aggregation methods
   - Forest plot visualization
   - Publication bias assessment tools

3. **Integration Enhancements**
   - Direct BEAST integration for Bayesian analyses
   - PhyloBayes support
   - RevBayes integration

### Research Applications
1. **Cross-Study Synthesis**
   - Phylogenomic meta-analyses using standardized effect sizes
   - Comparative support assessment across taxonomic groups
   - Method comparison studies using effect size framework

2. **Quality Assessment**
   - Alignment quality metrics based on effect size distributions
   - Model adequacy assessment using normalized support measures
   - Data sufficiency evaluation

## Technical Notes

### Effect Size Implementation Details
- Site-specific likelihood differences calculated via PAUP* constraint analyses
- Signal statistics computed using NumPy statistical functions
- Fallback hierarchy: standard → robust → weighted effect sizes
- Cohen's d thresholds adapted for phylogenetic signal interpretation

### Statistical Foundation
- Based on established Cohen's d framework from psychology/statistics
- Adapted for phylogenetic likelihood differences
- Provides dimensionless, standardized support measures
- Enables meta-analytic synthesis across studies

## Contact & Maintenance
- Primary maintainer: James McInerney
- Language: Python 3.8+
- Dependencies: BioPython, NumPy, PAUP*, MrBayes (optional)
- License: MIT

This context provides the foundation for continued development and application of panDecay's effect size capabilities in phylogenetic research.