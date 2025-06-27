# ML/Bayesian Decay Branch Support Analysis Report (v1.0.3)

Date: 2025-06-27 18:50:47

## Analysis Parameters

- Alignment file: `Primate10.nex`
- Data type: `dna`
- Analysis mode: `both`
- ML Model string: `JC`
- PAUP* `lset` command: `lset nst=1 basefreq=empirical rates=equal pinvar=0;`
- Bayesian software: `mrbayes`
- Bayesian model: `JC`
- MCMC generations: `1000000`
- Burnin fraction: `0.25`
- Marginal likelihood method: `harmonic mean`
- Bootstrap analysis: Performed

## Summary Statistics

- ML tree log-likelihood: **6303.660911**
- Number of internal branches tested: 9
- Avg ML log-likelihood difference (constrained vs ML): 29.8614
- Min ML log-likelihood difference: 0.7954
- Max ML log-likelihood difference: 91.1946
- Branches with significant AU support (p < 0.05): 6 / 9 evaluated
- Avg Bayesian decay (marginal lnL difference): -3.8593
- Min Bayesian decay: -5.0350
- Max Bayesian decay: -1.7540
- Branches with strong Bayes Factor support (BF > 10): 0 / 9 evaluated

## Detailed Branch Support Results

| Clade ID | Taxa Count | Constrained lnL | LnL Diff from ML | AU p-value | Significant (AU) | Bayes Decay | Bayes Factor | BF Support | Bootstrap | Included Taxa (sample) |
|----------|------------ |-----------------|------------------|------------|-------------------- |-------------|--------------|------------ |----------- |--------------------------|
| Clade_10 | 3 | 6319.7655 | 16.1046 | 0.0490 | **Yes** | -4.6180 | 0.01 | None | 96 | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta |
| Clade_11 | 2 | 6338.1118 | 34.4509 | 0.0082 | **Yes** | -3.6800 | 0.03 | None | 96 | Macaca_fuscata, Macaca_mulatta |
| Clade_3 | 10 | 6337.5337 | 33.8728 | 0.0001 | **Yes** | -1.7540 | 0.17 | None | 100 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_4 | 9 | 6315.6955 | 12.0346 | 0.0751 | No | -3.5580 | 0.03 | None | 100 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_5 | 5 | 6327.7805 | 24.1195 | 0.0065 | **Yes** | -4.4570 | 0.01 | None | 100 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_6 | 4 | 6314.3441 | 10.6832 | 0.1105 | No | -5.0350 | 0.01 | None | 88 | Gorilla_gorilla, Homo_sapiens, Orangutan... |
| Clade_7 | 3 | 6349.1582 | 45.4973 | 0.0001 | **Yes** | -3.4590 | 0.03 | None | 100 | Gorilla_gorilla, Homo_sapiens, Pan_troglodytes |
| Clade_8 | 2 | 6304.4563 | 0.7954 | 0.6061 | No | -4.1940 | 0.02 | None | 68 | Homo_sapiens, Pan_troglodytes |
| Clade_9 | 4 | 6394.8556 | 91.1946 | 0.0001 | **Yes** | -3.9790 | 0.02 | None | 100 | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta... |

## Interpretation Guide

### ML Analysis
- **LnL Diff from ML**: Log-likelihood of the best tree *without* the clade minus ML tree's log-likelihood. More negative (larger absolute difference) implies stronger support for the clade's presence in the ML tree.
- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.

### Bayesian Analysis
- **Bayes Decay**: Marginal log-likelihood difference between unconstrained and constrained analyses. More negative values indicate stronger support for the clade.
- **Bayes Factor**: Exponential of the Bayes Decay value. Represents the ratio of marginal likelihoods. Common interpretations:
  - BF > 100: Very strong support
  - BF > 10: Strong support
  - BF > 3: Moderate support
  - BF > 1: Weak support
  - BF ≤ 1: No support

- **Bootstrap**: Bootstrap support value (percentage of bootstrap replicates in which the clade appears). Higher values (e.g., > 70) suggest stronger support for the clade.
