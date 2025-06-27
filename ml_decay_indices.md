# Bayesian Decay Branch Support Analysis Report (v1.0.3)

Date: 2025-06-27 21:34:48

## Analysis Parameters

- Alignment file: `Primate10.nex`
- Data type: `dna`
- Analysis mode: `bayesian`
- Bayesian software: `mrbayes`
- Bayesian model: `JC`
- MCMC generations: `100000`
- Burnin fraction: `0.25`
- Marginal likelihood method: `stepping-stone`

## Summary Statistics

- Number of internal branches tested: 10
- Avg Bayesian decay (marginal lnL difference): -4.0647
- Min Bayesian decay: -6.2013
- Max Bayesian decay: -0.9072

**⚠️ WARNING**: 10/10 branches have negative Bayes Decay values.
This suggests potential issues with:
- MCMC convergence (consider increasing --bayes-ngen)
- Harmonic mean estimation reliability
- Model specification

- Branches with strong Bayes Factor support (BF > 10): 0 / 10 evaluated

## Detailed Branch Support Results

| Clade ID | Taxa Count | Bayes Decay | Bayes Factor | BF Support | Included Taxa (sample) |
|----------|------------ |-------------|--------------|------------ |--------------------------|
| Clade_10 | 4 | -5.1806 | 0.01 | None | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta... |
| Clade_11 | 3 | -6.1194 | 0.00 | None | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta |
| Clade_12 | 2 | -4.6907 | 0.01 | None | Macaca_fuscata, Macaca_mulatta |
| Clade_3 | 11 | -0.9072 | 0.40 | None | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_4 | 10 | -3.1414 | 0.04 | None | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_5 | 9 | -5.5651 | 0.00 | None | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_6 | 5 | -6.2013 | 0.00 | None | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_7 | 4 | -3.4735 | 0.03 | None | Gorilla_gorilla, Homo_sapiens, Orangutan... |
| Clade_8 | 3 | -4.0080 | 0.02 | None | Gorilla_gorilla, Homo_sapiens, Pan_troglodytes |
| Clade_9 | 2 | -1.3603 | 0.26 | None | Homo_sapiens, Pan_troglodytes |

## Interpretation Guide

### Bayesian Analysis
- **Bayes Decay**: Marginal log-likelihood difference (unconstrained - constrained). **Positive values** indicate support for the clade, with larger positive values indicating stronger support. **Negative values** suggest the constrained analysis had higher marginal likelihood, which may indicate:
  - Poor chain convergence or insufficient MCMC sampling
  - Complex posterior distribution requiring more steps
  - Genuine lack of support for the clade
- **Bayes Factor**: Exponential of the Bayes Decay value. Represents the ratio of marginal likelihoods. Common interpretations:
  - BF > 100: Very strong support
  - BF > 10: Strong support
  - BF > 3: Moderate support
  - BF > 1: Weak support
  - BF ≤ 1: No support

