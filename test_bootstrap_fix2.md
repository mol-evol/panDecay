# ML-Decay Branch Support Analysis Report (v1.0.3)

Date: 2025-06-27 18:24:24

## Analysis Parameters

- Alignment file: `Primate10.nex`
- Data type: `dna`
- Analysis mode: `ml`
- ML Model string: `JC`
- PAUP* `lset` command: `lset nst=1 basefreq=empirical rates=equal pinvar=0;`

## Summary Statistics

- ML tree log-likelihood: **6303.660911**
- Number of internal branches tested: 9
- Avg ML log-likelihood difference (constrained vs ML): 29.8614
- Min ML log-likelihood difference: 0.7954
- Max ML log-likelihood difference: 91.1946
- Branches with significant AU support (p < 0.05): 6 / 9 evaluated

## Detailed Branch Support Results

| Clade ID | Taxa Count | Constrained lnL | LnL Diff from ML | AU p-value | Significant (AU) | Included Taxa (sample) |
|----------|------------ |-----------------|------------------|------------|-------------------- |--------------------------|
| Clade_10 | 3 | 6319.7655 | 16.1046 | 0.0446 | **Yes** | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta |
| Clade_11 | 2 | 6338.1118 | 34.4509 | 0.0106 | **Yes** | Macaca_fuscata, Macaca_mulatta |
| Clade_3 | 10 | 6337.5337 | 33.8728 | 0.0001 | **Yes** | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_4 | 9 | 6315.6955 | 12.0346 | 0.0639 | No | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_5 | 5 | 6327.7805 | 24.1195 | 0.0086 | **Yes** | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_6 | 4 | 6314.3441 | 10.6832 | 0.1098 | No | Gorilla_gorilla, Homo_sapiens, Orangutan... |
| Clade_7 | 3 | 6349.1582 | 45.4973 | 0.0001 | **Yes** | Gorilla_gorilla, Homo_sapiens, Pan_troglodytes |
| Clade_8 | 2 | 6304.4563 | 0.7954 | 0.6080 | No | Homo_sapiens, Pan_troglodytes |
| Clade_9 | 4 | 6394.8556 | 91.1946 | 0.0001 | **Yes** | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta... |

## Interpretation Guide

### ML Analysis
- **LnL Diff from ML**: Log-likelihood of the best tree *without* the clade minus ML tree's log-likelihood. More negative (larger absolute difference) implies stronger support for the clade's presence in the ML tree.
- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.

