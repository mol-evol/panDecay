# panDecay ML Branch Support Analysis Report (v1.0.3)

Date: 2025-05-21 15:42:35

## Analysis Parameters

- Alignment file: `Primate10.nex`
- Data type: `dna`
- Model string: `JC`
- PAUP* `lset` command: `lset nst=1 basefreq=empirical rates=equal pinvar=0;`
- Bootstrap analysis: Performed

## Summary Statistics

- ML tree log-likelihood: **6303.660911**
- Number of internal branches tested: 9
- Avg log-likelihood difference (constrained vs ML): 29.8614
- Min log-likelihood difference: 0.7954
- Max log-likelihood difference: 91.1946
- Branches with significant AU support (p < 0.05): 5 / 9 evaluated

## Detailed Branch Support Results

| Clade ID | Taxa Count | Constrained lnL | LnL Diff from ML | AU p-value | Significant (AU) | Bootstrap | Included Taxa (sample) |
|----------|------------|-----------------|------------------|------------|-------------------- |----------- |--------------------------|
| Clade_10 | 3 | 6319.7655 | 16.1046 | 0.0502 | No | 97 | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta |
| Clade_11 | 2 | 6338.1118 | 34.4509 | 0.0070 | **Yes** | 98 | Macaca_fuscata, Macaca_mulatta |
| Clade_3 | 10 | 6337.5337 | 33.8728 | 0.0001 | **Yes** | 100 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_4 | 9 | 6315.6955 | 12.0346 | 0.0694 | No | 95 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_5 | 5 | 6327.7805 | 24.1195 | 0.0090 | **Yes** | 99 | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_6 | 4 | 6314.3441 | 10.6832 | 0.1099 | No | 93 | Gorilla_gorilla, Homo_sapiens, Orangutan... |
| Clade_7 | 3 | 6349.1582 | 45.4973 | 0.0001 | **Yes** | 100 | Gorilla_gorilla, Homo_sapiens, Pan_troglodytes |
| Clade_8 | 2 | 6304.4563 | 0.7954 | 0.5991 | No | 56 | Homo_sapiens, Pan_troglodytes |
| Clade_9 | 4 | 6394.8556 | 91.1946 | 0.0001 | **Yes** | 100 | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta... |

## Interpretation Guide

- **LnL Diff from ML**: Log-likelihood of the best tree *without* the clade minus ML tree's log-likelihood. More negative (larger absolute difference) implies stronger support for the clade's presence in the ML tree.
- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.
- **Bootstrap**: Bootstrap support value (percentage of bootstrap replicates in which the clade appears). Higher values (e.g., > 70) suggest stronger support for the clade.
