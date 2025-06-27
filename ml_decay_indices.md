# ML-Decay Branch Support Analysis Report (v1.0.3)

Date: 2025-06-26 21:33:16

## Analysis Parameters

- Alignment file: `Primate10.nex`
- Data type: `dna`
- Model string: `JC`
- PAUP* `lset` command: `lset nst=1 basefreq=equal rates=equal pinvar=0;`

## Summary Statistics

- ML tree log-likelihood: **6424.202404**
- Number of internal branches tested: 9
- Avg log-likelihood difference (constrained vs ML): 32.3840
- Min log-likelihood difference: 2.1541
- Max log-likelihood difference: 97.3797
- Branches with significant AU support (p < 0.05): 6 / 9 evaluated

## Detailed Branch Support Results

| Clade ID | Taxa Count | Constrained lnL | LnL Diff from ML | AU p-value | Significant (AU) | Included Taxa (sample) |
|----------|------------|-----------------|------------------|------------|-------------------- |--------------------------|
| Clade_10 | 3 | 6442.0396 | 17.8372 | 0.0320 | **Yes** | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta |
| Clade_11 | 2 | 6462.1394 | 37.9370 | 0.0050 | **Yes** | Macaca_fuscata, Macaca_mulatta |
| Clade_3 | 10 | 6460.0404 | 35.8380 | 0.0001 | **Yes** | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_4 | 9 | 6437.0260 | 12.8236 | 0.0720 | No | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_5 | 5 | 6452.7837 | 28.5813 | 0.0034 | **Yes** | Gibbon, Gorilla_gorilla, Homo_sapiens... |
| Clade_6 | 4 | 6436.2025 | 12.0001 | 0.0981 | No | Gorilla_gorilla, Homo_sapiens, Orangutan... |
| Clade_7 | 3 | 6471.1075 | 46.9051 | 0.0001 | **Yes** | Gorilla_gorilla, Homo_sapiens, Pan_troglodytes |
| Clade_8 | 2 | 6426.3565 | 2.1541 | 0.4913 | No | Homo_sapiens, Pan_troglodytes |
| Clade_9 | 4 | 6521.5821 | 97.3797 | 0.0001 | **Yes** | Macaca_fascicularis, Macaca_fuscata, Macaca_mulatta... |

## Interpretation Guide

- **LnL Diff from ML**: Log-likelihood of the best tree *without* the clade minus ML tree's log-likelihood. More negative (larger absolute difference) implies stronger support for the clade's presence in the ML tree.
- **AU p-value**: P-value from the Approximately Unbiased test comparing the ML tree against the alternative (constrained) tree. Lower p-values (e.g., < 0.05) suggest the alternative tree (where the clade is broken) is significantly worse than the ML tree, thus supporting the clade.
