# admixcov_data

Variance decomposition analysis of the UK ([Patterson et al. 2022](https://doi.org/10.1038/s41586-021-04287-4)) and Bohemia ([Papac et al. 2021](https://doi.org/10.1126/sciadv.abi6941)) human ancient DNA datasets associated with the paper "The contribution of admixture, selection, and genetic drift to four thousand years of human allele frequency change" by Alexis Simon and Graham Coop.

Depends on:
- `snakemake`
- `conda`

To run this workflow:
```
snakemake -c {threads} --use-conda
```
