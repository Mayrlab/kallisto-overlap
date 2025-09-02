[![DOI](https://zenodo.org/badge/160115858.svg)](https://zenodo.org/doi/10.5281/zenodo.10895236)

# Purpose
This project benchmarks the performance of `kallisto quant` on simulated 3'-end tag-based RNA-seq data. In particular, it focuses on simulating reads from multi-UTR genes and assesses the effectiveness of using truncated transcripts to represent and distinguish 3' UTR isoforms.

The accompanying manuscript is openly available at:

> Fansler, M.M., Mitschka, S. & Mayr, C. Quantifying 3′UTR length from scRNA-seq data reveals changes independent of gene expression. *Nat Commun* **15**, 4050 (2024). [https://doi.org/10.1038/s41467-024-48254-9](https://doi.org/10.1038/s41467-024-48254-9)

# Requirements
 - [conda](https://conda.io/miniconda.html)
 - [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

# Running

    snakemake --use-conda
