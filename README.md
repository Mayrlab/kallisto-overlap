[![DOI](https://zenodo.org/badge/160115858.svg)](https://zenodo.org/doi/10.5281/zenodo.10895236)

# Purpose
This project benchmarks the performance of `kallisto quant` on simulated 3'-end tag-based RNA-seq data. In particular, it focuses on simulating reads from multi-UTR genes and assesses the effectiveness of using truncated transcripts to represent and distinguish 3' UTR isoforms.

# Requirements
 - [conda](https://conda.io/miniconda.html)
 - [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

# Running

    snakemake --use-conda
