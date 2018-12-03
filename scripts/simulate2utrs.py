#!/usr/bin/env python
from numpy.random import seed
from scimulator.fastq import *

wcs = snakemake.wildcards

position = [0, int(wcs.distance)]
counts = [int(wcs.countsDistal), int(wcs.countsProximal)]
mu = [float(snakemake.params.mu)]*2
sd = [float(snakemake.params.sd)]*2
readLength = int(snakemake.params.readLength)

if not snakemake.params.seq:
    seed(9102)
    reads, txs = generate_reads(position, counts, mu, sd, readLength)
else:
    reads, txs = generate_reads(position, counts, mu, sd, readLength, seq=snakemake.params.seq)

txs_to_fasta(txs, snakemake.output.fasta, truncate=int(wcs.txWidth))
reads_to_fastq(reads, snakemake.output.fastq)
