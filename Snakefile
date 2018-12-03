#!/usr/bin/env snakemake

configfile: "config.yaml"


rule all:
    input:
        expand("data/kallisto/Rac1.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}/abundance.tsv",
               distance=[50,100,150,200,250,300,350,400,450,500,550,600,650,700],
               countsProximal=[50,100], countsDistal=[50,100], txWidth=[350,400,450,500,550,600])

rule simulate_two_utr:
    output:
        fastq="data/fastq/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.fq",
        fasta="data/fasta/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.fa"
    params:
        seq=lambda wcs: None if wcs.gene == "random" else config['seq'][wcs.gene],
        mu=config['peaks']['mu'],
        sd=config['peaks']['sd'],
        readLength=config['readLength']
    conda:
        "envs/kallisto.yaml"
    script:
        "scripts/simulate2utrs.py"

rule kallisto_index:
    input:
        fasta="data/fasta/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.fa"
    output:
        kdx="data/kdx/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.kdx"
    conda:
        "envs/kallisto.yaml"
    shell:
        """
        kallisto index -i {output.kdx} {input.fasta}
        """

rule kallisto_quant:
    input:
        kdx="data/kdx/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.kdx",
        fastq="data/fastq/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}.fq"
    output:
        tsv="data/kallisto/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}/abundance.tsv",
        bam="data/kallisto/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}/pseudoalignments.bam"
    params:
        outDir="data/kallisto/{gene}.d{distance}.ctp{countsProximal}.ctd{countsDistal}.txw{txWidth}"
    conda:
        "envs/kallisto.yaml"
    shell:
        """
        kallisto quant -i {input.kdx} -o {params.outDir} --single -l1 -s1 --fr-stranded --pseudobam --plaintext {input.fastq}
        """

