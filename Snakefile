#!/usr/bin/env snakemake

configfile: "config.yaml"


rule all:
    input:
        expand("data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/Rac1.r{replicate}/confusion.txt",
               replicate=list(range(10)),
               distance=[50,100,150,200,250,300,350,400,450,500,550,600,650,700],
               countsProximal=[50,100], countsDistal=[50,100], txWidth=[350,400,450,500,550,600])

rule simulate_two_utr:
    output:
        fastq="data/fastq/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.fq",
        fasta="data/fasta/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.fa"
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
        fasta="data/fasta/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.fa"
    output:
        kdx="data/kdx/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.kdx"
    conda:
        "envs/kallisto.yaml"
    shell:
        """
        kallisto index -i {output.kdx} {input.fasta}
        """

rule kallisto_quant:
    input:
        kdx="data/kdx/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.kdx",
        fastq="data/fastq/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}.fq"
    output:
        tsv="data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}/abundance.tsv",
        bam="data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}/pseudoalignments.bam"
    params:
        outDir="data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}"
    conda:
        "envs/kallisto.yaml"
    shell:
        """
        kallisto quant -i {input.kdx} -o {params.outDir} --single -l1 -s1 --fr-stranded --pseudobam --plaintext {input.fastq}
        """

rule summarize_bam:
    input:
        bam="data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}/pseudoalignments.bam"
    output:
        "data/kallisto/txw{txWidth}/d{distance}/ctp{countsProximal}.ctd{countsDistal}/{gene}.r{replicate}/confusion.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view {input.bam} | cut -f1,3 | \
        awk -F'[.\t]' '{{ print $1, $3 }}' | sort | uniq -c > {output}
        """
