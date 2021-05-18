################################################################################
# Calculate GC content and read length of each read in the fasta files that went into MALT
#
# Irina Velsko, 12/02/2021
################################################################################

from glob import glob
import os
import re

workdir: "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/output"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.malt_mapped_readlength.tsv", sample=SAMPLES.keys()),
        expand("{sample}.malt_mapped_gc.tsv", sample=SAMPLES.keys()),
        expand("FastQC/{sample}_MALT_unmapped_fastqc.html", sample=SAMPLES.keys())

rule malt_map_rl:
    output:
        "{sample}.malt_mapped_readlength.tsv"
    message: "Run samtools stats RL on malt mapped {wildcards.sample}"
    group: "samtools"
    params: 
        sam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        samtools stats {params.sam} | grep ^RL | cut -f 2- > {output}
        """

rule malt_map_gc:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{sample}.unmapped.blastn.sam.gz"
    output:
        "{sample}.malt_mapped_gc.tsv"
    message: "Run samtools stats GC on {wildcards.sample}"
    group: "samtools"
    params: 
        sam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        samtools stats {input} | grep ^GCF | cut -f 2- > {output}
        """

rule malt_ummap_rl:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{sample}.unmapped.blastn.sam.gz"
    output:
        "unmapped_fastq/{sample}_MALT_unmapped.fastq.gz"
    message: "Run samtools stats RL on {wildcards.sample}"
    group: "gc_rl"
    params: 
        fastq = "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{sample}.unmapped.fastq.gz"
    shell:
        """
        zcat {input} | grep -v "^@" | awk -F"\t" '{{print $1}}' | sort | uniq > {wildcards.sample}.malt_mapped.tsv
        seqkit grep -v -f {wildcards.sample}.malt_mapped.tsv {params.fastq} > unmapped_fastq/{wildcards.sample}_MALT_unmapped.fastq
        pigz -p 8 unmapped_fastq/{wildcards.sample}_MALT_unmapped.fastq
        """

rule malt_ummap_fastqc:
    input:
        "unmapped_fastq/{sample}_MALT_unmapped.fastq.gz"
    output:
        "FastQC/{sample}_MALT_unmapped_fastqc.html"
    message: "Run FastQC on unmapped {wildcards.sample}"
    group: "gc_rl"
    params: 
    threads: 4
    shell:
        """
        fastqc {input} -t {threads} -o FastQC
        """    
        
        