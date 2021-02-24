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
        expand("{sample}.gc_rl.tsv.gz", sample=SAMPLES.keys())

rule malt_map_rl:
    output:
        "{sample}.malt_mapped_readlength.tsv"
    message: "Run samtools stats RL on malt mapped {wildcards.sample}"
    group: "samtools"
    params: 
        sam = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        samtools stats {params.sam}.unmapped.blastn.sam.gz | grep ^RL | cut -f 2- > {output}
        """

rule malt_map_gc:
    input:
        "{sample}.unmapped.blastn.sam.gz"
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
        "{sample}.malt_mapped_gc.tsv"
    output:
        "./unmapped_fasta/{sample}_MALT_unmapped.fasta.gz"
    message: "Run emboss infoseq on {wildcards.sample}"
    group: "infoseq"
    params: 
        rl = lambda wildcards: SAMPLES[wildcards.sample]
        fastq = /projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{sample}.unmapped.fastq.gz
    shell:
        """
        zcat {input} | grep -v "^@" | awk -F"\t" '{print $1}' | sort | uniq > {wildcards.sample}.malt_mapped.tsv
        zcat {params.fastq} | seqkit fq2fa | seqkit grep -v -f {wildcards.sample}.malt_mapped.tsv > {wildcards.sample}_MALT_unmapped.fasta
        pigz -p 4 {wildcards.sample}_MALT_unmapped.fasta
        """

rule malt_ummap_fastqc:
    input:
        "./unmapped_fasta/{sample}_MALT_unmapped.fasta.gz"
    output:
        "./FastQC/{sample}.html"
    message: "Run emboss infoseq on {wildcards.sample}"
    group: "infoseq"
    params: 
        rl = lambda wildcards: SAMPLES[wildcards.sample]
        fastq = /projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{sample}.unmapped.fastq.gz
    shell:
        """
        fastqc {input} -o ./FastQC
        """
      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        