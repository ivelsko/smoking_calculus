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
for sample in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/*.fasta.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

#### NOTSUBDELR ###################################################################
NOTSUBDELR = {}
for notsubdelr in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/malt_checks/sub10M/ELR*unmapped.blastn.sam.gz"):
	NOTSUBDELR[os.path.basename(notsubdelr).split(".u")[0]] = notsubdelr
################################################################################

#### NOTSUBDIVE ###################################################################
NOTSUBDIVE = {}
for notsubdive in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/malt_checks/sub10M/IVE*unmapped.blastn.sam.gz"):
	NOTSUBDIVE[os.path.basename(notsubdive).split(".u")[0]] = notsubdive
################################################################################


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.malt_mapped_readlength.tsv", sample=SAMPLES.keys()),
        expand("{sample}.malt_mapped_gc.tsv", sample=SAMPLES.keys()),
        expand("FastQC/{sample}_MALT_unmapped_fastqc.html", sample=SAMPLES.keys()),
        expand("{notsubdelr}.malt_mapped_readlength.tsv", notsubdelr=NOTSUBDELR.keys()),
        expand("{notsubdelr}.malt_mapped_gc.tsv", notsubdelr=NOTSUBDELR.keys()),
        expand("FastQC/{notsubdelr}_MALT_unmapped_fastqc.html", notsubdelr=NOTSUBDELR.keys()),
        expand("{notsubdive}.malt_mapped_readlength.tsv", notsubdive=NOTSUBDIVE.keys()),
        expand("{notsubdive}.malt_mapped_gc.tsv", notsubdive=NOTSUBDIVE.keys()),
        expand("FastQC/{notsubdive}_MALT_unmapped_fastqc.html", notsubdive=NOTSUBDIVE.keys())

rule malt_map_rl_elr:
    output:
        "{notsubdelr}.malt_mapped_readlength.tsv"
    message: "Run samtools stats RL on malt mapped {wildcards.notsubdelr}"
    group: "samtools"
    params: 
        sam = lambda wildcards: NOTSUBDELR[wildcards.notsubdelr]
    shell:
        """
        samtools stats {params.sam} | grep ^RL | cut -f 2- > {output}
        """

rule malt_map_gc_elr:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{notsubdelr}.unmapped.blastn.sam.gz"
    output:
        "{notsubdelr}.malt_mapped_gc.tsv"
    message: "Run samtools stats GC on {wildcards.notsubdelr}"
    group: "samtools"
    params: 
        sam = lambda wildcards: NOTSUBDELR[wildcards.notsubdelr]
    shell:
        """
        samtools stats {input} | grep ^GCF | cut -f 2- > {output}
        """

rule malt_ummap_rl_elr:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{notsubdelr}.unmapped.blastn.sam.gz"
    output:
        "unmapped_fastq/{notsubdelr}_MALT_unmapped.fastq.gz"
    message: "Run samtools stats RL on {wildcards.notsubdelr}"
    group: "gc_rl"
    params: 
        fastq = "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{notsubdelr}.unmapped.fastq.gz"
    shell:
        """
        zcat {input} | grep -v "^@" | awk -F"\t" '{{print $1}}' | sort | uniq > {wildcards.notsubdelr}.malt_mapped.tsv
        seqkit grep -v -f {wildcards.notsubdelr}.malt_mapped.tsv {params.fastq} > unmapped_fastq/{wildcards.notsubdelr}_MALT_unmapped.fastq
        pigz -p 8 unmapped_fastq/{wildcards.notsubdelr}_MALT_unmapped.fastq
        """

rule malt_ummap_fastqc_elr:
    input:
        "unmapped_fastq/{notsubdelr}_MALT_unmapped.fastq.gz"
    output:
        "FastQC/{notsubdelr}_MALT_unmapped_fastqc.html"
    message: "Run FastQC on unmapped {wildcards.notsubdelr}"
    group: "gc_rl"
    params: 
    threads: 4
    shell:
        """
        fastqc {input} -t {threads} -o FastQC
        """    
        
rule malt_map_rl_ive:
    output:
        "{notsubdive}.malt_mapped_readlength.tsv"
    message: "Run samtools stats RL on malt mapped {wildcards.notsubdive}"
    group: "samtools"
    params: 
        sam = lambda wildcards: NOTSUBDELR[wildcards.notsubdive]
    shell:
        """
        samtools stats {params.sam} | grep ^RL | cut -f 2- > {output}
        """

rule malt_map_gc_ive:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{notsubdive}.unmapped.blastn.sam.gz"
    output:
        "{notsubdive}.malt_mapped_gc.tsv"
    message: "Run samtools stats GC on {wildcards.notsubdive}"
    group: "samtools"
    params: 
        sam = lambda wildcards: NOTSUBDELR[wildcards.notsubdive]
    shell:
        """
        samtools stats {input} | grep ^GCF | cut -f 2- > {output}
        """

rule malt_ummap_rl_ive:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{notsubdive}.unmapped.blastn.sam.gz"
    output:
        "unmapped_fastq/{notsubdive}_MALT_unmapped.fastq.gz"
    message: "Run samtools stats RL on {wildcards.notsubdive}"
    group: "gc_rl"
    params: 
        fastq = "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{notsubdive}.unmapped.fastq.gz"
    shell:
        """
        zcat {input} | grep -v "^@" | awk -F"\t" '{{print $1}}' | sort | uniq > {wildcards.notsubdive}.malt_mapped.tsv
        seqkit grep -v -f {wildcards.notsubdive}.malt_mapped.tsv {params.fastq} > unmapped_fastq/{wildcards.notsubdive}_MALT_unmapped.fastq
        pigz -p 8 unmapped_fastq/{wildcards.notsubdive}_MALT_unmapped.fastq
        """

rule malt_ummap_fastqc_ive:
    input:
        "unmapped_fastq/{notsubdive}_MALT_unmapped.fastq.gz"
    output:
        "FastQC/{notsubdive}_MALT_unmapped_fastqc.html"
    message: "Run FastQC on unmapped {wildcards.notsubdive}"
    group: "gc_rl"
    params: 
    threads: 4
    shell:
        """
        fastqc {input} -t {threads} -o FastQC
        """  

rule malt_ummap_rl_subd:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{sample}.fasta.gz"
    output:
        "unmapped_fastq/{sample}_MALT_unmapped.fastq.gz"
    message: "Run samtools stats RL on {wildcards.sample}"
    group: "gc_rl"
    params: 
        fastq = "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{sample}.unmapped.fastq.gz"
    shell:
        """
        zcat {sample}_MALT_unmapped.fasta.gz | bioawk -c fastx '{histo[length($seq)]++} END {for (l in histo) print l,histo[l]}' | sort -n > {sample}.malt_unmapped_readlength.tsv
        """

rule malt_ummap_rl_subd:
    input:
        "/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/{sample}.fasta.gz"
    output:
        "unmapped_fastq/{sample}_MALT_unmapped.fastq.gz"
    message: "Run samtools stats RL on {wildcards.sample}"
    group: "gc_rl"
    params: 
        fastq = "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{sample}.unmapped.fastq.gz"
    shell:
        """
        zcat {input} | grep "^>" | awk -F"\t" '{{print $1}}' | sort | uniq | sed 's/\>//g' > {wildcards.sample}.malt_mapped.tsv
        seqkit grep -v -f {wildcards.sample}.malt_mapped.tsv {params.fastq} > unmapped_fastq/{wildcards.sample}_MALT_unmapped.fastq
        pigz -p 8 unmapped_fastq/{wildcards.sample}_MALT_unmapped.fastq
        """
      
rule malt_ummap_fastqc_subd:
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
