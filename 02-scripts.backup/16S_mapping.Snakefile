################################################################################
# Run SourceTracker on MID/CMB/Radcliffe/Kilteasheen/IVE/ELR calculus sequence data
#
# Irina Velsko, 18/17/2021
################################################################################

from glob import glob
import os
import re

workdir: "/projects1/microbiome_calculus/smoking_calculus/04-anaysis/sourceTracker/16S_mapping"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/projects1/microbiome_calculus/smoking_calculus/04-anaysis/sourceTracker/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".f")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.16S.gz", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("{sample}.sai")
    message: "Align sample {wildcards.sample} against the SILVA database using BWA aln"
    group: "bwa"
    params: 
        reffa = "/projects1/microbiome_sciences/reference_databases/source/silva/release_128_DNA/silva128_thymines.fasta",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 4
    shell:
        """
        bwa aln -n 0.02 -l 1024 -t {threads} \ # -n 0.01 -l 32 use these parameters b/c they're in Zandra's file from James and we want to use her sources
            {params.reffa} \
            {params.fastq} > {output}
        """

rule bwa_samse:
    input:
        "{sample}.sai"
    output:
        temp("{sample}.16Smapped.bam")
    message: "Generate alignment file for sample {wildcards.sample} against the SILVA database"
    group: "bwa"
    params:
        reffa = "/projects1/microbiome_sciences/reference_databases/source/silva/release_128_DNA/silva128_thymines.fasta",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bwa samse \
            {params.reffa} \
            {input} \
            {params.fastq} | \
        samtools view -Sb -F 4 - | samtools sort - -o {output}
        """

rule samtools_fasta:
    input:
        "{sample}.16Smapped.bam"
    output:
        "{sample}.16Smapped.fa"
    message: "Convert the 16S-mapped reads back into a fasta file"
    group: "bwa"
    params:
    shell:
        """
        samtools fasta -s {input} {output}
        """
        
rule samtools_header_fix:
    input:
        "{sample}.16Smapped.fa"
    output:
        "{sample}.16S.fa"
    message: "Convert the 16S-mapped reads back into a fasta file"
    group: "bwa"
    params:
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        awk -v name="$HEADER" '/^>/{print ">"name "_" ++i; next}{print}' {input} > {output}
        pigz -p 4 {output}
        """
        
