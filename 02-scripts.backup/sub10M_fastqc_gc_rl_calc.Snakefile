################################################################################
# Calculate GC content and read length of each read in the fasta files that went into MALT
#
# Irina Velsko, 12/02/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/FastQC/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/VLC*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.unmapped.10M_fastqc.zip", sample=SAMPLES.keys())

rule sub10M:
    output:
        "{sample}.unmapped.10M_fastqc.zip"
    message: "Run FASTQC on sub10M {wildcards.sample}"
    group: "infoseq"
    params: 
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        fastqc {params.fasta} -o /mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/FastQC
        """
                