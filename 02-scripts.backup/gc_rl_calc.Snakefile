################################################################################
# Calculate GC content and read length of each read in the fasta files that went into MALT
#
# Irina Velsko, 12/02/2021
################################################################################

from glob import glob
import os
import re

workdir: "/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/*.fastq"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.gc_rl.tsv.gz", sample=SAMPLES.keys())

rule gc_rl:
    output:
        "{sample}.gc_rl.tsv.gz"
    message: "Run emboss infoseq on {wildcards.sample}"
    group: "infoseq"
    params: 
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        infoseq -auto -outfile {wildcards.sample}.gc_rl.tsv -only -name -length -pgc {params.fasta}
        pigz -p 8 {wildcards.sample}.gc_rl.tsv
        pigz -p 8 /projects1/microbiome_calculus/smoking_calculus/04-analysis/GC_RL/input/{wildcards.sample}.unmapped.fastq
        """
        