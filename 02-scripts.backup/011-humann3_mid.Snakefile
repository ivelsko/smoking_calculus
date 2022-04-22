################################################################################
# Run HUMAnN3 on MID and CMB data 
#
# Irina Velsko, 11/03/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/humann3/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/humann3/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("output/{sample}.unmapped_genefamilies.tsv", sample=SAMPLES.keys())


rule humann2:
    output:
        "output/{sample}.unmapped_genefamilies.tsv"
    message: "Run humann3 on {wildcards.sample}"
    params: 
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 16
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann3 --input {params.fastq} --output output --threads {threads}
        """

