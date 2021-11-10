################################################################################
# Run MetaPhlAn3 on MID, CMB, Kilteasheen, Radcliffe, ELR, IVE, JAE samples
#
# Irina Velsko, 23/06/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/metaphlan3"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/metaphlan3/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.metaphlan3.tsv", sample=SAMPLES.keys())

rule metaphlan3:
    output:
        "{sample}.metaphlan3.tsv"
    message: "Run metaphlan3 {wildcards.sample}"
    params: 
#         fastq = "/mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/eager2/eager2_run/samtools/filter/*.gz"
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate mpa3
        set -u

        metaphlan {params.fastq} --input_type fastq --nproc 12 > {output}
        """

        
        