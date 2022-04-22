################################################################################
# Run MetaPhlAn3 on MID, CMB, Kilteasheen, Radcliffe, ELR, IVE, JAE samples
# subsampled to no more than 10M reads
# Irina Velsko, 01/04/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/metaphlan3"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("sub10M/{sample}.sub10M.metaphlan3.tsv", sample=SAMPLES.keys())

rule metaphlan3:
    output:
        "sub10M/{sample}.sub10M.metaphlan3.tsv"
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

        
        