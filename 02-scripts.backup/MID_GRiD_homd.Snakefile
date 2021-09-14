################################################################################
# Determine differences in replication activity between species in the HOMD database
# using GRiD to check for differences between heavy & light/non-smokers
#
# Irina Velsko 04/02/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/GRiD/HOMD"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/eager2/eager2_run/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
         "output/done_file.txt"

rule grid:
    output: touch("output/done_file.txt")
    message: "Run GRiD with the HOMD database"
    params: 
        reads = "mnt/archgen/microbiome_calculus/smoking_calculus/03-preprocessing/eager2/eager2_run/samtools/filter",
        out_folder = "/mnt/archgen/microbiome_calculus/smoking_calculus/04-analysis/GRiD/HOMD/output",
        dir = "/mnt/archgen/microbiome_sciences/reference_databases/built/HOMD_2021_07_29"
    threads: 24
    shell:
        """
	set +u
	source $HOME/miniconda3/etc/profile.d/conda.sh
	conda activate GRiD
	set -u
        grid multiplex -r /{params.reads} \
        -e fastq.gz \
        -o {params.out_folder} \
        -d {params.dir} \
        -m \
        -n {threads}
        """

