#!/usr/bin/env bash

#SBATCH -c 4
#SBATCH --mem=32000
#SBATCH --partition=short
#SBATCH -o /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.out
#SBATCH -e /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.err
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-use=velsko@shh.mpg.de
#SBATCH --array=0-71%4
#SBATCH -J "10M_fastqc"

SAMPLES=( $(find /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/ -name '*.gz' -type f) )
SAMPLENAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

fastqc "${SAMPLENAME}" -o /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/FastQC

