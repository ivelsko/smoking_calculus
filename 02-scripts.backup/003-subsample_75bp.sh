#!/bin/bash

#SBATCH -n 4
#SBATCH --mem 12G
#SBATCH --partition=medium
#SBATCH -o /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.sub75bp.out
#SBATCH -e /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.sub75bp.err
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-use=velsko@shh.mpg.de
#SBATCH -J "sub75bp"

for f in /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/read_len_check/full_files/*.gz
do
	bioawk -c fastx '{if (length($seq) < 76){print "@"$name" "$comment"\n"$seq"\n+\n"$qual}}' $f > /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/read_len_check/sub75bp/$(basename $f .gz).75bp.gz
done

