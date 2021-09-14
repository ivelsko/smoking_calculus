#!/bin/bash

#SBATCH -n 4
#SBATCH --mem 12G
#SBATCH --partition=short
#SBATCH -o /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.sub10M.out
#SBATCH -e /projects1/clusterhomes/velsko/slurm_logs/slurm.%j.sub10M.err
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-use=velsko@shh.mpg.de
#SBATCH -J "sub10M"

for f in /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/full_files/*.gz
do
	seqtk sample -s10000 $f 10000000 > /projects1/microbiome_calculus/smoking_calculus/03-preprocessing/depth_check/sub10M/$(basename $f .gz).10M.gz
done

