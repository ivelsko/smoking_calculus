#!/usr/bin/env bash
sbatch \
-c 112 \
--mem 1850000 \
--partition=supercruncher \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=velsko@shh.mpg.de \
-J "maltnt" \
--wrap="/projects1/microbiome_calculus/smoking_calculus/02-scripts.backup/007-malt-genbank-nt_2017_2step_85pc_supp_0.01 \
/projects1/microbiome_calculus/smoking_calculus/03-preprocessing/eager2/eager2_run/samtools/filter/*.gz \
/projects1/microbiome_calculus/smoking_calculus/04-analysis/malt_checks/nt"

