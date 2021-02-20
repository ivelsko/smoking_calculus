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
-J "Malt75bp" \
--wrap="/projects1/microbiome_calculus/evolution/02-scripts.backup/008-malt-genbank-refseq_bacarchhomo_gcs_20181122_4step_85pc_supp_0.01 \
/projects1/microbiome_calculus/smoking_calculus/03-preprocessing/read_len_check/sub75bp/*.gz \
/projects1/microbiome_calculus/smoking_calculus/04-analysis/malt_checks/sub75bp"

