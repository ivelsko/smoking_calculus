#!/usr/bin/env bash

OUTDIR=/projects1/microbiome_calculus/smoking_calculus/04-analysis/sourceTracker/shotgun/rarefied_table

sbatch \
	-c 12 \
	--mem=24000 \
	-o ~/slurm_logs/slurm.%j.out \
	-e ~/slurm_logs/slurm.%j.err \
	--partition=long \
	--mail-type=fail \
	--mail-type=time_limit \
	--mail-user=velsko@shh.mpg.de \
	--export=ALL \
	-J "Sourcetracker" \
	--wrap="Rscript \
	/projects1/users/velsko/bin/sourcetracker-1.0.1/sourcetracker_for_qiime.r \
	-i /projects1/microbiome_calculus/smoking_calculus/04-analysis/sourceTracker/shotgun/shotgun_sourcetracker_table_noblanks_rarefied.txt \
	-m /projects1/microbiome_calculus/smoking_calculus/00-documentation.backup/source_tracker_mappingfile_new.tsv \
	-o "$OUTDIR"/shotgun_sourcetracker_table_rarefied_reduced_sources_log_"$(date +"%Y%m%d")" \
	--train_rarefaction 10000 \
	-r 51500 \
	-v"

