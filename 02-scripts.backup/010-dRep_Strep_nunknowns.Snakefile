################################################################################
# Run dRep on Strep genomes of unknown and known Strep groups
# to place the unknowns in groups based on ANI, if possible
#
# Irina Velsko, 23/06/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/abpCapture/01-data.backup/probe_design/reference_genomes/dRep"


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "dRep_ANImf_out/data_tables/Ndb.csv"
#         "dRep_gANI_out/data_tables/Ndb.csv"

rule dRep:
    output:
        "dRep_ANImf_out/data_tables/Ndb.csv"
    message: "Run dRep to cluster Strep genomes"
    params: 
        genomes = "/mnt/archgen/microbiome_calculus/abpCapture/01-data.backup/probe_design/reference_genomes/dRep/Strep_genomes_for_dRep.tsv"
    threads: 24
    shell:
        """
        dRep compare -p {threads} dRep_ANImf_out/ -g {params.genomes} --S_algorithm ANImf -pa 0.95 -sa 0.99
        """

# rule dRep:
#     output:
#         "dRep_gANI_out/data_tables/Ndb.csv"
#     message: "Run dRep to cluster Strep genomes"
#     params: 
#         genomes = "/mnt/archgen/microbiome_calculus/abpCapture/01-data.backup/probe_design/reference_genomes/dRep/Strep_genomes_for_dRep.tsv"
#     threads: 24
#     shell:
#         """
#         dRep compare -p {threads} dRep_gANI_out/ -g {params.genomes} --S_algorithm gANI -pa 0.95 -sa 0.99
#         """
     
        