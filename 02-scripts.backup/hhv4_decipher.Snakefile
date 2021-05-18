################################################################################
# Project: Smoking calculus
# Part: HHV4 presence
# Step: Sequence aigning and clustering
#
# Irina Velsko, 23/02/2020
################################################################################

from snakemake.utils import R

workdir: "/projects1/microbiome_calculus/smoking_calculus/01-data/hhv4"

# /projects1/clusterhomes/velsko/R/x86_64-pc-linux-gnu-library/3.5/

rule all:
    input: 
        "hhv4_extraction_DECIPHER_aligned.fasta",
        "hhv4_clustersML_08.nwk",
        "hhv4_clustersML_05.nwk",
        "hhv4_clustersML_02.nwk",
        "hhv4_clustersML_08.tsv",
        "hhv4_clustersML_05.tsv",
        "hhv4_clustersML_02.tsv"

rule hhv4_align_cluster:
    output:
        hhv4_aln = "hhv4_extraction_DECIPHER_aligned.fasta",
        hhv4_clustersML_08 = "hhv4_clustersML_08.tsv",
        hhv4_clustersML_05 = "hhv4_clustersML_05.tsv",
        hhv4_clustersML_02 = "hhv4_clustersML_02.tsv",
        hhv4_nwk_08 = "hhv4_clustersML_08.nwk",
        hhv4_nwk_05 = "hhv4_clustersML_05.nwk",
        hhv4_nwk_02 = "hhv4_clustersML_02.nwk"
    message: "Align hhv4 genomes and cluster them"
    params:
        hhv4_fas = "hhv4_ncbi-genomes-2021-03-02/all_hhv4_genomes.fasta"
    threads:16

    run:
        R("""
          library(data.table)
          library(DECIPHER)
          library(tidyverse)
          
          hhv4_fas <- "{params.hhv4_fas}"
          
          # specify the path to the FASTA file (in quotes)
          hhv4_fas <- "{params.hhv4_fas}"

          # load the sequences from the file
          hhv4_seqs <- readDNAStringSet(hhv4_fas)
          
          # perform the alignment
          hhv4_aligned <- AlignSeqs(hhv4_seqs, processors = 16)

          # save the alignment
          writeXStringSet(hhv4_aligned, "{output.hhv4_aln}")

          # make a distance matrix from the alignment
          hhv4_dm <- DistanceMatrix(hhv4_aligned, correction = "Jukes-Cantor")

          # IdClusters returning both a data frame with the cluster and a dendrogram
          hhv4_clustersML_08 <- IdClusters(hhv4_dm, method = "ML", myXStringSet = hhv4_aligned, type="both", cutoff = 0.8, showPlot = T, processors = 4)
          hhv4_clustersML_05 <- IdClusters(hhv4_dm, method = "ML", myXStringSet = hhv4_aligned, type="both", cutoff = 0.5, showPlot = T, processors = 4)
          hhv4_clustersML_02 <- IdClusters(hhv4_dm, method = "ML", myXStringSet = hhv4_aligned, type="both", cutoff = 0.2, showPlot = T, processors = 4)
          
          # save the clusters as a table
          hhv4_clustersML_08_df <- hhv4_clustersML_08[1] %>%
            as.data.frame(.) %>%
            rownames_to_column("Genome")
          fwrite(hhv4_clustersML_08_df, "{output.hhv4_clustersML_08}", quote = F, sep = "\t")
          
          hhv4_clustersML_05_df <- hhv4_clustersML_05[1] %>%
            as.data.frame(.) %>%
            rownames_to_column("Genome")
          fwrite(hhv4_clustersML_05_df, "{output.hhv4_clustersML_05}", quote = F, sep = "\t")
          
          hhv4_clustersML_02_df <- hhv4_clustersML_02[1] %>%
            as.data.frame(.) %>%
            rownames_to_column("Genome")
          fwrite(hhv4_clustersML_02_df, "{output.hhv4_clustersML_02}", quote = F, sep = "\t")
          
          # save the dendrograms for plotting
          hhv4_clustersML_08_dend <- hhv4_clustersML_08[[2]]
          WriteDendrogram(hhv4_clustersML_08_dend, file = "{output.hhv4_nwk_08}")
          
          hhv4_clustersML_05_dend <- hhv4_clustersML_05[[2]]
          WriteDendrogram(hhv4_clustersML_05_dend, file = "{output.hhv4_nwk_05}")
          
          hhv4_clustersML_02_dend <- hhv4_clustersML_02[[2]]
          WriteDendrogram(hhv4_clustersML_02_dend, file = "{output.hhv4_nwk_02}")

        """)
        
