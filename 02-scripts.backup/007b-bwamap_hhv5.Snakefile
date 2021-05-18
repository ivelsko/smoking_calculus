################################################################################
# Align sequencing data against multiple different reference genes
#
# Irina Velsko 01/04/2021
################################################################################

from glob import glob
import os
import re

workdir: "/projects1/microbiome_calculus/smoking_calculus/04-analysis/hhv_mapping/hhv5"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/projects1/microbiome_calculus/smoking_calculus/04-analysis/gc_rl_sub_mapped/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}_hhv5.mapped_rmdup.cov", sample=SAMPLES.keys()),
        expand("{sample}_hhv5.mapped_rmdup.bam", sample=SAMPLES.keys()),
        expand("{sample}_hhv5.mapped_rmdup.bam.bai", sample=SAMPLES.keys())

rule bwa_mem:
    output:
        temp("{sample}_hhv5.bam")
    message: "Align sample {wildcards.sample} against hhv5 using BWA aln"
    group: "bwa"
    params:
        reffa = "/projects1/microbiome_calculus/smoking_calculus/01-data/hhv5/ncbi-genomes-2021-04-01/GCF_000845245.1_ViralProj14559_genomic.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 4
    shell:
        """
        bwa mem -t {threads} \
            {params.reffa} \
            {params.fastq} > {output}
        """

rule bwa_mem_sort:
    input:
        "{sample}_hhv5.bam"
    output:
        bam = "{sample}_hhv5.mapped.bam"
    message: "Generate alignment file for sample {wildcards.sample} against hhv5"
    group: "bwa"
    params:
    shell:
        """
        samtools view -Sb -F 4 {input} | samtools sort - -o {output.bam}
        """

rule samtools_mem_rmdup:
    input:
        "{sample}_hhv5.mapped.bam"
    output:
        bam = "{sample}_hhv5.mapped_rmdup.bam"
    message: "Remove duplicate mapped reads for sample {wildcards.sample} against hhv5"
    group: "bwa"
    params:
    shell:
        """
        samtools rmdup -s {input} {output.bam}
        """

# rule samtools_index_mem:
#     input:
#         "{sample}_hhv5.mapped.bam"
#     output:
#         "{sample}_hhv5.mapped.bam.bai"
#     message: "Index bam file for sample {wildcards.sample} mapped against hhv5"
#     group: "bwa"
#     params:
#     shell:
#         """
#         samtools index {input}
#         """

rule samtools_mem_rmdup_index:
    input:
        "{sample}_hhv5.mapped_rmdup.bam"
    output:
        "{sample}_hhv5.mapped_rmdup.bam.bai"
    message: "Index bam file for sample {wildcards.sample} mapped against hhv5"
    group: "bwa"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_rmdup_coverage:
    input:
        "{sample}_hhv5.mapped_rmdup.bam"
    output:
        "{sample}_hhv5.mapped_rmdup.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against hhv5"
    group: "bwa"
    params:
    shell:
        """
        bedtools genomecov -ibam {input} > {output}
        """

# rule bedtools_coverage:
#     input:
#         "{sample}_hhv5.mapped_rmdup.bam"
#     output:
#         "{sample}_hhv5.mapped_rmdup.cov"
#     message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against hhv5"
#     group: "bwa"
#     params:
#         refbed = "/projects1/microbiome_calculus/smoking_calculus/01-data/hhv5/ncbi-genomes-2021-03-02/hhv5.bed",
#     shell:
#         """
#         bedtools coverage -a {params.refbed} -b {input} -hist | grep -v all > {output}
#         """
