#!/usr/bin/env bash

# =============================================================================
# Script: Preprocessing Ribo-seq reads of single-cell dataset
# =============================================================================
# Sources:
# - VanInsberghe 2021 single-cell Ribo-seq dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481)
# - Gene annotation and transcript FASTA: [GENCODE v47](https://www.gencodegenes.org/human)
# - Reference genome: hg38 ([UCSC Genome Browser](https://genome.ucsc.edu/))
#
# Tools:
# - SRA Toolkit
# - STAR
# - cutadapt
# - Nextflow
#
# Notes:
# - All tools are assumed to be installed and available in your PATH
# - All FASTQ files from VanInsberghe 2021 single-cell Ribo-seq dataset are assumed to be downloaded via SRA Toolkit
# - Paths to input files may need to be adjusted for your system
# =============================================================================

# Step 1: Identify contaminating sequences
join -1 4 -2 1 -t$'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12 \
<(sort -k4,4 ../ref/gencode.v47.annotation.bed) \
<(awk '$3~/transcript/ || /tRNA/ || /rRNA/ || /snRNA/ || /snoRNA/ {print $12}' ../ref/gencode.v47.annotation.gtf | sed 's/"//g;s/;//' | sort) \
| awk '$7==$8' > ../ref/gencode.v47.annotation.contaminations.bed

# Step 2: Extract and deduplicate contaminant sequences
bedtools getfasta -bed ../ref/gencode.v47.annotation.contaminations.bed -fi ../ref/hg38/hg38.fa -name -split -s \
> ../ref/contaminants/hg38.contaminations.fa
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ../ref/contaminants/hg38.contaminations.fa \
| awk '!seen[$0]++ {print $1 "\n" $NF}' > ../ref/contaminants/hg38.contaminations.fas
mv ../ref/contaminants/hg38.contaminations.fas ../ref/contaminants/hg38.contaminations.fa

# Step 3: Set Nextflow environment and run the pipeline
export NXF_HOME=../../sc/scRiboSeq_manuscript/data_processing/.nextflow
export NXF_WORK=../../sc/scRiboSeq_manuscript/data_processing/work
export NXF_TEMP=../../sc/scRiboSeq_manuscript/data_processing/tmp

NXF_VER=22.10.1 nextflow run main.nf \
  -with-conda \
  --reads '../../sc/scRiboSeq_manuscript/data_processing/fastq/*_{1,2}.fastq.gz' \
  --genomeFasta '../ref/hg38/hg38.fa' \
  --annotations '../ref/gencode.v47.annotation.gtf' \
  --contaminationFasta '../ref/contaminants/hg38.contaminations.fa' \
  -c ../../sc/scRiboSeq_manuscript/data_processing/local_override.config -resume
