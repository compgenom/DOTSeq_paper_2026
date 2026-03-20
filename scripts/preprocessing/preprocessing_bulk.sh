#!/usr/bin/env bash

# =============================================================================
# Script: Preprocessing Ribo-seq and RNA-seq reads of bulk dataset
# =============================================================================
# Sources:
# - Ly 2024 HeLa cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - Gene annotation and transcript FASTA: [GENCODE v47](https://www.gencodegenes.org/human)
# - Reference genome: hg38 ([UCSC Genome Browser](https://genome.ucsc.edu/))
#
# Tools:
# - SRA Toolkit
# - STAR
# - cutadapt
#
# Notes:
# - All tools are assumed to be installed and available in your PATH
# - All FASTQ files from Ly 2024 HeLa cell cycle dataset are assumed to be downloaded via SRA Toolkit
# - Paths to input files may need to be adjusted for your system
# =============================================================================

# Trim and align RNA-seq reads
for i in ../data/ly_2024/rna/chx/*.fastq.gz; do
  cutadapt -j 16 -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz "$i"

  STAR --runMode alignReads \
    --runThreadN 32 \
    --outFilterType BySJout --outFilterMismatchNmax 2 \
    --genomeDir ../ref/hg38_star_index \
    --readFilesIn $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix data/bulk/$(basename "$i" .fastq.gz) \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd --outSAMattributes All
done

# Trim and align Ribo-seq reads
for i in ../data/ly_2024/ribo/chx/*.fastq.gz; do
  cutadapt -j 16 -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz "$i"

  STAR --runMode alignReads \
    --runThreadN 32 \
    --outFilterType BySJout --outFilterMismatchNmax 2 \
    --genomeDir ../ref/hg38_star_index \
    --readFilesIn $(dirname "$i")/$(basename "$i" .fastq.gz).trimmed.fasta.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix data/bulk/$(basename "$i" .fastq.gz) \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterMultimapNmax 1 --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd --outSAMattributes All
done
