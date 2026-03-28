#!/usr/bin/env bash
set -euo pipefail

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

mkdir -p data/bulk/alignment

# RNA processing
for i in data/bulk/rna/chx/*.fastq.gz; do
  [ -e "$i" ] || continue

  base=$(basename "$i" .fastq.gz)

  cutadapt -j 16 -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o data/bulk/rna/${base}.trimmed.fasta.gz "$i"

  STAR --runMode alignReads \
    --runThreadN 32 \
    --genomeDir ref/hg38_star_index \
    --readFilesIn data/bulk/rna/${base}.trimmed.fasta.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix data/bulk/alignment/${base}_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterType BySJout \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All
done

# Ribo processing
for i in data/bulk/ribo/chx/*.fastq.gz; do
  [ -e "$i" ] || continue # if file doesn't exist it will skip.

  base=$(basename "$i" .fastq.gz)

  cutadapt -j 16 -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o data/bulk/ribo/${base}.trimmed.fasta.gz "$i"

  STAR --runMode alignReads \
    --runThreadN 32 \
    --genomeDir ref/hg38_star_index \
    --readFilesIn data/bulk/ribo/${base}.trimmed.fasta.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix data/bulk/alignment/${base}_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterType BySJout \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All
done


