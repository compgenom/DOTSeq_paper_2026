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
# - STAR
# - cutadapt
#
# Usage:
# bash preprocessing_bulk.sh \
#   --rna data/bulk/rna/chx \
#   --ribo data/bulk/ribo/chx \
#   --out data/bulk \
#   --ref ref/hg38_star_index \
#   --threads 32
# =============================================================================

# ---------------------------
# Default parameters
# ---------------------------
RNA_DIR="data/bulk/rna/chx"
RIBO_DIR="data/bulk/ribo/chx"
OUT_DIR="data/bulk"
REF_DIR="ref/hg38_star_index"
THREADS=32

# ---------------------------
# Parse arguments
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rna) RNA_DIR="$2"; shift 2 ;;
    --ribo) RIBO_DIR="$2"; shift 2 ;;
    --out) OUT_DIR="$2"; shift 2 ;;
    --ref) REF_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ---------------------------
# RNA processing
# ---------------------------
echo "Processing RNA-seq..."

for i in "${RNA_DIR}"/*.fastq.gz; do
  [ -e "$i" ] || continue

  base=$(basename "$i" .fastq.gz)

  cutadapt -j "$THREADS" -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o "${RNA_DIR}/${base}.trimmed.fasta.gz" "$i"

  STAR --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$REF_DIR" \
    --readFilesIn "${RNA_DIR}/${base}.trimmed.fasta.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${RNA_DIR}/${base}_" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterType BySJout \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All
done

echo "RNA-seq processing completed."

# ---------------------------
# Ribo processing
# ---------------------------
echo "Processing Ribo-seq..."

for i in "${RIBO_DIR}"/*.fastq.gz; do
  [ -e "$i" ] || continue

  base=$(basename "$i" .fastq.gz)

  cutadapt -j "$THREADS" -m 15 -u 8 -e 0.1 --match-read-wildcards \
    -a TCGTATGCCGTCTTCTGCTTG -O 1 \
    -o "${RIBO_DIR}/${base}.trimmed.fasta.gz" "$i"

  STAR --runMode alignReads \
    --runThreadN "$THREADS" \
    --genomeDir "$REF_DIR" \
    --readFilesIn "${RIBO_DIR}/${base}.trimmed.fasta.gz" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${RIBO_DIR}/${base}_" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outFilterType BySJout \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 1 \
    --outFilterMatchNmin 16 \
    --alignEndsType EndToEnd \
    --outSAMattributes All
done

echo "Ribo-seq processing completed."