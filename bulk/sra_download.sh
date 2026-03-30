#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script: Retrieve raw sequencing reads for DOTSeq analysis
# =============================================================================
# Sources:
# - Ly 2024 HeLa cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - VanInsberghe 2021 single-cell Ribo-seq dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481)
#
# Tools:
# - SRA Toolkit
#
# Notes:
# - Adjust output directories as needed the default is data/bulk/rna data/bulk/ribo and data/sc/fastq
#
# Usage:
# bash sra_downloads.sh \
#   --bulk_rna data/bulk/rna \
#   --bulk_ribo data/bulk/ribo \
#   --threads 16 \
# =============================================================================


# ---------------------------
# Default parameters
# ---------------------------
mkdir -p "data"
mkdir -p "data/bulk"
mkdir -p "data/bulk/rna"
mkdir -p "data/bulk/ribo"

BULK_RNA_DIR="data/bulk/rna"
BULK_RIBO_DIR="data/bulk/ribo"
THREADS=16
JOBS=2

# ---------------------------
# Parse arguments
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bulk_rna) BULK_RNA_DIR="$2"; shift 2 ;;
    --bulk_ribo) BULK_RIBO_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ---------------------------
# Create directories
# ---------------------------
mkdir -p "${BULK_RNA_DIR}/chx"
mkdir -p "${BULK_RIBO_DIR}/chx"

# ---------------------------
# Bulk RNA-seq (CHX-treated)
# ---------------------------
RNA_SRRS=(
SRR24230462 SRR24230466 SRR24230472
SRR24230474 SRR24230477 SRR24230479
)

echo "Downloading bulk RNA-seq..."

for srr in "${RNA_SRRS[@]}"; do
    fasterq-dump "$srr" --threads "$THREADS" -O "${BULK_RNA_DIR}/chx"
    pigz -p "$THREADS" "${BULK_RNA_DIR}/chx/${srr}.fastq"
done

echo "Bulk RNA-seq downloads completed."

# ---------------------------
# Bulk Ribo-seq (CHX-treated)
# ---------------------------
RIBO_SRRS=(
SRR24230465 SRR24230467 SRR24230469
SRR24230471 SRR24230480 SRR24230482
)

echo "Downloading bulk Ribo-seq..."

for srr in "${RIBO_SRRS[@]}"; do
    fasterq-dump "$srr" --threads "$THREADS" -O "${BULK_RIBO_DIR}/chx"
    pigz -p "$THREADS" "${BULK_RIBO_DIR}/chx/${srr}.fastq"
done

echo "Bulk Ribo-seq downloads completed."
echo "All bulk downloads completed."