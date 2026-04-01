#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script: Retrieve raw sequencing reads for DOTSeq analysis
# =============================================================================
# Sources:
# - VanInsberghe 2021 single-cell Ribo-seq dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481)
#
# Tools:
# - SRA Toolkit
#
# Notes:
# - Please make sure to clone the forked repository into sc/ prior to downloading the raw sequencing reads
# - Adjust output directories as needed the default is sc/scRiboSeq_manuscript/data_processing
#
# Usage:
# bash sc/sra_download_sc.sh \
#   --sc_dir sc/scRiboSeq_manuscript/data_processing \
#   --threads 16 \
#   --jobs 2
# =============================================================================

# ---------------------------
# Default parameters
# ---------------------------
SC_DIR="sc/scRiboSeq_manuscript/data_processing"
THREADS=16
JOBS=2

# ---------------------------
# Parse arguments
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sc_dir) SC_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ---------------------------
# Create directories
# ---------------------------
FASTQ_DIR="${SC_DIR}/fastq"
mkdir -p "$FASTQ_DIR"

# ---------------------------
# SRR list
# ---------------------------
SC_SRRS=(
SRR13125087 SRR13125088 SRR13125089 SRR13125090 SRR13125091 SRR13125092
SRR13125093 SRR13125094 SRR13125095 SRR13125096 SRR13125097 SRR13125098
SRR13125099 SRR13125100 SRR13125101 SRR13125102 SRR13125103 SRR13125104 
SRR13125105 SRR14530594 SRR14530595 SRR14530596 SRR14530597 SRR14530598
SRR14530599 SRR14530600 SRR14530601
)

echo "Downloading single-cell RNA-seq data..."

download_srr() {
    local srr="$1"

    echo "[START] $srr"

    fasterq-dump "$srr" \
        --split-files \
        --threads "$THREADS" \
        -O "$FASTQ_DIR"

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] Download failed for $srr"
        return 1
    fi

    pigz -p "$THREADS" "$FASTQ_DIR/${srr}"_*.fastq

    echo "[DONE] $srr"
}

export FASTQ_DIR THREADS
export -f download_srr

# ---------------------------
# Run in parallel
# ---------------------------
parallel -j "$JOBS" download_srr ::: "${SC_SRRS[@]}"

echo "All single-cell downloads completed."