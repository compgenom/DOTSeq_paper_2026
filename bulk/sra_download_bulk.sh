#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script: Retrieve Ribo-seq and RNA-seq reads of bulk dataset
# =============================================================================
# Sources:
# - Ly 2024 HeLa cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
#
# Tools:
# - SRA Toolkit
#
# Notes:
# - All tools are assumed to be installed and available in your PATH
# - Paths to input files may need to be adjusted for your system
# =============================================================================

# RNA download
for srr in SRR24230462 SRR24230466 SRR24230472 SRR24230474 SRR24230477 SRR24230479; do
    fasterq-dump "$srr" --threads 16 -O ../../../dotseq/data/ly_2024/rna/chx
    pigz -p 16 ../../../dotseq/data/ly_2024/rna/chx/"$srr".fastq
done

echo "All RNA-seq downloads completed."

# Ribo download
for srr in SRR24230465 SRR24230467 SRR24230469 SRR24230471 SRR24230480 SRR24230482; do
    fasterq-dump "$srr" --threads 16 -O ../../../dotseq/data/ly_2024/ribo/chx
    pigz -p 16 ../../../dotseq/data/ly_2024/ribo/chx/"$srr".fastq
done

echo "All Ribo-seq downloads completed."