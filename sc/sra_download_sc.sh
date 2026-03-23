#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Script: Retrieve Ribo-seq reads for DOTSeq analysis
# =============================================================================
# Sources:
# - VanInsberghe 2021 single-cell Ribo-seq dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481)
#
# Tools:
# - SRA Toolkit
#
# Notes:
# - Adjust output directories as needed
# =============================================================================

# List of SRRs to download
SRRS=(
SRR13125087 SRR13125088 SRR13125089 SRR13125090 SRR13125091 SRR13125092 SRR13125093 SRR13125094 SRR13125095
SRR13125096 SRR13125097 SRR13125098 SRR13125099 SRR13125100 SRR13125101 SRR13125102 SRR13125103 SRR13125104 
SRR13125105 SRR14530594 SRR14530595 SRR14530596 SRR14530597 SRR14530598 SRR14530599 SRR14530600 SRR14530601)

# Output directory
OUTDIR=scripts/preprocessing/sc/fastq
mkdir -p "$OUTDIR"

# Export the output dir for GNU parallel
export OUTDIR

# Function to download and compress a single SRR
download_srr() {
    local srr="$1"
    echo "Downloading $srr ..."
    fasterq-dump "$srr" --split-files --threads 16 -O "$OUTDIR"
    echo "Compressing $srr ..."
    pigz -p 16 "$OUTDIR/${srr}"_*.fastq
    echo "$srr completed."
}

export -f download_srr

# Run two SRRs in parallel
parallel -j 2 download_srr ::: "${SRRS[@]}"

echo "All downloads completed."