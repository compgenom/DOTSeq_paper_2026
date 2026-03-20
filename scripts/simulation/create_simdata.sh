#!/usr/bin/env bash

# =============================================================================
# Script: Generate simulation data
# =============================================================================
# Sources:
# - Ly 2024 cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - Gene annotations: ORF-level annotation using DOTSeq's getORF() function
# - Count table: The read counts summarised from DOTSeq. 
# Notes:
# - All FASTQ files from Ly 2024 cell cycle dataset are assumed to be downloaded via SRA Toolkit
# - Paths to input files may need to be adjusted for your system
# =============================================================================

mkdir -p data/simdata/

for r in uORF_up_mORF_down uORF_down_mORF_up; do
  for sf in 0.5 1 2 4; do
    for ms in 0.5 1 2 4; do
      for g in $(seq 0.5 0.5 2); do
        Rscript scripts/simulation/create_simdata.R \
          -c ../data/ly_2024/getExonicReads_countReads.rds \
          -a ../ref/gr_orfs.rds \
          -o data/simdata \
          -g $g -r $r \
          --size_factor $sf \
          --min_size $ms \
          -s 1 -b 3
      done
    done
  done
done