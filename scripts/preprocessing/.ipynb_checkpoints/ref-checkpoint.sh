#!/usr/bin/env bash

# =============================================================================
# Script: Download reference files from databases
# =============================================================================
# Sources:
# - Gene annotation and transcript FASTA: [GENCODE v47](https://www.gencodegenes.org/human)
# - Reference genome: hg38 ([UCSC Genome Browser](https://genome.ucsc.edu/))
#
# Tools:
# - SRA Toolkit
# - STAR
#
# Notes:
# - All tools are assumed to be installed and available in your PATH
# - Paths to input files may need to be adjusted for your system
# =============================================================================

# Download reference files from databases

wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \ 
  -O ../ref/hg38/hg38.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz \
  -O ../ref/gencode.v47.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz \
  -O ../ref/gencode.v47.pc_transcripts.fa.gz

# Generate STAR index
STAR --runMode genomeGenerate \
  --runThreadN 32 \
  --genomeFastaFiles ../ref/hg38/hg38.fa.gz \
  --sjdbGTFfile ../ref/gencode.v47.annotation.gtf.gz \
  --genomeDir ../ref/hg38_star_index
