#!/usr/bin/env bash

# =============================================================================
# Script: Download reference files from databases and generate STAR index
# =============================================================================
#
# Tools:
# - SRA Toolkit
# - STAR
#
# Notes:
# - Paths to input files may need to be adjusted for your system
# =============================================================================

# Download reference files from databases

wget https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \ 
  -O ref/hg38/hg38.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz \
  -O ref/gencode.v47.annotation.gtf.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.pc_transcripts.fa.gz \
  -O ref/gencode.v47.pc_transcripts.fa.gz

wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  -O ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz
  
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_rna.fna.gz \
  -O ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz
  
# Generate STAR index
STAR --runMode genomeGenerate \
  --runThreadN 32 \
  --genomeFastaFiles ref/hg38/hg38.fa.gz \
  --sjdbGTFfile ref/gencode.v47.annotation.gtf.gz \
  --genomeDir ref/hg38_star_index
