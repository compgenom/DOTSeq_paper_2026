#!/usr/bin/env Rscript

# =============================================================================
# Script: Generate the ORF-level annotation using DOTSeq's `getORFs()` function
# =============================================================================
# Sources:
# Gene annotation and transcript FASTA: [MANE](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4)
#
# Notes:
# - Paths to input files may need to be adjusted for your system
# =============================================================================

suppressPackageStartupMessages({
    library(DOTSeq)
    library(Matrix)
    library(S4Vectors)
    library(GenomicRanges)
    library(argparse)
})

parser <- ArgumentParser(description = "Generate the ORF-level annotation")

parser$add_argument("-a", "--annotation", type = "character",
                    default = "ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz")

parser$add_argument("-s", "--sequences", type = "character",
                    default = "ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz")

parser$add_argument("-o", "--out-dir", type = "character",
                    default = "ref")

opt <- parser$parse_args()

gr <- getORFs(
    sequences = opt$sequences,
    annotation = opt$annotation,
    organism = "Homo sapiens",
    circ_seqs = NULL,
    start_codons = "ATG",
    stop_codons = "TAA|TAG|TGA",
    min_len = 0,
    longest_orf = TRUE,
    verbose = TRUE
  )
    
saveRDS(gr, file.path(opt$out_dir,"gr_orfs.rds"))

message("ORF annotation done. ORF GRanges saved to ", opt$out_dir)