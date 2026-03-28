#!/usr/bin/env Rscript

# =============================================================================
# Script: Bulk cell-cycle dataset analysis
# =============================================================================
# Sources:
# - Ly 2024 cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - Gene annotation and transcript FASTA: [MANE](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4)
#
# Notes:
# - All FASTQ files from Ly 2024 cell cycle dataset are assumed to be downloaded via SRA Toolkit
# - Paths to input files may need to be adjusted for your system
# =============================================================================

suppressPackageStartupMessages({
    library(DOTSeq)
    library(Matrix)
    library(S4Vectors)
    library(GenomicRanges)
    library(SingleCellExperiment)
    library(ggplot2)
    library(RColorBrewer)
    library(viridis)
    library(dplyr)
    library(tidyr)
    library(FNN)
    library(glmmTMB)
    library(emmeans)
    library(DHARMa)
    library(splines)
    library(purrr)
    library(scales)
    library(argparse)
})

parser <- ArgumentParser(description = "DOTSeq bulk analysis")

parser$add_argument("-ss", "--start", type = "integer", default = 1,
                    help = "Step to start the pipeline from")

parser$add_argument("-a", "--annotation", type = "character",
                    default = "ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz")

parser$add_argument("-s", "--sequences", type = "character",
                    default = "ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz")

parser$add_argument("-bam", "--bam-dir", type = "character",
                    default = "data/bulk/alignment")

parser$add_argument("-mat", "--mat-dir", type = "character",
                    default = "data/bulk/quantification")

parser$add_argument("-o", "--out-dir", type = "character",
                    default = "results/bulk")

opt <- parser$parse_args()

# Output directories relative to DOTSeq_paper_2026
out_dir     <- opt$out_dir
fig_dir     <- file.path(out_dir, "figures")

# Create output directories if they don’t exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (opt$start == 1) {
  if (is.null(opt$annotation) || is.null(opt$sequences) || is.null(opt$bam_dir) || is.null(opt$mat_dir) || 
      is.null(opt$gr_dir) || is.null(opt$out_dir)) {
    stop("To start from alignment files, you must provide --outdir, --annotation, --sequences, --bam-dir, --mat-dir and --gr-dir")
  }
  
  # Step 1: Generate the ORF-level annotation using DOTSeq's getORFs() function
  gr <- getORFs(
    sequences = opt$sequences,
    annotation = opt$annotation,
    organism = "Homo sapiens",
    start_codons = "ATG",
    stop_codons = "TAA|TAG|TGA",
    min_len = 0,
    longest_orf = TRUE,
    verbose = TRUE
  )
    
  saveRDS(gr, opt$gr_dir)
  message("ORF annotation done. ORF GRanges saved to", opt$gr_dir)

  # Step 2: Clean up BAM files; Prepare the counts and condition tables
  bam_list <- list.files(
    path = opt$bam_dir,
    pattern = "Aligned.sortedByCoord.out.bam$",
    recursive = TRUE,
    full.names = TRUE
    )
    
  exonic_dir <- file.path(opt$bam_dir, "exonic_bam")
    
  dir.create(exonic_dir, recursive = TRUE, showWarnings = FALSE)
    
  gr <- readRDS(opt$gr_dir)
  bam_output_dir <- exonic_dir
  getExonicReads(gr = gr, bam_files = bam_list, bam_output_dir = bam_output_dir, coding_genes_only = TRUE)
    
  meta <- read.table("ref/metadata_PRJNA957808.txt.gz")
  names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
  cond <- meta[meta$treatment == "chx", ]
  cond$treatment <- NULL
    
  bam_files <- list.files(exonic_dir, pattern = ".bam$", full.names = TRUE)
  bam_files <- bam_files[basename(bam_files) %in% paste0(cond$run, "_Aligned.sortedByCoord.out.exonic.sorted.bam")]
    
  cnt <- countReads(gr = gr, bam_files = bam_files)
  names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))
    
  # Step 3: Create a DOTSeqDataSets object and run the DOTSeq workflow 
  # Since the ORF-level annotation was prepared using the getORF() function from DOTSeq, we used DOTSeqDataSetsFromSummarizeOverlaps() instead of the DOTSeqDataSetsFromFeatureCounts() function
  d <- DOTSeqDataSetsFromSummarizeOverlaps(
        count_table = cnt, 
        condition_table = cond, 
        annotation = gr
    )
  d <- DOTSeq(datasets = d)
  saveRDS(d, file.path(opt$mat_dir, "bulk.rds"))
}

if (opt$start <=2) {
  if (is.null(opt$mat_dir) || is.null(opt$out_dir) || is.null(opt$gr_dir)) {
    stop("To start with precomputed per-ORF matrices, you must provide --gr-dir, --mat-dir and --out-dir")
  }
    
    # Step 4: Extract and inspect results from DOTSeq using the getContrasts() function
    d <- readRDS(file.path(opt$mat_dir, "bulk.rds"))
    
    results <- getContrasts(d, type = "interaction")
    ou <- results$DOU[results$DOU$contrast == "Mitotic_Cycling - Interphase", ]
    te <- results$DTE[results$DTE$contrast == "Mitotic_Cycling - Interphase", ]
    results <- merge(ou, te, by = c("orf_id", "contrast"), all = TRUE)
    
    # Step 5: Visualisation using the plotDOT() function
    
    pdf(file.path(fig_dir, "venn.pdf"),  width = 3, height = 2.5)
    plotDOT(plot_type = "venn", results = results, force_new_device = FALSE)
    dev.off()
    
    pdf(file.path(fig_dir, "composite.pdf"), width = 4, height = 4)
    plotDOT(
        plot_type = "composite", 
        results = results, 
        plot_params = list(color_by = "significance", legend_position = "bottomright"),
        force_new_device = FALSE
    )
    dev.off()
    
    pdf(file.path(fig_dir, "composite_orf_type.pdf"), width = 5, height = 5)
    plotDOT(
        plot_type = "composite", 
        results = results, 
        data = getDOU(d), 
        plot_params = list(color_by = "orf_type", legend_position = "bottomright"),
        force_new_device = FALSE
    )
    dev.off()
    
    pdf(file.path(fig_dir, "volcano.pdf"), width = 5, height = 5)
    plotDOT(
        plot_type = "volcano", 
        results = results,
        id_mapping = TRUE,
        plot_params = list(color_by = "significance", top_hits = 3, legend_position = "topright"),
        force_new_device = FALSE
    )
    dev.off()
    
    pdf(file.path(fig_dir, "volcano_orf_type.pdf"), width = 5, height = 5)
    plotDOT(
        plot_type = "volcano", 
        results = results,
        data = getDOU(d),
        id_mapping = TRUE,
        plot_params = list(color_by = "orf_type", top_hits = 3, legend_position = "top"),
        force_new_device = FALSE
    )
    dev.off()
    
    pdf(file.path(fig_dir, "heatmap.pdf"), width = 5, height = 5)
    plotDOT(
        plot_type = "heatmap", 
        results = results, 
        data = getDOU(d), 
        id_mapping = TRUE, 
        plot_params = list(rank_by = "significance", top_hits = 50),
        force_new_device = FALSE
    )
    dev.off()
    
    pdf(file.path(fig_dir, "usage.pdf"), width = 5, height = 5)
    orderby <- c("Mitotic_Cycling", "Mitotic_Arrest", "Interphase")
    id <- "CSDE1"
    plotDOT(
        plot_type = "usage",
        data = getDOU(d), 
        gene_id = id, 
        id_mapping = TRUE, 
        plot_params = list(order_by = orderby),
        force_new_device = FALSE
    )
    dev.off()
}