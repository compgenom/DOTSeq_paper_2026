#!/usr/bin/env Rscript

# =============================================================================
# Script: Preprocessing Ribo-seq and RNA-seq reads of bulk dataset
# =============================================================================
# Sources:
# - Ly 2024 HeLa cell cycle dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA957808)
# - Gene annotation and transcript FASTA: [GENCODE v47](https://www.gencodegenes.org/human)
# - Reference genome: hg38 ([UCSC Genome Browser](https://genome.ucsc.edu/))
#
# Notes:
# - Paths to input files may need to be adjusted for your system
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DOTSeq)
  library(SummarizedExperiment)
})

# Define options
option_list <- list(
  make_option(c("-c", "--count_table"), type = "character", help = "count table rds file", metavar = "string"),
  make_option(c("-a", "--annotation"), type = "character", help = "GRanges annotation rds file", metavar = "string"),
  make_option(c("-o", "--out_dir"), type = "character", default = ".", help = "Output directory [default %default]"),
  make_option(c("-g", "--gcoeff"), type = "double", help = "gcoeff value", metavar = "numeric"),
  make_option(c("-r", "--regulation_type"), type = "character", help = "Regulation type", metavar = "string"),
  make_option(c("--size_factor"), type = "double", default = NA, help = "Size factor [optional]"),
  make_option(c("--min_size"), type = "double", default = NA, help = "Minimum size [optional]"),
  make_option(c("-s", "--num_samples"), type = "integer", default = 1, help = "Number of samples [default %default]"),
  make_option(c("-b", "--num_batches"), type = "integer", default = 2, help = "Number of batches [default %default]"),
  make_option(c("--diagplot_ribo"), action = "store_true", default = FALSE, help = "Enable ribo diagnostic plots"),
  make_option(c("--diagplot_rna"), action = "store_true", default = FALSE, help = "Enable RNA diagnostic plots")
)

# Parse arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
if (is.null(opt$count_table) || is.null(opt$annotation) || is.null(opt$gcoeff) || is.null(opt$regulation_type)) {
  stop(
    "--count_table, --annotation, --gcoeff, and --regulation_type are required.\n\n", 
    "Usage example:\n", 
    "Rscript create_simdata.R \\ \n", 
    "-c ../data/ly_2024/getExonicReads_countReads.rds -a ../ref/gr_orfs.rds \\ \n",
    "-o ~/benchmarking/data/ \\ \n", #need to check with Lim which directory is this
    "-g 1.5 -r uORF_up_mORF_down \\ \n", 
    "--size_factor 1 --min_size 2"
  )
}

gr <- readRDS(opt$annotation)
cnt <- readRDS(opt$count_table)
cnt <- as.data.frame(cnt)
names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))

meta <- read.table("ref/metadata.txt.gz")
names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
cond <- meta[meta$treatment == "chx", ]
cond$treatment <- NULL

# Construct DOTSeqDataSets
d <- DOTSeqDataSetsFromSummarizeOverlaps(
  count_table = cnt,
  condition_table = cond,
  annotation = gr
)
raw_counts <- assay(getDOU(d))
raw_counts <- raw_counts[, grep("Cycling|Interphase", colnames(raw_counts))]
ribo <- raw_counts[, grep("ribo", colnames(raw_counts))]
rna <- raw_counts[, grep("rna", colnames(raw_counts))]
rowranges <- rowRanges(getDOU(d))

# Run simulation
set.seed(42)

# Conditional simulation based on optional args
if (!is.na(opt$size_factor) && !is.na(opt$min_size)) {
  message("using regulation_type: ", opt$regulation_type, "; log2FC: ", opt$gcoeff, "; size_factor: ", opt$size_factor, "; min_size: ", opt$min_size)
  simData <- simDOT(
    ribo = ribo,
    rna = rna,
    annotation = rowranges,
    regulation_type = opt$regulation_type,
    gcoeff = opt$gcoeff * log(2),
    size_factor = opt$size_factor,
    min_size = opt$min_size,
    num_samples = opt$num_samples,
    num_batches = opt$num_batches,
    diagplot_ribo = opt$diagplot_ribo,
    diagplot_rna = opt$diagplot_rna
  )

} else {
  message("using regulation_type: ", opt$regulation_type, "; log2FC: ", opt$gcoeff)
  simData <- simDOT(
    ribo = ribo,
    rna = rna,
    annotation = rowranges,
    regulation_type = opt$regulation_type,
    gcoeff = opt$gcoeff * log(2),
    num_samples = opt$num_samples,
    num_batches = opt$num_batches,
    diagplot_ribo = opt$diagplot_ribo,
    diagplot_rna = opt$diagplot_rna
  )
}

# Ensure output directory exists
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

# Construct output filename
outfile <- file.path(
  opt$out_dir,
  if (!is.na(opt$size_factor) && !is.na(opt$min_size)) {
    paste0(
      opt$regulation_type, 
      "_g", opt$gcoeff,
      "_sf", opt$size_factor, 
      "_ms", opt$min_size, 
      ".rds"
    )
  } else if (!is.na(opt$size_factor)) {
    paste0(
      opt$regulation_type, 
      "_g", opt$gcoeff, 
      "_sf", opt$size_factor, 
      ".rds"
    )
  }
)

# Save and checksum
saveRDS(simData, file = outfile)

