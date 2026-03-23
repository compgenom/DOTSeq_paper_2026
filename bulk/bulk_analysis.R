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
})

# Step 1: Generate the ORF-level annotation using DOTSeq's getORF() function

annotation <- "ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz"
sequences <- "ref/MANE.GRCh38.v1.4.ensembl_rna.fna.gz"

gr <- getORFs(
    sequences,
    annotation,
    organism = "Homo sapiens",
    circ_seqs = NULL,
    start_codons = "ATG",
    stop_codons = "TAA|TAG|TGA",
    min_len = 0,
    longest_orf = TRUE,
    verbose = TRUE
)

saveRDS(gr, "ref/gr_orfs.rds")

# Step 2: Clean up BAM files; Prepare the counts and condition tables


bam_list <- list.files(
    path = "data/bulk",
    pattern = "Aligned.sortedByCoord.out.bam$",
    recursive = TRUE,
    full.names = TRUE
)

gr <- readRDS("ref/gr_orfs.rds")
bam_output_dir <- "data/bulk/exonic_bam"
getExonicReads(gr = gr, bam_files = bam_list, bam_output_dir = bam_output_dir, coding_genes_only = TRUE)

meta <- read.table("ref/metadata.txt.gz")
names(meta) <- c("run", "strategy", "replicate", "treatment", "condition")
cond <- meta[meta$treatment == "chx", ]
cond$treatment <- NULL

bam_files <- list.files("data/bulk/exonic_bam", pattern = ".bam$", full.names = TRUE) # get it from the previous runs by Lim because I can't install STAR and cutadapt at the moment
bam_files <- bam_files[basename(bam_files) %in% paste0(cond$run, "_Aligned.sortedByCoord.out.exonic.sorted.bam")]

cnt <- countReads(gr = gr, bam_files = bam_files)
names(cnt) <- gsub(".*(SRR[0-9]+).*", "\\1", names(cnt))

# Step 3: Create a DOTSeqDataSets object and run the DOTSeq workflow 
## Since the ORF-level annotation was prepared using the getORF() function from DOTSeq, we used DOTSeqDataSetsFromSummarizeOverlaps() instead of the DOTSeqDataSetsFromFeatureCounts() function
d <- DOTSeqDataSetsFromSummarizeOverlaps(
    count_table = cnt, 
    condition_table = cond, 
    annotation = gr
)
d <- DOTSeq(datasets = d)
saveRDS(d, "results/bulk/d.rds")

# Step 4: Extract and inspect results from DOTSeq using the getContrasts() function
results <- getContrasts(d, type = "interaction")
ou <- results$DOU[results$DOU$contrast == "Mitotic_Cycling - Interphase", ]
te <- results$DTE[results$DTE$contrast == "Mitotic_Cycling - Interphase", ]
results <- merge(ou, te, by = c("orf_id", "contrast"), all = TRUE)

# Step 5: Visualisation using the plotDOT() function

# Figure 2(A)
pdf("results/bulk/figures/heatmap.pdf", width = 4, height = 9)
plotDOT(
    plot_type = "heatmap", 
    results = results, 
    data = getDOU(d), 
    id_mapping = mapping, 
    plot_params = list(rank_by = "significance", top_hits = 50),
    force_new_device = FALSE
)
dev.off()

# Figure 2(B)
pdf("results/bulk/figures/volcano_orf_type.pdf", width = 5, height = 5)
plotDOT(
    plot_type = "volcano", 
    results = results,
    data = getDOU(d),
    id_mapping = mapping,
    plot_params = list(color_by = "orf_type", top_hits = 3, legend_position = "top"),
    force_new_device = FALSE
)
dev.off()

# Figure 2(C)
pdf("results/bulk/figures/usage.pdf", width = 3.5, height = 8)
orderby <- c("Mitotic_Cycling", "Mitotic_Arrest", "Interphase")
id <- "CSDE1"
plotDOT(
    plot_type = "usage",
    data = getDOU(d), 
    gene_id = id, 
    id_mapping = mapping, 
    plot_params = list(order_by = orderby),
    force_new_device = FALSE
)
dev.off()

# Figure 2(D)
pdf("results/bulk/figures/venn.pdf", width = 3, height = 2.5)
plotDOT(plot_type = "venn", results = results, force_new_device = FALSE)
dev.off()

# Figure 2(E)
pdf("results/bulk/figures/composite_orf_type.pdf", width = 5, height = 5)
plotDOT(
    plot_type = "composite", 
    results = results, 
    data = getDOU(d), 
    plot_params = list(color_by = "orf_type", legend_position = "bottomright"),
    force_new_device = FALSE
)
dev.off()