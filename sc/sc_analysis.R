#!/usr/bin/env Rscript

# =============================================================================
# Script: Single-cell analysis
# =============================================================================
# Sources:
# - VanInsberghe 2021 single-cell Ribo-seq dataset: [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA680481)
# - Gene annotation and transcript FASTA: [MANE](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4)
#
# Notes:
# - All FASTQ files from Vanlnsberghe 2021 single-cell Ribo-seq dataset are assumed to be downloaded via SRA Toolkit
# - Paths to input files may need to be adjusted for your system
# =============================================================================

suppressPackageStartupMessages({
    library(DOTSeq)
    library(Matrix)
    library(S4Vectors)
    library(GenomicRanges)
    library(SingleCellExperiment)
    library(zellkonverter)
    library(Seurat)
    library(SeuratObject) 
    library(aricode)
    library(ggplot2)
    library(RColorBrewer)
    library(viridis)
    library(ggExtra)
    library(dplyr)
    library(tidyr)
    library(FNN)
    library(glmmTMB)
    library(emmeans)
    library(DHARMa)
    library(splines)
    library(purrr)
    library(igraph)
    library(scales)
})

# Step 1: Generate the ORF-level annotation using DOTSeq's getORFs() function
# Skip this step if you have already generated the ORF-level annotation in bulk_analysis.R

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

# Step 2: Count single-cell reads/UMIs from BAM over genomic features using DOTSeq's countReadsSingleCell() function

bam_list <- list.files(
    path = "sc/data_processing/work_addCB",
    pattern = "_Aligned.sortedByCoord.out_CB.bam$",
    recursive = TRUE,
    full.names = TRUE
)

bam_list <- bam_list[7:37]
outdir <- "data/sc" 
for (i in seq_along(bam_list)) {
    mat <- countReadsSingleCell(gr, bam_list[i])
    fname <- tools::file_path_sans_ext(basename(bam_list[i]))
    saveRDS(mat, file = file.path(outdir, paste0(fname, ".mat.rds")))
}

# Step 3: Match per-ORF count matrices with sample metadata

`%||%` <- function(a, b) if (!is.null(a)) a else b
printf <- function(...) cat(sprintf(...), "\n")

# Where per-ORF matrices live
matrix_dir <- "sc/data_processing/RPF_update/scRibo/quantification"

# Output directories relative to DOTSeq_paper_2026
out_dir     <- "results/sc"
fig_dir     <- file.path(out_dir, "figures")

# Create output directories if they don’t exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## SRA meta (expects: Run, cell_type, treatment, ...)
sra <- read.csv("ref/SraRunTable_PRJNA680481.csv")
sra <- sra[sra$cell_type == "hTERT RPE-1", , drop = FALSE]

## Per-ORF matrices named "<Run>_Aligned.sortedByCoord.out_CB.mat.rds"
rds_files <- list.files(matrix_dir, pattern = "_Aligned\\.sortedByCoord\\.out_CB\\.mat\\.rds$", full.names = TRUE)
if (!length(rds_files)) stop("No per-ORF .rds matrices found.")
mats_acc  <- vapply(strsplit(basename(rds_files), "_"), `[`, character(1), 1)
sra_sub <- sra[sra$Run %in% mats_acc, , drop = FALSE]
sra_sub$matrix <- file.path(matrix_dir, paste0(sra_sub$Run, "_Aligned.sortedByCoord.out_CB.mat.rds"))
stopifnot(nrow(sra_sub) > 0, all(file.exists(sra_sub$matrix)))
printf("Included runs: %d", nrow(sra_sub))

# Step 4: Load and merge per-ORF matrices with unified feature space

read_one <- function(path, run) {
    mat <- readRDS(path)
    stopifnot(inherits(mat, "dgCMatrix"))
    colnames(mat) <- paste0(run, "_", colnames(mat))
    mat
}
mats <- setNames(vector("list", nrow(sra_sub)), sra_sub$Run)
for (i in seq_len(nrow(sra_sub))) mats[[i]] <- read_one(sra_sub$matrix[i], sra_sub$Run[i])
printf("Read %d per-ORF matrices.", length(mats))

# Count number of cells
num_cells <- function(m, nonempty_only = TRUE) {
    stopifnot(inherits(m, "dgCMatrix"))
    if (nonempty_only) sum(Matrix::colSums(m) > 0) else ncol(m)
}

for (i in seq_len(nrow(sra_sub))) sra_sub$num_cells[i] <- num_cells(mats[[i]])
aggregate(num_cells ~ treatment, data = sra_sub, FUN = function(x) sum(x, na.rm = TRUE))


all_orfs <- Reduce(union, lapply(mats, rownames)) |> sort()

pad_rows_to <- function(M, genes) {
    miss <- setdiff(genes, rownames(M))
    if (length(miss)) {
        Z <- Matrix(0, nrow = length(miss), ncol = ncol(M), sparse = TRUE)
        rownames(Z) <- miss
        M <- rbind(M, Z)
    }
    M[genes, , drop = FALSE]
}
mats <- lapply(mats, pad_rows_to, genes = all_orfs)
combined <- do.call(cbind, mats)  # ORFs x cells
stopifnot(inherits(combined, "dgCMatrix"))
printf("Combined ORF matrix: %d ORFs x %d cells", nrow(combined), ncol(combined))

## Drop all-zero ORFs
keep_rows <- Matrix::rowSums(combined) > 0
M <- combined[keep_rows, , drop = FALSE]
printf("After removing all-zero ORFs: %d ORFs x %d cells", nrow(M), ncol(M))


# Step 5: Map ORFs to genes; split uORF/mORF; aggregate to gene level

gr <- readRDS("ref/gr_orfs.rds"); stopifnot(is(gr, "GRanges"))
idx <- match(rownames(M), names(gr))
if (anyNA(idx)) {
    n_miss <- sum(is.na(idx))
    counts_lost <- sum(Matrix::rowSums(M[is.na(idx), , drop = FALSE]))
    stop(sprintf("Annotation missing for %d/%d ORFs; lost counts=%g", n_miss, nrow(M), counts_lost))
}
gr_aligned <- gr[idx]
gene <- as.character(mcols(gr_aligned)$gene_id)
orf  <- as.character(mcols(gr_aligned)$orf_type)
stopifnot(length(gene) == nrow(M), length(orf) == nrow(M))

rows_morf <- which(orf == "mORF")
rows_uorf <- which(orf == "uORF")
if (!length(rows_morf)) stop("No 'mORF' rows found.")
if (!length(rows_uorf)) stop("No 'uORF' rows found.")

## Totals sanity
total_all  <- sum(M)
total_mORF <- sum(M[rows_morf, , drop = FALSE])
total_uORF <- sum(M[rows_uorf, , drop = FALSE])
printf("Counts total=%g, mORF=%g, uORF=%g, used frac=%.3f",
       total_all, total_mORF, total_uORF, (total_mORF + total_uORF) / total_all)

## Aggregate to genes
genes_factor <- factor(gene)
i_m <- as.integer(genes_factor[rows_morf]); j_m <- rows_morf
i_u <- as.integer(genes_factor[rows_uorf]); j_u <- rows_uorf

G_m <- sparseMatrix(i = i_m, j = j_m, x = 1, dims = c(nlevels(genes_factor), nrow(M)))
G_u <- sparseMatrix(i = i_u, j = j_u, x = 1, dims = c(nlevels(genes_factor), nrow(M)))

M_morf <- G_m %*% M; rownames(M_morf) <- levels(genes_factor)
M_uorf <- G_u %*% M; rownames(M_uorf) <- levels(genes_factor)
stopifnot(inherits(M_morf, "dgCMatrix"), inherits(M_uorf, "dgCMatrix"), identical(dim(M_morf), dim(M_uorf)))
printf("Gene-level matrices: %d genes x %d cells", nrow(M_morf), ncol(M_morf))

# Step 6: Create a SCE and Seurat v5 layered object

cell_ids <- colnames(M_morf)
gene_ids <- rownames(M_morf)
sce <- SingleCellExperiment(
    assays  = list(uORF = M_uorf, mORF = M_morf),
    colData = S4Vectors::DataFrame(cell_id = cell_ids),
    rowData = S4Vectors::DataFrame(gene_id = gene_ids)
)
srr_from_cells <- sub("_.*", "", cell_ids)
sce$treatment  <- sra$treatment[match(srr_from_cells, sra$Run)]

sobj <- CreateSeuratObject(counts = M_morf, assay = "mORF", project = "uORF_mORF")
sobj[["uORF"]] <- CreateAssayObject(counts = M_uorf)

sobj$SRR <- factor(srr_from_cells)
sobj$treatment <- sra$treatment[match(srr_from_cells, sra$Run)]
sobj$treatment[is.na(sobj$treatment)] <- "Unknown"

## Helper: totals from v5 layers
get_counts_totals <- function(obj, assay) {
    lyr <- tryCatch(Layers(obj[[assay]]), error = function(e) character(0))
    if (!length(lyr)) stop(sprintf("Assay '%s' has no layers.", assay))
    target <- if ("counts" %in% lyr) "counts" else lyr[1]
    Matrix::colSums(LayerData(obj, assay = assay, layer = target))
}

## Totals and fraction
mORF_tot_counts <- get_counts_totals(sobj, "mORF")
uORF_tot_counts <- get_counts_totals(sobj, "uORF")
sobj$uORF_mORF_ratio <- uORF_tot_counts / (mORF_tot_counts + 1)

# Step 7: Create uORF- and mORF-only UMAP 

## mORF-only UMAP (default assay mORF)
DefaultAssay(sobj) <- "mORF"
sobj <- NormalizeData(sobj, layer = "counts")      # writes 'data'
sobj <- FindVariableFeatures(sobj)
sobj <- ScaleData(sobj, layer = "data")
sobj <- RunPCA(sobj)
sobj <- FindNeighbors(sobj, dims = 1:20)
sobj <- RunUMAP(sobj, dims = 1:20)
sobj <- FindClusters(sobj, resolution = 0.5)

## palette
stage_colors <- c(
    "Interphase treatment"              = "#619CFF",
    "Interphase and Mitotic treatments" = "#00BA38",
    "Interphase and G0 treatments"      = "#F8766D",
    "Unknown"                           = "grey70"
)

## UMAPs (mORF) Figure S5(A, C)
p_umap_ratio <- FeaturePlot(
    sobj, features = "uORF_mORF_ratio", reduction = "umap"
) + scale_color_viridis(option = "D") +
    labs(title = "UMAP coloured by uORF proportion (mORF data)") +
    theme_classic()

p_umap_trt <- DimPlot(sobj, reduction = "umap", group.by = "treatment") +
    scale_color_manual(values = stage_colors) +
    labs(title = "UMAP by coloured by treatment (mORF data)") +
    theme_classic() # + theme(legend.position = "bottom")

print(p_umap_ratio); print(p_umap_trt)


## uORF-only UMAP (sparse-aware; keep cells with uORF UMIs > 50) Figure S5(B, D)
DefaultAssay(sobj) <- "uORF"
u_counts <- LayerData(sobj, assay = "uORF", layer = "counts")
keep_u   <- Matrix::colSums(u_counts) > 50
s_u      <- subset(sobj, cells = colnames(sobj)[keep_u])

s_u <- NormalizeData(s_u, layer = "counts")
s_u <- FindVariableFeatures(s_u, nfeatures = 2000)
s_u <- ScaleData(s_u, layer = "data")
s_u <- RunPCA(s_u, features = VariableFeatures(s_u), npcs = 30, verbose = FALSE)
s_u <- FindNeighbors(s_u, reduction = "pca", dims = 1:20)
s_u <- RunUMAP(s_u, reduction = "pca", dims = 1:20)

p_umap_uorf_ratio <- FeaturePlot(
    s_u, features = "uORF_mORF_ratio", reduction = "umap"
) + scale_color_viridis(option = "D") + 
    labs(title = "UMAP coloured by uORF proportion (uORF data)") + theme_classic()

p_umap_uorf_trt <- DimPlot(s_u, reduction = "umap", group.by = "treatment") +
    scale_color_manual(values = stage_colors) +
    labs(title = "UMAP colored by treatments (uORF data)") +
    theme_classic() # + theme(legend.position = "bottom")

print(p_umap_uorf_ratio); print(p_umap_uorf_trt)

# Step 8: Create a WNN joint manifold Figure 3C, 3D)
##     - mORF: SCTransform + PCA  (SCT_mORF)
##     - uORF: LogNormalize + PCA (assay = "uORF")

library(future)
options(future.globals.maxSize = 1 * 1024^3)  # 1 GB limit

## Subset to cells with sufficient uORF depth
u_counts <- LayerData(sobj, assay = "uORF", layer = "counts")
u_n      <- Matrix::colSums(u_counts)

thr <- 50L
cells_ok <- colnames(sobj)[u_n > thr]
message(sprintf("WNN subset: keeping %d / %d cells (%.1f%%) with uORF UMIs > %d",
                length(cells_ok), ncol(sobj), 100*length(cells_ok)/ncol(sobj), thr))

s_wnn <- subset(sobj, cells = cells_ok)   # metadata (treatment, SRR, uORF_mORF_ratio) is retained

## SCT per assay on the subset (uORF totals > 0) 
set.seed(123)  # reproducibility
s_wnn <- SCTransform(s_wnn, assay = "mORF", new.assay.name = "SCT_mORF", verbose = FALSE)
s_wnn <- SCTransform(s_wnn, assay = "uORF", new.assay.name = "SCT_uORF", verbose = FALSE)

## PCA per assay with distinct reduction names 
s_wnn <- RunPCA(s_wnn, assay = "SCT_mORF", reduction.name = "pca_mORF", npcs = 30)
s_wnn <- RunPCA(s_wnn, assay = "SCT_uORF", reduction.name = "pca_uORF", npcs = 30)

## WNN neighbors + UMAP across both PC spaces 
s_wnn <- FindMultiModalNeighbors(
    s_wnn,
    reduction.list = list("pca_mORF", "pca_uORF"),
    dims.list      = list(1:20, 1:20)
)
s_wnn <- RunUMAP(s_wnn, nn.name = "weighted.nn", reduction.name = "wnn.umap")

## Plots (orientation only) 
p_wnn_trt <- DimPlot(s_wnn, reduction = "wnn.umap", group.by = "treatment") +
    scale_color_manual(values = stage_colors) +
    labs(title = sprintf("WNN-UMAP (mORF + uORF); uORF UMIs > %d", thr)) +
    theme_classic() # + theme(legend.position = "bottom")

p_wnn_ratio <- FeaturePlot(s_wnn, reduction = "wnn.umap", features = "uORF_mORF_ratio") +
    labs(title = "WNN-UMAP colored by uORF fraction") + theme_classic()

print(p_wnn_trt); print(p_wnn_ratio)

## quick sanity check if the fraction structured on this graph
Glist <- tryCatch(Seurat::Graphs(s_wnn), error = function(e) list())
if ("weighted.nn" %in% names(Glist)) {
    A <- as(Glist[["weighted.nn"]], "dgCMatrix")
    y <- s_wnn$uORF_mORF_ratio
    Wy <- as.vector((A %*% y) / pmax(1, Matrix::rowSums(A)))   # neighbor-averaged fraction
    morans_I <- suppressWarnings(cor(y, Wy, use = "pairwise.complete.obs"))
    message(sprintf("Graph Moran's I on WNN (subset): %.4f", morans_I))
}

# Step 9: GLMM on uORF fraction (depth-neutral), SRR random effect

DefaultAssay(sobj) <- "mORF"  # not required for GLMM inputs below

df <- data.frame(
    u = as.numeric(uORF_tot_counts),
    m = as.numeric(mORF_tot_counts),
    n = as.numeric(uORF_tot_counts + mORF_tot_counts),
    Treatment = factor(sobj$treatment),
    SRR = factor(sobj$SRR)
)
df <- subset(df, n > 0)

fit_bb <- glmmTMB(
    cbind(u, n - u) ~ Treatment + (1 | SRR),
    family = betabinomial(), data = df
)
print(summary(fit_bb))

# Step 10: EMMs (probability scale) + pairwise differences (Figure 3A & 3B)

emm_trt <- emmeans(fit_bb, ~ Treatment, type = "response")
emm_df  <- as.data.frame(emm_trt)
y_col   <- if ("prob" %in% names(emm_df)) "prob" else
    if ("emmean" %in% names(emm_df)) "emmean" else
        if ("response" %in% names(emm_df)) "response" else stop("No prob/emmean/response column.")
lcl_col <- if ("lower.CL" %in% names(emm_df)) "lower.CL" else "asymp.LCL"
ucl_col <- if ("upper.CL" %in% names(emm_df)) "upper.CL" else "asymp.UCL"

p_emm <- ggplot(emm_df, aes(x = Treatment, y = .data[[y_col]])) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = .data[[lcl_col]], ymax = .data[[ucl_col]]), width = 0.15) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "uORF fraction (u/(u+m))", x = NULL,
         title = "Model-adjusted uORF fraction by treatment (beta-binomial GLMM)") +
    theme_classic() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
print(p_emm)

## Pairwise differences on probability scale (percentage points)
pw <- pairs(emmeans(fit_bb, ~ Treatment), adjust = "tukey")
pw_or <- as.data.frame(summary(pw, infer = TRUE, type = "response"))   # has 'odds.ratio'

## Robust column names
y_col_pw <- "odds.ratio"
lcl2 <- if ("lower.CL" %in% names(pw_or)) "lower.CL" else "asymp.LCL"
ucl2 <- if ("upper.CL" %in% names(pw_or)) "upper.CL" else "asymp.UCL"

## Plot ORs (log scale) 
p_pw_or <- ggplot(pw_or, aes(x = contrast, y = .data[[y_col_pw]])) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = .data[[lcl2]], ymax = .data[[ucl2]]), width = 0.15) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_y_log10() +
    labs(x = NULL, y = "Odds ratio (uORF fraction; log scale)",
         title = "Pairwise treatment contrasts (odds ratios)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

print(p_pw_or)

# Step 15: Save figures

ggsave(file.path(fig_dir, "umap_morf_ratio.pdf"), p_umap_ratio, width = 4.75, height = 4)
ggsave(file.path(fig_dir, "umap_morf_trt.pdf"),   p_umap_trt,   width = 6.45, height = 4)
ggsave(file.path(fig_dir, "umap_uorf_ratio.pdf"), p_umap_uorf_ratio, width = 4.75, height = 4)
ggsave(file.path(fig_dir, "umap_uorf_trt.pdf"),   p_umap_uorf_trt,   width = 6.45, height = 4)
ggsave(file.path(fig_dir, "wnn_trt.pdf"),         p_wnn_trt,       width = 6.45, height = 4)
ggsave(file.path(fig_dir, "wnn_ratio.pdf"),       p_wnn_ratio,     width = 4.75, height = 4)
ggsave(file.path(fig_dir, "emm_by_trt.pdf"),      p_emm,           width = 4, height = 4)
ggsave(file.path(fig_dir, "pw_diffs.pdf"),        p_pw_or,         width = 4, height = 4.75)

