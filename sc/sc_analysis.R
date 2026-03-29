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
  library(argparse)
})

# Define options
parser <- ArgumentParser(description = "DOTSeq single cell analysis")

parser$add_argument("-ss", "--start", type = "integer", default = 1,
                    help = "Step to start the pipeline from")

parser$add_argument("-gr", "--gr-dir", type = "character",
                    default = "ref/gr_orfs.rds")

parser$add_argument("-bam", "--bam-dir", type = "character",
                    default = "data/sc/alignment")

parser$add_argument("-mat", "--mat-dir", type = "character",
                    default = "data/sc/quantification")

parser$add_argument("-umi", "--umi-thr", type = "integer",
                    default = 50)

parser$add_argument("-o", "--out-dir", type = "character",
                    default = "results/sc")

opt <- parser$parse_args()

# Output directories relative to DOTSeq_paper_2026
out_dir     <- opt$out_dir
fig_dir     <- file.path(out_dir, "figures")

# Create output directories if they don’t exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (opt$start == 1) {
  if (is.null(opt$bam_dir) || is.null(opt$mat_dir) || is.null(opt$gr_dir) || is.null(opt$out_dir)) {
    stop("To start from read counting step, you must provide --out-dir, --bam-dir, --gr-dir, and --mat-dir")
  }
  
  # Step 1: Count single-cell reads/UMIs from BAM over genomic features using DOTSeq's countReadsSingleCell() function
  gr <- readRDS(opt$gr_dir)  
  bam_list <- list.files(opt$bam_dir, pattern = "_Aligned.sortedByCoord.out_CB.bam$", recursive = TRUE, full.names = TRUE)
  dir.create(opt$mat_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (i in seq_along(bam_list)) {
    mat <- countReadsSingleCell(gr, bam_list[i])
    fname <- tools::file_path_sans_ext(basename(bam_list[i]))
    saveRDS(mat, file = file.path(opt$mat_dir, paste0(fname, ".mat.rds")))
  }
  message("Read counting done. Per-ORF matrices saved to ", opt$mat_dir)
}

if (opt$start <=2) {
  if (is.null(opt$mat_dir) || is.null(opt$out_dir) || is.null(opt$gr_dir)) {
    stop("To start with precomputed per-ORF matrices, you must provide --gr-dir --rds-dir and --out-dir")
  }
  
  # Step 3: Match per-ORF count matrices with sample metadata
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  printf <- function(...) cat(sprintf(...), "\n")
  
  selected_runs <- c(
    "SRR13125084", "SRR13125088", "SRR13125092", "SRR13125094",
    "SRR13125096", "SRR13125097", "SRR13125102", "SRR13125103",
    "SRR13125104", "SRR14530593", "SRR14530595", "SRR14530596",
    "SRR14530602", "SRR14530605", "SRR14530606")
  
  ## SRA meta (expects: Run, cell_type, treatment, ...)
  sra <- read.csv("ref/SraRunTable_PRJNA680481.csv")
  sra <- sra[sra$cell_type == "hTERT RPE-1", , drop = FALSE]

  ## Per-ORF matrices named "<Run>_Aligned.sortedByCoord.out_CB.mat.rds"
  rds_files <- file.path(matrix_dir, paste0(selected_runs, "_Aligned.sortedByCoord.out_CB.mat.rds"))
  if (!length(rds_files)) stop("No per-ORF .rds matrices found.")
  # Only keep runs present in sra
  sra_sub <- sra[sra$Run %in% selected_runs, , drop = FALSE]
  sra_sub$matrix <- file.path(matrix_dir, paste0(sra_sub$Run, "_Aligned.sortedByCoord.out_CB.mat.rds"))
  stopifnot(nrow(sra_sub) > 0, all(file.exists(sra_sub$matrix)))
  printf("Included runs: %d", nrow(sra_sub))
    
  # ## Per-ORF matrices named "<Run>_Aligned.sortedByCoord.out_CB.mat.rds"
  # rds_files <- list.files(opt$mat_dir, pattern = "_Aligned\\.sortedByCoord\\.out_CB\\.mat\\.rds$", full.names = TRUE)
  # if (!length(rds_files)) stop("No per-ORF .rds matrices found.")
  # mats_acc  <- vapply(strsplit(basename(rds_files), "_"), `[`, character(1), 1)
  # sra_sub <- sra[sra$Run %in% mats_acc, , drop = FALSE]
  # sra_sub$matrix <- file.path(opt$mat_dir, paste0(sra_sub$Run, "_Aligned.sortedByCoord.out_CB.mat.rds"))
  # stopifnot(nrow(sra_sub) > 0, all(file.exists(sra_sub$matrix)))
  # printf("Included runs: %d", nrow(sra_sub))
  
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
  
  gr <- readRDS(opt$gr_dir); stopifnot(is(gr, "GRanges"))
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
  m_counts <- LayerData(sobj, assay = "mORF", layer = "counts")
  keep_m   <- Matrix::colSums(m_counts) > opt$umi_thr
  s_m      <- subset(sobj, cells = colnames(sobj)[keep_m])
  
  s_m <- NormalizeData(s_m, layer = "counts")      # writes 'data'
  s_m <- FindVariableFeatures(s_m)
  s_m <- ScaleData(s_m, layer = "data")
  s_m <- RunPCA(s_m)
  s_m <- FindNeighbors(s_m, dims = 1:20)
  s_m <- RunUMAP(s_m, dims = 1:20)
  s_m <- FindClusters(s_m, resolution = 0.5)
  
  ## palette
  stage_colors <- c(
    "Interphase treatment"              = "#619CFF",
    "Interphase and Mitotic treatments" = "#00BA38",
    "Interphase and G0 treatments"      = "#F8766D",
    "Unknown"                           = "grey70"
  )
  
  ## UMAPs (mORF)
  p_umap_ratio <- FeaturePlot(
    s_m, features = "uORF_mORF_ratio", reduction = "umap"
  ) + scale_color_viridis(option = "D") +
    labs(title = "UMAP coloured by uORF proportion (mORF data)") +
    theme_classic()
  
  p_umap_trt <- DimPlot(s_m, reduction = "umap", group.by = "treatment") +
    scale_color_manual(values = stage_colors) +
    labs(title = "UMAP coloured by treatment (mORF data)") +
    theme_classic()
  
  print(p_umap_ratio); print(p_umap_trt)
  
  
  ## uORF-only UMAP (sparse-aware; keep cells with uORF UMIs > opt$umi_thr)
  DefaultAssay(sobj) <- "uORF"
  u_counts <- LayerData(sobj, assay = "uORF", layer = "counts")
  keep_u   <- Matrix::colSums(u_counts) > opt$umi_thr
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
  
  # Step 8: Create a WNN joint manifold
  ##     - mORF: SCTransform + PCA  (SCT_mORF)
  ##     - uORF: LogNormalize + PCA (assay = "uORF")
  
  library(future)
  options(future.globals.maxSize = 1 * 1024^3)  # 1 GB limit
  
  ## Subset to cells with sufficient uORF depth
  u_counts <- LayerData(sobj, assay = "uORF", layer = "counts")
  u_n      <- Matrix::colSums(u_counts)
  
  thr <- opt$umi_thr
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
  
  # Step 10: EMMs (probability scale) + pairwise differences
  
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
  
  # Step 11: DHARMa diagnostics
  
  set.seed(1)
  res <- simulateResiduals(fit_bb, n = 1000)
  plot(res)
  print(testUniformity(res))
  print(testDispersion(res))
  
  ## DHARMa violins by treatment
  res_df <- data.frame(Treatment = df$Treatment, dharma = res$scaledResiduals)
  p_res_dharma <- ggplot(res_df, aes(Treatment, dharma, fill = Treatment)) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.1, color = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey30") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() + theme(legend.position = "none",
                            axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(title = "DHARMa scaled residuals by Treatment", y = "Scaled residual", x = NULL)
  print(p_res_dharma)
  
  ## Pearson residuals violin
  res_pear <- residuals(fit_bb, type = "pearson")
  res_df2  <- data.frame(Treatment = df$Treatment, pearson = res_pear)
  p_res_pear <- ggplot(res_df2, aes(Treatment, pearson, fill = Treatment)) +
    geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_boxplot(width = 0.1, color = "white", outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    theme_classic() + theme(legend.position = "none",
                            axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(title = "Pearson residuals by Treatment", y = "Pearson residual", x = NULL)
  print(p_res_pear)
  
  # Step 12: Leave-one-SRR-out robustness
  
  runs <- levels(df$SRR)
  emm_full <- as.data.frame(emmeans(fit_bb, ~ Treatment, type = "response")) %>% mutate(kind = "full")
  y_col_loo <- y_col  # reuse selected column name from EMMs
  
  emm_loo <- map_dfr(runs, function(r) {
    fit <- glmmTMB(cbind(u, n - u) ~ Treatment + (1 | SRR),
                   family = betabinomial(), data = subset(df, SRR != r))
    as.data.frame(emmeans(fit, ~ Treatment, type = "response")) %>% mutate(kind = paste0("LOO_", r))
  })
  emm_all <- bind_rows(emm_full, emm_loo)
  
  p_loo <- ggplot(emm_all, aes(Treatment, .data[[y_col_loo]], color = (kind == "full"))) +
    geom_point(position = position_jitter(width = .1), alpha = .8) +
    scale_color_manual(values = c("gray60", "black"), guide = "none") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(y = "uORF fraction (u/(u+m))", x = NULL,
         title = "Leave-one-SRR-out robustness: EMMs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  print(p_loo)
  
  # Step 13: Graph Moran’s I for uORF fraction on the current graph
  
  Glist <- tryCatch(Seurat::Graphs(sobj), error = function(e) list())
  gname <- if ("weighted.nn" %in% names(Glist)) "weighted.nn" else names(Glist)[1] %||% NA_character_
  if (!is.na(gname)) {
    A <- as(Glist[[gname]], "dgCMatrix")
    y <- sobj$uORF_mORF_ratio
    Wy <- as.vector((A %*% y) / pmax(1, Matrix::rowSums(A)))
    morans_I <- suppressWarnings(cor(y, Wy, use = "pairwise.complete.obs"))
    printf("Graph Moran's I (on %s): %.4f", gname, morans_I)
  } else {
    printf("No Seurat graph found for Moran's I computation.")
  }
  
  # Step 14: Rebuild WNN on s_wnn and compare clusters vs Treatment with NMI/AMI/ARI
  
  s_wnn <- FindMultiModalNeighbors(
    s_wnn,
    reduction.list = list("pca_mORF", "pca_uORF"),
    dims.list      = list(1:20, 1:20),
    weighted.nn.name = "weighted.nn",   # Neighbor
    snn.graph.name   = "wsnn",          # Graph
    knn.graph.name   = "wknn"           # Graph
  )
  ## Check:
  names(SeuratObject::Neighbors(s_wnn))
  names(SeuratObject::Graphs(s_wnn))
  
  ## Cluster on the WNN SNN graph (Leiden) and recompute UMAP from the WNN neighbor object
  s_wnn <- FindClusters(s_wnn, graph.name = "wsnn", algorithm = 4,
                        resolution = 0.5, cluster.name = "wnn_leiden")
  s_wnn <- RunUMAP(s_wnn, nn.name = "weighted.nn", reduction.name = "wnn.umap")
  
  ## External validation: cluster-vs-treatment
  cl <- s_wnn$wnn_leiden
  tr <- s_wnn$treatment
  
  nmi <- aricode::NMI(cl, tr, variant = "max")  # normalized mutual information
  ami <- aricode::AMI(cl, tr)                   # adjusted MI (corrected for chance)
  ari <- aricode::ARI(cl, tr)                   # adjusted Rand Index (corrected for chance)
  
  cat(sprintf("NMI(cluster, Treatment) = %.3f\n", nmi))
  cat(sprintf("AMI(cluster, Treatment) = %.3f\n", ami))
  cat(sprintf("ARI(cluster, Treatment) = %.3f\n", ari))
  
  
  ## Permutation baseline for NMI/AMI/ARI (keeps cluster sizes fixed; shuffles Treatment)
  set.seed(1)
  B <- 500
  nmib <- numeric(B); amib <- numeric(B); arib <- numeric(B)
  cl <- s_wnn$wnn_leiden
  tr <- s_wnn$treatment
  for (b in seq_len(B)) {
    tr_b <- sample(tr)
    nmib[b] <- aricode::NMI(cl, tr_b, variant = "max")
    amib[b] <- aricode::AMI(cl, tr_b)
    arib[b] <- aricode::ARI(cl, tr_b)
  }
  quantile(amib, c(.5,.95,.99))
  
  p_emp <- (sum(amib >= ami) + 1) / (length(amib) + 1)  # one-sided, clustering >= observed
  p_two <- 2 * min(p_emp, 1 - p_emp)                  
  p_emp
  
  # Step 15: Save figures
  
  ggsave(file.path(fig_dir, "umap_morf_ratio.pdf"), p_umap_ratio, width = 4.75, height = 4)
  ggsave(file.path(fig_dir, "umap_morf_trt.pdf"),   p_umap_trt,   width = 6.45, height = 4)
  ggsave(file.path(fig_dir, "umap_uorf_ratio.pdf"), p_umap_uorf_ratio, width = 4.75, height = 4)
  ggsave(file.path(fig_dir, "umap_uorf_trt.pdf"),   p_umap_uorf_trt,   width = 6.45, height = 4)
  ggsave(file.path(fig_dir, "wnn_trt.pdf"),         p_wnn_trt,       width = 6.45, height = 4)
  ggsave(file.path(fig_dir, "wnn_ratio.pdf"),       p_wnn_ratio,     width = 4.75, height = 4)
  ggsave(file.path(fig_dir, "emm_by_trt.pdf"),      p_emm,           width = 4, height = 4)
  ggsave(file.path(fig_dir, "pw_diffs.pdf"),        p_pw_or,         width = 4, height = 4.75)
  ggsave(file.path(fig_dir, "dharma_violin.pdf"),   p_res_dharma,    width = 4, height = 4)
  ggsave(file.path(fig_dir, "pearson_violin.pdf"),  p_res_pear,      width = 4, height = 4)
  ggsave(file.path(fig_dir, "loo_srr_emms.pdf"),    p_loo,           width = 4, height = 4)
}

sessionInfo()