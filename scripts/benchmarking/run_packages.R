#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(iCOBRA)
library(xtail)
library(riborex)
library(anota2seq)
library(DOTSeq)

input_dir <- "data/simdata/"
output_dir <- "results/benchmarking/"

dir.create(output_dir, recursive = TRUE)

regulation_type <- c("uORF_up_mORF_down", "uORF_down_mORF_up")
sf_vals <- c(0.5, 1, 2, 4)
ms_vals <- c(0.5, 1, 2, 4)
gcoeff_vals <- seq(0.5, 2, by = 0.5)

param_grid <- expand.grid(gcoeff = gcoeff_vals, size_factor = sf_vals, min_size = ms_vals, scenario = regulation_type)

for (i in seq_len(nrow(param_grid))) {
    # Read simulated data
    g <- param_grid$gcoeff[i]
    r <- param_grid$scenario[i]
    sf <- param_grid$size_factor[i]
    ms <- param_grid$min_size[i]
    simData <- readRDS(
      file.path(
        input_dir, 
        paste0(r, "_g", g, "_sf", sf, "_ms", ms, ".rds")
      )
    )
    
    # Get counts and split Ribo and RNA
    cnt <- assay(getDOU(simData))
    ribo_cols <- grep("\\.ribo$", colnames(cnt))
    rna_cols <- grep("\\.rna$", colnames(cnt))
    ribocnt <- cnt[, ribo_cols]
    rnacnt <- cnt[, rna_cols]
    
    # Define condition vector
    ribocond <- sub(".*\\.(condition[01])\\..*", "\\1", colnames(ribocnt))
    rnacond <- sub(".*\\.(condition[01])\\..*", "\\1", colnames(rnacnt))
    
    # Prepare experiment outline for RiboDiff
    exp_outline <- colData(getDOU(simData))
    exp_outline$replicate <- NULL
    exp_outline$batch <- NULL
    exp_outline$Sample <- rownames(exp_outline)
    names(exp_outline) <- c("Conditions", "Data_Type", "Samples")
    exp_outline$Data_Type <- as.character(exp_outline$Data_Type)
    exp_outline$Data_Type[exp_outline$Data_Type == "0"] <- "Ribo-Seq"
    exp_outline$Data_Type[exp_outline$Data_Type == "1"] <- "RNA-Seq"
    exp_outline$Conditions <- as.character(exp_outline$Conditions)
    exp_outline$Conditions[exp_outline$Conditions == "0"] <- "Control"
    exp_outline$Conditions[exp_outline$Conditions == "1"] <- "Treated"
    exp_outline <- exp_outline[c("Samples", "Data_Type", "Conditions")]
    # exp_outline_path <- paste0("data/", r, "_g", g, "_sf", sf, "_exp_outline.txt")
    exp_outline_path <- "data/simdata/exp_outline.txt"
    write.csv(exp_outline, row.names = FALSE, quote = FALSE, file = exp_outline_path)
    
    # Prepare count data
    cnt <- cbind(Entry = rownames(cnt), cnt)
    write.table(cnt, file = "data/simdata/counts.txt", row.names = FALSE, quote = FALSE, sep = "\t")
    
    if (identical(ribocond, rnacond)) {
        # Run RiboDiff
        system("/opt/miniconda2/bin/python /opt/RiboDiff/scripts/TE.py -e data/simdata/exp_outline.txt -c data/simdata/counts.txt -o results/benchmarking/ribodiff.out")
        res_rd <- read.csv("results/benchmarking/ribodiff.out", sep = "\t")
      
        # Run Riborex
        start_rrx <- Sys.time()
        res_rrx <- riborex(rnaCntTable = rnacnt, riboCntTable = ribocnt, rnaCond = rnacond, riboCond = ribocond)
        end_rrx <- Sys.time()
        elapsed_rrx <- DOTSeq:::runtime(end_rrx, start_rrx)
        
        if (!is.null(elapsed_rrx$mins)) {
          message(sprintf("Riborex runtime: %d mins %.3f secs", elapsed_rrx$mins, elapsed_rrx$secs))
        } else {
          message(sprintf("Riborex runtime: %.3f secs", elapsed_rrx$secs))
        }
        
        # Run xtail
        start_xt <- Sys.time()
        res_xt <- xtail(mrna = rnacnt, rpf = ribocnt, condition = ribocond, bins = 1000, threads = 1)
        res_xt <- resultsTable(res_xt)
        end_xt <- Sys.time()
        elapsed_xt <- DOTSeq:::runtime(end_xt, start_xt)
        
        if (!is.null(elapsed_xt$mins)) {
          message(sprintf("xtail runtime: %d mins %.3f secs", elapsed_xt$mins, elapsed_xt$secs))
        } else {
          message(sprintf("xtail runtime: %.3f secs", elapsed_xt$secs))
        }
        
        # Run anota2seq
        start_apv <- Sys.time()
        ads <- anota2seqDataSetFromMatrix(dataP = ribocnt, dataT = rnacnt, phenoVec = ribocond, dataType = "RNAseq", normalize = TRUE)
        ads <- anota2seqAnalyze(Anota2seqDataSet = ads, analysis = "translation")
        res_apv <- anota2seqGetOutput(
            ads, analysis = "translation",
            output = "full",
            selContrast = 1,
            getRVM = TRUE
        )
        end_apv <- Sys.time()
        elapsed_apv <- DOTSeq:::runtime(end_apv, start_apv)
        
        if (!is.null(elapsed_apv$mins)) {
          message(sprintf("anota2seq runtime: %d mins %.3f secs", elapsed_apv$mins, elapsed_apv$secs))
        } else {
          message(sprintf("anota2seq runtime: %.3f secs", elapsed_apv$secs))
        }
        
        # Run DOTSeq
        start_d <- Sys.time()
        fmla <- ~ condition * strategy + batch
        d <- DOTSeq(datasets = simData, formula = fmla)
        res_d <- getContrasts(d, type = "interaction")
        end_d <- Sys.time()
        elapsed_d <- DOTSeq:::runtime(end_d, start_d)
        
        if (!is.null(elapsed_d$mins)) {
          message(sprintf("DOTSeq runtime: %d mins %.3f secs", elapsed_d$mins, elapsed_d$secs))
        } else {
          message(sprintf("DOTSeq runtime: %.3f secs", elapsed_d$secs))
        }
        
    } else {
      stop("Ribo-seq and RNA-seq sample indices didn't match.")
    }
    
    # Create a padj dataframe with the full list of features
    all_features <- rownames(rowData(getDOU(simData)))
    all_status <- rowData(getDOU(simData))$status
    
    # Number of features
    n <- length(all_features)
    
    # Create an empty matrix first
    padj_matrix <- matrix(NA_real_, nrow = n, ncol = 6)
    colnames(padj_matrix) <- c("Xtail:padj", "Riborex:padj", "anota2seq:padj", "RiboDiff:padj", "DTE:padj", "DOU:padj")
    
    # Convert to data frame and assign row names
    padj_df <- as.data.frame(padj_matrix)
    rownames(padj_df) <- all_features
    
    # Fill in values where available
    padj_df[rownames(res_xt), "Xtail:padj"] <- res_xt$pvalue.adjust
    padj_df[rownames(res_rrx), "Riborex:padj"] <- res_rrx$padj
    padj_df[rownames(res_apv), "anota2seq:padj"] <- res_apv[, "apvRvmPAdj"]
    padj_df[res_rd$geneIDs, "RiboDiff:padj"] <- res_rd$padj
    padj_df[res_d$DTE$orf_id, "DTE:padj"] <- res_d$DTE$padj
    padj_df[res_d$DOU$orf_id, "DOU:padj"] <- res_d$DOU$lfsr
    
    # Create COBRAData
    cobra <- COBRAData(
        padj = padj_df,
        truth = data.frame(row.names = all_features, status = all_status)
    )
    
    saveRDS(cobra, file = file.path(output_dir, paste0("cobradata_", r, "_g", g, "_sf", sf, "_ms", ms, ".rds")))
}

invisible(file.remove(c("data/simdata/exp_outline.txt", "data/simdata/counts.txt", "results/benchmarking/ribodiff.out")))
