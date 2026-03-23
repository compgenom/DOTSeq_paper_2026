#!/usr/bin/env Rscript

# =============================================================================
# Script: Plotting the results from benchmarking analysis
# =============================================================================
# Input files:
# - iCOBRA RDS generated from the run_packages.R
#
# Notes:
# - Paths to input files may need to be adjusted for your system
# - There may be extra figures generated. 
# =============================================================================

library(iCOBRA)
library(ggplot2)
library(dplyr)

input_dir <- "results/benchmarking/"
output_dir <- "results/benchmarking/figures/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# FDR-TPR plots (Figure 4)
# # Please note the scale might be different from the manuscript
# # There will be extra plots instead of just 4 shown in Figure 4 in the manuscript. To generate the exact 4 plots in Figure 4, please set the following parameters instead.
# # The 4 exact combinations for Figure 4:
# param_grid <- data.frame(
#   scenario    = rep("uORF_up_mORF_down", 4),
#   gcoeff      = c(0.5, 1, 1.5, 2),
#   size_factor = c(0.5, 1, 1, 2),
#   min_size    = c(1, 1, 1, 1)
# )

regulation_type <- c("uORF_up_mORF_down")
gcoeff_vals <- seq(0.5, 2, by = 0.5)
sf_vals <- c(0.5, 1, 2)
ms_vals <- 1
param_grid <- expand.grid(gcoeff = gcoeff_vals, size_factor = sf_vals, min_size = ms_vals, scenario = regulation_type)

for (i in seq_len(nrow(param_grid))) {
  g <- param_grid$gcoeff[i]
  r <- param_grid$scenario[i]
  sf <- param_grid$size_factor[i]
  ms <- param_grid$min_size[i]
  color_scheme <- c(
    "#a58f52", "#dbb512", "#3bcdff", "#0072e1", "#a4badd", "#43474e"
  )
  cobra <- readRDS(file = paste0(input_dir, "cobradata_", r, "_g", g, "_sf", sf, "_ms", ms, ".rds"))
  names(padj(cobra)) <- c("1 Xtail", "3 Riborex", "2 anota2seq", "4 RiboDiff", "5 DOTSeq DTE", "6 DOTSeq DOU")
  cobra_perf <- calculate_performance(cobra,
                                      binary_truth = "status",
                                      aspects = c("fdrtpr", "fdrtprcurve"),
                                      splv = "none", thrs = c(0.01, 0.05, 0.1),
                                      maxsplit = 4)
  cobra_plot <- prepare_data_for_plot(cobra_perf, colorscheme = color_scheme)
  # plot_fdrtprcurve(cobra_plot, plottype = "curve")
  p <- plot_fdrtprcurve(cobra_plot, plottype = "points", title = paste0("r: ", r, " g: ", g, " sf: ", sf, " ms: ", ms))
  print(p)
  ggsave(
  filename = paste0(output_dir, "fdr_tpr_", r, "_g", g, "_sf", sf, "_ms", ms, ".png"),
  plot = p,
  width = 6,
  height = 5
)
}


# Boxplot (Figure 5)

regulation_type <- c("uORF_up_mORF_down")
gcoeff_vals <- seq(0.5, 2, by = 0.5)
sf_vals <- c(0.5, 1, 2)
ms_vals <- 1
param_grid <- expand.grid(gcoeff = gcoeff_vals, size_factor = sf_vals, min_size = ms_vals, scenario = regulation_type)

perf <- data.frame()
for (i in seq_len(nrow(param_grid))) {
  g <- param_grid$gcoeff[i]
  r <- param_grid$scenario[i]
  sf <- param_grid$size_factor[i]
  ms <- param_grid$min_size[i]
  method_colors <- c("#43474e", "#a4badd", "#0072e1", "#3bcdff", "#dbb512", "#a58f52")
  cobra <- readRDS(file = paste0(input_dir, "cobradata_", r, "_g", g, "_sf", sf, "_ms", ms, ".rds"))
  names(padj(cobra)) <- c("Xtail", "Riborex", "anota2seq", "RiboDiff", "DOTSeq DTE", "DOTSeq DOU")
  cobra_perf <- calculate_performance(cobra,
                                      binary_truth = "status",
                                      aspects = c("fdrtpr", "fdrtprcurve"),
                                      splv = "none", thrs = c(0.01, 0.05, 0.1),
                                      maxsplit = 4)
  perf <- rbind(perf, fdrtpr(cobra_perf))
}

perf$method <- factor(perf$method,
                      levels = c("DOTSeq DOU", "DOTSeq DTE", "RiboDiff", "Riborex", "anota2seq", "Xtail"))

p_box <- ggplot(perf, aes(x = method, y = FDR)) +
  geom_boxplot(aes(color = method), fill = "white", outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(aes(color = method), shape = 21, width = 0.2, fill = NA, show.legend = FALSE) +
  geom_hline(data = data.frame(thr = c("thr0.01", "thr0.05", "thr0.1"), y = c(0.01, 0.05, 0.1)),
             aes(yintercept = y), linetype = "dashed", color = "black") +
  scale_color_manual(values = method_colors, name = "Method") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11)
  ) +
  facet_wrap(~ thr, labeller = labeller(
    thr = c("thr0.01" = "Threshold 0.01", "thr0.05" = "Threshold 0.05", "thr0.1" = "Threshold 0.1")
  )) +
  ylab("False Discovery Rate (FDR)")

# Save it
ggsave(
  filename = paste0(output_dir, "boxplot_fdr.png"),
  plot = p_box,
  width = 7,
  height = 5
)

# p-value hist plot (Figure 6)

regulation_type <- c("uORF_up_mORF_down")
gcoeff_vals <- seq(0.5, 2, by = 0.5)
sf_vals <- c(0.5, 1, 2)
ms_vals <- 1
param_grid <- expand.grid(gcoeff = gcoeff_vals, size_factor = sf_vals, min_size = ms_vals, scenario = regulation_type)

padj_dfs <- data.frame()
for (i in seq_len(nrow(param_grid))) {
  g <- param_grid$gcoeff[i]
  r <- param_grid$scenario[i]
  sf <- param_grid$size_factor[i]
  ms <- param_grid$min_size[i]

  cobra <- readRDS(file = paste0(input_dir, "cobradata_", r, "_g", g, "_sf", sf, "_ms", ms, ".rds"))
  padj_df <- padj(cobra)[, c("DOU:padj", "DTE:padj")]
  padj_df$scenario <- paste0(r, "_g", g, "_sf", sf, "_ms", ms)
  padj_dfs <- rbind(padj_dfs, padj_df)
}


dou_df <- na.omit(padj_dfs[, c("DOU:padj", "scenario")])
dou_df <- dou_df %>%
  group_by(scenario) %>%
  mutate(median_padj = median(`DOU:padj`, na.rm = TRUE)) %>%
  ungroup()
# Compute order based on median_padj
scenario_order <- dou_df %>%
  distinct(scenario, median_padj) %>%
  arrange(median_padj) %>%
  pull(scenario)
# Apply reordering
dou_df$scenario <- factor(dou_df$scenario, levels = scenario_order)

dte_df <- na.omit(padj_dfs[, c("DTE:padj", "scenario")])
dte_df <- dte_df %>%
  group_by(scenario) %>%
  mutate(median_padj = median(`DTE:padj`, na.rm = TRUE)) %>%
  ungroup()
# Compute order based on median_padj
scenario_order <- dte_df %>%
  distinct(scenario, median_padj) %>%
  arrange(median_padj) %>%
  pull(scenario)
# Apply reordering
dte_df$scenario <- factor(dte_df$scenario, levels = scenario_order)


# Add a 'type' column to each
dou_df$type <- "DOU"
dte_df$type <- "DTE"

# Rename columns to a common name for plotting
dou_df <- dou_df %>% rename(padj = `DOU:padj`)
dte_df <- dte_df %>% rename(padj = `DTE:padj`)

# Combine into one long-format data frame
combined_df <- bind_rows(dou_df, dte_df)

# Facet grid plot
p_hist <- ggplot(combined_df, aes(x = padj, fill = scenario)) + 
  geom_histogram(position = "stack", bins = 30) + 
  facet_wrap(~type, scales = "free", labeller = as_labeller(c(DOU = "DOU (Local False Sign Rate)", DTE = "DTE (FDR-adjusted p-value)"))) +
  theme_bw() + 
  theme(legend.position = "none") +
  scale_fill_viridis_d(option = "G", direction = -1) +
  labs(x = "Statistical significance measures", y = "Number of ORFs")

# Save it
ggsave(
  filename = paste0(output_dir, "hist_padj.png"),
  plot = p_hist,
  width = 8,
  height = 5
)
