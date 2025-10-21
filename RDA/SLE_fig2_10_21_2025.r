# fig 2: volcano plots
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
# volcano plot function
get_volcano_plot = function(df, fc_thresh = 0.5, pval_thresh = 0.05, p_title = "Volcano Plot") {
  df$neg_log10_pval = -log10(df$p_val_adj)
  ggplot(df, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(color = "cornflowerblue", alpha = 0.8, size = 2) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      title = p_title
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}
# load results
sef_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/filter_before_sef_cM_fixed_pvals.csv", row.names = 1)
sef_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/sef_cd8_fixed_pvals.csv", row.names = 1)
sef_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/sef_cd4_fixed_pvals.csv", row.names = 1)
SCT_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/SCT_fair_comparison_pb_cM_10_15_2025.csv", row.names = 1)
SCT_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/SCT_fair_comparison_pb_cd8_10_15_2025.csv", row.names = 1)
SCT_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/cell_type_fixing_results_09_18_2025/SCT_fair_comparison_pb_cd4_10_15_2025.csv", row.names = 1)

# begin plotting volcano plots for uniquely identified pb DEGs
cell_types = c("cM","CD8+","CD4+")
pseudobulk_specific = list(cM = setdiff(SCT_cM$gene, rownames(sef_cM)),
                           cd4 = setdiff(SCT_cd4$gene, rownames(sef_cd4)),
                           cd8 = setdiff(SCT_cd8$gene, rownames(sef_cd8)))
cM_pb_unique = subset(SCT_cM, gene %in% pseudobulk_specific[["cM"]])
cd8_pb_unique = subset(SCT_cd8, gene %in% pseudobulk_specific[["cd8"]])
cd4_pb_unique = subset(SCT_cd4, gene %in% pseudobulk_specific[["cd4"]])

get_volcano_plot(cd8_pb_unique, fc_thresh = 0.2, pval_thresh = 0.05, p_title = "CD8+ PB-specific Volcano Plot")
get_volcano_plot(cM_pb_unique, fc_thresh = 0.2, pval_thresh = 0.05, p_title = "Monocyte PB-specific Volcano Plot")
get_volcano_plot(cd4_pb_unique, fc_thresh = 0.2, pval_thresh = 0.05, p_title = "CD4+ PB-specific Volcano Plot")

