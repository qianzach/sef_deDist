# comparison with SCT, pseudobulk methods
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pbmcapply)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript fair_pb.R <cell_type>\n  <cell_type> must be one of: cM, cd8, cd4", call. = FALSE)
}
cell_type <- tolower(args[[1]])
valid_ct <- c("cM", "cd8", "cd4") # in case i didn't use a capital letter
if (cell_type == "cm") cell_type <- "cM"
if (!cell_type %in% c("cM", "cd8", "cd4")) {
  stop("cell_type must be one of: cM, cd8, cd4", call. = FALSE)
}

covariate_df = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")

if (cell_type == "cd8") {
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_all_genes_list_no_ribo.RDS") #renamed to cd8 for consistency, original was "updated_cd8"
  min_donors = 50
  donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/cd8_sObj_sct.RDS") #must run
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj) 
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd8_fixed_pvals.csv", row.names = 1)
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
  
  
}else if(cell_type == "cd4"){
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_all_genes_list_no_ribo.RDS") #renamed to cd4 for consistency, original was "TEST_cd4"
  min_donors = 50
  donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/cd4_sObj_sct.RDS") #must run
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj) 
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd4_fixed_pvals.csv", row.names = 1)
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
  
}else if(cell_type == "cM"){
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_all_genes_list_no_ribo.RDS") #renamed to cM for consistency, original was "monocytes"
  min_donors = 50
  donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
  genes_of_interest = names(donors_to_use_per_gene)
  print(length(genes_of_interest))
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/cM_sObj_sct.RDS") #filtered mitochondrial counts before, as should; original was "classical_monocyte"
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj)
  
  sef_res = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cM_fixed_pvals.csv", row.names = 1) # renamed to just "sef_cM", original was "filter_before_" 
  sig = subset(sef_res, bonferroni_pval < 0.05)
  dim(sig)
}

# revised pseudobulk
gene_list = names(donors_to_use_per_gene)
pb_sObj <- AggregateExpression(sObj, assays = "SCT", return.seurat = T, group.by = c("donor_id","disease"))
pb_sObj$donor_id <- sub("^g", "", pb_sObj$donor_id)
Idents(pb_sObj) <- pb_sObj$disease

fair_pseudo_wilcoxon = function(pb_sObj, donors_i, gene_i){
  single_gene_sObj = subset(pb_sObj,
                            subset  = donor_id %in% donors_i,
                            features = gene_i)
  design_matrix = cbind(single_gene_sObj@meta.data[,2:3],
                        expression = as.vector(single_gene_sObj@assays$SCT@layers$data)) #use on pseudobulk values
  
  x = design_matrix$expression[design_matrix$disease == "normal"]
  y = design_matrix$expression[design_matrix$disease == "systemic lupus erythematosus"]
  wilcox_res <- wilcox.test( # same thing
    x = x,
    y = y,
    alternative = "two.sided"
  )
  logFC = mean(log2(y + 1)) - mean(log2(x + 1))
  return(
    data.frame(
      gene = gene_i,
      p_value = wilcox_res$p.value,
      statistic = unname(wilcox_res$statistic),
      logFC = logFC,
      n_normal = sum(design_matrix$disease == "normal"),
      n_sle = sum(design_matrix$disease == "systemic lupus erythematosus"),
      stringsAsFactors = FALSE
    )
  )
}

results_list <- pbmclapply(
  genes_of_interest,
  mc.cores = 3,
  function(g) {
    donors_i <- donors_to_use_per_gene[[g]]
    try(
      fair_pseudo_wilcoxon(pb_sObj, donors_i, g),
      silent = TRUE
    )
  }
)
results_df = do.call(rbind, Filter(Negate(is.null), results_list))
results_df$p_val_adj = p.adjust(results_df$p_value, method = "bonferroni")
revised_pseudobulk = results_df %>%
  arrange(p_val_adj)
print(dim(subset(revised_pseudobulk, p_val_adj < 0.05)))

if (cell_type == "cd8") {
  write.csv(subset(revised_pseudobulk, p_val_adj < 0.05), "/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd8_10_15_2025.csv")
  # sig_pb = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd8_10_15_2025.csv", row.names = 1)
  dim(sig_pb)
} else if(cell_type == "cd4"){
  write.csv(subset(revised_pseudobulk, p_val_adj < 0.05), "/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd4_10_15_2025.csv")
  # sig_pb = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd4_10_15_2025.csv", row.names = 1)
  dim(sig_pb)
} else if(cell_type == "cM"){
  write.csv(subset(revised_pseudobulk, p_val_adj < 0.05), "/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cM_10_15_2025.csv")
  # sig_pb = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cM_10_15_2025.csv", row.names = 1)
  dim(sig_pb)
}
