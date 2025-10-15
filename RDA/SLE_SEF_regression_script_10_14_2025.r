rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(parallel)
library(doSNOW)
library(foreach) 
library(doParallel) 
library(pbmcapply)  
library(Ake)
source("main_github_10_14_2025.R") 

####### load data ####### 
cell_type = "CD4"
covariate_df = read.csv("covariate_data.csv")
if(cell_type == "CD4"){
  # sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_cd4_sObj_sct.RDS")
  # DefaultAssay(sObj) = "SCT" #must run
  # sObj@assays$SCT$scale.data = NULL
  # sObj = PrepSCTFindMarkers(sObj)
  # sObj_metadata = sObj@meta.data
  # exprMatReal = sObj@assays$SCT$counts
  # rm(sObj)
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_github_cd4_seurat_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_github_cd4_seurat_exprMat.RDS")
} else if (cell_type == "CD8"){
  # sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/cd8_sObj_sct.RDS")
  # DefaultAssay(sObj) = "SCT" #must run
  # sObj@assays$SCT$scale.data = NULL
  # sObj = PrepSCTFindMarkers(sObj)
  # sObj_metadata = sObj@meta.data
  # exprMatReal = sObj@assays$SCT$counts
  # rm(sObj)
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/updated_github_cd8_seurat_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/updated_github_cd8_seurat_exprMat.RDS")
} else if (cell_type == "cM"){
  # sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/classical_monocyte_sObj_sct.RDS")
  # DefaultAssay(sObj) = "SCT" #must run
  # sObj@assays$SCT$scale.data = NULL
  # sObj = PrepSCTFindMarkers(sObj)
  # sObj_metadata = sObj@meta.data
  # exprMatReal = sObj@assays$SCT$counts
  # rm(sObj)
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/monocytes_seurat_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/monocytes_seurat_exprMat.RDS")
  
}
donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/monocytes_norm_results/monocytes_donors_for_all_genes_list_no_ribo.RDS")
min_donors = 50
donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
genes_of_interest = names(donors_to_use_per_gene)
print(length(genes_of_interest))



####### sef regression ####### 
p = 2
res_list <- setNames(pbmclapply(genes_of_interest, function(name) {
  tryCatch({
    exprMat1 = exprMatReal[name, , drop = F]
    stablerunSeuratCounts_Concise(exprMat = exprMat1, donor_list = donors_to_use_per_gene[[name]], sObj_meta = sObj_metadata, p = p, plot_flag = F) 
  }, error = function(e) {
    message(paste("Error in processing:", name, "-", e$message))
    return(NULL)
  })
}, mc.cores = 3, ignore.interactive = TRUE), genes_of_interest)
res_list = Filter(Negate(is.null), res_list)
pvals <- unlist(sapply(res_list, function(x) x$pval_1))
pval_df <- data.frame(p_value = pvals, row.names = names(pvals))
pval_df$bonferroni_pval = p.adjust(pval_df$p_value, method = "bonferroni")
dim(subset(pval_df, bonferroni_pval < 0.05)) 
pval_df <- pval_df %>%
  arrange(bonferroni_pval)
dim(pval_df)
head(pval_df)
