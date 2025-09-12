####### load relevant packages ####### 
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(anndata)
library(parallel)
library(doSNOW)
library(foreach) #NCT
library(doParallel) #NCT 
library(Hmisc) # pathway analysis
library(pbmcapply)  # Provides progress bar support
library(Ake)
library(presto)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/rdaMain072025.R")  

####### load data #######
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/covariate_data.csv")
donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/cd8_donors_for_all_genes_list_no_ribo.RDS")
donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = 20)
genes_of_interest = names(donors_to_use_per_gene)
sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/cd8_seurat_metadata.RDS")
exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/cd8_seurat_exprMat.RDS")

#  ---- parallel computing  ----
print(length(names(donors_to_use_per_gene)))
print(length(genes_of_interest))
p = 2
cd8_res_hli_list <- setNames(pbmclapply(genes_of_interest, function(name) {
  tryCatch({
    exprMat1 = exprMatReal[name, , drop = F]
    stablerunSeuratCounts_Concise(exprMat = exprMat1, donor_list = donors_to_use_per_gene[[name]], sObj_meta = sObj_metadata, p = 2, plot_flag = F) 
  }, error = function(e) {
    message(paste("Error in processing:", name, "-", e$message))
    return(NULL)
  })
}, mc.cores = 3, ignore.interactive = TRUE), genes_of_interest)


#  ---- p-value computations  ----
cd8_res_hli_list_filtered = Filter(Negate(is.null), cd8_res_hli_list)
cd8_hli_pvals <- unlist(sapply(cd8_res_hli_list_filtered, function(x) x$pval_1))
cd8_hli_pval_df <- data.frame(p_value = cd8_hli_pvals, row.names = names(cd8_hli_pvals))
cd8_hli_pval_df$bonferroni_pval = p.adjust(cd8_hli_pval_df$p_value, method = "bonferroni")
dim(subset(cd8_hli_pval_df, bonferroni_pval < 0.05)) 
cd8_hli_pval_df <- cd8_hli_pval_df %>%
  arrange(bonferroni_pval)



