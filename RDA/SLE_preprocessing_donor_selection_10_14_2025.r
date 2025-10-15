####### load relevant packages ####### 
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(parallel)
library(doSNOW)
library(foreach) 
library(doParallel) 
library(pbmcapply)  # Provides progress bar support
library(Ake)
source("main_github_10_14_2025.R") 

####### load data ####### 
cell_type = "CD4"
covariate_df = read.csv("covariate_data.csv")
if(cell_type == "CD4"){
  sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_cd4_sObj_sct.RDS")
} else if (cell_type == "CD8"){
  sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/cd8_sObj_sct.RDS")
} else if (cell_type == "cM"){
  sObj = readRDS("./density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/classical_monocyte_sObj_sct.RDS")
}
DefaultAssay(sObj) = "SCT" #must run
sObj@assays$SCT$scale.data = NULL
sObj = PrepSCTFindMarkers(sObj)
sObj_metadata = sObj@meta.data
exprMatReal = sObj@assays$SCT$counts

####### preprocessing ####### 
#  ---- median expression (use means for ties at 0)  ----
mean_median_expr <- compute_mean_median_expr(exprMatReal, sort.by = "median")
mean_median_expr <- mean_median_expr %>% arrange(desc(mean))
write.csv(mean_median_expr, file = paste0(cell_type,"_mean_median_expr_all_genes.csv"))
donors_to_use_per_gene = hongzhe_filtration_method(sObj = sObj, gene_vector = mean_median_expr$gene,
                                                   exprMat = exprMatReal) # make function name anonymous
donors_to_use_per_gene <- donors_to_use_per_gene[!grepl("^RPL|^RPS", names(donors_to_use_per_gene)) & names(donors_to_use_per_gene) != "MALAT1"]
saveRDS(donors_to_use_per_gene, paste0(cell_type,"_donors_for_all_genes_list_no_ribo.RDS"))

