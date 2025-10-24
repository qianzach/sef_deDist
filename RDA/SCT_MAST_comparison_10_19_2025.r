# comparison with MAST
rm(list = ls())  
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pbmcapply)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/rdaMain072025.R") 
print("Using SCT pseudobulk counts...")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript fair_pb.R <cell_type>\n  <cell_type> must be one of: cM, cd8, cd4", call. = FALSE)
}
cell_type <- tolower(args[[1]])
valid_ct <- c("cM", "cd8", "cd4") #in case i didn't use a capital letter
if (cell_type == "cm") cell_type <- "cM"
if (!cell_type %in% c("cM", "cd8", "cd4")) {
  stop("cell_type must be one of: cM, cd8, cd4", call. = FALSE)
}

covariate_df = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/covariate_data.csv")

if (cell_type == "cd8") {
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/cd8_sObj_sct.RDS") #must run
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj) 
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/updated_cd8_donors_for_all_genes_list_no_ribo.RDS")
  min_donors = 50
  donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
  

  
}else if(cell_type == "cd4"){
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_cd4_sObj_sct.RDS") #must run
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj) 
  donors_to_use_per_gene = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/TEST/TEST_cd4_donors_for_all_genes_list_no_ribo.RDS")
  min_donors = 50
  donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
  
  
}else if(cell_type == "cM"){
  sObj = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/classical_monocyte_sObj_sct.RDS") #filtered mitochondrial counts before, as should
  sObj@assays$SCT$scale.data = NULL
  sObj = PrepSCTFindMarkers(sObj)
}


Idents(sObj) <- sObj$disease
gene_list = names(donors_to_use_per_gene)
DefaultAssay(sObj) = "SCT"
mast_de_results <- FindMarkers(
  object = sObj,
  logfc.threshold = 0,
  ident.1 = "normal",
  ident.2 = "systemic lupus erythematosus",
  test.use = "MAST",
  features = gene_list
)
mast_de_results = subset(mast_de_results, p_val_adj < 0.05)



