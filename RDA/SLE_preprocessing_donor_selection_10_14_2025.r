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
library(pbmcapply)
source("main_github_10_14_2025.R") 

####### load data ####### 
cell_type = "cd8"
covariate_df = read.csv("covariate_data.csv")
if(cell_type == "cd4"){
  sObj = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/cd4_sObj_sct.RDS") # changed from TEST_cd4
} else if (cell_type == "cd8"){
  sObj = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/cd8_sObj_sct.RDS")
} else if (cell_type == "cM"){
  sObj = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/sObjs/cM_sObj_sct.RDS") # changed from "classical_monocyte"
}
DefaultAssay(sObj) = "SCT" #must run
sObj@assays$SCT$scale.data = NULL
sObj = PrepSCTFindMarkers(sObj)
sObj_metadata = sObj@meta.data
exprMatReal = sObj@assays$SCT$counts

# saveRDS(exprMatReal, paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/",cell_type,"_exprMat.RDS"))
# saveRDS(sObj_metadata, paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/",cell_type,"_metadata.RDS"))

####### preprocessing ####### 
donors_to_use_per_gene = filtration_method(sObj = sObj, gene_vector = rownames(exprMatReal),
                                                   exprMat = exprMatReal) # make function name anonymous
donors_to_use_per_gene = donors_to_use_per_gene[!grepl("^RPL|^RPS", names(donors_to_use_per_gene)) & names(donors_to_use_per_gene) != "MALAT1"]
# saveRDS(donors_to_use_per_gene, paste0(cell_type,"_donors_for_all_genes_list_no_ribo.RDS"))

min_donors = 50
donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
saveRDS(donors_to_use_per_gene, paste0(cell_type,"_donors_for_genes_processed.RDS"))


