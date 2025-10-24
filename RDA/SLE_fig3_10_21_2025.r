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
source("./density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 

####### load data ####### 
# cell_type = "cd4"
cell_type = "cM"
# cell_type = "cd8"
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")

if(cell_type == "cd4"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_metadata.RDS") #changed from "TEST_cd4_seurat"
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_exprMat.RDS") #changed from "TEST_cd4_seurat"
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_all_genes_list_no_ribo.RDS")
} else if (cell_type == "cd8"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_metadata.RDS") #changed from "updated_cd8_seurat"
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_exprMat.RDS") #changed from "updated_cd8_seurat"
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_all_genes_list_no_ribo.RDS") #renamed to cd8 for consistency, original was "updated_cd8"
} else if (cell_type == "cM"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_sObj_metadata.RDS")
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_exprMatReal.RDS")
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_all_genes_list_no_ribo.RDS") #renamed to cM for consistency, original was "monocytes"
}

min_donors = 50
donors_to_use_per_gene = filter_genes_by_disease_donor_counts(donors_to_use_per_gene, covariate_df, min_donors = min_donors)
genes_of_interest = names(donors_to_use_per_gene)
print(length(genes_of_interest))

# specific_genes = c("S100A9", "ISG15","CYBA", S100A4) # cM; S100A4 S.7
# specific_genes = c("HLA-A","NKG7", "LGALS1", MT2A, LTB) # cd8 MT2A, LTB S.7
# specific_genes = c("B2M", "IL32","KLF6", TSC22D3) # cd4 TSC22D3 S.7

genes_of_interest = names(donors_to_use_per_gene)
don = "1101"
print(length(genes_of_interest))
genes_of_interest = specific_genes
p = 2
res_hli_list <- setNames(pbmclapply(genes_of_interest, function(name) { #run sef regression for the genes of interest
  tryCatch({
    exprMat1 = exprMatReal[name, , drop = F]
    stablerunSeuratCounts_Concise(exprMat = exprMat1, donor_list = donors_to_use_per_gene[[name]], sObj_meta = sObj_metadata, p = p, plot_flag = F) 
  }, error = function(e) {
    message(paste("Error in processing:", name, "-", e$message))
    return(NULL)
  })
}, mc.cores = 3, ignore.interactive = TRUE), genes_of_interest)


for (g in specific_genes) {
  print(g)
  carrier_test = plot_hli_carrier_with_hist(exprMat = exprMatReal, gene = g, sObj_meta = sObj_metadata, donors_to_use_list = donors_to_use_per_gene)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-carrier-",g, ".pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(carrier_test)
  dev.off()
  
  comparison_test = group_model_comparison_plot(sef_object = res_hli_list, gene = g)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-model-comparison-",g, ".pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(comparison_test)
  dev.off()
  
  ind_test = get_individual_sef_counts_hli_model(donor_id = don, gene = g, exprMat = exprMatReal, sObj_meta = sObj_metadata, donor_list = donors_to_use_per_gene)
  GraphName = paste0("./density_estimation/sef_lupus_results/exploration/to_github_repo/test_PDFs/",cell_type,"-dke-",g,"-donor-",don,"-ind-model.pdf" )
  pdf(GraphName, width = 5, height = 5)
  print(ind_test)
  dev.off()
}

