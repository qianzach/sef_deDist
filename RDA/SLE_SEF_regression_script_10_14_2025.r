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
cell_type = "CD4"
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")
if(cell_type == "CD4"){
  sObj_metadata = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_metadata.RDS") #changed from "TEST_cd4_seurat"
  exprMatReal = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_exprMat.RDS") #changed from "TEST_cd4_seurat"
  donors_to_use_per_gene = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_all_genes_list_no_ribo.RDS")
} else if (cell_type == "CD8"){
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

sig = subset(pval_df, bonferroni_pval < 0.05)
er1 = run_go_enrichment_plot(gene_symbols = row.names(sig), maxGSSize = 500, cutoff = 0.6)


# WARNING: GO IDs have updated since original downstream analysis (GO pathways have slightly changed in size)
# results are, for the most part, the same except for more granular distinctions of provided pathways (ex. cytokine prod in cM, regulation of protein complex assembly in cd4)
# interesting pathways denote the ones identified by results
#saved only the results, reconstruct using:
if(cell_type == "cM"){
  interesting_pathways = c("positive regulation of cytokine production","regulation of innate immune response", "regulation of endocytosis",
                           "antigen processing and presentation", "regulation of leukocyte proliferation")
  res = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/results/enrichment_res/all_cM_interesting_pathways.csv", row.names = 1)
  res$Count <- sapply(strsplit(res$geneID, "/"), length)
  ego_reconstructed <- new("enrichResult")
  ego_reconstructed@result <- res
  ego_reconstructed@organism <- "Homo sapiens"
  ego_reconstructed@keytype <- "ENTREZID"
  ego_reconstructed@ontology <- "BP"
  ego_reconstructed@readable <- TRUE
  barplot(height = ego_reconstructed) + ggtitle("Monocytes") + guides(size = "none")
  
} else if (cell_type == "cd8"){
  interesting_pathways = c("regulation of lymphocyte apoptotic process","regulation of T cell activation", "leukocyte mediated cytotoxicity",
                           "peptide antigen assembly with MHC protein complex", "positive regulation of leukocyte cell-cell adhesion")
  res = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/results/enrichment_res/all_cd8_interesting_pathways.csv", row.names = 1)
  res$Count <- sapply(strsplit(res$geneID, "/"), length)
  ego_reconstructed <- new("enrichResult")
  ego_reconstructed@result <- res
  ego_reconstructed@organism <- "Homo sapiens"
  ego_reconstructed@keytype <- "ENTREZID"
  ego_reconstructed@ontology <- "BP"
  ego_reconstructed@readable <- TRUE
  barplot(height = ego_reconstructed) + ggtitle("CD8+") + guides(size = "none")
} else if(cell_type == "cd4"){
  interesting_pathways = c("regulation of T cell activation","antigen processing and presentation", "lymphocyte differentiation",
                           "regulation of protein-containing complex assembly", "positive regulation of neutrophil migration")
  res = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/results/enrichment_res/all_cd4_interesting_pathways.csv", row.names = 1)
  res$Count <- sapply(strsplit(res$geneID, "/"), length)
  ego_reconstructed <- new("enrichResult")
  ego_reconstructed@result <- res
  ego_reconstructed@organism <- "Homo sapiens"
  ego_reconstructed@keytype <- "ENTREZID"
  ego_reconstructed@ontology <- "BP"
  ego_reconstructed@readable <- TRUE
  barplot(height = ego_reconstructed) + ggtitle("CD4+") + guides(size = "none")
}

