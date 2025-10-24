# validation with perez et al supplementary genes
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
source("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/main_github_10_14_2025.R") 

perez = read.csv("./density_estimation/sef_lupus_results/exploration/perez_validation_DEGs.csv", row.names = 1)
covariate_df = read.csv("./density_estimation/sef_lupus_results/exploration/to_github_repo/DATA/covariate_data.csv")
donors_cd8 = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd8_donors_for_all_genes_list_no_ribo.RDS") #renamed to cd8 for consistency, original was "updated_cd8"
donors_cd4 = readRDS("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cd4_donors_for_all_genes_list_no_ribo.RDS") #renamed to cd4 for consistency, original was "TEST_cd4"
donors_cM = readRDS("./density_estimation/sef_lupus_results/exploration/to_github_repo/objects/cM_donors_for_all_genes_list_no_ribo.RDS") #renamed to cM for consistency, original was "monocytes"
min_donors = 50
donors_cd8 = filter_genes_by_disease_donor_counts(donors_cd8, covariate_df, min_donors = min_donors)
donors_cd4 = filter_genes_by_disease_donor_counts(donors_cd4, covariate_df, min_donors = min_donors)
donors_cM = filter_genes_by_disease_donor_counts(donors_cM, covariate_df, min_donors = min_donors)
cd8_validated_genes = intersect(rownames(perez)[perez$T8 == 1], names(donors_cd8))
cd4_validated_genes= intersect(rownames(perez)[perez$T4 == 1], names(donors_cd4))
cM_validated_genes = intersect(rownames(perez)[perez$cM == 1], names(donors_cM))

# load results
sef_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cM_fixed_pvals.csv", row.names = 1)
sef_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd8_fixed_pvals.csv", row.names = 1)
sef_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/sef_cd4_fixed_pvals.csv", row.names = 1)
SCT_cM = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cM_10_15_2025.csv", row.names = 1)
SCT_cd8 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd8_10_15_2025.csv", row.names = 1)
SCT_cd4 = read.csv("/Users/zaqian/Desktop/density_estimation/sef_lupus_results/exploration/to_github_repo/results/SCT_fair_comparison_pb_cd4_10_15_2025.csv", row.names = 1)

# cd8
sum(row.names(sef_cd8) %in% cd8_validated_genes)/length(cd8_validated_genes)
sum(SCT_cd8$gene %in% cd8_validated_genes)/length(cd8_validated_genes)

#cd4
sum(row.names(sef_cd4) %in% cd4_validated_genes)/length(cd4_validated_genes) 
sum(SCT_cd4$gene %in% cd4_validated_genes)/length(cd4_validated_genes)

# cM
sum(row.names(sef_cM) %in% cM_validated_genes)/length(cM_validated_genes)
sum(SCT_cM$gene %in% cM_validated_genes)/length(cM_validated_genes)

