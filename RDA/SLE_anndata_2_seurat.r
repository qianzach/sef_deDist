# conversion from anndata to seurat script
library(dplyr) 
library(ggplot2)
library(tidyverse)
library(Seurat)
library(anndata)
library(reticulate)
library(Matrix)
filepath = "/Users/zaqian/Desktop/OLD_LUPUS/raw_lupus.h5ad" # extracts the X as the raw counts matrix
data <- anndata::read_h5ad(filename = filepath)
sc <- import("scanpy")
sc$pp$filter_cells(data, min_genes = 200)
sc$pp$filter_genes(data, min_cells = 30)
colnames(data) = data$var$feature_name
meta = data$obs
counts = data$X
rm(data)
counts <- as(counts, "CsparseMatrix")
counts = t(counts)
sObj =  CreateSeuratObject(counts = counts, meta.data = meta)
sObj = subset(sObj, subset = (self_reported_ethnicity %in% c("Asian", "European"))) # sanity check
saveRDS(sObj,"/Users/zaqian/Desktop/lupus/lupus_seurat.rds")

cell_type_of_interest ="classical monocyte" #"CD4-positive, alpha-beta T cell", "CD8-positive, alpha-beta T cell"
saveRDS(subset(sObj, subset = (cell_type == cell_type_of_interest)),  "/Users/zaqian/Desktop/density_estimation/sef_lupus_results/cell_type_specific_seurat_objects/classical_monocyte_seurat_SLE.rds")