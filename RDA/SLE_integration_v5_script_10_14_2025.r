library(Seurat)
options(future.globals.maxSize = 70 * 1024^3)  # Setting to 70 GB
cell_type = "CD8"
cell_type_specific_lupus_filepath = paste0(cell_type,"_SLE.rds") #load in cell-type specific R object
sObj = readRDS(cell_type_specific_lupus_filepath)
sObj[["percent.mt"]] <- PercentageFeatureSet(sObj, pattern = "^MT-")
sObj <- subset(sObj, subset = percent.mt < 10)

sObj[["RNA"]] <- split(sObj[["RNA"]], f = sObj$Processing_Cohort)
sObj <- SCTransform(sObj , vst.flavor = "v2", conserve.memory = T, variable.features.n = 3000, ncells = 3000) 
sObj <- RunPCA(sObj)
sObj <- RunUMAP(sObj, dims = 1:30)

sObj = IntegrateLayers(object = sObj, method = "RPCAIntegration", orig.reduction = "pca", normalization.method ="SCT")
sObj <- PrepSCTFindMarkers(object = sObj)
saveRDS(sObj, paste0(cell_type,"_sObj_sct.RDS"))
