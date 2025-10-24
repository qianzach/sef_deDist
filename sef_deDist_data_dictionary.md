## Data Dictionary for *Differential Density Analysis in Single-Cell Genomics Using Specially Designed Exponential Families*

This markdown document provides information on structure of downstream analysis and reproducibility of figures.

--- 

### 1. Data Overview
#### Single Cell Data
The original raw data from the original Perez et al (2022) study can be accessed on the Chan-Zuckerberg Initiative (CZI)
Cell Science's public domain in the `.h5ad` format (AnnData v0.10). It contains cell-specific metadata, including batch, sex, age, disease status, cell_type, ethnicity etc.

In our study, we conduct cell-type specific normalization using `Seurat v5` in `R`. As such, we subset the original raw data by cell type (classical monocytes, CD4+ T-cells, CD8+ T-cells)
and save the individual Seurat objects. Using these Seurat objects, we then undergo SCTransform-based integration using the aforementioned version `v5` of `Seurat`
to adjust for batch effects. Our script was run in an R environment for each specific cell type. 

---

### 2. File Overview
We briefly explain the file purpose in our data analysis and simulations.

#### Simulations

| File Name | Format     | Description |
|---------------|----------|------------|
| `simulations_main_10_13_2025`  | `.R`   | `main` file defining functions used in simulations |
| `refined_simulations_10_13_2025`  | `.R`   | Workflow for Figure 1: QQ-Plots and power analysis in main contents for Poisson-Gamma and ZINB |
| `refined_supp_simulations_10_13_2025`  | `.R`   | Supplementary tests and figure construction code (Figures S.1-S.4) |
---

#### Real Data Analysis
In the real data analysis section, we provide scripts and code for preprocessing, as well as downstream analysis (inference, enrichment, etc.) and method comparisons.

| File Name | Format     | Description |
|---------------|----------|------------|
| `SLE_get_raw_counts` | `.ipynb`     | Extracts raw count data and stores in CZI-downloaded AnnData object |
| `SLE_anndata_2_seurat` | `.R`     | Script converts from original AnnData to Seurat object and does basic preprocessing; creates cell-type specific Seurat objects; run on command line or `R` |
| `SLE_integration_v5_script_10_14_2025` | `.R`  | Performs SCTransform-based integration for each cell type; run on command line or `R` |
| `SLE_preprocessing_donor_selection_10_14_2925`| `.R`   | Conducts preprocessing and identifies, for each gene in a cell type, donors that satisfy the filtration criteria for downstream analysis |
| `SLE_SEF_regression_script_10_14_2025`      | `.R`  | Runs cell-type specific SEF regression modeling and testing; code to reproduce Figure 5 |
| `SLE_compare_methods_10_19_2025`  | `.R`   | Procedural pipeline for pseudobulk methods |
| `SLE_MAST_comparison_10_19_2025`  | `.R`   | Procedural pipeline for MAST cell-level method |
| `SLE_table1_fig2_10_21_2025`  | `.R`   | Workflow to reproduce Table 1 numbers and Figure 2 |
| `SLE_validation_comparison_perez_10_21_2025`  | `.R`   | Workflow to reproduce proportions in Table 2 |
| `SLE_fig3_10_21_2025`  | `.R`   | Workflow to reproduce Figure 3; Figures S.5, S.6, S.7 |
| `SLE_permutation_script_10_14_2025`         | `.R`   | Script to conduct permutation test and reproduce Figure 4  (*note: memory-intensive*); run on command line |
| `main_github_10_14_2025`  | `.R`   | `main` file defining functions used in real data analysis and downstream analysis |

---

### 3. Object Overview
In this section, we overview the formatted data objects used in our simulations and real data analysis. Without loss of generality, we mention cell-type specific objects once. 

| Object Name | Format     | Description |
|---------------|----------|------------|
| `covariate_data`           | `.csv`     | Covariate data with donor ID and disease status specification|
| `CELLTYPE_seurat_exprMat` | `.RDS`  | SCTransform, library-size corrected counts matrix (genes x cell) |
| `CELLTYPE_seurat_metadata`| `.RDS`   | Contains metadata from original data (cell-level rows) |
| `CELLTYPE_donors_for_all_genes_list_no_ribo`      | `.RDS`  | List of arrays, each element contains an array of donors that satisfy filtration criteria for a given gene |
| `sef_CELLTYPE_fixed_pvals`         | `.csv`   | Cell type-specific differentially distributed genes with adjusted p-values (Bonferroni) using SEF regression |
| `SCT_fair_comparison_pb_CELLTYPE_10_15_2025`         | `.csv`   | Cell type specific DEGs with adjusted p-values (Bonferroni) pseudobulked data and log fold change values |
| `all_CELLTYPE_interesting_pathways`         | `.csv`   | Includes `clusterProfiler` enrichment analysis results |
| `perez_validation_DEGs`         | `.csv`   | File including DEGs from original Perez et al (2022) study's Supplementary Table 3 |

---
