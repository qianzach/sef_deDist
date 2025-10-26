## Data Dictionary for *Differential Density Analysis in Single-Cell Genomics Using Specially Designed Exponential Families*

In this section, we provide information on the structure of downstream analysis and reproducibility of figures for the manuscript.

--- 

### 1. Data Overview
#### Single Cell Data
The original raw data from the original Perez et al (2022) study can be accessed on the [Chan-Zuckerberg Initiative (CZI)
Cell Science's public domain](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2) in the `.h5ad` format (AnnData).

It contains cell-specific metadata, including batch, sex, age, disease status, cell_type, ethnicity etc. 
We provide our processed real data for each cell type in Zenodo, which can be accessed using the link in the `README` contained in the Github repository.
In our study, we conduct cell-type specific normalization using standard settings with `Seurat v5` in `R`.
This includes SCTransform-based integration to adjust for batch effects, and using `PrepSCTFindMarkers()` to provide corrected counts.

---

### 2. File Overview
We briefly explain the file purpose in our data analysis and simulations.

#### Simulations
The simulations section is split into three `R` files. We have a main file that defines all functions used in simulations. The remaining two files contain the main simulation results and that from the supplementary materials.

| File Name | Format     | Description |
|---------------|----------|------------|
| `simulations_main_10_13_2025`  | `.R`   | `main` file defining functions used in simulations |
| `refined_simulations_10_13_2025`  | `.R`   | Workflow for Figure 1: QQ-Plots and power analysis in main contents for Poisson-Gamma and ZINB |
| `refined_supp_simulations_10_13_2025`  | `.R`   | Supplementary tests and figure construction code (Figures S.1-S.4) |
---

#### Real Data Analysis
In the real data analysis section, we provide code to demonstrate the main real data contributions using the proposed SEF methodology as well as replicate the key figures with respect to these data.
This includes the regression pipeline and enrichment figures. 

Within the `walkthrough` directory of the repository, we include the main function file, denoted as `main_github_10_14_2025.R`, which defines all the functions needed to conduct our analysis.
We apply these in practice in the walkthrough titled `SLE_SEF_regression_walkthrough_10_24_2025.R`, where we test for distributional differences between SLE and control populations with the CD8+ T-cell type.
This walkthrough also includes code snippets to reproduce enrichment analysis. This can be replicated for the two other cell types. 

In `SLE_fig3_10_21_2025.R`, we provide code snippets that reproduce figures 3, S.5, S.6, and S.7,
which illustrate carrier densities, individual-level densities, and group-wise density comparisons for uniquely-identified genes (S100A4, TSC22D3, MT2A, and LTB.)

In the remaining files, we provide scripts run on the command line for the pseudobulk data and related comparisons. 

*Note on enrichment results: The Gene Ontology resource has updated its term counts and annotations since the original downstream analysis was conducted.
As a result, a few of the pathways of interest have more general or specific distinctions of our original pathways.
This includes two pathways: cytokine production in monocytes and regulation of protein complex assembly in CD4+ T-cells. Despite this, our results retain consistency in interpretation.*


| File Name | Format     | Description |
|---------------|----------|------------|
| `SLE_SEF_regression_walkthrough_10_24_2025`      | `.R`  | Runs cell-type specific SEF regression modeling and testing; code to reproduce Figure 5 |
| `SLE_compare_methods_10_19_2025`  | `.R`   | Procedural pipeline for pseudobulk methods |
| `SLE_table1_fig2_10_21_2025`  | `.R`   | Workflow to reproduce Table 1 numbers and Figure 2 |
| `SLE_validation_comparison_perez_10_21_2025`  | `.R`   | Workflow to reproduce proportions in Table 2 |
| `SLE_fig3_figS7_10_21_2025`  | `.R`   | Workflow to reproduce Figure 3; Figures S.5, S.6, S.7 |
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
| `CELLTYPE_donors_for_genes_processed`      | `.RDS`  | List of arrays, each element contains an array of donors that have at least 100 cells and 20% non-zero values for all genes tested per cell type |
| `pb_CELLTYPE_sObj`         | `.RDS`   | Pseudobulked Seurat object after using standard `AggregateExpression()` |
| `sef_CELLTYPE_fixed_pvals`         | `.csv`   | Cell type-specific differentially distributed genes with adjusted p-values (Bonferroni) using SEF regression |
| `SCT_fair_comparison_pb_CELLTYPE_10_15_2025`         | `.csv`   | Cell type specific DEGs with adjusted p-values (Bonferroni) pseudobulked data and log fold change values |
| `all_CELLTYPE_interesting_pathways`         | `.csv`   | Includes `clusterProfiler` enrichment analysis results |
| `perez_validation_DEGs`         | `.csv`   | File including DEGs from original Perez et al (2022) study's Supplementary Table 3 |

---
