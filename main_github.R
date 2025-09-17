library(dplyr) 
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(MASS) 
library(patchwork)
library(grid)
library(parallel)
library(truncnorm)
library(Seurat)
library(pbmcapply)
library(anndata)
library(parallel)
library(doSNOW)
library(foreach) #NCT
library(doParallel) #NCT 
library(Hmisc)
library(Ake) #new KDE approach


#  ---- BIN COUNT HELPER (USED IN RDA AS WELL) ----
get_bin_counts = function(Y, bin_edges, K) {
  # Function to get bin counts for a whole matrix 
  # Flatten the matrix to use cut and table efficiently
  Y_flat = as.vector(t(Y)) 
  # Use cut() to classify each element of Y into bins defined by bin_edges
  bin_indices <- cut(Y_flat, breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
  # Count the occurrences in each bin per row
  counts_matrix <- matrix(0, nrow = nrow(Y), ncol = K) 
  for (i in 1:nrow(Y)) {
    row_indices = bin_indices[ ((i - 1) * ncol(Y) + 1):(i * ncol(Y)) ]
    counts_matrix[i, ] = tabulate(row_indices, nbins = K) 
  }
  return(counts_matrix)
}

get_countdata_bin_counts = function(Y, grid_points){
  # Used for single cell data
  
  n = length(Y)
  K = length(grid_points)
  
  # Initialize an n x K matrix to store bin counts
  counts_matrix = matrix(0, nrow = n, ncol = K)
  
  for (i in seq_along(Y)) {
    # Factor ensures tabulation only over grid_points
    tab <- table(factor(Y[[i]], levels = grid_points))
    counts_matrix[i, ] <- as.numeric(tab)
  }
  
  return(counts_matrix)
}

get_sc_bin_counts = function(Y, bin_edges, K){
  # Used for single cell data or list of arrays (varying cell counts)
  n = length(Y)
  counts_matrix = matrix(0, nrow = n, ncol = K) # Initialize an n x K matrix to store bin counts
  
  # Loop over each element in the list
  for (i in seq_along(Y)) {
    # Use cut() to classify the elements of Y[[i]] into bins
    bin_indices <- cut(Y[[i]], breaks = bin_edges, include.lowest = TRUE, labels = FALSE)
    
    # Count occurrences in each bin
    counts_matrix[i, ] = tabulate(bin_indices, nbins = K)
  }
  return(counts_matrix)
}

#  ---- PROCESSING HELPER FUNCTIONS  ----
subset_exprMat_by_donors <- function(exprMat, metadata, donor_vector) {
  ##EXPECTS A VECTOR/ARRAY NOT MATRIX FOR exprMat
  selected_cells <- rownames(metadata)[metadata$donor_id %in% donor_vector]
  exprMat_subset <- exprMat[, selected_cells, drop = FALSE]
  return(exprMat_subset)
}

subset_metadata_by_donors <- function(sObj, donor_vector) {
  # Identify the matching cell indices
  cell_indices <- which(sObj$donor_id %in% donor_vector)
  
  # Subset exprMat (genes x cells)
  sObj_meta <- sObj@meta.data[cell_indices, , drop = F]
  
  return(sObj_meta)
}

filter_genes_by_disease_donor_counts <- function(donors_to_use_per_gene, covariate_df, min_donors = 20) {
  if (!all(c("donor_id", "disease_status") %in% colnames(covariate_df))) {
    stop("covariate_df must contain 'donor_id' and 'disease_status' columns.")
  }
  disease_levels <- sort(unique(covariate_df$disease_status))  # Ensure consistent level order
  
  if (length(disease_levels) != 2) {
    stop("disease_status must be binary.")
  }
  
  filtered_donors <- Filter(function(donors) {
    matched_df <- covariate_df[covariate_df$donor_id %in% donors, ]
    
    # Force factor to include both disease levels
    matched_df$disease_status <- factor(matched_df$disease_status, levels = disease_levels)
    
    donor_counts <- table(matched_df$disease_status)
    
    all(donor_counts >= min_donors)
  }, donors_to_use_per_gene)
  
  return(filtered_donors)
}

#  ---- RUN MODEL FUNCTIONS  ----
sef_regression = function(exprMat, sObj_meta, donor_list, p = 2, plot_flag = F, h = NULL){
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(exprMat = exprMat,metadata = sObj_meta, donor_vector = donor_list))
  #as.vector(subset_exprMat_by_donors(sObj = sObj, exprMat = exprMat, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  #apply poisson VST based on sqrt + 1/2 transformation
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation
  #Gaussian KDE evaluated at integers
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ #uneven bins
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
      print(K)
    } else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    dev.off()
    carrier_est = carrier$est.fn
    # est_grid = carrier$eval.points #will be equivalent
    # est_midpoints = est_grid
    # K = unique_pts
    # K = length(est_midpoints)
    hist_plot_breaks = carrier$hist$breaks
  } else{
    K = 100
    print(paste0("greater than or equal to 100 unique points: using equidistant bins as points for Gaussian KDE; carrier density estimate (", unique_pts," unique points), setting K bin count to: ", K))
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  if (plot_flag == T) {
    hist(y_agg,
         probability = TRUE,       # Scale histogram to density
         breaks = hist_plot_breaks,              # Adjust bin count as needed
         col = "lightgray",
         border = "white",
         main = "Histogram + Carrier Density",
         xlab = "Transformed Expression",
         ylab = "Density")
    lines(est_midpoints, carrier_est, col = "dodgerblue", lwd = 2)
  }
  
  print("computing group 1 bin matrices")
  # counts in each bin for both groups
  control_donors =unique(sObj_meta[sObj_meta$disease == "normal", "donor_id"]) #group 1
  y1 = lapply(control_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  print("computing group 2 bin matrices")
  disease_donors = unique(sObj_meta[sObj_meta$disease != "normal", "donor_id"]) #group 2
  y2 = lapply(disease_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  if (unique_pts < 100) {
    Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K)
    Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) 
  }
  else{
    Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K)
    Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) 
  }
  binwidth = diff(est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  print("running SEF regression model")
  tk = est_midpoints
  X = rep(1, length(est_midpoints)) # Construct design matrix
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  # For now we consider sum first to avoid numerical issue
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  #setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) 
  scale_factor2 = cellSum2/sum(carrier_est)
  carrier_scale1 = carrier_est*scale_factor1 
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  df1 <- df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 <- df2[df2$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  # Fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  # Some numerical issue
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  # Some numerical issue
  convergence_issue <- (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for one fewer 
    varb <- paste0("X_", 1:p)
    X <- rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X <- cbind(X, tk^dim)
    }
    formula <- as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + "))) # rebuild formula
    df1_reduced <- df1[, -ncol(df1), drop = FALSE] # drop last columns
    df2_reduced <- df2[, -ncol(df2), drop = FALSE]
    
    sef_gp1 <- tryCatch(glm(formula, family = poisson(link = "log"), data = df1_reduced),
                        error = function(e) {
                          message("Retry for sef_gp1 failed: ", conditionMessage(e))
                          return(NULL)
                        })
    sef_gp2 <- tryCatch(glm(formula, family = poisson(link = "log"), data = df2_reduced),
                        error = function(e) {
                          message("Retry for sef_gp2 failed: ", conditionMessage(e))
                          return(NULL)
                        })
  }
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  
  
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  df_plot <- data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  
  if (plot_flag == T) {
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  return(list(y1 = y1, y2 = y2, sef_df1 = sef_df1, sef_df2 = sef_df2, df1 = df1, df2 = df2, 
              carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid, hist_plot_breaks = hist_plot_breaks,
              X =X, K = K, p = p, n_total_cells = n_total_cells, binwidth = binwidth,
              Smat1 = Smat1, Smat2 = Smat2, Smat = Smat, cellSum1 = cellSum1, cellSum2 = cellSum2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2))
}
sef_inference <- function(X = NULL, K = NULL, sef_df1 = NULL, sef_df2 = NULL, cellSum1 = NULL, cellSum2 = NULL, est_midpoints = NULL, binwidth = NULL, beta_est1 = NULL, beta_est2 = NULL, 
                          Smat1 = NULL, Smat2 = NULL, p = NULL, sef_regression_obj = NULL)
{
  if (!is.null(sef_regression_obj)) {   # If sef_regression_obj is provided, extract components from it
    X <- sef_regression_obj$X
    K <- sef_regression_obj$K
    n_total_cells <- sef_regression_obj$n_total_cells
    sef_df1 <- sef_regression_obj$sef_df1
    sef_df2 <- sef_regression_obj$sef_df2
    cellSum1 <- sef_regression_obj$cellSum1
    cellSum2 <- sef_regression_obj$cellSum2
    est_midpoints <- sef_regression_obj$est_midpoints
    binwidth <- sef_regression_obj$binwidth
    beta_est1 <- sef_regression_obj$beta_est1
    beta_est2 <- sef_regression_obj$beta_est2
    Smat1 <- sef_regression_obj$Smat1
    Smat2 <- sef_regression_obj$Smat2
    p <- sef_regression_obj$p
  }
  
  if (any(sapply(list(X, K, sef_df1, sef_df2, cellSum1, cellSum2, n_total_cells,
                      est_midpoints, binwidth, beta_est1, beta_est2, 
                      Smat1, Smat2, p), is.null))) {
    stop("Some required parameters are missing. Either provide them individually or via sef_regression_obj.")
  }
  beta_diff = beta_est2 - beta_est1
  S1sum <- colSums(Smat1)
  S2sum <- colSums(Smat2)
  
  
  #1. marginal distribution (no condition on intercept)
  G_1 <- t(X) %*% (sef_df1 * X) * cellSum1 / sum(sef_df1)
  G_2 <- t(X) %*% (sef_df2 * X) * cellSum2 / sum(sef_df2)
  
  pairwise_diff <- outer(est_midpoints, est_midpoints, "-")
  n_total_cells <- cellSum1 + cellSum2
  M <- dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  
  Z11 <- t(X) %*% (diag(K) - cellSum1 * as.vector(exp(X %*% beta_est1)) * M / sum(sef_df1))
  Z12 <- t(X) %*% (-cellSum1 * as.vector(exp(X %*% beta_est1)) * M / sum(sef_df1))
  Z21 <- t(X) %*% (-cellSum2 * as.vector(exp(X %*% beta_est2)) * M / sum(sef_df2))
  Z22 <- t(X) %*% (diag(K) - cellSum2 * as.vector(exp(X %*% beta_est2)) * M / sum(sef_df2))
  
  L1 <- solve(G_1, Z11) - solve(G_2, Z21)
  L2 <- solve(G_1, Z12) - solve(G_2, Z22)
  
  m_vec1 <- rowSums(Smat1)
  m_vec2 <- rowSums(Smat2)
  Smat1_ave <- Smat1 / m_vec1
  Smat2_ave <- Smat2 / m_vec2
  
  Smat1_centered <- t(Smat1_ave) - colSums(Smat1) / sum(m_vec1)
  Smat2_centered <- t(Smat2_ave) - colSums(Smat2) / sum(m_vec2)
  
  tmp1 <- diag(S1sum) - 2 * t(Smat1) %*% Smat1_ave + t(Smat1_ave) %*% Smat1_ave + 
    Smat1_centered %*% (t(Smat1_centered) * m_vec1^2)
  tmp2 <- diag(S2sum) - 2 * t(Smat2) %*% Smat2_ave + t(Smat2_ave) %*% Smat2_ave + 
    Smat2_centered %*% (t(Smat2_centered) * m_vec2^2)
  
  Cov_bar <- L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  
  chi_stat1 <- as.numeric(beta_diff[-1] %*% solve(Cov_bar[-1, -1], beta_diff[-1]))
  pval_1 <- 1 - pchisq(chi_stat1, df = p)
  print(paste0("marginal distribution pval: ", pval_1))
  
  return(list(beta_diff = beta_diff, pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_bar = Cov_bar))
}

stablerunSeuratCounts_Concise = function(exprMat, sObj_meta, donor_list, p = 2, plot_flag = T, h = NULL){
  regression_object = sef_regression(exprMat = exprMat, sObj_meta = sObj_meta, donor_list = donor_list, p = p, plot_flag = plot_flag, h = h) # , transformation = transformation)
  test_object = sef_inference(sef_regression_obj = regression_object)
  to_remove <- c("X", "Smat1", "Smat2", "Smat", "sef_gp1", "sef_gp2") # too large and unnecessary at the end
  regression_statistics <- regression_object[!(names(regression_object) %in% to_remove)]
  return(c(regression_statistics, test_object))
}

#  ---- PERMUTATION FUNCTIONS (changed to match sef_regression)  ----
sef_permuted_regression = function(exprMat, sObj_meta, donor_list, covariate_df, p = 2, plot_flag = F, h = NULL, seed = 100){
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(exprMat = exprMat, metadata = sObj_meta, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  #apply poisson VST based on sqrt + 1/2 transformation
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation
  #Gaussian KDE evaluated at integers
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ #uneven bins
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
      print(K)
    } else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    dev.off()
    carrier_est = carrier$est.fn
    hist_plot_breaks = carrier$hist$breaks
  } else{
    K = 100
    print(paste0("greater than or equal to 100 unique points: using equidistant bins as points for Gaussian KDE; carrier density estimate (", unique_pts," unique points), setting K bin count to: ", K))
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  if (plot_flag == T) {
    hist(y_agg,
         probability = TRUE,       # Scale histogram to density
         breaks = hist_plot_breaks,              # Adjust bin count as needed
         col = "lightgray",
         border = "white",
         main = "Histogram + Carrier Density",
         xlab = "Transformed Expression",
         ylab = "Density")
    lines(est_midpoints, carrier_est, col = "dodgerblue", lwd = 2)
  }
  
  print(paste0("permuting samples with SLE and healthy controls with seed ", seed, "..."))
  set.seed(seed) # for replicability
  covariate_df = subset(covariate_df, donor_id %in% donor_list)
  permuted = covariate_df %>% mutate(disease_status = sample(disease_status))
  
  print("done")
  
  print("computing group 1 bin matrices")
  # counts in each bin for both groups
  control_donors = subset(permuted, disease_status == 0)$donor_id #group 1
  y1 = lapply(control_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  print("computing group 2 bin matrices")
  disease_donors = subset(permuted, disease_status == 1)$donor_id #group 2
  y2 = lapply(disease_donors, function(d){ #list of arrays
    rows = which(as.character(sObj_meta$donor_id) == d)
    y_agg[rows]
  })
  
  if (unique_pts < 100) {
    Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K)
    Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) 
  }
  else{
    Smat1 = get_sc_bin_counts(Y = y1, bin_edges = est_grid, K = K)
    Smat2 = get_sc_bin_counts(Y = y2, bin_edges = est_grid, K = K) 
  }
  binwidth = diff(est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  print("running SEF regression model")
  tk = est_midpoints
  X = rep(1, length(est_midpoints)) # Construct design matrix
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  # For now we consider sum first to avoid numerical issue
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  #setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) 
  scale_factor2 = cellSum2/sum(carrier_est)
  carrier_scale1 = carrier_est*scale_factor1 
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  df1 <- df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 <- df2[df2$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  # Fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  # Some numerical issue
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  # Some numerical issue
  convergence_issue <- (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for one fewer 
    varb <- paste0("X_", 1:p)
    X <- rep(1, length(est_midpoints))
    for (dim in 1:p) {
      X <- cbind(X, tk^dim)
    }
    formula <- as.formula(paste0("sum_cts ~ offset(log(carrier_scale)) + ", paste(varb, collapse = " + "))) # rebuild formula
    df1_reduced <- df1[, -ncol(df1), drop = FALSE] # drop last columns
    df2_reduced <- df2[, -ncol(df2), drop = FALSE]
    
    sef_gp1 <- tryCatch(glm(formula, family = poisson(link = "log"), data = df1_reduced),
                        error = function(e) {
                          message("Retry for sef_gp1 failed: ", conditionMessage(e))
                          return(NULL)
                        })
    sef_gp2 <- tryCatch(glm(formula, family = poisson(link = "log"), data = df2_reduced),
                        error = function(e) {
                          message("Retry for sef_gp2 failed: ", conditionMessage(e))
                          return(NULL)
                        })
  }
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  
  
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  df_plot <- data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  
  if (plot_flag == T) {
    plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
      scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
      theme_minimal()
    print(plt)
  }
  
  return(list(control_donors = control_donors, disease_donors = disease_donors, y1 = y1, y2 = y2, 
              sef_df1 = sef_df1, sef_df2 = sef_df2, df1 = df1, df2 = df2, 
              carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid,
              X =X, K = K, p = p, n_total_cells = n_total_cells, binwidth = binwidth,
              Smat1 = Smat1, Smat2 = Smat2, Smat = Smat, cellSum1 = cellSum1, cellSum2 = cellSum2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2))
}

stablePermutationTest = function(exprMat, sObj_meta, donor_list, covariate_df, p = 2, plot_flag = T, h = NULL, seed = 100){
  regression_object = sef_permuted_regression(exprMat = exprMat, sObj_meta = sObj_meta, donor_list = donor_list, covariate_df = covariate_df, p = p, plot_flag = plot_flag, h = h, seed = seed) 
  test_object = sef_inference(sef_regression_obj = regression_object)
  to_remove <- c("y1","y2","sef_gp1","sef_gp2", "df1","df2", "X", "Smat1", "Smat2", "Smat", "sef_gp1", "sef_gp2") # too large and unnecessary at the end
  regression_statistics <- regression_object[!(names(regression_object) %in% to_remove)]
  
  return(c(regression_statistics, test_object))
}

####### FIGURE FUNCTIONS ####### 
gg_qqplot <- function(ps, ci = 0.95) {
  #plotting qq plot with 95% confidence interval band 
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  max_lim <- max(df$expected, df$observed, df$clower, df$cupper, na.rm = TRUE) # df$clower, df$cupper,
  
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1, fill = "steelblue"
    ) +
    geom_point(aes(expected, observed), shape = 21, size = 3, color = "black", fill = "steelblue" ) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
    theme_bw(base_size = 22) +
    xlab(log10Pe) + ylab(log10Po) 
} 
plot_hli_carrier_with_hist = function(exprMat, gene, sObj_meta, donors_to_use_list, h = NULL){
  exprMat = exprMat[gene, , drop = F]
  donor_list = donors_to_use_list[[gene]]
  
  print("subsetting seurat object")
  sObj_meta = sObj_meta[sObj_meta$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(metadata = sObj_meta, exprMat = exprMat, donor_vector = donor_list))
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print("applying poisson VST to data")
  #apply poisson VST based on sqrt + 1/2 transformation
  y_agg = sqrt(y_agg + 1/2) #bartlett poisson transformation
  #Gaussian KDE evaluated at integers
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  print(paste0("unique points: ", unique_pts))
  
  if(unique_pts < 100){ #uneven bins
    print(paste0("under 100 unique points: using discrete KDE with a Gamma kernel as carrier density estimate"))
    if (is.null(h)){
      if (unique_pts < 15){
        K = 30
        h = 0.005
      } else if (unique_pts < 20){
        K = 40
        h = 0.015
      } else if (unique_pts < 35){
        K = 70
        h = 0.02
      } else if(unique_pts < 40){
        K = 80
        h = 0.05
      } else if (unique_pts < 50){
        K = 100
        h = 0.1
      } else{
        K = 100
        h = 0.2
      }
      print(paste0("bandwidth not specified. h = ", h))
    }
    else{
      K = unique_pts*2
    }
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    pdf(NULL)
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    dev.off()
    carrier_est = carrier$est.fn
    hist_plot_breaks = carrier$hist$breaks
  } else{
    K = 100
    print(paste0("greater than or equal to 100 unique points: using equidistant bins as points for Gaussian KDE; carrier density estimate (", unique_pts," unique points), setting K bin count to: ", K))
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    hist_plot_breaks = est_grid
  }
  
  p <- ggplot() +
    geom_histogram(aes(x = y_agg, y = ..density..),
                   breaks = hist_plot_breaks,
                   fill = "skyblue1", alpha = 0.75,color = "white") +
    geom_line(aes(x = est_midpoints, y = carrier_est),
              color = "dodgerblue", size = 1) +
    labs(
      title = paste0("Carrier Density for ", gene),
      x = "Expression",
      y = "Density"
    ) + theme_minimal()
  
  return(p)
}
group_model_comparison_plot <- function(sef_object, gene) {
  if (is.null(sef_object)) {
    return("empty object")
  }else{
    sef_1 = sef_object[[gene]]$sef_df1
    sef_2 = sef_object[[gene]]$sef_df2
    est_grid = sef_object[[gene]]$est_grid
    est_midpoints = sef_object[[gene]]$est_midpoints
    df_plot <- data.frame(
      x = rep(est_midpoints, 2),
      y = c(sef_1, sef_2),
      group = rep(c("Control", "SLE"), each = length(est_midpoints))
    )
    
    p <- ggplot(df_plot, aes(x = x, y = y, color = group, fill = group)) +
      geom_line(size = 1.2) +
      geom_area(alpha = 0.2, position = "identity") +
      scale_color_manual(values = c("Control" = "dodgerblue", "SLE" = "firebrick")) +
      scale_fill_manual(values = c("Control" = "dodgerblue", "SLE" = "firebrick")) +
      labs(
        title = paste("Group-Level Model Density:", gene),
        x = "Expression",
        y = "Density",
        color = "Group",
        fill = "Group"
      ) +
      theme_minimal(base_size = 14)
    
    return(p)
  }
}
get_individual_sef_counts_hli_model = function(donor_id, gene, exprMat, sObj_meta, donor_list, p = 2, h = NULL, plot_flag = F){
  #run global model first, receive the carrier density as well
  res = stablerunSeuratCounts_Concise(exprMat = exprMat[gene, , drop = F], sObj_meta = sObj_meta, donor_list = donor_list[[gene]], p = p, h = h, plot_flag = plot_flag) #group-specific model execution
  print(paste0("now running individual-specific model for donor ", donor_id))
  
  #get midpoints, carrier density estimate, gridpoints, status
  ind_carrier_est = res$carrier_est
  ind_est_midpoints = res$est_midpoints
  ind_est_grid = res$est_grid
  ind_status = as.character(sObj_meta[sObj_meta$donor_id == donor_id,]$disease[1])
  y_agg = c(unlist(res$y1), unlist(res$y2))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  n_total_cells = res$n_total_cells
  K = res$K
  h = res$h
  p = res$p
  hist_plot_breaks = res$hist_plot_breaks
  
  if (ind_status == "normal") {
    grp_model = res$sef_df1 #Normal
  }else{
    grp_model = res$sef_df2 #SLE
  }
  
  ind_sObj_meta = sObj_meta[sObj_meta$donor_id == donor_id, , drop = FALSE]
  ind_cell_barcodes <- rownames(ind_sObj_meta) #to subset individual-specific vector of expression counts
  y_ind <- as.vector(exprMat[gene, ind_cell_barcodes, drop = F])
  # y_ind = as.vector(ind_sObj@assays$SCT$counts[gene,]) #individual-specific vector of expression counts
  y_ind = sqrt(y_ind + 1/2)
  n_indiv_cells = length(y_ind)
  
  tk_ind = ind_est_midpoints
  X_ind = rep(1, length(ind_est_midpoints))
  for (dim in 1:p){
    X_ind = cbind(X_ind, tk_ind^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X_ind) = c("Intercept", varb)
  
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  
  Smat_ind = get_sc_bin_counts(Y = list(y_ind), bin_edges = ind_est_grid, K = K)
  ind_SSum = colSums(Smat_ind)
  binwidth = diff(ind_est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  
  # setup for modeling individual-wise
  cellSum_ind = length(unlist(y_ind))
  scale_factor_ind = cellSum_ind/sum(ind_carrier_est)
  ind_carrier_scale = ind_carrier_est*scale_factor_ind
  df1 = as.data.frame(cbind( ind_SSum, ind_carrier_scale, X_ind))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  df1 <- df1[df1$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+')))
  sef_gp_ind = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  # Some numerical issue
  beta_est_ind = as.vector(sef_gp_ind$coefficients)
  sef_df1 = as.vector( ind_carrier_est * exp(X_ind %*% beta_est_ind) )
  
  #plotting machinery
  p_df <- data.frame(
    x = ind_est_midpoints,
    y = sef_df1,
    group = "Donor-Specific"
  )
  p_df2 <- data.frame(
    x = ind_est_midpoints,
    y = grp_model,
    group = "Group Model"
  )
  p_df_combined <- rbind(p_df, p_df2)
  
  df_hist <- data.frame(y_ind = y_ind)
  
  plt <- ggplot() +
    geom_histogram(data = df_hist, aes(x = y_ind, y = ..density..), #histogram
                   breaks = hist_plot_breaks,
                   fill = "skyblue1", alpha = 0.75, color = "white") +
    geom_line(data = p_df_combined, aes(x = x, y = y, color = group), size = 1.2) + #group versus donor specific density estimates
    scale_color_manual(values = c("Donor-Specific" = "dodgerblue",
                                  "Group Model" = "firebrick")) +
    theme_minimal() +
    labs(
      x = "Expression",
      y = "Density",
      title = paste0("Donor ", donor_id, " - Specific Density: ", gene),
      color = "Model"
    )
  return(plt)
  print("done")
  
}
run_go_enrichment_plot <- function(gene_symbols, orgDb = org.Hs.eg.db, ontology = "BP", minGSSize = 10, maxGSSize = 200) {
  require(clusterProfiler)
  require(enrichplot)
  require(org.Hs.eg.db)  # load if default used
  require(ggplot2)
  gene_df <- clusterProfiler::bitr(gene_symbols,
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = orgDb)
  if (nrow(gene_df) == 0) {
    stop("No gene symbols were converted to Entrez IDs. Check your gene list.")
  }
  # GO enrichment
  ego <- enrichGO(gene = gene_df$ENTREZID,
                  OrgDb = orgDb,
                  minGSSize = minGSSize,
                  maxGSSize = maxGSSize,
                  keyType = "ENTREZID",
                  ont = ontology,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05,
                  readable = TRUE)
  ego = clusterProfiler::simplify(ego, cutoff = 0.6, by = "p.adjust")
}



