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
library(foreach)
library(doParallel) 
library(Hmisc)
library(Ake)
library(truncnorm)

########## 1. DATA GENERATION FUNCTIONS ##########
simu_hier_vary_cellcts <- function(seed = 100, lower = 300, upper = 1000, alpha1, alpha2, tau1, tau2, sigma1, sigma2, n1, n2) {
  Y1_list <- vector("list", n1) 
  Y2_list <- vector("list", n2)
  set.seed(seed = seed)
  for (i in seq_len(n1)) {  # GROUP 1
    n_cells_i <- sample(lower:upper, 1)
    alpha_1i <- rnorm(1, alpha1, sqrt(tau1))
    Y1_list[[i]] <- rnorm(n_cells_i, mean = alpha_1i, sd = sqrt(sigma1))
  }
  
  for (i in seq_len(n2)) {  # GROUP 2
    n_cells_i <- sample(lower:upper, 1)
    alpha_2i <- rnorm(1, alpha2, sqrt(tau2))
    Y2_list[[i]] <- rnorm(n_cells_i, mean = alpha_2i, sd = sqrt(sigma2))
  }
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}
simu_poigamma_vary_cellcts = function(seed = 100, lower = 300, upper = 1000, alpha1, alpha2, beta1, beta2, n1, n2){
  # Simulate Poisson-Gamma mixture model data with alpha and beta are the shape and size parameters
  Y1_list <- vector("list", n1) 
  Y2_list <- vector("list", n2)
  set.seed(seed = seed)
  for (i in 1:n1) { #GROUP 1
    n_cells_i <- sample(lower:upper, 1)
    lambda1i = rgamma(1, shape = alpha1, rate = beta1)
    Y1_list[[i]] <- rpois(n_cells_i, lambda1i)
  }
  for (i in 1:n2) { #GROUP 2 
    n_cells_i <- sample(lower:upper, 1)
    lambda2i = rgamma(1, shape = alpha2, rate = beta2)
    Y2_list[[i]] <- rpois(n_cells_i, lambda2i)
  }  
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}
simu_zinb_vary_cellcts <- function(n1, n2, mu, mu2, theta, sigma_sq = 0.5, pi = 0.5, pi2 = 0.5,
                                   b_pi = 1, b_theta = 1, de_type = 'dispersion', seed = 100, lower = 300, upper = 1000) {
  set.seed(seed)
  
  # Output lists
  Y1_list <- vector("list", n1)
  Y2_list <- vector("list", n2)
  
  # Individual-specific means
  mu_i_group1 <- rtruncnorm(n1, a = 0, mean = mu, sd = sqrt(sigma_sq))
  
  if (de_type == 'dispersion') {
    mu2 <- mu
    theta2 <- theta * b_theta
    pi2 <- pi
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
    
  } else if (de_type == 'mixture') {
    if (b_pi < 0 || b_pi > 1) stop("b_pi must be between 0 and 1 for mixture")
    
    theta2 <- theta * b_theta
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
    
    for (i in seq_len(n1)) {
      n_cells_i <- sample(lower:upper, 1)
      Y1_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
        if (runif(1) < pi) {
          0
        } else {
          if (runif(1) < b_pi) {
            rnbinom(1, size = theta, mu = mu_i_group1[i])
          } else {
            rnbinom(1, size = theta2, mu = mu_i_group2[i])
          }
        }
      })
    }
    
    for (i in seq_len(n2)) {
      n_cells_i <- sample(lower:upper, 1)
      Y2_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
        if (runif(1) < pi2) {
          0
        } else {
          if (runif(1) < (1 - b_pi)) {
            rnbinom(1, size = theta, mu = mu_i_group1[i])
          } else {
            rnbinom(1, size = theta2, mu = mu_i_group2[i])
          }
        }
      })
    }
    
    return(list(Y1 = Y1_list, Y2 = Y2_list))
    
  } else if (de_type == "mean") {
    theta2 <- theta
    pi2 <- pi
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
  } else {
    stop('Invalid de_type specified. Choose either "mean", "dispersion", or "mixture".')
  }
  
  # Default ZINB simulation for mean/dispersion
  for (i in seq_len(n1)) {
    n_cells_i <- sample(lower:upper, 1)
    Y1_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
      if (runif(1) < pi) 0 else rnbinom(1, size = theta, mu = mu_i_group1[i])
    })
  }
  
  for (i in seq_len(n2)) {
    n_cells_i <- sample(lower:upper, 1)
    Y2_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
      if (runif(1) < pi2) 0 else rnbinom(1, size = theta2, mu = mu_i_group2[i])
    })
  }
  
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}

########## 2. SEF HELPER FUNCTIONS ##########
#  ---- DISCRETE GAUSSIAN KERNEL FUNCTION ----
continuous_discrete_kde <- function(data, eval_points, bw = 0.25) {
  dx <- diff(eval_points)[1]  # Step size between grid points; assumes equally spaced grid
  raw_density <- sapply(eval_points, function(x) {
    mean(dnorm((x - data) / bw)) / bw
  })
  normalized_density <- raw_density / sum(raw_density * dx)
  return(normalized_density)
}
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

#  ---- POWER ANALYSIS SPECIFIC  ----

gamma_params <- function(mu, V) {
  if(V <= mu) stop("Variance must be greater than the mean for overdispersion.")
  alpha <- mu^2/(V - mu)
  beta  <- mu/(V - mu)
  return(list(alpha = alpha, beta = beta))
}

compute_moments = function(y, p) {
  #computes pth MOMENTS and stores as a matrix
  return(sapply(1:p, function(k) y^k))
} 

mom_cov_diff_cell_cts = function(Y, p) {
  # Compute covariance matrix for Moments (for MoM estimator) with different cell counts
  n = length(Y)  # number of individuals (length of the list)
  total_samples <- 0
  overall_sum <- rep(0, p)  # Sum for overall mean
  within_var_sum <- matrix(0, nrow = p, ncol = p)  # Within-individual variance
  between_var_sum <- matrix(0, nrow = p, ncol = p)  # Between-individual variance 
  
  mean_list = matrix(0, nrow = n, ncol = p)  # matrix of p moments for each individual
  m_vector <- numeric(n)  # Store the number of samples for each individual
  
  # Within-individual variance
  for (i in 1:n) {  # iterate for each individual
    samples_i = Y[[i]]  # Access the ith individual's data
    m_i = length(samples_i) 
    m_vector[i] = m_i
    mean_samples_i = mean(samples_i)
    
    # Compute the moments for each sample
    moments_i = t(sapply(samples_i, compute_moments, p = p))  # ncell_i x p moment
    
    mean_t_i = colMeans(moments_i)  # Compute the mean of p moments for individual
    mean_list[i, ] = mean_t_i
    
    # Update overall sum for computing the overall mean
    overall_sum = overall_sum + m_i * mean_t_i
    total_samples = total_samples + m_i
    
    # Compute within-individual variance
    diffs = t(moments_i) - mean_t_i  # (T_i - Tbar_i)
    within_var_sum = within_var_sum + diffs %*% t(diffs)  # sum (T_i - Tbar_i)^2
  }
  
  # Compute the overall mean of the moments
  overall_mean <- overall_sum / total_samples
  
  # Between-individual variance
  for (i in 1:n) {
    m_i = m_vector[i]
    mean_t_i = mean_list[i, ]
    diff <- mean_t_i - overall_mean
    between_var_sum <- between_var_sum + m_i^2 * (diff %*% t(diff))  # m_i^2, not m_i
  }
  
  # Covariance calculation
  Cov = within_var_sum / total_samples^2 + between_var_sum / total_samples^2
  return(Cov)
}

compute_empirical_power = function(pval_matrix, alpha = 0.05){
  #compute empirical power using the #of rejections we make divided by the total number of simulations
  pmat = as.matrix(pval_matrix)
  return(colMeans(pmat < alpha))
}

extract_and_combine <- function(results, attr_name) {
  #helper function: function to extract and combine each attribute into a data frame
  attr_list <- lapply(results, `[[`, attr_name) # extract specific attribute from each simulation result
  
  # convert list of vectors into a data frame with rows for each simulation
  attr_df <- do.call(rbind, lapply(attr_list, function(x) as.data.frame(t(x))))
  
  return(attr_df)
}

extract_and_combine_scalar <- function(results, attr_name) {
  # extract specific scalar attribute from each simulation result
  attr_list <- lapply(results, function(res) res[[attr_name]])
  attr_vector <- unlist(attr_list, use.names = FALSE) # convert list into a data frame
  
  return(attr_vector)
}

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
  max_lim <- max(df$expected, df$observed, df$clower, df$cupper, na.rm = TRUE) 
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

########## 3. SEF FUNCTIONS ##########
#  ---- GAUSSIAN  ----
simulateNormalRealistic = function(repID = NULL, de_type = "mean", p = 2,  K = 75,  lower = 300, upper = 1000, n1 = 100, n2 = 100, idx, alpha1, alpha2, tau1, tau2, sigma1, sigma2){
  #save values: effect size, binwidth, betas, covariances
  if (de_type == "mean") {
    effect_size = alpha2 - alpha1
  } else if(de_type == "variance"){
    effect_size = sigma2 - sigma1
  } else{
    warning("Shift type not found. Please use valid shift type.")
    return(NULL)
  }
  
  seed = idx + 62
  simData = simu_hier_vary_cellcts(seed = seed, lower = lower, upper = upper, alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
  n_total_cells = length(y_agg)
  print(paste0("total cells: ", n_total_cells))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u )
  carrier_est = carrier$y
  est_grid = seq(l, u, length.out = K + 1) 
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2  
  
  hist_df <- data.frame(y_agg = y_agg)
  breaks <- seq(l, u, length.out = K + 2)
  line_df <- data.frame(x = est_midpoints, y = carrier_est)
  carrier_plot <- ggplot(hist_df, aes(x = y_agg)) +
    geom_histogram(aes(y = ..density..),
                   # breaks = est_grid,
                   fill = "skyblue1", alpha = 0.75,
                   color = "white") +
    geom_line(data = line_df, aes(x = x, y = y),
              color = "dodgerblue", linewidth = 1) +
    labs(
      title = "",
      x = "Expression",
      y = "Density"
    ) +
    theme_minimal()
  print(carrier_plot)
  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  Smat1 = get_sc_bin_counts(Y = Y1, bin_edges = est_grid, K = K)
  Smat2 = get_sc_bin_counts(Y = Y2, bin_edges = est_grid, K = K)
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  
  print("running SEF regression model")
  #Construct design matrix
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  # For now we consider sum first to avoid numerical issue
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
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
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  
  # Plot differences
  df_plot <- data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
    scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
    theme_minimal()
  print(plt)
  
  #1. marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1) # cell sum times the bindwidth or divide the sum of density
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M <- dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  Z11 = t(X) %*% ( diag(K) -  cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # Z_{t1, j} for j in T1
  Z12 = t(X) %*% ( -cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # For j in T2
  Z21 = t(X) %*% ( -cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2)  ) # For j in T1
  Z22 = t(X) %*% ( diag(K) -  cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2) ) # For j in T2
  L1 = solve(G_1, Z11) - solve(G_2, Z21) # For S1 data 
  L2 = solve(G_1, Z12) - solve(G_2, Z22)# For S2 data  
  m_vec1 = rowSums(Smat1) 
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1/m_vec1
  Smat2_ave = Smat2/m_vec2
  Smat1_centered = t(Smat1_ave) - colSums(Smat1)/sum(m_vec1)  
  Smat2_centered = t(Smat2_ave) - colSums(Smat2)/sum(m_vec2)
  tmp1 = diag(S1sum) - 2* t(Smat1)%*%Smat1_ave + t(Smat1_ave)%*%Smat1_ave + 
    Smat1_centered%*%( t(Smat1_centered) * m_vec1^2 )
  tmp2 = diag(S2sum) - 2* t(Smat2)%*%Smat2_ave + t(Smat2_ave)%*%Smat2_ave + 
    Smat2_centered%*%( t(Smat2_centered) * m_vec2^2 )
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  
  chi_stat1 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1], beta_diff[-1]) ) #marginal distribution test
  pval_1 = 1 - pchisq( chi_stat1, df = p)
  print(pval_1)
  
  # Model without intercept effects
  tc_1 = NULL
  tc_2 = NULL
  Y1_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  Y2_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  for (dim in 1:p){
    tc_1 = cbind(tc_1, tk^dim)
    tc_2 = cbind(tc_2, tk^dim)
    Y1_suffi = cbind(Y1_suffi, y1^dim) 
    Y2_suffi = cbind(Y2_suffi, y2^dim)
  }  
  T_c1 = colMeans(Y1_suffi)
  T_c2 = colMeans(Y2_suffi)
  X_t1 = t( t(tc_1) - T_c1 ) #centered design matrix based on mean of group 1
  X_t2 = t( t(tc_2) - T_c2 ) #centered design matrix based on mean of group 2
  X_t1f = cbind(1, X_t1)
  X_t2f = cbind(1, X_t2) 
  colnames(X_t1f) = c("Intercept", varb)  
  colnames(X_t2f) = c("Intercept", varb) 
  df1_c = as.data.frame(cbind( S1sum, carrier_scale1, X_t1f))
  df2_c = as.data.frame(cbind( S2sum, carrier_scale2, X_t2f))
  colnames(df1_c)[1] = "sum_cts"
  colnames(df1_c)[2] = "carrier_scale"
  colnames(df2_c)[1] = "sum_cts" 
  colnames(df2_c)[2] = "carrier_scale"
  df1_c <- df1_c[df1_c$carrier_scale > 0, , drop = FALSE]
  df2_c <- df2_c[df2_c$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  # Fit glm model for each group
  sef_gp1_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df1_c))  # Some numerical issue
  sef_gp2_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df2_c))  # Some numerical issue
  beta_est1_c = as.vector(sef_gp1_ct$coefficients)
  beta_est2_c = as.vector(sef_gp2_ct$coefficients)
  beta_diff_c = beta_est1_c - beta_est2_c
  beta_diff_c
  sef_df1_c = as.vector( carrier_est * exp(X_t1f %*% beta_est1_c) )
  sef_df2_c = as.vector( carrier_est * exp(X_t2f %*% beta_est2_c) )
  
  # 3. Inference without intercept effects
  G_t1 = t(X_t1) %*%  ( sef_df1_c * X_t1)*cellSum1/sum(sef_df1_c)
  G_t1
  G_t2 = t(X_t2) %*% ( sef_df2_c * X_t2)*cellSum2/sum(sef_df2_c) 
  G_t2  
  Z11c = t(X_t1)%*% (diag(K) - cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  Z12c = t(X_t1)%*%( -cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  Z21c = t(X_t2)%*%( -cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  Z22c = t(X_t2)%*%( diag(K) - cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  L1c = solve(G_t1, Z11c) - solve(G_t2, Z21c)
  L2c = solve(G_t1, Z12c) - solve(G_t2, Z22c)
  Cov_bar_sub = L1c%*%tmp1%*%t(L1c) + L2c%*%tmp2%*%t(L2c)
  chi_stat3 = beta_diff_c[-1]%*%solve(Cov_bar_sub, beta_diff_c[-1])
  pval_3 = 1 - pchisq( chi_stat3, df = p) 
  
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, n1 = n1, n2 = n2, binwidth = binwidth, X = X, y_agg = y_agg, y1 = y1, y2 = y2,
              carrier_est = carrier_est, est_midpoints = est_midpoints, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2,
              sef_df1 = sef_df1, sef_df2 = sef_df2, sef_df1_c = sef_df1_c, sef_df2_c = sef_df2_c,
              pval_1 = pval_1, pval_3 = pval_3[1,1], chi_stat1 = chi_stat1,chi_stat3 = chi_stat3, Cov_1 = Cov_bar, Cov_3 = Cov_bar_sub,
              beta_diff = beta_diff, beta_diff_c = beta_diff_c, carrier_plt = carrier_plot, comparison_plt = plt))
}
#  ---- POISSON GAMMA  ----
simulatePGRealistic = function(repID = NULL, de_type = "mean", p = 2,  K = NULL,  lower = 300, upper = 1000, n1 = 100, n2 = 100, idx, alpha1, alpha2, beta1, beta2, regularize = F){
  if (de_type == "mean") {
    effect_size = alpha2/beta2 - alpha1/beta1
  } else if(de_type == "variance"){
    var1 = alpha1/beta1 + alpha1/(beta1^2)
    var2 = alpha2/beta2 + alpha2/(beta2^2)
    effect_size = var2 - var1 #variance test
  } else{
    warning("Shift type not found. Please use valid shift type.")
    return(NULL)
  }
  seed = idx + 62
  set.seed(seed)
  simData = simu_poigamma_vary_cellcts(seed =seed, lower = lower, upper = upper, alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2,n1 = n1, n2 = n2)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
  n_total_cells = length(y_agg)
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  print(paste0("total cells: ", n_total_cells))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  
  if (is.null(K)) {
    est_grid = seq(l,u)
    est_midpoints = est_grid
    K = length(est_midpoints)
    print(paste0("Number of midpoints K not specified. Using total number of boundary points in pooled data as new K: ", K))
    
    # option 2: discrete Gaussian KDE
    carrier_est = continuous_discrete_kde(data = y_agg, eval_points = est_midpoints, bw = 0.05)
    
  } else{
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
  }
  
  hist_df <- data.frame(y_agg = y_agg)
  breaks <- seq(l, u, length.out = K + 2)
  line_df <- data.frame(x = est_midpoints, y = carrier_est)
  carrier_plot <- ggplot(hist_df, aes(x = y_agg)) +
    geom_histogram(aes(y = ..density..),
                   breaks = breaks,
                   fill = "skyblue1", alpha = 0.75,
                   color = "white") +
    geom_line(data = line_df, aes(x = x, y = y),
              color = "dodgerblue", linewidth = 1) +
    labs(
      title = "",
      x = "Expression",
      y = "Density"
    ) +
    theme_minimal()
  print(carrier_plot)
  
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  Smat1 = get_countdata_bin_counts(Y = Y1, grid_points = est_grid)
  Smat2 = get_countdata_bin_counts(Y = Y2, grid_points = est_grid)
  # Smat1 = get_sc_bin_counts(Y = Y1, bin_edges = est_grid, K = K)
  # Smat2 = get_sc_bin_counts(Y = Y2, bin_edges = est_grid, K = K)
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  
  print("running SEF regression model")
  #Construct design matrix
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
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
    print(paste0("new test with p = ", p))
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
  
  
  # Plot differences
  df_plot <- data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
    scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
    theme_minimal()
  print(plt)
  
  #1. marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1) # cell sum times the bindwidth or divide the sum of density
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M <- dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  Z11 = t(X) %*% ( diag(K) -  cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # Z_{t1, j} for j in T1
  Z12 = t(X) %*% ( -cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # For j in T2
  Z21 = t(X) %*% ( -cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2)  ) # For j in T1
  Z22 = t(X) %*% ( diag(K) -  cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2) ) # For j in T2
  L1 = solve(G_1, Z11) - solve(G_2, Z21) # For S1 data 
  L2 = solve(G_1, Z12) - solve(G_2, Z22)# For S2 data  
  m_vec1 = rowSums(Smat1) 
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1/m_vec1
  Smat2_ave = Smat2/m_vec2
  Smat1_centered = t(Smat1_ave) - colSums(Smat1)/sum(m_vec1)  
  Smat2_centered = t(Smat2_ave) - colSums(Smat2)/sum(m_vec2)
  tmp1 = diag(S1sum) - 2* t(Smat1)%*%Smat1_ave + t(Smat1_ave)%*%Smat1_ave + 
    Smat1_centered%*%( t(Smat1_centered) * m_vec1^2 )
  tmp2 = diag(S2sum) - 2* t(Smat2)%*%Smat2_ave + t(Smat2_ave)%*%Smat2_ave + 
    Smat2_centered%*%( t(Smat2_centered) * m_vec2^2 )
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  if(regularize == T){
    #regularize Covariance by eigenvalue
    eig_vals = eigen(Cov_bar)$values
    max_eigenvalue = median(eig_vals)
    Cov_bar  = Cov_bar + diag(rep(1, dim(Cov_bar)[1]))*max_eigenvalue
  }
  
  chi_stat1 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1], beta_diff[-1]) ) #marginal distribution test
  pval_1 = 1 - pchisq( chi_stat1, df = p)
  print(pval_1)
  
  # Model without intercept effects
  tc_1 = NULL
  tc_2 = NULL
  Y1_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  Y2_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  for (dim in 1:p){
    tc_1 = cbind(tc_1, tk^dim)
    tc_2 = cbind(tc_2, tk^dim)
    Y1_suffi = cbind(Y1_suffi, y1^dim) 
    Y2_suffi = cbind(Y2_suffi, y2^dim)
  }  
  T_c1 = colMeans(Y1_suffi)
  T_c2 = colMeans(Y2_suffi)
  X_t1 = t( t(tc_1) - T_c1 ) #centered design matrix based on mean of group 1
  X_t2 = t( t(tc_2) - T_c2 ) #centered design matrix based on mean of group 2
  X_t1f = cbind(1, X_t1)
  X_t2f = cbind(1, X_t2) 
  colnames(X_t1f) = c("Intercept", varb)  
  colnames(X_t2f) = c("Intercept", varb) 
  df1_c = as.data.frame(cbind( S1sum, carrier_scale1, X_t1f))
  df2_c = as.data.frame(cbind( S2sum, carrier_scale2, X_t2f))
  colnames(df1_c)[1] = "sum_cts"
  colnames(df1_c)[2] = "carrier_scale"
  colnames(df2_c)[1] = "sum_cts" 
  colnames(df2_c)[2] = "carrier_scale"
  df1_c <- df1_c[df1_c$carrier_scale > 0, , drop = FALSE]
  df2_c <- df2_c[df2_c$carrier_scale > 0, , drop = FALSE]
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  # Fit glm model for each group
  sef_gp1_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df1_c))  # Some numerical issue
  sef_gp2_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df2_c))  # Some numerical issue
  beta_est1_c = as.vector(sef_gp1_ct$coefficients)
  beta_est2_c = as.vector(sef_gp2_ct$coefficients)
  beta_diff_c = beta_est1_c - beta_est2_c
  beta_diff_c
  sef_df1_c = as.vector( carrier_est * exp(X_t1f %*% beta_est1_c) )
  sef_df2_c = as.vector( carrier_est * exp(X_t2f %*% beta_est2_c) )
  
  # 3. Inference without intercept effects
  G_t1 = t(X_t1) %*%  ( sef_df1_c * X_t1)*cellSum1/sum(sef_df1_c)
  G_t1
  G_t2 = t(X_t2) %*% ( sef_df2_c * X_t2)*cellSum2/sum(sef_df2_c) 
  G_t2  
  Z11c = t(X_t1)%*% (diag(K) - cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  Z12c = t(X_t1)%*%( -cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  Z21c = t(X_t2)%*%( -cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  Z22c = t(X_t2)%*%( diag(K) - cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  L1c = solve(G_t1, Z11c) - solve(G_t2, Z21c)
  L2c = solve(G_t1, Z12c) - solve(G_t2, Z22c)
  Cov_bar_sub = L1c%*%tmp1%*%t(L1c) + L2c%*%tmp2%*%t(L2c)
  chi_stat3 = beta_diff_c[-1]%*%solve(Cov_bar_sub, beta_diff_c[-1])
  pval_3 = 1 - pchisq( chi_stat3, df = p) 
  
  return(list(de_type = de_type, effect_size = effect_size, y_agg = y_agg, est_midpoints = est_midpoints, carrier_est = carrier_est, binwidth = binwidth, K = K, y1 = y1, y2 = y2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, beta_diff = beta_diff, beta_est1_c = beta_est1_c, beta_est2_c = beta_est2_c, beta_diff_c = beta_diff_c,
              sef_df1 = sef_df1, sef_df2 = sef_df2, sef_df1_c = sef_df1_c, sef_df2_c = sef_df2_c, pval_1 = pval_1, pval_3 = pval_3, Cov_1 = Cov_bar, Cov_3 = Cov_bar_sub,
              carrier_plt = carrier_plot, comparison_plt = plt, breaks = breaks))
}
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL  ----
simulateZINB = function(n1, n2, de_type = "dispersion", idx, mu, mu2, sigma_sq = 0.5, theta, b_theta = 1, b_pi = 1, pi = 0.5, pi2 = 0.5, lower = 300, upper = 1000, K = NULL, p = 2, repID = NULL){
  #save values: effect size, binwidth, betas, covariances
  if (de_type == "dispersion") {
    theta2 = b_theta*theta
    effect_size = theta2 - theta
  } else if(de_type == "mixture"){
    effect_size = b_pi
  } else{
    effect_size = mu2 - mu
  }
  seed = idx + 62
  simData = simu_zinb_vary_cellcts(n1 = n1, n2 = n2, lower = lower, upper = upper, mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, pi = pi, pi2 = pi2, 
                                   b_pi = b_pi, b_theta = b_theta, de_type = de_type, seed = seed)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
  unique_pts = length(unique(y_agg))
  print(paste0("unique points: ", unique_pts))
  
  # Carrier density 
  n_total_cells = length(y_agg)
  print(paste0("total cells: ", n_total_cells))
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  
  if (is.null(K)) {
    est_grid = seq(l,u)
    est_midpoints = est_grid
    K = length(est_midpoints)
    print(paste0("Number of midpoints K not specified. Using total number of boundary points in pooled data as new K: ", K))
    
    
    # option 2: discrete Gaussian KDE
    carrier_est = continuous_discrete_kde(data = y_agg, eval_points = est_midpoints, bw = 0.05)
    
    
  } else{
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
  }
  
  hist_df <- data.frame(y_agg = y_agg)
  breaks <- seq(l, u, length.out = K + 2)
  line_df <- data.frame(x = est_midpoints, y = carrier_est)
  carrier_plot <- ggplot(hist_df, aes(x = y_agg)) +
    geom_histogram(aes(y = ..density..),
                   breaks = breaks,
                   fill = "skyblue1", alpha = 0.75,
                   color = "white") +
    geom_line(data = line_df, aes(x = x, y = y),
              color = "dodgerblue", linewidth = 1) +
    labs(
      title = "",
      x = "Expression",
      y = "Density"
    ) +
    theme_minimal()
  print(carrier_plot)
  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  Smat1 = get_countdata_bin_counts(Y = Y1, grid_points = est_grid)
  Smat2 = get_countdata_bin_counts(Y = Y2, grid_points = est_grid)
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  
  print("running SEF regression model")
  #Construct design matrix
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  # For now we consider sum first to avoid numerical issue
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  
  if (unique_pts < 75) {
    binwidth = median(diff(est_midpoints))
  }else{
    binwidth = diff(est_midpoints)[1] #is identical to (u - l)/K in large enough settings
  }
  
  #setup for modeling
  cellSum1 = length(unlist(y1))
  cellSum2 = length(unlist(y2))
  scale_factor1 = cellSum1/sum(carrier_est) 
  scale_factor2 = cellSum2/sum(carrier_est)
  carrier_scale1 = carrier_est*scale_factor1 
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  df1 <- df1[df1$carrier_scale > 0, , drop = FALSE]
  df2 <- df2[df2$carrier_scale > 0, , drop = FALSE]
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  
  # Fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  # Some numerical issue
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  # Some numerical issue
  convergence_issue <- (!is.null(sef_gp1) && !sef_gp1$converged) || (!is.null(sef_gp2) && !sef_gp2$converged)
  if (convergence_issue) {
    message("Convergence issue detected in at least one model. Retrying with reduced columns...")
    p = p - 1 # we test for one fewer 
    print(paste0("new test with p = ", p))
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
  
  
  # Plot differences
  df_plot <- data.frame(
    est_midpoints = rep(est_midpoints, 2),
    sef_value = c(sef_df1, sef_df2),
    group = factor(rep(c("Group 1", "Group 2"), each = length(est_midpoints)))
  )
  plt = ggplot(df_plot, aes(x = est_midpoints, y = sef_value, fill = group, color = group)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = 0, ymax = sef_value), alpha = 0.3, color = NA) +
    scale_color_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
    labs(x = "Expression", y = "Density", color = "Group", fill = "Group") +
    theme_minimal()
  print(plt)
  
  #1. marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1) # cell sum times the bindwidth or divide the sum of density
  G_1
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2)
  G_2
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M <- dnorm(pairwise_diff, sd = binwidth) / n_total_cells
  Z11 = t(X) %*% ( diag(K) -  cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # Z_{t1, j} for j in T1
  Z12 = t(X) %*% ( -cellSum1*as.vector( exp(X %*% beta_est1))*M/sum(sef_df1) ) # For j in T2
  Z21 = t(X) %*% ( -cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2)  ) # For j in T1
  Z22 = t(X) %*% ( diag(K) -  cellSum2*as.vector(exp(X %*% beta_est2))*M/sum(sef_df2) ) # For j in T2
  L1 = solve(G_1, Z11) - solve(G_2, Z21) # For S1 data
  L2 = solve(G_1, Z12) - solve(G_2, Z22)# For S2 data
  m_vec1 = rowSums(Smat1)
  m_vec2 = rowSums(Smat2)
  Smat1_ave = Smat1/m_vec1
  Smat2_ave = Smat2/m_vec2
  Smat1_centered = t(Smat1_ave) - colSums(Smat1)/sum(m_vec1)
  Smat2_centered = t(Smat2_ave) - colSums(Smat2)/sum(m_vec2)
  tmp1 = diag(S1sum) - 2* t(Smat1)%*%Smat1_ave + t(Smat1_ave)%*%Smat1_ave +
    Smat1_centered%*%( t(Smat1_centered) * m_vec1^2 )
  tmp2 = diag(S2sum) - 2* t(Smat2)%*%Smat2_ave + t(Smat2_ave)%*%Smat2_ave +
    Smat2_centered%*%( t(Smat2_centered) * m_vec2^2 )
  Cov_bar = L1 %*% tmp1 %*% t(L1) + L2 %*% tmp2 %*% t(L2)
  
  chi_stat1 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1], beta_diff[-1]) ) #marginal distribution test
  pval_1 = 1 - pchisq( chi_stat1, df = p)
  print(pval_1)
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, n1 = n1, n2 = n2, binwidth = binwidth, X = X, y_agg = y_agg, y1 = y1, y2 = y2,
              carrier_est = carrier_est, est_midpoints = est_midpoints, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2, sef_gp1 = sef_gp1, sef_gp2 = sef_gp2,
              sef_df1 = sef_df1, sef_df2 = sef_df2, 
              pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_1 = Cov_bar, 
              beta_diff = beta_diff,
              carrier_plt = carrier_plot, comparison_plt = plt))
}
########## 4. COMPETING METHOD FUNCTIONS ##########
#  ---- METHOD OF MOMENTS  ----
simulateMoMRealistic = function(distribution = c("gaussian","nb", "pg", "zinb"), de_type = "mean", p = 2, repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
  seed = idx + 62
  set.seed(seed)
  
  effect_size  = NULL
  extra_args <- list(...)
  
  if (distribution == "gaussian") {
    # simu_hier uses alpha1, alpha2, tau1, tau2, sigma1, sigma2
    simData <- simu_hier_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper, ...)
    effect_size <- switch(de_type,
                          mean = extra_args$alpha2 - extra_args$alpha1,
                          variance = extra_args$sigma2 - extra_args$sigma1,
                          stop("Invalid de_type for gaussian"))
  } else if (distribution == "nb") {
    # simu_NB uses mu1, mu2, theta1, theta2, sigma_sq1, sigma_sq2
    effect_size <- switch(de_type,
                          mean = extra_args$mu2 - extra_args$mu1,
                          variance = (extra_args$mu2 + extra_args$mu2^2/extra_args$theta2),
                          stop("Invalid de_type for nb"))
    simData <- simu_NB_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper, ...)
  } else if (distribution == "pg") {
    # simu_poigamma uses alpha1, alpha2, beta1, beta2
    simData <- simu_poigamma_vary_cellcts(seed = seed, lower = 300, upper = 1000, n1 = n1, n2 = n2,...)
  } else if (distribution == "zinb"){
    # simu_zinb uses mu, mu2, theta, b_theta, b_pi, pi2 for mixture
    effect_size <- switch(de_type,
                          dispersion = {
                            theta2 <- b_theta * theta
                            theta2 - theta
                          },
                          mixture = extra_args$b_pi,
                          mean = {
                            warning("Mean shift inherently includes variance shift")
                            mu2 - mu
                          },
                          {
                            warning("Invalid de_type for zinb; returning NULL")
                            NULL
                          })
    simData <- simu_zinb_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper,  de_type = de_type, ...)
  } else{
    print("Invalid distribution type")
    return(NULL)
  }
  
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) 
  y2 = as.vector(unlist(Y2))
  
  # Moment estimator
  start.time.mom = Sys.time()
  Cov1 = mom_cov_diff_cell_cts(Y1, p)
  Cov2 = mom_cov_diff_cell_cts(Y2, p)
  T1 = y1
  for (dim in 2:p){
    T1 = cbind(T1, y1^dim)
  }
  T2 = y2
  for (dim in 2:p){
    T2 = cbind(T2, y2^dim)
  }
  T1_bar = colMeans(T1)
  T2_bar = colMeans(T2)
  
  chisq_mom = as.numeric(t((T1_bar - T2_bar)) %*% diag((1 / diag(Cov1 + Cov2))) %*% (T1_bar - T2_bar))
  pv_mom = 1 - pchisq(chisq_mom, df = p) #we use p
  end.time.mom = Sys.time()
  time.taken.mom = end.time.mom - start.time.mom
  return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, pval = pv_mom, chisq_mom = chisq_mom, Cov1 = Cov1, Cov2 = Cov2, p = p, time_taken = time.taken.mom))
}
#  ---- PSEUDOBULK  ----
old_simulatePseudobulkRealistic = function(distribution = c("gaussian","nb", "pg", "zinb"), test_type = c("ttest","ks"), de_type = "mean", repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
  seed = idx + 62
  set.seed(seed)
  
  effect_size  = NULL
  extra_args <- list(...)
  for (argname in names(extra_args)) {
    assign(argname, extra_args[[argname]])
  }
  
  if (distribution == "gaussian") {
    # simu_hier uses alpha1, alpha2, tau1, tau2, sigma1, sigma2
    simData <- simu_hier_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper, ...)
    effect_size <- switch(de_type,
                          mean = extra_args$alpha2 - extra_args$alpha1,
                          variance = extra_args$sigma2 - extra_args$sigma1,
                          stop("Invalid de_type for gaussian"))
  } else if (distribution == "nb") {
    # simu_NB uses mu1, mu2, theta1, theta2, sigma_sq1, sigma_sq2
    effect_size <- switch(de_type,
                          mean = extra_args$mu2 - extra_args$mu1,
                          variance = (extra_args$mu2 + extra_args$mu2^2/extra_args$theta2),
                          stop("Invalid de_type for nb"))
    simData <- simu_NB_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper, ...)
  } else if (distribution == "pg") {
    # simu_poigamma uses alpha1, alpha2, beta1, beta2
    simData <- simu_poigamma_vary_cellcts(seed = seed, lower = lower, upper = upper, n1 = n1, n2 = n2,...)
  } else if (distribution == "zinb"){
    # simu_zinb uses mu, mu2, theta, b_theta, b_pi, pi2 for mixture
    effect_size <- switch(de_type,
                          dispersion = {
                            theta2 <- b_theta * theta
                            theta2 - theta
                          },
                          mixture = b_pi,
                          mean = {
                            warning("Mean shift inherently includes variance shift")
                            mu2 - mu
                          },
                          {
                            warning("Invalid de_type for zinb; returning NULL")
                            NULL
                          })
    simData <- simu_zinb_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper,  de_type = de_type,...)
  } else{
    print("Invalid distribution type")
    return(NULL)
  }
  
  Y1 = simData$Y1
  Y2 = simData$Y2
  
  pbY1 = vapply(Y1, mean, numeric(1))
  pbY2 = vapply(Y2, mean, numeric(1))
  
  pb.df1 = data.frame(pbY1, rep(0,length(n1)))
  colnames(pb.df1) = c("pseudobulk", "group")
  pb.df2 = data.frame(pbY2, rep(1,length(n2)))
  colnames(pb.df2) = c("pseudobulk", "group")
  pb.total = rbind(pb.df1, pb.df2)
  
  if (test_type == "ttest") {
    pb.test = t.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    
    pv_pb = pb.test$p.value
    
    mean1_pb = pb.test$estimate[1]
    mean2_pb = pb.test$estimate[2]
    
    test_stat_pb = pb.test$statistic
    sd_pb = pb.test$stderr
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, mean1_pb = mean1_pb, mean2_pb = mean2_pb, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, sd = sd_pb, test_type = test_type))
    
  } else if(test_type == "ks"){
    
    pb.test = ks.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    pv_pb = pb.test$p.value
    test_stat_pb = pb.test$statistic
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
    
  } else{
    warning("Please use valid test")
    return(NULL)
  }
}

simulatePseudobulkRealistic = function(distribution = c("gaussian", "pg", "zinb"), test_type = c("ttest","ks","ftest", "wilcoxon"), de_type = "mean", repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
  seed = idx + 62
  set.seed(seed)
  
  effect_size  = NULL
  extra_args <- list(...)
  for (argname in names(extra_args)) {
    assign(argname, extra_args[[argname]])
  }
  
  if (distribution == "gaussian") {
    # simu_hier uses alpha1, alpha2, tau1, tau2, sigma1, sigma2
    simData <- simu_hier_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper, ...)
    effect_size <- switch(de_type,
                          mean = extra_args$alpha2 - extra_args$alpha1,
                          variance = extra_args$sigma2 - extra_args$sigma1,
                          stop("Invalid de_type for gaussian"))
  } else if (distribution == "pg") {
    # simu_poigamma uses alpha1, alpha2, beta1, beta2
    simData <- simu_poigamma_vary_cellcts(seed = seed, lower = lower, upper = upper, n1 = n1, n2 = n2,...)
  } else if (distribution == "zinb"){
    # simu_zinb uses mu, mu2, theta, b_theta, b_pi, pi2 for mixture
    effect_size <- switch(de_type,
                          dispersion = {
                            theta2 <- b_theta * theta
                            theta2 - theta
                          },
                          mixture = b_pi,
                          mean = {
                            warning("Mean shift inherently includes variance shift")
                            mu2 - mu
                          },
                          {
                            warning("Invalid de_type for zinb; returning NULL")
                            NULL
                          })
    simData <- simu_zinb_vary_cellcts(seed = seed, n1 = n1, n2 = n2, lower = lower, upper = upper,  de_type = de_type, ...)
  } else{
    print("Invalid distribution type")
    return(NULL)
  }
  
  Y1 = simData$Y1
  Y2 = simData$Y2
  
  pbY1 = vapply(Y1, mean, numeric(1))
  pbY2 = vapply(Y2, mean, numeric(1))
  
  pb.df1 = data.frame(pbY1, rep(0,length(n1)))
  colnames(pb.df1) = c("pseudobulk", "group")
  pb.df2 = data.frame(pbY2, rep(1,length(n2)))
  colnames(pb.df2) = c("pseudobulk", "group")
  pb.total = rbind(pb.df1, pb.df2)
  
  if (test_type == "ttest") {
    pb.test = t.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    
    pv_pb = pb.test$p.value
    
    mean1_pb = pb.test$estimate[1]
    mean2_pb = pb.test$estimate[2]
    
    test_stat_pb = pb.test$statistic
    sd_pb = pb.test$stderr
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, mean1_pb = mean1_pb, mean2_pb = mean2_pb, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, sd = sd_pb, test_type = test_type))
    
  } else if(test_type == "ks"){
    
    pb.test = ks.test(pb.df1$pseudobulk , pb.df2$pseudobulk)
    pv_pb = pb.test$p.value
    test_stat_pb = pb.test$statistic
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
    
  } else if(test_type == "ftest"){
    pb.test = var.test(pb.df1$pseudobulk, pb.df2$pseudobulk, ratio = 1, alternative = "two.sided")
    test_stat_pb = pb.test$statistic
    pv_pb = pb.test$p.value
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
    
  } else if(test_type == "wilcoxon"){
    pb.test = wilcox.test(pb.df1$pseudobulk, pb.df2$pseudobulk, alternative = "two.sided", exact = NULL, conf.int = FALSE)
    pv_pb = pb.test$p.value
    test_stat_pb = pb.test$statistic
    return(list(distribution = distribution, effect_size = effect_size, de_type = de_type, Y1 = Y1, Y2 = Y2,
                pval = pv_pb, test_stat = test_stat_pb, test_type = test_type))
  } else{
    warning("Please use valid test")
    return(NULL)
  }
}

