rm(list = ls())  
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
library(truncnorm)

########## 1. DATA GENERATION FUNCTIONS ##########
#  ---- SAME CELLS ACROSS ALL PEOPLE  ----
simu_hier = function(alpha1, alpha2, tau1, tau2, sigma1, sigma2, n1, n2, ncells){
  #simulates hierarchical data with y_ij ~ N (alpha_i, sigma^2) where alpha_i ~ N(alpha, tau^2) where y_ij ~ f_i and E[f_i] ~ N(alpha, sigma^2 + tau^2)
  #equivalent to SEF Hierarchical Model for fixed beta1, beta2 = 0
  Y1 = matrix(data = NA, nrow = n1, ncol = ncells)
  Y2 = matrix(data = NA, nrow = n2, ncol = ncells)
  for (i in 1:n1) { #GROUP 1
    alpha_1i = rnorm(1, alpha1, sqrt(tau1))
    Y1[i,] = rnorm(ncells, mean = alpha_1i, sd = sqrt(sigma1))
  }
  for (i in 1:n2) { #GROUP 2 
    alpha_2i = rnorm(1, alpha2, sqrt(tau2))
    Y2[i,] = rnorm(ncells, mean = alpha_2i, sd = sqrt(sigma2))
  }
  return(list(Y1 = Y1, Y2 = Y2))
} 

simu_poigamma = function(alpha1, alpha2, beta1, beta2, n1, n2, ncells){
  # Simulate Poisson-Gamma mixture model data
  # Empirical Bayes formulation with y_ij ~ Poi(Lambda), Lambda ~ gamma(alpha, beta)
  #alpha and beta are the shape and size parameters
  Y1 = matrix(data = NA, nrow = n1, ncol = ncells)
  Y2 = matrix(data = NA, nrow = n2, ncol = ncells)
  for (i in 1:n1) { #GROUP 1
    lambda1i = rgamma(1, shape = alpha1, rate = beta1)
    Y1[i,] = rpois(ncells, lambda1i)
  }
  for (i in 1:n2) { #GROUP 2 
    lambda2i = rgamma(1, shape = alpha2, rate = beta2)
    Y2[i,] = rpois(ncells, lambda2i)
  }  
  return(list(Y1 = Y1, Y2 = Y2))
}

simu_zinb <- function(n1, n2, ncells, mu, mu2, sigma_sq = 0.5, theta, pi = 0.5, pi2 = 0.5,
                      b_pi = 1, b_theta = 1, de_type = 'dispersion', seed = 123) {
  # Set seed for reproducibility
  set.seed(seed)
  # Initialize matrices for group 1 and group 2
  Y1 <- matrix(data = NA, nrow = n1, ncol = ncells)
  Y2 <- matrix(data = NA, nrow = n2, ncol = ncells)
  
  mu_i_group1 <- rtruncnorm(n1, a = 0, mean = mu, sd = sqrt(sigma_sq)) # ind-specific means for Group 1
  
  # Adjust parameters for Group 2 based on DE type
  if (de_type == 'dispersion') {
    mu2 <- mu
    theta2 <- theta*b_theta #scale
    pi2 = pi
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq)) #generate individual specific means
  } else if (de_type == 'mixture') {
    if(b_pi < 0 || b_pi > 1 ){
      stop("mixture must be between 0 and 1")
    }
    theta2 <- theta*b_theta
    pi2 = pi2
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
    for (i in 1:n1) { # Simulate ZINB counts for Group 1
      for (j in 1:ncells) {
        if (runif(1) < pi) {
          Y1[i, j] <- 0
        } else {
          Y1[i, j] <- b_pi*rnbinom(1, size = theta, mu = mu_i_group1[i]) + (1 - b_pi)*rnbinom(1, size = theta2, mu = mu_i_group2[i])
        }
      }
    }
    for (i in 1:n2) { # Simulate ZINB counts for Group 2
      for (j in 1:ncells) {
        if (runif(1) < pi2) {
          Y2[i, j] <- 0
        } else {
          Y2[i, j] <- (1 - b_pi)*rnbinom(1, size = theta, mu = mu_i_group1[i]) + b_pi*rnbinom(1, size = theta2, mu = mu_i_group2[i])
        }
      }
    }
    return(list(Y1 = Y1, Y2 = Y2))
    
  } else if (de_type == "mean"){
    theta2 = theta
    pi2 = pi
  } else {
    stop('Invalid de_type specified. Choose either "mean" or "variance".')
  }
  
  # Simulate ZINB counts for Group 1
  for (i in 1:n1) {
    for (j in 1:ncells) {
      if (runif(1) < pi) {
        Y1[i, j] <- 0
      } else {
        Y1[i, j] <- rnbinom(1, size = theta, mu = mu_i_group1[i])
      }
    }
  }
  
  # Simulate ZINB counts for Group 2
  for (i in 1:n2) {
    for (j in 1:ncells) {
      if (runif(1) < pi2) {
        Y2[i, j] <- 0
      } else {
        Y2[i, j] <- rnbinom(1, size = theta2, mu = mu_i_group2[i])
      }
    }
  }
  
  return(list(Y1 = Y1, Y2 = Y2))
}
simu_general <- function(distribution_choice, n1 = 100, n2 = 100, ncells = 500, ...) {
  dispatch_list <- list(
    normal = simu_hier,
    poigamma = simu_poigamma,
    zinb = simu_zinb
  )
  
  if (!(distribution_choice %in% names(dispatch_list))) {
    stop("Invalid distribution_choice. Choose from: ", paste(names(dispatch_list), collapse = ", "))
  }
  
  chosen_fun <- dispatch_list[[distribution_choice]]
  
  # Call with common + distribution-specific args
  result <- chosen_fun(n1 = n1, n2 = n2, ncells = ncells, ...)
  return(result)
}

#  ---- DIFFERING CELLS ACROSS ALL PEOPLE  ----
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
    # Y1[i,] = rpois(ncells, lambda1i)
  }
  for (i in 1:n2) { #GROUP 2 
    n_cells_i <- sample(lower:upper, 1)
    lambda2i = rgamma(1, shape = alpha2, rate = beta2)
    Y2_list[[i]] <- rpois(n_cells_i, lambda2i)
    # Y2[i,] = rpois(ncells, lambda2i)
  }  
  return(list(Y1 = Y1_list, Y2 = Y2_list))
}
simu_NB_vary_cellcts <- function(seed = 100, mu1 = 5, mu2 = 5, theta1 = 1, theta2 = 1, sigma_sq1 = 0.1, sigma_sq2 = 0.1, lower = 300, upper = 1000, n1, n2) {
  Y1_list <- vector("list", n1) 
  Y2_list <- vector("list", n2)
  if (mu1 < 0 || mu2 < 0) {
    print("Must be positive mean parameter.")
    return(NULL)
  } else if(theta1 < 0 || theta2 < 0){
    print("Must be positive dispersion parameter.")
    return(NULL)
  }
  set.seed(seed = seed)
  for (i in seq_len(n1)) {  # GROUP 1
    n_cells_i <- sample(lower:upper, 1)
    mu_1i <- rtruncnorm(1, a = 0, mu1, sqrt(sigma_sq1))
    # var1 = mu1 + mu1^2/theta1
    Y1_list[[i]] <- rnbinom(n = n_cells_i, mu = mu1, size = theta1)
  }
  
  for (i in seq_len(n2)) {  # GROUP 2
    n_cells_i <- sample(lower:upper, 1)
    mu_2i <- rtruncnorm(1, a = 0, mu2, sqrt(sigma_sq2))
    # var2 = mu2 + mu2^2/theta2
    Y2_list[[i]] <- rnbinom(n = n_cells_i, mu = mu2, size = theta2)
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
    if (b_pi < 0 || b_pi > 1) stop("mixture must be between 0 and 1")
    theta2 <- theta * b_theta
    mu_i_group2 <- rtruncnorm(n2, a = 0, mean = mu2, sd = sqrt(sigma_sq))
    
    for (i in seq_len(n1)) {
      n_cells_i <- sample(lower:upper, 1)
      Y1_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
        if (runif(1) < pi) {
          0
        } else {
          b_pi * rnbinom(1, size = theta, mu = mu_i_group1[i]) +
            (1 - b_pi) * rnbinom(1, size = theta2, mu = mu_i_group2[i])
        }
      })
    }
    for (i in seq_len(n2)) {
      n_cells_i <- sample(lower:upper, 1)
      Y2_list[[i]] <- sapply(seq_len(n_cells_i), function(j) {
        if (runif(1) < pi2) {
          0
        } else {
          (1 - b_pi) * rnbinom(1, size = theta, mu = mu_i_group1[i]) +
            b_pi * rnbinom(1, size = theta2, mu = mu_i_group2[i])
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

edited_simu_zinb_vary_cellcts <- function(n1, n2, mu, mu2, theta, sigma_sq = 0.5, pi = 0.5, pi2 = 0.5,
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
#  ---- DISCRETE GAUSSIAN KERNEL FUNCTION (USED IN RDA AS WELL) ----
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
  #COMPUTES pth MOMENTS and stores as a matrix
  return(sapply(1:p, function(k) y^k))
} 

mom_cov = function(Y, p) {
  # Compute the Covariance matrix for Moments (for MoM estimator)
  n = nrow(Y) # individuals per group
  total_samples <- 0
  overall_sum <- rep(0, p)  # Sum for overall mean
  within_var_sum <- matrix(0, nrow = p, ncol = p)  # Within-individual variance
  between_var_sum <- matrix(0, nrow = p, ncol = p)  # Between-individual variance 
  
  mean_list = matrix(0, nrow = nrow(Y), ncol = p) # matrix of p moments for each individual
  m_vector <- numeric(n)  # Store the number of samples for each individual
  
  #within individual variance
  for (i in 1:nrow(Y)){ #iterate for each individual
    samples_i = Y[i, ]
    m_i = length(samples_i) 
    m_vector[i] = m_i
    mean_samples_i = mean(samples_i)
    moments_i = t(sapply(samples_i, compute_moments, p = p)) #ncell_i x 1 moment 
    
    mean_t_i = colMeans(moments_i) # Compute the mean of p moments for individual
    mean_list[i, ] = mean_t_i
    # Update overall sum for computing the overall mean
    overall_sum = overall_sum + m_i * mean_t_i
    total_samples = total_samples + m_i 
    # Compute within-individual variance
    diffs = t(moments_i) - mean_t_i #(T_i - Tbar_i)
    within_var_sum = within_var_sum + diffs%*%t(diffs) #sum (T_i - Tbar_i)^2
  }
  # Compute the overall mean of the moments
  overall_mean <- overall_sum / total_samples
  for (i in 1:nrow(Y)){
    m_i = m_vector[i]
    mean_t_i = mean_list[i, ]
    diff <- mean_t_i - overall_mean
    between_var_sum <- between_var_sum + m_i^2 * (diff %*% t(diff)) # Should be m_i^2 not m_i
  }
  Cov = within_var_sum/total_samples^2 + between_var_sum/total_samples^2
  return(Cov)
}

mom_cov_diff_cell_cts = function(Y, p) {
  # Compute the Covariance matrix for Moments (for MoM estimator) with different cell counts
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
  
  # Extract the specific attribute from each simulation result
  attr_list <- lapply(results, `[[`, attr_name)
  
  # Convert the list of vectors into a data frame with rows for each simulation
  attr_df <- do.call(rbind, lapply(attr_list, function(x) as.data.frame(t(x))))
  
  # Add a simulation ID column
  #attr_df <- cbind(sim_id = 1:n_sim, attr_df)
  
  return(attr_df)
}

extract_and_combine_scalar <- function(results, attr_name) {
  # Extract the specific scalar attribute from each simulation result
  attr_list <- lapply(results, function(res) res[[attr_name]])
  
  # Convert the list into a data frame
  attr_vector <- unlist(attr_list, use.names = FALSE)
  
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
simulateNormal = function(de_type = "mean", p = 2, repID = NULL ,K = 75, n1 = 100, n2 = 100, ncells = 500, idx, alpha1, alpha2, tau1, tau2, sigma1, sigma2){
  if (de_type == "mean") {
    effect_size = alpha2 - alpha1
  } else if(de_type == "variance"){
    effect_size = sigma2 - sigma1
  } else{
    warning("Shift type not found. Please use valid shift type.")
    return(NULL)
  }
  seed = idx + 62
  set.seed(seed)
  simData = simu_hier(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2, ncells = ncells)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(c(Y1)) 
  y2 = as.vector(c(Y2))
  y_agg = c(y1, y2)
  # Carrier density 
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u ) 
  est_grid = seq(l, u, length.out = K + 1) 
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  # Counts in each bin 
  Smat1 = get_bin_counts(Y1, bin_edges = est_grid, K)
  Smat2 = get_bin_counts(Y2, bin_edges = est_grid, K)
  Smat = rbind(Smat1, Smat2) 
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
  carrier_est = carrier$y
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
  #setup for modeling
  cellSum1 = ncells*n1 
  cellSum2 = ncells*n2
  scale_factor1 = cellSum1/sum(carrier$y) 
  scale_factor2 = cellSum2/sum(carrier$y)
  carrier_scale1 = carrier_est*scale_factor1
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  # Fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  # Marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1)# times the bindwidth or divide the sum of density
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M = dnorm(pairwise_diff, sd = carrier$bw)/(ncells*(n1+n2)) # Divided by n = n1+n2
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
  
  # Inference without intercept effects
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
  
  #save values: effect size, binwidth, betas, covariances
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, ncells = ncells, n1 = n1, n2 = n2, binwidth = binwidth,
              carrier_est = carrier_est, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2, est_midpoints = est_midpoints,
              pval_1 = pval_1, pval_3 = pval_3, chi_stat1 = chi_stat1, chi_stat3 = chi_stat3, Cov_1 = Cov_bar, Cov_3 =Cov_bar_sub,
              beta_diff = beta_diff, beta_diff_c = beta_diff_c))
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

simulatePG = function(de_type = "mean", p = 2, K = 75,  repID = NULL, n1 = 100, n2 = 100, ncells = 500, idx, alpha1, alpha2, beta1, beta2){
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
  simData = simu_poigamma(alpha1 = alpha1, alpha2 = alpha2, beta = beta1, beta2 = beta2,n1 = n1, n2 = n2, ncells = ncells)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(c(Y1)) 
  y2 = as.vector(c(Y2))
  y_agg = c(y1, y2)
  # Carrier density 
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u ) 
  est_grid = seq(l, u, length.out = K + 1) 
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  # Counts in each bin 
  Smat1 = get_bin_counts(Y1, bin_edges = est_grid, K)
  Smat2 = get_bin_counts(Y2, bin_edges = est_grid, K)
  Smat = rbind(Smat1, Smat2) 
  tk = est_midpoints
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, tk^dim)
  } 
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb) 
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
  scale_factor = (ncells*nrow(Smat))/sum(carrier$y) 
  carrier_est = carrier$y*scale_factor # Better to sum  
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
  #setup for modeling
  cellSum1 = ncells*n1 
  cellSum2 = ncells*n2
  scale_factor1 = cellSum1/sum(carrier$y) 
  scale_factor2 = cellSum2/sum(carrier$y)
  carrier_scale1 = carrier_est*scale_factor1
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
  colnames(df1)[1] = "sum_cts"
  colnames(df1)[2] = "carrier_scale"
  colnames(df2)[1] = "sum_cts" 
  colnames(df2)[2] = "carrier_scale"
  formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  # Fit glm model for each group
  sef_gp1 = tryCatch(glm( formula, family = poisson(link="log"), data = df1))  # Some numerical issue
  sef_gp2 = tryCatch(glm( formula, family = poisson(link="log"), data = df2))  # Some numerical issue
  beta_est1 = as.vector(sef_gp1$coefficients)
  beta_est2 = as.vector(sef_gp2$coefficients)
  beta_diff = beta_est1 - beta_est2 
  sef_df1 = as.vector( carrier_est * exp(X %*% beta_est1) )
  sef_df2 = as.vector( carrier_est * exp(X %*% beta_est2) )
  
  #1. marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1)# times the bindwidth or divide the sum of density
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M = dnorm(pairwise_diff, sd = carrier$bw)/(ncells*(n1+n2)) # Divided by n = n1+n2
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
  
  
  #3. Inference without intercept effects
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
  
  
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, ncells = ncells, n1 = n1, n2 = n2, binwidth = binwidth,
              carrier_est = carrier_est, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2, est_midpoints = est_midpoints,
              pval_1 = pval_1, pval_3 = pval_3, chi_stat1 = chi_stat1, chi_stat3 = chi_stat3, Cov_1 = Cov_bar, Cov_3 = Cov_bar_sub,
              beta_diff = beta_diff, beta_diff_c = beta_diff_c))
}
#  ---- NEGATIVE BINOMIAL  ----
simulateNBRealistic = function(repID = NULL, de_type = "mean", p = 2,  K = NULL, lower = 300, upper = 1000, n1 = 100, n2 = 100, idx = 1, 
                               mu1 = 5, mu2 = 5, theta1 = 1, theta2 = 1, sigma_sq1 = 0.1, sigma_sq2 = 0.1){
  if (de_type == "mean") {
    effect_size = mu2 - mu1
  } else if(de_type == "variance"){
    var1 = mu1 + mu1^2/theta1
    var2 = mu2 + mu2^2/theta2
    effect_size = var2 - var1 #variance test
  } else{
    warning("Shift type not found. Please use valid shift type.")
    return(NULL)
  }
  seed = idx + 62
  set.seed(seed)
  simData = simu_NB_vary_cellcts(seed = seed, lower = lower, upper = upper, n1 = n1, n2 = n2, mu1 = mu1, mu2 = mu2, theta1 = theta1, theta2 = theta2, sigma_sq1 = sigma_sq1, sigma_sq2 = sigma2)
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  y_agg = c(y1, y2)
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
                   fill = "lightgray",
                   color = "white") +
    geom_line(data = line_df, aes(x = x, y = y),
              color = "dodgerblue", linewidth = 1) +
    labs(
      title = "Histogram + Carrier",
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
  
  return(list(de_type = de_type, effect_size = effect_size, y_agg = y_agg, est_midpoints = est_midpoints, carrier_est = carrier_est, binwidth = binwidth, K = K, y1 = y1, y2 = y2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, beta_diff = beta_diff, beta_est1_c = beta_est1_c, beta_est2_c = beta_est2_c, beta_diff_c = beta_diff_c,
              sef_df1 = sef_df1, sef_df2 = sef_df2, sef_df1_c = sef_df1_c, sef_df2_c = sef_df2_c, pval_1 = pval_1, pval_3 = pval_3, Cov_1 = Cov_bar, Cov_3 = Cov_bar_sub,
              carrier_plt = carrier_plot, comparison_plt = plt))
  
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
simulateZINBfixed = function(repID, idx, mu, mu2, sigma_sq = 0.5, theta, b_theta = 1, b_pi = 1, pi = 0.5, pi2 = 0.5, n1, n2, ncells, de_type = "dispersion", K, p){
  seed = idx + 62
  simData = simu_zinb(n1, n2, ncells, mu, mu2, sigma_sq = sigma_sq, theta = theta, pi = pi, pi2 = pi2, 
                      b_pi = b_pi, b_theta = b_theta, de_type = de_type, seed = seed)
  
  Y1 = simData$Y1
  Y2 = simData$Y2
  y1 = as.vector(c(Y1)) 
  y2 = as.vector(c(Y2))
  y_agg = c(y1, y2)
  # Carrier density 
  l = min(y_agg) #lower bound
  u = max(y_agg) #upper bound
  carrier = density(x = y_agg, kernel = c("gaussian"), n = K, from = l, to = u ) 
  est_grid = seq(l, u, length.out = K + 1) 
  est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2  
  Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
  Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
  # Counts in each bin 
  Smat1 = get_bin_counts(Y1, bin_edges = est_grid, K)
  Smat2 = get_bin_counts(Y2, bin_edges = est_grid, K)
  Smat = rbind(Smat1, Smat2) 
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
  carrier_est = carrier$y
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  sum(carrier$y*binwidth)
  
  #setup for modeling
  cellSum1 = ncells*n1 
  cellSum2 = ncells*n2
  scale_factor1 = cellSum1/sum(carrier$y) 
  scale_factor2 = cellSum2/sum(carrier$y)
  carrier_scale1 = carrier_est*scale_factor1
  carrier_scale2 = carrier_est*scale_factor2
  df1 = as.data.frame(cbind( S1sum, carrier_scale1, X))
  df2 = as.data.frame(cbind( S2sum, carrier_scale2, X))
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
  
  #1. marginal distribution (no condition on intercept)
  G_1 = t(X) %*% ( sef_df1 * X )*cellSum1/sum(sef_df1)# times the bindwidth or divide the sum of density
  G_1 
  G_2 = t(X) %*% ( sef_df2 * X)*cellSum2/sum(sef_df2) 
  G_2 
  pairwise_diff = outer( est_midpoints, est_midpoints, "-")
  M = dnorm(pairwise_diff, sd = carrier$bw)/(ncells*(n1+n2)) # Divided by n = n1+n2
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
  
  #2. Schur complement to remove effect of intercept by conditioning on estimated value
  chi_stat2 = as.numeric( beta_diff[-1]%*%solve(Cov_bar[-1, -1] - Cov_bar[-1, 1]%*%t(Cov_bar[1, -1])/Cov_bar[1, 1], beta_diff[-1]) ) # Using Schur decompostion 
  pval_2 = 1 - pchisq(chi_stat2, df = p)
  
  #3. Inference without intercept effects
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
  
  #save values: effect size, binwidth, betas, covariances
  if (de_type == "dispersion") {
    theta2 = b_theta*theta
    effect_size = theta2 - theta
  }
  else if(de_type == "mixture"){
    effect_size = b_pi
  }
  else{
    effect_size = mu2 - mu
  }
  
  return(list(effect_size = effect_size, p = p, K = K, ncells = ncells, n1 = n1, n2 = n2, binwidth = binwidth,
              pval_1 = pval_1, pval_2 = pval_2, pval_3 = pval_3[1,1], sef_gp1 = sef_gp1, sef_gp2 = sef_gp2,
              chi_stat1 = chi_stat1, chi_stat2 = chi_stat2, chi_stat3 = chi_stat3,
              beta_diff = beta_diff, beta_diff_c = beta_diff_c,
              Cov_1 = Cov_bar, Cov_2 = Cov_bar[-1, -1] - Cov_bar[-1, 1]%*%t(Cov_bar[1, -1])/Cov_bar[1, 1], Cov_3 =Cov_bar_sub
  ))
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
simulatePseudobulkRealistic = function(distribution = c("gaussian","nb", "pg", "zinb"), test_type = c("ttest","ks"), de_type = "mean", repID = NULL, n1 = 100, n2 = 100, idx = 1, lower = 300, upper = 1000, ...){
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


########## 5. SEF METHOD FUNCTIONS ##########
#  ---- HELPER FUNCTIONS  ----
compute_mean_median_expr <- function(exprMat, sort.by = c("median", "mean")) {
  sort.by <- match.arg(sort.by)
  gene_indices <- seq_len(nrow(exprMat))
  gene_names <- rownames(exprMat)
  
  # Parallel computation of mean and median
  stats_list <- pbmclapply(
    gene_indices,
    function(i) {
      expr_vals <- exprMat[i, ]
      list(
        gene = gene_names[i],
        mean = mean(expr_vals),
        median = median(expr_vals)
      )
    },
    mc.cores = parallel::detectCores() - 1
  )
  stats_df <- do.call(rbind, lapply(stats_list, as.data.frame))
  stats_df$mean <- as.numeric(stats_df$mean)
  stats_df$median <- as.numeric(stats_df$median)
  stats_df <- stats_df[order(-stats_df[[sort.by]]), ]
  return(stats_df)
}


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

filter_genes_by_disease_donor_counts <- function(donors_to_use_per_gene, covariate_df, min_donors = 40) {
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

hongzhe_filtration_method <- function(sObj, gene_vector, exprMat, n_cores = parallel::detectCores() - 1) {
  donors <- unique(sObj$donor_id)
  
  # Ensure alignment
  stopifnot(ncol(exprMat) == ncol(sObj))
  
  # Convert to character to avoid factor issues
  donor_ids <- as.character(sObj$donor_id)
  
  # Parallel over genes
  results <- pbmclapply(gene_vector, function(gene) {
    if (!(gene %in% rownames(exprMat))) return(NULL)
    
    gene_expr <- exprMat[gene, ]
    passing_donors <- c()
    
    for (donor in donors) {
      cell_idx <- which(donor_ids == donor)
      if (length(cell_idx) < 100) next
      
      donor_expr <- gene_expr[cell_idx]
      nonzero_prop <- sum(donor_expr != 0) / length(donor_expr)
      
      if (nonzero_prop >= 0.2) {
        passing_donors <- c(passing_donors, donor)
      }
    }
    
    return(passing_donors)
  }, mc.cores = n_cores)
  
  names(results) <- gene_vector
  return(results) #should be a list with names for each gene
}
#  ---- RUN MODEL FUNCTIONS  ----

# library(Ake)
# testgene = "HLA-A"
# yt = as.vector(sqrt(exprMatReal[testgene,] + 1/2))
# l = min(yt)
# u = max(yt)
# print(length(unique(yt)))
# h = 0.02
# K = 60
# grid1 =  seq(l, u, length.out = K + 1)
# est<-dke.fun(yt,h,"discrete","GA", x = grid1,a0 = l, a1 = u)
# plot(est)
# est$est.fn

get_carrier_density_norm = function(exprMat, sObj, donor_list, plot_flag = T, h = NULL){
  print("subsetting seurat object")
  sObj_meta = sObj@meta.data[sObj$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(sObj = sObj, exprMat = exprMat, donor_vector = donor_list))
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
    carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
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
    p <- ggplot() +
      geom_histogram(aes(x = y_agg, y = ..density..),
                     breaks = hist_plot_breaks,
                     fill = "lightgray", color = "white") +
      geom_line(aes(x = est_midpoints, y = carrier_est),
                color = "dodgerblue", size = 1) +
      labs(
        title = "Histogram + Carrier Density",
        x = "Expression",
        y = "Density"
      ) + theme_minimal()
  }
  return(list(y_agg = y_agg, carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid, carrier =carrier, carrier_plt = p))
}

get_carrier_density_gauss = function(exprMat, sObj, donor_list, plot_flag = T, transformation = T, h = NULL){
  print("subsetting seurat object")
  sObj_meta = sObj@meta.data[sObj$donor_id %in% donor_list, , drop = FALSE]
  print("subset exprMat vec")
  y_agg = as.vector(subset_exprMat_by_donors(sObj = sObj, exprMat = exprMat, donor_vector = donor_list))
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
    print(paste0("under 100 unique points: using discrete KDE with a Gaussian kernel as carrier density estimate"))
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
    } else{
      K = 2*unique_pts
    }
    print(paste0("K = ", K))
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    # carrier = dke.fun(y_agg, h, "discrete", "GA", x = est_midpoints, a0 = l, a1 = u)
    carrier_est = continuous_discrete_kde(data = y_agg, eval_points = est_midpoints)
    hist_plot_breaks = hist(y_agg)$breaks
    
    # hist_plot_breaks = carrier$hist$breaks
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
         xlab = "Expression",
         ylab = "Density")
    lines(est_midpoints, carrier_est, col = "dodgerblue", lwd = 2)
  }
  return(list(y_agg = y_agg, carrier_est = carrier_est, est_midpoints = est_midpoints, est_grid = est_grid))
}

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


#  ---- PERMUTATION FUNCTIONS  ----
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
        K = 60
        h = 0.02
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
  ego = simplify(ego, cutoff = 0.6, by = "p.adjust")
}








