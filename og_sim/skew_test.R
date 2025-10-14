# normal vs gamma; allow for fixed individual errors? thought about not using a hierarchical model...

# note: can consider G1 ~N(1,1) G2 ~ Exp(1) alternatively
sim_skew_comparison <- function(seed = 100, lower = 300, upper = 1000, n1 = 100, n2 = 100,
                                shape2 = 9, scale2 = 2/3, mu1 = 6, var = 4, tau1 = 0.1, tau2 = 0.1) {
  set.seed(seed)
  Y1_list <- vector("list", n1)
  Y2_list <- vector("list", n2)
  for (i in seq_len(n1)) {
    eps = rnorm(1, mean = 0, sd = sqrt(tau1))
    n_cells_i <- sample(lower:upper, 1)
    mu_i = mu1 + eps #comment out after
    Y1_list[[i]] <- rnorm(n_cells_i, mean = mu_i, sd = sqrt(var)) #+ eps  #Group 1 (normal)
  }

  for (i in seq_len(n2)) {
    eps = rnorm(1, mean = 0, sd = sqrt(tau2))
    n_cells_i <- sample(lower:upper, 1)
    shape2_i = shape2 + eps
    Y2_list[[i]] <- rgamma(n_cells_i, shape = shape2_i, scale = scale2) #+ eps  # Group 2 (gamma)
  }

  return(list(Y_gaussian = Y1_list, Y_gamma = Y2_list))
}


out <- sim_skew_comparison(seed = 1, n1 = 1000, n2 = 1000)
gauss_data <- unlist(out$Y_gaussian)
skew_data <- unlist(out$Y_gamma)
summary_stats <- data.frame(
  model = c("Gaussian", "Skewed (Gamma)"),
  mean = c(mean(gauss_data), mean(skew_data)),
  variance = c(var(gauss_data), var(skew_data)),
  skewness = c(skewness(gauss_data), skewness(skew_data)))
print(summary_stats)

simulateSkewRealstic = function(repID = NULL, p = 2,  K = NULL, lower = 300, upper = 1000, n1 = 100, n2 = 100, 
                                idx =1, plot_flag = F,
                                tau1 = 0.01, tau2 = 0.01, mu1 = 6, var = 4, shape2 = 9, scale2 = 2/3){
  seed = idx + 62
  set.seed(seed)
  
  simData = sim_skew_comparison(seed = seed, lower =lower, upper = upper, n1 = n1, n2 = n2,
                               tau1 = tau1, tau2 = tau2, mu1 = mu1, var =var, shape2 = shape2, scale2 = scale2)
  Y1 = simData$Y_gaussian
  Y2 = simData$Y_gamma
  y1 = as.vector(unlist(Y1)) #vector of all cells in grp 1
  y2 = as.vector(unlist(Y2)) #vector of all cells in grp 2
  
  if(plot_flag == T){
    df <- bind_rows(
      data.frame(value = y1, group = "Group 1"),
      data.frame(value = y2, group = "Group 2")
    )
    
    hist_data <- df %>%
      group_by(group) %>%
      summarise(hist = list(hist(value, breaks = K, plot = FALSE)), .groups = "drop") %>%
      mutate(hist_df = lapply(hist, \(h) data.frame(x = h$mids, density = h$density))) %>%
      unnest(hist_df)
    
    raw_plt = ggplot(hist_data, aes(x = x, fill = group)) +
      geom_ribbon(aes(ymin = 0, ymax = density), alpha = 0.3, color = NA) +
      scale_fill_manual(values = c("Group 1" = "firebrick", "Group 2" = "dodgerblue")) +
      labs(x = "Expression", y = "Density") +
      theme_minimal() + title("raw histograms")
    print(raw_plt)
  }

  true_effect_size = 2/sqrt(shape2)
  empirical_effect_size = abs(skewness(y2) - skewness(y1))
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
    K = length(est_midpoints)*2
    print(paste0("Number of midpoints K not specified. Using total number of boundary points in pooled data as new K: ", K))
    carrier_est = continuous_discrete_kde(data = y_agg, eval_points = est_midpoints, bw = 0.05)
    
    Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
    Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
    Smat1 = get_countdata_bin_counts(Y = Y1, grid_points = est_grid)
    Smat2 = get_countdata_bin_counts(Y = Y2, grid_points = est_grid)
  } else{
    est_grid = seq(l, u, length.out = K + 1)
    est_midpoints  = (est_grid[-1] + est_grid[-length(est_grid)])/2
    carrier = density(x = y_agg, kernel = c("gaussian"), bw = "nrd0", n = K, from = l, to = u )
    carrier_est = carrier$y
    
    Smat1 = matrix(NA, nrow = n1, ncol = K) #group 1 COUNTS 
    Smat2 = matrix(NA, nrow = n2, ncol = K) #group 2 COUNTS 
    Smat1 = get_sc_bin_counts(Y = Y1, bin_edges = est_grid, K = K)
    Smat2 = get_sc_bin_counts(Y = Y2, bin_edges = est_grid, K = K)
  }
  binwidth = (u - l)/K #in theory should be equal to diff(est_midpoints)
  
  carrier_plot = NULL
  breaks = NULL
  if(plot_flag == T){
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
  }
  
  Smat = rbind(Smat1, Smat2) #combine two bin matrices, used in SEF regression
  print("running SEF regression model")
  #Construct design matrix
  tk = est_midpoints
  mean_tk = mean(tk)
  X = rep(1, length(est_midpoints))
  for (dim in 1:p){
    X = cbind(X, (tk - mean_tk)^dim/sd(tk))
  }
  varb = paste0("X_", 1:p)
  colnames(X) = c("Intercept", varb)
  
  S1sum = colSums(Smat1)
  S2sum = colSums(Smat2)
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
      X = cbind(X, (tk - mean_tk)^dim/sd(tk))
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
  
  plt = NULL
  if(plot_flag == T){
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
  }
  
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
  return(list(true_effect_size = true_effect_size, empirical_effect_size = empirical_effect_size, y_agg = y_agg, 
              est_midpoints = est_midpoints, carrier_est = carrier_est, binwidth = binwidth, K = K, y1 = y1, y2 = y2,
              beta_est1 = beta_est1, beta_est2 = beta_est2, beta_diff = beta_diff,
              sef_df1 = sef_df1, sef_df2 = sef_df2, pval_1 = pval_1, Cov_1 = Cov_bar,
              carrier_plt = carrier_plot, comparison_plt = plt, breaks = breaks))
}

output2 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = 10, plot_flag = T)
output3 = simulateSkewRealstic(repID = NULL, p = 3, K = 50, idx = 10, plot_flag = T)

library(parallel) #depends on the individual specific error
library(doSNOW)
library(foreach)
library(doParallel)
n1 = 100; n2 = 100
lower = 300; upper = 1000
seed = 100
maxIter = 300
pv_p2 = NULL
pv_p3 = NULL
nCPUS = 5
tau1 = 0.1; tau2 = 0.1
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}
output <- foreach(i = 1:maxIter, .packages = c("ggplot2","e1071","dplyr","tidyr")) %dopar% {
  # p2 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = i, tau1 = tau1, tau2 = tau2)
  # p3 = simulateSkewRealstic(repID = NULL, p = 3, K = 50, idx = i, tau1 = tau1, tau2 = tau2)
  p2 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = 10, plot_flag = T,
                            mu1 = 2, var = 1, shape2 = 4, scale2 = 0.5)
  p3 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = 10, plot_flag = T,
                            mu1 = 2, var = 1, shape2 = 4, scale2 = 0.5)
  list(p2 = p2, p3 = p3)
}
pv_p2 = cbind(pv_p2, sapply(output, function(x) x$p2$pval_1))
pv_p3 = cbind(pv_p3, sapply(output, function(x) x$p3$pval_1))
stopCluster(cl)
q = 0.05 
print(mean(pv_p2 < 0.05))
print(mean(pv_p3 < 0.05))

GraphName = "./density_estimation/simulations/supplementary2025/updated_gaussian_vs_gamma_p2_test.pdf"
pdf(GraphName)
gg_qqplot(pv_p2) + ggtitle("Skewness Test (p = 2)")
dev.off()

# ----- POWER ANALYSIS -----
grids = read.csv("./density_estimation/simulations/supplementary2025/skewness_grid.csv")
#grids = grids[3:11,]
grids <- grids[nrow(grids):1, ]

mus = grids$mean
vars = grids$variance
shapes = grids$shape_k
scales = grids$scale_theta
skewness_effect_size = grids$skewness
n1 = 100; n2 = 100
lower = 300; upper = 1000
seed = 100
maxIter = 200
pv_p2 = NULL
pv_p3 = NULL
nCPUS = 5
tau1 = 0.1; tau2 = 0.1
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}
for (i in c(1:length(mus))) {
  mu_i = mus[i]; var_i = vars[i]
  shape_i = shapes[i]; scale_i = scales[i]
  print(c(mu_i, var_i, shape_i, scale_i))
  print(paste0("effect size for skewness: ", skewness_effect_size[i]))

  output <- foreach(i = 1:maxIter, .packages = c("ggplot2","e1071","dplyr","tidyr")) %dopar% {
    p2 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = i, tau1 = tau1, tau2 = tau2,
                              mu1 = mu_i, var = var_i, shape2 = shape_i, scale2 = scale_i)
    p3 = simulateSkewRealstic(repID = NULL, p = 3, K = 50, idx = i, tau1 = tau1, tau2 = tau2,
                              mu1 = mu_i, var = var_i, shape2 = shape_i, scale2 = scale_i)
    list(p2 = p2, p3 = p3)
  }
  pv_p2 = cbind(pv_p2, sapply(output, function(x) x$p2$pval_1))
  pv_p3 = cbind(pv_p3, sapply(output, function(x) x$p3$pval_1))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("p = 2", "p = 3")
df_pw = data.frame(effect_size = skewness_effect_size, p_2 = colMeans( pv_p2 < q ),
                   p_3 = colMeans(pv_p3 < q))

matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-skewness-n1-", n1, "-n2-", n2, "-p-2-3.csv" )
write_csv(df_pw, matrixName)
df_pw = read.csv(matrixName)
df_pw = df_pw[2:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("p_2" = "orangered2", "p_3" = "royalblue2")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  #theme_grey(base_size = 22) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Skewness (" ~ 2 / sqrt(k) ~ ")"), 
       y = "Power", 
       title = "Empirical Power (Skewness)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/Powerplot-skewness-n1-", n1, "-n2-", n2, "-p-2-3.pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

out2 = simulateSkewRealstic(repID = NULL, p = 2, K = 50, idx = 100, tau1 = tau1, tau2 = tau2,
                          mu1 = 6, var = 9, shape2 = 4, scale2 = 1.5,plot_flag = T)
out3 = simulateSkewRealstic(repID = NULL, p = 3, K = 50, idx = i, tau1 = tau1, tau2 = tau2,
                            mu1 = 2, var = 1, shape2 = 4, scale2 = 0.5,plot_flag = T)

