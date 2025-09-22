# supplementary simulations
rm(list = ls())  
library(dplyr) 
library(truncnorm)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(mvtnorm)
library(MASS) 
library(pbmcapply)
library(patchwork)
library(grid)
library(parallel)
library(doSNOW)
library(foreach)
library(reshape)
source("/Users/zaqian/Desktop/density_estimation/simulations/main062025.R") 
########## QQ Plots ##########
#  ---- GAUSSIAN p = 3  ----
maxIter = 100
nCPUS = 8
tau1 = 2; tau2 = 2
sigma1 = 10; sigma2 = 10
n1 = 100; n2 = 100
p = 3
K = 75
alpha2 = 0
alpha1 = 0
repID = 77
pvVec_margin = NULL 
# helper to run one block using your original code
run_block <- function(idx_range) {
  pbmclapply(idx_range, function(i) {
    simulateNormalRealistic(
      repID = repID, de_type = "mean",
      p = p, idx = i,
      alpha1 = alpha1, alpha2 = alpha2,
      tau1 = tau1, tau2 = tau2,
      sigma1 = sigma1, sigma2 = sigma2,
      n1 = n1, n2 = n2, K = K
    )
  }, mc.cores = nCPUS)
}
# run the three chunks: 1–100, 101–200, 201–300
blocks <- list(1:100, 101:200, 201:300)
outputs_list <- lapply(blocks, run_block)
output <- do.call(c, outputs_list)
pvVec_margin <- sapply(output, function(x) x$pval_1)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, "-K-", K, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Gaussian Model (p = 3)")
print(realistic_qqplot)
dev.off() 
#  ---- GAUSSIAN n1 = 100 n2 = 200  ----
maxIter = 300
nCPUS = 8
tau1 = 2; tau2 = 2
sigma1 = 10; sigma2 = 10
n1 = 100; n2 = 200
p = 2
K = 75
alpha2 = 0
alpha1 = 0
repID = 77


cores <- 5

pv_list <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res <- simulateNormalRealistic(
      repID = repID, de_type = "mean",
      p = 2, idx = i,
      alpha1 = alpha1, alpha2 = alpha2,
      tau1 = tau1, tau2 = tau2,
      sigma1 = sigma1, sigma2 = sigma2,
      n1 = n1, n2 = n2, K = K
    )
    pv <- res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = cores, mc.set.seed = TRUE)

pvVec_margin <- unlist(pv_list, use.names = FALSE)

cat("Kept", length(pvVec_margin), "p-values out of", maxIter, "\n")
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, "-K-", K, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Gaussian Model")
print(realistic_qqplot)
dev.off() 




#  ---- POISSON GAMMA p = 3  ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "PG"
maxIter = 300 
nCPUS = 10 
n1 = 100; n2 = 100
p = 3
alpha2 = 20#5
alpha1 = 20#5
beta1 = 2
beta2 = 2
repID = 77  
pvVec_margin = NULL 
pv_list <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res <- simulatePGRealistic(
      repID = repID,
      idx = i,
      alpha1 = alpha1,
      alpha2 = alpha2,
      beta1 = beta1,
      beta2 = beta2,
      n1 = n1,
      n2 = n2,
      de_type = "mean",
      p = p
    )
    pv <- res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = cores, mc.set.seed = TRUE)

pvVec_margin <- unlist(pv_list, use.names = FALSE)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-15-10-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Poisson Gamma (p = 3)")
print(realistic_qqplot)
dev.off() 

#  ---- POISSON GAMMA n1 = 100 n2 = 200  ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "PG"
maxIter = 300 
nCPUS = 10 
n1 = 100; n2 = 200
p = 2
alpha2 = 20#5
alpha1 = 20#5
beta1 = 2
beta2 = 2
repID = 77  
pvVec_margin = NULL 
pv_list <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res <- simulatePGRealistic(
      repID = repID,
      idx = i,
      alpha1 = alpha1,
      alpha2 = alpha2,
      beta1 = beta1,
      beta2 = beta2,
      n1 = n1,
      n2 = n2,
      de_type = "mean",
      p = p
    )
    pv <- res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = cores, mc.set.seed = TRUE)

pvVec_margin <- unlist(pv_list, use.names = FALSE)
gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-15-10", DataType, "-n1-", n1, "-n2-", n2, "-p-", p,".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot = gg_qqplot(unlist(pvVec_margin)) + ggtitle("Poisson Gamma Model")
print(realistic_qqplot)
dev.off() 

#  ---- ZINB p = 3  ----
args <- commandArgs(trailingOnly = TRUE)
DataType = "zinb"
maxIter = 300 
nCPUS = 10 
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_theta = 1
pi = 0.5 ; pi2 = 0.5
b_pi = 1
de_type = "dispersion"
repID = 77
p = 3
cores = 5
pvVec_margin = NULL 
pv_list <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res = simulateZINB(
      repID = repID, idx = i,
      mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
      theta = theta, b_theta = 1,
      pi = pi, pi2 = pi2,
      n1 = n1, n2 = n2,
      de_type = de_type, p = p
    )
    pv <- res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = cores, mc.set.seed = TRUE)
pvVec_margin <- unlist(pv_list, use.names = FALSE)

gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB (p = 3)")
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot =gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB (p = 3)")
print(realistic_qqplot)
dev.off() 


#  ---- ZINB n1=100,n2=200  ----
args <- commandArgs(trailingOnly = TRUE)
DataType = "zinb"
maxIter = 300 
nCPUS = 10 
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 200
sigma_sq = 5
b_theta = 1
pi = 0.5 ; pi2 = 0.5
b_pi = 1
de_type = "dispersion"
repID = 77
p = 2
cores = 5
pvVec_margin = NULL 
pv_list <- pbmclapply(1:maxIter, function(i) {
  tryCatch({
    res = simulateZINB(
      repID = repID, idx = i,
      mu = mu, mu2 = mu2, sigma_sq = sigma_sq,
      theta = theta, b_theta = 1,
      pi = pi, pi2 = pi2,
      n1 = n1, n2 = n2,
      de_type = de_type, p = p
    )
    pv <- res$pval_1
    if (is.null(pv) || !is.numeric(pv) || !is.finite(pv)) return(NULL)
    as.numeric(pv)
  }, error = function(e) {
    NULL  # skip failures
  })
}, mc.cores = cores, mc.set.seed = TRUE)
pvVec_margin <- unlist(pv_list, use.names = FALSE)

gg_qqplot(unlist(pvVec_margin))
GraphName = paste0("./density_estimation/simulations/logs2025/figures/supplemental/varying_cell_cts_qqplot-", DataType, "-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
realistic_qqplot =gg_qqplot(unlist(pvVec_margin)) + ggtitle("ZINB Model")
print(realistic_qqplot)
dev.off() 

########## Power Analysis ##########
#  ---- GAUSSIAN  MEAN p = 3 ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "gaussian"
de_type = "mean"
maxIter = 300 
nCPUS = 5
tau1 = 2; tau2 = 2
sigma1 = 10; sigma2 = 10
n1 = 100; n2 = 100
p = 3
K = 75
alpha2s = seq(0, 1, by = 0.1)
alpha1 = 0 
repID = 77  

cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

for (idx_alpha in 1:length(alpha2s)) {
  alpha2 <- alpha2s[idx_alpha]
  print(paste0("alpha2=", alpha2))
  
  
  # Run maxIter simulations in parallel.
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateNormalRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = K,idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
print(output$binwidth)
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (alpha2s - alpha1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))

matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
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
  labs(x = expression("Effect Size (" ~ alpha^{(2)} - alpha^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (Gaussian)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- GAUSSIAN VARIANCE p = 3 ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "gaussian"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
sigma1 = 10
sigma2s = seq(10, 14,0.25) 
n1 = 100; n2 = 100
p = 3
K = 75
alpha2 = 0
alpha1 = 0 
repID = 77  
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

for (idx_sigma in 1:length(sigma2s)) {
  sigma2 <- sigma2s[idx_sigma]
  print(paste0("sigma2=", sigma2))
  
  
  # Run maxIter simulations in parallel.
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateNormalRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = K,idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)

df_pw = data.frame(effect_size = (sigma2s - sigma1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:10,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
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
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (Gaussian)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- POISSON GAMMA MEAN p = 3 ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "pg"
de_type = "mean"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 3
repID = 77 

mu1  <- 10
V1   <- 15
params1 <- gamma_params(mu1, V1)
alpha1 <- params1$alpha
beta1  <- params1$beta
target_means <- seq(10, 12, length.out = 15)
target_variance <- 15 # must be > any target mean
if(any(target_means >= target_variance)) {
  stop("Each target mean must be less than the target variance.")
}
alpha2s <- target_means^2 / (target_variance - target_means)
beta2s  <- target_means / (target_variance - target_means)
result_df <- data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 <- alpha2s[idx]
  beta2 <- beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom  <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_means - mu1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))

matrixName = paste0("./density_estimation/simulations/logs2025/figures/varying_cells/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:10,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ mu^{(2)} - mu^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- POISSON GAMMA VARIANCE p = 3----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "pg"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 100
p = 3
repID = 77 
mu1  <- 10
V1   <- 15
params1 <- gamma_params(mu1, V1)
alpha1 <- params1$alpha
beta1  <- params1$beta
target_means <- 10
target_variance <- seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s <- target_means^2 / (target_variance - target_means)
beta2s  <- target_means / (target_variance - target_means)
result_df <- data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 <- alpha2s[idx]
  beta2 <- beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom  <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_variance - V1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw

matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)

df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 
#  ---- ZERO-INFLATED NEGATIVE BINOMIAL VARIANCE  p = 3----
args <- commandArgs(trailingOnly = TRUE)
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 100
sigma_sq = 5
b_thetas = seq(1, 1.70, by =  0.07)
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "dispersion"
repID = 77
K = 100
p = 3
DataType = "zinb"
maxIter = 300 
nCPUS = 8
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(b_thetas))
for (idx_b in 1:length(b_thetas)) {
  b_theta <- b_thetas[idx_b]
  print(paste0("scale factor b_theta =", b_theta, " or (variance) = ", b_theta*theta))
  
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    tryCatch({
      simulateZINB(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                   b_pi = b_pi, pi = pi, pi2 = pi2,
                   repID, idx = i, n1 = n1, n2 = n2, de_type = de_type, K = NULL, p = p)
    }, error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(NULL)
    })
  }
  output1 <- Filter(Negate(is.null), output) # Filter out the NULLs (i.e., the errors)
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    simulateMoMRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                         b_pi = b_pi, pi = pi, pi2 = pi2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID, n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                           b_pi = b_pi, pi = pi, pi2 = pi2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                       b_pi = b_pi, pi = pi, pi2 = pi2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  #find number of numerical issues in our test cases per simulation study
  bad_test_cases[idx_b] = length(output) - length(output1)
  
  # Extract the p-value vectors and bind them as columns.
  if (is.null(pvMat_margin)) {
    pvMat_margin <- list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] <- sapply(output1, function(x) x$pval_1)
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing

Allmethods = c("margin", "mom", "ttest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef <- sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw
matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size Factor (" ~ theta^{(2)}/theta^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)

GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 


#  ---- GAUSSIAN VARIANCE n1 = 100; n2 = 200 ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "gaussian"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
sigma1 = 10
sigma2s = seq(10, 14,0.25) 
n1 = 100; n2 = 200
p = 2
K = 75
alpha2 = 0
alpha1 = 0 
repID = 77  
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

for (idx_sigma in 1:length(sigma2s)) {
  sigma2 <- sigma2s[idx_sigma]
  print(paste0("sigma2=", sigma2))
  
  
  # Run maxIter simulations in parallel.
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateNormalRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = K,idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)

df_pw = data.frame(effect_size = (sigma2s - sigma1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw
matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:9,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
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
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (Gaussian)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- POISSON GAMMA VARIANCE n1 = 100; n2 = 200 ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "pg"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
n1 = 100; n2 = 200
p = 2
repID = 77 
mu1  <- 10
V1   <- 15
params1 <- gamma_params(mu1, V1)
alpha1 <- params1$alpha
beta1  <- params1$beta
target_means <- 10
target_variance <- seq(15, 23, length.out = 10)
if(any(target_means >= target_variance)) { # means must be less than the constant variance
  stop("Each target mean must be less than the target variance.")
}
alpha2s <- target_means^2 / (target_variance - target_means)
beta2s  <- target_means / (target_variance - target_means)
result_df <- data.frame(
  target_mean = target_means,
  target_variance = target_variance,
  alpha2 = alpha2s,
  beta2 = beta2s
)
print(result_df)

pvMat_margin = NULL 
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL

cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
for (idx in 1:(dim(result_df)[1])) {
  alpha2 <- alpha2s[idx]
  beta2 <- beta2s[idx]
  print(paste0("mean =", result_df$target_mean[idx], ", variance = ", result_df$target_variance[idx]))
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulatePGRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                        n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = NULL, idx = i)
  }
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1, beta2 = beta2, 
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom  <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)
q = 0.05 
Allmethods = c("margin", "mom", "ttest", "ks")
df_pw = data.frame(effect_size = (target_variance - V1), margin = colMeans( pvMat_margin < q ),
                   mom = colMeans(pvMat_mom < q),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw

matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)

df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (PG)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

#  ---- ZERO-INFLATED NEGATIVE BINOMIAL VARIANCE n1 = 100; n2 = 200----
args <- commandArgs(trailingOnly = TRUE)
mu = 10; mu2 = 10
theta = 5
n1 = 100; n2 = 200
sigma_sq = 5
b_thetas = seq(1, 1.70, by =  0.07)
b_pi = 1
pi = 0.5 ; pi2 = 0.5
de_type = "dispersion"
repID = 77
K = 100
p = 2
DataType = "zinb"
maxIter = 300 
nCPUS = 8
pvMat_margin = NULL
pvMat_mom = NULL
pvMat_ttest = NULL
pvMat_ks = NULL
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoParallel(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 
bad_test_cases = array(NA, dim = length(b_thetas))
for (idx_b in 1:length(b_thetas)) {
  b_theta <- b_thetas[idx_b]
  print(paste0("scale factor b_theta =", b_theta, " or (variance) = ", b_theta*theta))
  
  
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    tryCatch({
      simulateZINB(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                   b_pi = b_pi, pi = pi, pi2 = pi2,
                   repID, idx = i, n1 = n1, n2 = n2, de_type = de_type, K = NULL, p = p)
    }, error = function(e) {
      message(sprintf("Error in iteration %d: %s", i, e$message))
      return(NULL)
    })
  }
  output1 <- Filter(Negate(is.null), output) # Filter out the NULLs (i.e., the errors)
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    simulateMoMRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                         b_pi = b_pi, pi = pi, pi2 = pi2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID, n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2","truncnorm")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                           b_pi = b_pi, pi = pi, pi2 = pi2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(mu = mu, mu2 = mu2, sigma_sq = sigma_sq, theta = theta, b_theta = b_theta,
                                       b_pi = b_pi, pi = pi, pi2 = pi2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  #find number of numerical issues in our test cases per simulation study
  bad_test_cases[idx_b] = length(output) - length(output1)
  
  # Extract the p-value vectors and bind them as columns.
  if (is.null(pvMat_margin)) {
    pvMat_margin <- list()
  }
  pvMat_margin[[length(pvMat_margin) + 1]] <- sapply(output1, function(x) x$pval_1)
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval))
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl) #stop parallel computing

Allmethods = c("margin", "mom", "ttest", "ks")
init_b_theta = b_thetas[1] #used purely for dataframe
q = 0.05 
margin_sef <- sapply(pvMat_margin, function(vec) mean(vec < q, na.rm = TRUE))
df_pw = data.frame(effect_size = (b_thetas*theta/theta), margin = margin_sef,
                   mom = colMeans(pvMat_mom < q, na.rm = T),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
df_pw
matrixName = paste0("./density_estimation/simulations/supplementary2025/PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:8,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3", "ttest" = "royalblue2", "ks" = "goldenrod1")
obj_pw = ggplot(df_pw2, aes(x = effect_size, y = Power, group = Methods, col=Methods, shape=Methods, linetype = Methods)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.5) + 
  theme_bw(base_size = 16) + theme(legend.position="bottom", legend.title = element_blank())+
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black", size = 1) +  
  geom_hline(yintercept = 0.05, linetype = "twodash", color = "red", size = 1) +  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  ) + 
  labs(x = expression("Effect Size Factor (" ~ theta^{(2)}/theta^{(1)} ~ ")"), 
       y = "Power", 
       title = "Empirical Power (ZINB)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)

GraphName = paste0("./density_estimation/simulations/supplementary2025/newPowerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 



########## Gaussian Sensitivity Analysis ##########
#  ---- Gaussian Variance-Only functions ----
simulateNormalSensitivityTest = function(repID = NULL, de_type = "mean", p = 2,  K = 75,  lower = 300, upper = 1000, n1 = 100, n2 = 100, idx, alpha1, alpha2, tau1, tau2, sigma1, sigma2){
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
  simData = simu_hier_vary_cellcts(seed = seed, lower = lower, upper = upper,alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2)
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

  if (p != 1) {
    X <- subset(X, select = -X_1)
    
    varb <- varb[varb != "X_1"]
    print(head(X))
    print(varb)
  }
  
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
  print(head(df1))
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
  
  # # Model without intercept effects
  # tc_1 = NULL
  # tc_2 = NULL
  # Y1_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  # Y2_suffi = NULL #used for computing the sample mean of raw counts (Tbar)
  # for (dim in 1:p){
  #   tc_1 = cbind(tc_1, tk^dim)
  #   tc_2 = cbind(tc_2, tk^dim)
  #   Y1_suffi = cbind(Y1_suffi, y1^dim) 
  #   Y2_suffi = cbind(Y2_suffi, y2^dim)
  # }  
  # T_c1 = colMeans(Y1_suffi)
  # T_c2 = colMeans(Y2_suffi)
  # X_t1 = t( t(tc_1) - T_c1 ) #centered design matrix based on mean of group 1
  # X_t2 = t( t(tc_2) - T_c2 ) #centered design matrix based on mean of group 2
  # X_t1f = cbind(1, X_t1)
  # X_t2f = cbind(1, X_t2) 
  # colnames(X_t1f) = c("Intercept", varb)  
  # colnames(X_t2f) = c("Intercept", varb) 
  # df1_c = as.data.frame(cbind( S1sum, carrier_scale1, X_t1f))
  # df2_c = as.data.frame(cbind( S2sum, carrier_scale2, X_t2f))
  # colnames(df1_c)[1] = "sum_cts"
  # colnames(df1_c)[2] = "carrier_scale"
  # colnames(df2_c)[1] = "sum_cts" 
  # colnames(df2_c)[2] = "carrier_scale"
  # df1_c <- df1_c[df1_c$carrier_scale > 0, , drop = FALSE]
  # df2_c <- df2_c[df2_c$carrier_scale > 0, , drop = FALSE]
  # formula = as.formula( paste0('sum_cts~offset(log(carrier_scale))+', paste(varb, collapse = '+'))) 
  # 
  # # Fit glm model for each group
  # sef_gp1_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df1_c))  # Some numerical issue
  # sef_gp2_ct = tryCatch(glm( formula, family = poisson(link="log"), data = df2_c))  # Some numerical issue
  # beta_est1_c = as.vector(sef_gp1_ct$coefficients)
  # beta_est2_c = as.vector(sef_gp2_ct$coefficients)
  # beta_diff_c = beta_est1_c - beta_est2_c
  # beta_diff_c
  # sef_df1_c = as.vector( carrier_est * exp(X_t1f %*% beta_est1_c) )
  # sef_df2_c = as.vector( carrier_est * exp(X_t2f %*% beta_est2_c) )
  # 
  # # 3. Inference without intercept effects
  # G_t1 = t(X_t1) %*%  ( sef_df1_c * X_t1)*cellSum1/sum(sef_df1_c)
  # G_t1
  # G_t2 = t(X_t2) %*% ( sef_df2_c * X_t2)*cellSum2/sum(sef_df2_c) 
  # G_t2  
  # Z11c = t(X_t1)%*% (diag(K) - cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  # Z12c = t(X_t1)%*%( -cellSum1*as.vector( exp(X_t1f %*% beta_est1_c))*M/sum(sef_df1_c) )
  # Z21c = t(X_t2)%*%( -cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  # Z22c = t(X_t2)%*%( diag(K) - cellSum2*as.vector( exp(X_t2f %*% beta_est2_c))*M/sum(sef_df2_c) )
  # L1c = solve(G_t1, Z11c) - solve(G_t2, Z21c)
  # L2c = solve(G_t1, Z12c) - solve(G_t2, Z22c)
  # Cov_bar_sub = L1c%*%tmp1%*%t(L1c) + L2c%*%tmp2%*%t(L2c)
  # chi_stat3 = beta_diff_c[-1]%*%solve(Cov_bar_sub, beta_diff_c[-1])
  # pval_3 = 1 - pchisq( chi_stat3, df = p) 
  
  return(list(effect_size = effect_size, de_type = de_type, p = p, K = K, n1 = n1, n2 = n2, binwidth = binwidth, X = X, y_agg = y_agg, y1 = y1, y2 = y2,
              carrier_est = carrier_est, est_midpoints = est_midpoints, carrier_scale1 = carrier_scale1, carrier_scale2 = carrier_scale2,
              sef_df1 = sef_df1, sef_df2 = sef_df2, pval_1 = pval_1, chi_stat1 = chi_stat1, Cov_1 = Cov_bar, 
              beta_diff = beta_diff, carrier_plt = carrier_plot, comparison_plt = plt))
}

#  ---- Gaussian Variance-Only test ----
args <- commandArgs(trailingOnly = TRUE)  
DataType = "gaussian"
de_type = "variance"
maxIter = 300 
nCPUS = 8 
tau1 = 2; tau2 = 2
sigma1 = 10
sigma2s = seq(10, 14,0.25) 
n1 = 100; n2 = 100
p = 2 #use p = 1 or 2 for comparison
K = 75
alpha2 = 0
alpha1 = 0 
repID = 77  
cl = makeCluster(nCPUS, type = 'SOCK' ) # outfile = logfile 
registerDoSNOW(cl)
combine = function(x, ...){
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
} 

pvMat_margin = NULL 
pvMat_mom = NULL 
pvMat_ttest = NULL
pvMat_ks = NULL

for (idx_sigma in 1:length(sigma2s)) {
  sigma2 <- sigma2s[idx_sigma]
  print(paste0("sigma2=", sigma2))
  
  
  # Run maxIter simulations in parallel.
  print("SEF:")
  output <- foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateNormalSensitivityTest(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2, n1 = n1, n2 = n2,repID = repID,de_type = de_type, p = p, K = K,idx = i)
  }
  
  print("MoM:")
  outputMoM = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    simulateMoMRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                         distribution = DataType, de_type = de_type, p = p,repID = repID,n1 = n1, n2 = n2, idx = i)
  }
  
  print("Pseudobulk:")
  outputPB = foreach(i = 1:maxIter, .packages = c("ggplot2")) %dopar% {
    t_testPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                           distribution = DataType,test_type = "ttest", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    ksPB = simulatePseudobulkRealistic(alpha1 = alpha1, alpha2 = alpha2, tau1 = tau1, tau2 = tau2, sigma1 = sigma1, sigma2 = sigma2,
                                       distribution = DataType,test_type = "ks", de_type = de_type, n1 = n1, n2 = n2, idx = i)
    
    list(t_testPB = t_testPB, ksPB = ksPB)
  }
  
  # Extract the p-value vectors and bind them as columns.
  pvMat_margin <- cbind(pvMat_margin, sapply(output, function(x) x$pval_1))
  pvMat_mom    <- cbind(pvMat_mom, sapply(outputMoM, function(x) x$pval)) # can comment out when p = 1
  pvMat_ttest    <- cbind(pvMat_ttest, sapply(outputPB, function(x) x$t_testPB$pval))
  pvMat_ks    <- cbind(pvMat_ks, sapply(outputPB, function(x) x$ksPB$pval))
}
stopCluster(cl)

df_pw = data.frame(effect_size = (sigma2s - sigma1), margin = colMeans( pvMat_margin < q ),
                   #mom = colMeans( pvMat_mom < q ),
                   ttest = colMeans(pvMat_ttest < q),
                   ks = colMeans(pvMat_ks < q))
matrixName = paste0("./density_estimation/simulations/supplementary2025/sensitivity-PowerMatrix-", DataType, "-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".csv" )
# write_csv(df_pw, matrixName)
# df_pw = read.csv(matrixName)
df_pw = df_pw[1:10,]
df_pw2 = reshape::melt( df_pw, id.vars = "effect_size")
colnames(df_pw2) = c("effect_size", "Methods", "Power")
custom_colors = c("margin" = "orangered2", "mom" = "darkolivegreen3","ttest" = "royalblue2", "ks" = "goldenrod1")
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
  labs(x = expression("Effect Size (" ~ sigma^{2 * "," ~ (2)} - sigma^{2 * "," ~ (1)} ~ ")"),
       y = "Power", 
       title = "Empirical Power (Gaussian p = 2)") + 
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) + 
  scale_color_manual(values = custom_colors)
print(obj_pw)
GraphName = paste0("./density_estimation/simulations/supplementary2025/sensitivity-Powerplot-", DataType,"-",de_type ,"-varying-cells-n1-", n1, "-n2-", n2, "-p-", p, ".pdf" )
pdf(GraphName, width = 5, height = 5)
print(obj_pw)
dev.off() 

