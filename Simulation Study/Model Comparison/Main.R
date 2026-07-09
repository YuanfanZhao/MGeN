# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, Mu, Sigma, Kurtosis, Ratio) {
  # 导入需要的package
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(mvtnorm)))
  suppressWarnings(suppressMessages(library(e1071)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressMessages(suppressWarnings(library(statmod)))
  
  # 导入生成随机数函数
  source("GS/rGaN.R")
  source("GS/rIGaN.R")
  source("GS/rIGauN.R")
  source("GS/rRIGauN.R")
  source("GS/rMix.R")
  
  # 导入密度函数
  source("Density/GaN_den.R")
  source("Density/IGaN_den.R")
  source("Density/IGauN_den.R")
  source("Density/RIGauN_den.R")
  source("Density/Mix_den.R")
  source("Density/NPMN_den.R")
  
  # 导入似然函数
  source("Loglikelihood/GaN_log.R")
  source("Loglikelihood/IGaN_log.R")
  source("Loglikelihood/IGauN_log.R")
  source("Loglikelihood/RIGauN_log.R")
  source("Loglikelihood/NPMN_log.R")
  
  # 导入MLE函数
  source("MLE/GaN_mle.R")
  source("MLE/IGaN_mle.R")
  source("MLE/IGauN_mle.R")
  source("MLE/RIGauN_mle.R")
  source("MLE/NPMN_mle.R")
  
  # 评价指标
  KL <- matrix(nrow = number, ncol = 6)
  ISE <- matrix(nrow = number, ncol = 6)
  MSE <- matrix(nrow = number, ncol = 6)
  Rank <- matrix(nrow = number, ncol = 6)
  
  iter <- 1
  while(iter <= number) {
    # 生成混合数据
    X <- rmix(size, Mu, Sigma, Kurtosis, Ratio)
    
    # 区分训练集和测试集
    ind <- sample(1:size, floor(0.8 * size))
    X.train <- X[ind, ]
    X.test <- X[-ind, ]
    Z.train <- Z[ind, ]
    Z.test <- Z[-ind, ]
    
    
    
    # 在训练集上拟合六个模型
    
    # 矩估计确定初值
    m <- colMeans(X.train)
    s <- var(X.train)
    kur <- apply(X.train, 2, kurtosis, type = 1) %>% mean()
    l1 <- (3 / kur) + 2
    l2 <- 3 / kur
    l3 <- (3 + sqrt(12 * kur + 9)) / (2 * kur)
    l4 <- 3 / kur
    
    GaN.fit <- GaN.MLE(X.train, m, s, l1)
    IGaN.fit <- IGaN.MLE(X.train, m, s, l2)
    IGauN.fit <- IGauN.MLE(X.train, m, s, l3)
    RIGauN.fit <- RIGauN.MLE(X.train, m, s, l4)
    
    # 在测试集上测试拟合效果
    true_parameters <- list(
      Mu = Mu, 
      Sigma = Sigma, 
      Kurtosis = Kurtosis, 
      Ratio = Ratio
    )
    
    # 把你的6个模型拟合得到的参数打包
    estimated_parameters <- list(
      GaN = list(Mu = GaN.fit$Mu, Sigma = GaN.fit$Sigma, Lambda = GaN.fit$Lambda),
      IGaN = list(Mu = IGaN.fit$Mu, Sigma = IGaN.fit$Sigma, Lambda = IGaN.fit$Lambda),
      IGauN = list(Mu = IGauN.fit$Mu, Sigma = IGauN.fit$Sigma, Lambda = IGauN.fit$Lambda),
      RIGauN = list(Mu = RIGauN.fit$Mu, Sigma = RIGauN.fit$Sigma, Lambda = RIGauN.fit$Lambda)
    )
    # KL散度和ISE
    evaluation1 <- KLandISE(X.test, Z.test, true_parameters, estimated_parameters)
    
    # MSE
    evaluation2 <- MSE(Z.test, true_parameters$Beta, estimated_parameters)
    
    # Rank
    Ran <- colMeans(matrix(c(rank(-as.numeric(evaluation1$KLD)), rank(-as.numeric(evaluation1$ISE)), rank(-as.numeric(evaluation2$MSE_mu))), nrow = 3))
    
    KL[iter, ] <- evaluation1$KLD
    ISE[iter, ] <- evaluation1$ISE
    MSE[iter, ] <- evaluation2$MSE_mu
    Rank[iter, ] <- Ran
    iter <- iter + 1
  }
  
  return(list(KL = KL, ISE = ISE, MSE = MSE, Rank = Rank))
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 1
size <- 500
alpha <- 2
Beta <- matrix(c(0.5, -1, 0, 0, 2, -0.5), ncol = 2)
Var.tau <- 1
Ratio <- c(0.3, 0.4, 0.1, 0, 0.2)

# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, alpha, Beta, Var.tau, Ratio)
stopCluster(cl)

# ==============================================================================================================================

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)


Simulation(1, size, alpha, Beta, Var.tau, Ratio)
