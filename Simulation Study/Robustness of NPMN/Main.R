# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, Mu, Sigma, nv, Ratio) {
  # 导入需要的package
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(mvtnorm)))
  suppressMessages(suppressWarnings(library(MASS)))
  suppressMessages(suppressWarnings(library(QRM)))
  suppressMessages(suppressWarnings(library(LaplacesDemon)))
  suppressWarnings(suppressMessages(library(e1071)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressMessages(suppressWarnings(library(statmod)))
  
  # 导入生成随机数函数
  source("GS/rMix.R")
  
  # 导入密度函数
  source("Density/GaN_den.R")
  source("Density/IGaN_den.R")
  source("Density/IGauN_den.R")
  source("Density/RIGauN_den.R")
  source("Density/Norm_den.R")
  source("Density/t_den.R")
  source("Density/Laplace_den.R")
  source("Density/Mix_den.R")
  source("Density/NPMN_den.R")
  
  # 导入似然函数
  source("Loglikelihood/GaN_log.R")
  source("Loglikelihood/IGaN_log.R")
  source("Loglikelihood/IGauN_log.R")
  source("Loglikelihood/RIGauN_log.R")
  source("Loglikelihood/Norm_log.R")
  source("Loglikelihood/t_log.R")
  source("Loglikelihood/Laplace_log.R")
  source("Loglikelihood/NPMN_log.R")
  
  # 导入MLE函数
  source("MLE/GaN_mle.R")
  source("MLE/IGaN_mle.R")
  source("MLE/IGauN_mle.R")
  source("MLE/RIGauN_mle.R")
  source("MLE/Norm_mle.R")
  source("MLE/t_mle.R")
  source("MLE/Laplace_mle.R")
  source("MLE/NPMN_mle.R")
  
  # 导入评价指标函数
  source("Criterion/KLandISE.R")
  
  # 导入辅助函数
  source("Auxiliary/Choose_M.R")
  
  # 评价指标
  KL.value <- matrix(nrow = number, ncol = 8)
  ISE.value <- matrix(nrow = number, ncol = 8)
  Rank.value <- array(dim = c(2, 8, number))
  
  iter <- 1
  while(iter <= number) {
    # 生成混合数据
    X <- rMix(size, Mu, Sigma, nv, Ratio)
    
    # 区分训练集和测试集
    ind <- sample(1:size, floor(0.8 * size))
    X.train <- X[ind, ]
    X.test <- X[-ind, ]
    
    # 在训练集上拟合五个模型
    
    # 矩估计确定初值
    m <- colMeans(X.train)
    s <- var(X.train)
    kur <- apply(X.train, 2, kurtosis, type = 1) %>% mean()
    l1 <- (3 / kur) + 2
    l2 <- 3 / kur
    l3 <- (3 + sqrt(12 * kur + 9)) / (2 * kur)
    l4 <- 3 / kur
    
    GaN.fit <- tryCatch({
      GaN.MLE(X.train, m, s, l1)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    IGaN.fit <- tryCatch({
      IGaN.MLE(X.train, m, s, l2)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    IGauN.fit <- tryCatch({
      IGauN.MLE(X.train, m, s, l3)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    RIGauN.fit <- tryCatch({
      RIGauN.MLE(X.train, m, s, l4)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    Norm.fit <- tryCatch({
      Norm.MLE(X.train)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    t.fit <- tryCatch({
      t.MLE(X.train)
    }, error = function(e) {
      NA
    })
    
    Laplace.fit <- tryCatch({
      Laplace.MLE(X.train)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    M.opt <- tryCatch({
      Choose.M(X.train, Mu, Sigma, nv, Ratio)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    NPMN.fit <- tryCatch({
      NPMN.MLE(X.train, min(M.opt))
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(GaN.fit)) & !all(is.na(IGaN.fit)) & !all(is.na(IGauN.fit)) & !all(is.na(RIGauN.fit)) & !all(is.na(NPMN.fit))) {
      if (GaN.fit$index == TRUE & IGaN.fit$index == TRUE & IGauN.fit$index == TRUE & RIGauN.fit$index == TRUE & NPMN.fit$index == TRUE) {
        # 在测试集上测试拟合效果
        true_parameters <- list(
          Mu = Mu,
          Sigma = Sigma, 
          nv = nv, 
          Ratio = Ratio
        )
        
        # 把你的8个模型拟合得到的参数打包
        estimated_parameters <- list(
          GaN = list(Mu = GaN.fit$Mu, Sigma = GaN.fit$Sigma, Lambda = GaN.fit$Lambda),
          IGaN = list(Mu = IGaN.fit$Mu, Sigma = IGaN.fit$Sigma, Lambda = IGaN.fit$Lambda),
          IGauN = list(Mu = IGauN.fit$Mu, Sigma = IGauN.fit$Sigma, Lambda = IGauN.fit$Lambda),
          RIGauN = list(Mu = RIGauN.fit$Mu, Sigma = RIGauN.fit$Sigma, Lambda = RIGauN.fit$Lambda), 
          Norm = list(Mu = Norm.fit$Mu, Sigma = Norm.fit$Sigma),
          t = list(Mu = t.fit$Mu, Sigma = t.fit$Sigma, nv = t.fit$nv),
          Laplace = list(Mu = Laplace.fit$Mu, Sigma = Laplace.fit$Sigma), 
          NPMN = list(Mu = NPMN.fit$Mu, Sigma = NPMN.fit$Sigma, p = NPMN.fit$p)
        )
        # KL散度和ISE
        evaluation <- KLandISE(X.test, true_parameters, estimated_parameters)
        
        # Rank
        Ran <- matrix(c(rank(-as.numeric(evaluation$KLD)), rank(-as.numeric(evaluation$ISE))), byrow = TRUE, nrow = 2)
        
        KL.value[iter, ] <- evaluation$KLD
        ISE.value[iter, ] <- evaluation$ISE
        Rank.value[, , iter] <- Ran
        iter <- iter + 1
      }
    }
  }
  
  return(list(KL = KL.value, ISE = ISE.value, Rank = Rank.value))
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 11
size <- 500
Mu <- c(-1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
nv <- 12
# Ratio <- c(1 / 3, 1 / 3, 1 / 3)
Ratio <- c(5 / 6, 1 / 12, 1 / 12)
# Ratio <- c(1 / 6, 2 / 3, 1 / 6)
# Ratio <- c(1 / 6, 1 / 6, 2 / 3)

# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, Mu, Sigma, nv, Ratio)
stopCluster(cl)

# ==============================================================================================================================

KL.value <- result_part[[1]]$KL
ISE.value <- result_part[[1]]$ISE
Rank.value <- result_part[[1]]$Rank

for (i in 2:num_cores) {
  KL.value <- rbind(KL.value, result_part[[i]]$KL)
  ISE.value <- rbind(ISE.value, result_part[[i]]$ISE)
  Rank.value <- abind(Rank.value, result_part[[i]]$Rank, along = 3)
}

colMeans(KL.value) %>% round(digits = 4)
colMeans(ISE.value) %>% round(digits = 4)
apply(Rank.value, c(1, 2), mean) %>% round(digits = 4)
apply(Rank.value, c(1, 2), mean) %>% colMeans() %>% round(digits = 4)

result <- data.frame(KL = colMeans(KL.value), 
                     ISE = colMeans(ISE.value), 
                     Rank1 = apply(Rank.value, c(1, 2), mean)[1, ], 
                     Rank2 = apply(Rank.value, c(1, 2), mean)[2, ])

write.csv(result, "t.csv")

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)