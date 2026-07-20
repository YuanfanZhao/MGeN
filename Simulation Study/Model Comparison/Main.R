# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, Beta, Sigma, Kurtosis, Ratio) {
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
  
  # 导入似然函数
  source("Loglikelihood/GaN_log.R")
  source("Loglikelihood/IGaN_log.R")
  source("Loglikelihood/IGauN_log.R")
  source("Loglikelihood/RIGauN_log.R")
  
  # 导入MLE函数
  source("MLE/GaN_mle.R")
  source("MLE/IGaN_mle.R")
  source("MLE/IGauN_mle.R")
  source("MLE/RIGauN_mle.R")
  
  # 导入评价指标函数
  source("Criterion/KLandISE.R")
  source("Criterion/MSE.R")
  
  # 评价指标
  KL.value <- matrix(nrow = number, ncol = 4)
  ISE.value <- matrix(nrow = number, ncol = 4)
  MSE.value <- matrix(nrow = number, ncol = 4)
  Rank.value <- array(dim = c(3, 4, number))
  
  iter <- 1
  while(iter <= number) {
    # 生成混合数据
    # Z <- cbind(rep(1, size))
    Z <- cbind(rep(1, size), rnorm(size), rpois(size, 2))
    X <- rMix(size, Z, Beta, Sigma, Kurtosis, Ratio)
    
    # 区分训练集和测试集
    ind <- sample(1:size, floor(0.8 * size))
    X.train <- X[ind, ]
    X.test <- X[-ind, ]
    Z.train <- Z[ind, ]
    Z.test <- Z[-ind, ]
    
    
    
    # 在训练集上拟合四个模型
    
    # 矩估计确定初值
    s <- var(X.train)
    kur <- apply(X.train, 2, kurtosis, type = 1) %>% mean()
    l1 <- (3 / kur) + 2
    l2 <- 3 / kur
    l3 <- (3 + sqrt(12 * kur + 9)) / (2 * kur)
    l4 <- 3 / kur
    
    GaN.fit <- tryCatch({
      GaN.MLE(X.train, Z.train, Beta, s, l1)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    IGaN.fit <- tryCatch({
      IGaN.MLE(X.train, Z.train, Beta, s, l2)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    IGauN.fit <- tryCatch({
      IGauN.MLE(X.train, Z.train, Beta, s, l3)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    RIGauN.fit <- tryCatch({
      RIGauN.MLE(X.train, Z.train, Beta, s, l4)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(GaN.fit)) & !all(is.na(IGaN.fit)) & !all(is.na(IGauN.fit)) & !all(is.na(RIGauN.fit))) {
      if (GaN.fit$index == TRUE & IGaN.fit$index == TRUE & IGauN.fit$index == TRUE & RIGauN.fit$index == TRUE) {
        # 在测试集上测试拟合效果
        true_parameters <- list(
          Beta = Beta,
          Sigma = Sigma, 
          Kurtosis = Kurtosis, 
          Ratio = Ratio
        )
        
        # 把你的4个模型拟合得到的参数打包
        estimated_parameters <- list(
          GaN = list(Beta = GaN.fit$Beta, Sigma = GaN.fit$Sigma, Lambda = GaN.fit$Lambda),
          IGaN = list(Beta = IGaN.fit$Beta, Sigma = IGaN.fit$Sigma, Lambda = IGaN.fit$Lambda),
          IGauN = list(Beta = IGauN.fit$Beta, Sigma = IGauN.fit$Sigma, Lambda = IGauN.fit$Lambda),
          RIGauN = list(Beta = RIGauN.fit$Beta, Sigma = RIGauN.fit$Sigma, Lambda = RIGauN.fit$Lambda)
        )
        # KL散度和ISE
        evaluation1 <- KLandISE(X.test, Z.test, true_parameters, estimated_parameters)
        
        # MSE
        evaluation2 <- MSE(Z.test, true_parameters$Beta, estimated_parameters, X.test)
        
        # Rank
        Ran <- matrix(c(rank(-as.numeric(evaluation1$KLD)), rank(-as.numeric(evaluation1$ISE)), rank(-as.numeric(evaluation2$MSE_pred))), byrow = TRUE, nrow = 3)
        
        KL.value[iter, ] <- evaluation1$KLD
        ISE.value[iter, ] <- evaluation1$ISE
        MSE.value[iter, ] <- evaluation2$MSE_pred
        Rank.value[, , iter] <- Ran
        iter <- iter + 1
      }
    }
  }
  
  return(list(KL = KL.value, ISE = ISE.value, MSE = MSE.value, Rank = Rank.value))
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 105
size <- 500
# Beta <- matrix(c(-1, 2), ncol = 2)
Beta <- matrix(c(1, -2.5, 0, 0, 1.5, -2), ncol = 2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
Kurtosis <- 2
# Ratio <- c(0.25, 0.25, 0.25, 0.25)
Ratio <- c(0.1, 0.1, 0.1, 0.7)

# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, Beta, Sigma, Kurtosis, Ratio)
stopCluster(cl)

# ==============================================================================================================================

KL.value <- result_part[[1]]$KL
ISE.value <- result_part[[1]]$ISE
MSE.value <- result_part[[1]]$MSE
Rank.value <- result_part[[1]]$Rank

for (i in 2:num_cores) {
  KL.value <- rbind(KL.value, result_part[[i]]$KL)
  ISE.value <- rbind(ISE.value, result_part[[i]]$ISE)
  MSE.value <- rbind(MSE.value, result_part[[i]]$MSE)
  Rank.value <- abind(Rank.value, result_part[[i]]$Rank, along = 3)
}

colMeans(KL.value) %>% round(digits = 4)
colMeans(ISE.value) %>% round(digits = 4)
colMeans(MSE.value) %>% round(digits = 4)
apply(Rank.value, c(1, 2), mean) %>% round(digits = 4)
apply(Rank.value, c(1, 2), mean) %>% colMeans() %>% round(digits = 4)

result <- data.frame(KL = colMeans(KL.value), 
                     ISE = colMeans(ISE.value), 
                     MSE = colMeans(MSE.value), 
                     Rank1 = apply(Rank.value, c(1, 2), mean)[1, ], 
                     Rank2 = apply(Rank.value, c(1, 2), mean)[2, ], 
                     Rank3 = apply(Rank.value, c(1, 2), mean)[3, ])

write.csv(result, "C4_RIGauN.csv")

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)


# Simulation(1, size, Beta, Sigma, Kurtosis, Ratio)
