library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

Simulation <- function(number, size, Mu, Sigma, Lambda, Prior = "Non-information", B = 2500, G = 10000) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(mvtnorm)))
  suppressWarnings(suppressMessages(library(e1071)))
  suppressWarnings(suppressMessages(library(Bessel)))
  suppressWarnings(suppressMessages(library(MCMCpack)))
  suppressWarnings(suppressMessages(library(abind)))
  suppressWarnings(suppressMessages(library(GIGrvg)))
  suppressWarnings(suppressMessages(library(statmod)))
  
  # 导入生成随机数函数、密度函数、似然函数以及MLE函数
  source("GS/rRIGauN.R")
  source("Density/RIGauN_den.R")
  source("Loglikelihood/RIGauN_log.R")
  source("MLE/RIGauN_mle.R")
  
  # 导入ECM算法、DA算法
  source("ECM/RIGauN_ecm.R")
  source("DA/RIGauN_da.R")
  
  # 开始Simulation
  k <- 1
  
  result <- list(Mu.Mode = matrix(nrow = number, ncol = length(Mu)), 
                 Sigma.Mode = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                 Mu.Mean = matrix(nrow = number, ncol = length(Mu)), 
                 Mu.std = matrix(nrow = number, ncol = length(Mu)), 
                 Mu.BCI = array(dim = c(length(Mu), 2, number)),
                 Sigma.Mean = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                 Sigma.std = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                 Sigma.BCI = array(dim = c(nrow(Sigma), ncol(Sigma), 2, number)))
  
  while (k <= number) {
    X <- rRIGauN(size, Mu, Sigma, Lambda)
    
    # 矩估计确定初值
    m <- colMeans(X)
    s <- var(X)
    kur <- apply(X, 2, kurtosis, type = 1) %>% mean()
    l <- 3 / kur
    
    # 设置先验信息
    if (Prior == "Non-information") {
      Mu0 = NULL; Kappa0 = NULL; Lambda0 = NULL; Nu0 = NULL
    } else if (Prior == "Conjugate") {
      Mu0 <- m
      Kappa0 <- 0.1
      Lambda0 <- s
      Nu0 <- ncol(X) + 2
    }
    
    
    Res <- tryCatch({
      RIGauN.MLE(X, m, s, l)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    if (!all(is.na(Res))) {
      if (Res$index == TRUE) {
        Mode.fit <- tryCatch({
          RIGauN.ECM(X, m, s, Res$Lambda, Mu0, Kappa0, Lambda0, Nu0)
        }, error = function(e) {
          NA
        }, warning = function(w) {
          NA
        })
        
        if (!all(is.na(Mode.fit))) {
          if (Mode.fit$index == TRUE) {
            Mean.fit <- RIGauN.DA(X, m, s, Res$Lambda, Mu0, Kappa0, Lambda0, Nu0, B, G)
            
            # 计算Beyesian可信区间
            indexx <- floor(c(0.025 * (G - B), 0.975 * (G - B)))
            Mu.sort <- apply(Mean.fit$Mu, 2, sort)
            Mu.L <- (Mu.sort[indexx[1], ] + Mu.sort[indexx[1] + 1, ]) / 2
            Mu.U <- (Mu.sort[indexx[2], ] + Mu.sort[indexx[2] + 1, ]) / 2
            Sigma.sort <- apply(Mean.fit$Sigma, c(1, 2), sort)
            Sigma.L <- (Sigma.sort[indexx[1], , ] + Sigma.sort[indexx[1] + 1, , ]) / 2
            Sigma.U <- (Sigma.sort[indexx[2], , ] + Sigma.sort[indexx[2] + 1, , ]) / 2
            
            # 储存结果
            result[[1]][k, ] <- Mode.fit$Mu
            result[[2]][, , k] <- Mode.fit$Sigma
            result[[3]][k, ] <- colMeans(Mean.fit$Mu)
            result[[4]][k, ] <- apply(Mean.fit$Mu, 2, var) |> sqrt()
            result[[5]][, , k] <- cbind(Mu.L, Mu.U)
            result[[6]][, , ] <- apply(Mean.fit$Sigma, c(1, 2), mean)
            result[[7]][, , k] <- apply(Mean.fit$Sigma, c(1, 2), var) |> sqrt()
            result[[8]][, , , k] <- abind(Sigma.L, Sigma.U, along = 3)
          }
        }
      }
    }
    k <- k + 1
  }
  
  return(result)
}

# 参数设置
num_cores <- detectCores()
all_times <- num_cores * 1
size <- 500
Mu <- c(-1, 2)
Sigma <- matrix(c(1, -0.5, -0.5, 2), nrow = 2)
Lambda <- 0.5

# 测试
# Simulation(1, size, Mu, Sigma, Lambda, "Conjugate", B = 200, G = 1000)

# 开始计时
start_time <- Sys.time()

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, Mu, Sigma, Lambda)
stopCluster(cl)

# 汇总并行运算结果
Mu.Mode <- result_part[[1]]$Mu.Mode
Sigma.Mode <- result_part[[1]]$Sigma.Mode
Mu.Mean <- result_part[[1]]$Mu.Mean
Mu.std <- result_part[[1]]$Mu.std
Mu.BCI <- result_part[[1]]$Mu.BCI
Sigma.Mean <- result_part[[1]]$Sigma.Mean
Sigma.std <- result_part[[1]]$Sigma.std
Sigma.BCI <- result_part[[1]]$Sigma.BCI

for (i in 2:num_cores) {
  Mu.Mode <- rbind(Mu.Mode, result_part[[i]]$Mu.Mode)
  Sigma.Mode <- abind(Sigma.Mode, result_part[[i]]$Sigma.Mode)
  Mu.Mean <- rbind(Mu.Mean, result_part[[i]]$Mu.Mean)
  Mu.std <- rbind(Mu.std, result_part[[i]]$Mu.std)
  Mu.BCI <- abind(Mu.BCI, result_part[[i]]$Mu.BCI)
  Sigma.Mean <- abind(Sigma.Mean, result_part[[i]]$Sigma.Mean)
  Sigma.std <- abind(Sigma.std, result_part[[i]]$Sigma.std)
  Sigma.BCI <- abind(Sigma.BCI, result_part[[i]]$Sigma.BCI)
}

# 输出结果

# 后验众数
Mu.Mode %>% colMeans() %>% round(digits = 4)
apply(Sigma.Mode, c(1, 2), mean) %>% round(digits = 4)

# μ的后验均值、标准差、置信区间和覆盖率
Mu.Mean %>% colMeans() %>% round(digits = 4)
Mu.std %>% colMeans() %>% round(digits = 4)
apply(Mu.BCI, c(1, 2), mean) %>% round(digits = 4)
rowMeans(Mu.BCI[, 1, ] <= Mu & Mu.BCI[, 2, ] >= Mu) %>% round(digits = 4)

# Σ的后验均值、标准差、置信区间和覆盖率
apply(Sigma.Mean, c(1, 2), mean) %>% round(digits = 4)
apply(Sigma.std, c(1, 2), mean) %>% round(digits = 4)
apply(Sigma.BCI, c(1, 2, 3), mean) %>% round(digits = 4)

index_BCI <- array(dim = c(nrow(Sigma), ncol(Sigma), all_times))
for (i in 1:all_times) {
  index_BCI[, , i] <- Sigma.BCI[, , 1, i] <= Sigma & Sigma.BCI[, , 2, i] >= Sigma
}
apply(index_BCI, c(1, 2), mean)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)