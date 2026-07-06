library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

Simulation <- function(number, size, Beta, Sigma, Lambda) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(statmod)))
  suppressWarnings(suppressMessages(library(mvtnorm)))
  suppressWarnings(suppressMessages(library(e1071)))
  suppressWarnings(suppressMessages(library(Bessel)))
  
  rRIGauN <- function(n, Z, Beta, Sigma, Lambda) {
    
    # 位置参数μ（n×d）
    Mu <- Z %*% Beta
    
    # 生成τ分布的随机数
    tau <- 1 / rinvgauss(n, mean = 1, shape = Lambda)
    
    # 随机表示方法生成多元IGauN分布随机数
    Z <- rmvnorm(n, rep(0, nrow(Sigma)), Sigma)
    X <- Mu + Z / sqrt(tau)
    return(X)
  }
  
  RIGauN.density <- function(x, z, Beta, Sigma, Lambda) {
    d <- ncol(Beta)
    
    # 位置参数Mu
    Mu <- t(Beta) %*% z
    
    # GIGau分布的参数a，b，p
    a <- as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) + Lambda
    b <- Lambda
    p <- (d + 1) / 2
    
    # 计算密度函数
    term1 <- exp(Lambda) * sqrt(Lambda) / (2 * pi)^((d + 1) / 2) / det(Sigma)^(0.5)
    term2 <- 2 * besselK(sqrt(a * b), p, expon.scaled = FALSE) / (a / b)^(p / 2)
    
    return(term1 * term2)
  }
  
  RIGauN.loglikeli <- function(X, Z, Beta, Sigma, Lambda) {
    d <- ncol(Beta)
    
    compute.density <- function(data) {
      x <- data[1:d]
      z <- data[-(1:d)]
      RIGauN.density(x, z, Beta, Sigma, Lambda)
    }
    
    # 返回密度函数的对数和
    return(sum(log(apply(cbind(X, Z), 1, compute.density))))
  }
  
  RIGauN.MLE <- function(X, Z, Beta_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-8) {
    index <- FALSE
    iter <- 1
    
    n <- nrow(X)
    d <- ncol(X)
    
    GIG.moment <- function(X, Z, Beta, Sigma, Lambda) {
      
      # 位置参数Mu
      Mu <- Z %*% Beta
      
      # GIGau分布的参数a，b，p
      compute.a <- function(data) {
        x <- data[1:d]
        mu <- data[-(1:d)]
        return(as.numeric(t(x - mu) %*% solve(Sigma) %*% (x - mu)) + Lambda)
      }
      a <- apply(cbind(X, Mu), 1, compute.a)
      b <- Lambda
      p <- (d + 1) / 2
      
      result1 <- BesselK(sqrt(a * b), p + 1) / BesselK(sqrt(a * b), p) * (b / a)^0.5
      result2 <- BesselK(sqrt(a * b), p - 1) / BesselK(sqrt(a * b), p) * (b / a)^(-0.5)
      
      return(list(E1 = result1, Ei1 = result2))
    }
    
    Loglikeli_new <- RIGauN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
    
    while (iter <= max_iter) {
      Beta <- Beta_new
      Sigma <- Sigma_new
      Lambda <- Lambda_new
      Loglikeli <- Loglikeli_new
      
      GIG.result <- GIG.moment(X, Z, Beta, Sigma, Lambda)
      
      a <- GIG.result$E1
      Bb <- mean(GIG.result$Ei1)
      Ba <- mean(a)
      A <- diag(a)
      
      Beta_new <- solve(t(Z) %*% A %*% Z) %*% (t(Z) %*% A %*% X)
      cat("Beta:", Beta_new, "\n")
      
      Q <- X - Z %*% Beta_new
      Sigma_new <- t(Q) %*% A %*% Q / n
      cat("Sigma:", Sigma_new, "\n")
      
      Lambda_new <- 1 / (Ba + Bb - 2)
      cat("Lambda:", Lambda_new, "\n")
      
      Loglikeli_new <- RIGauN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
      cat("Loglikelihood:", Loglikeli_new, "\n")
      
      if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
        index <- TRUE
        break
      }
      
      cat("Iteration number:", iter, "\n")
      iter <- iter + 1
    }
    
    result <- list(Beta = Beta_new, 
                   Sigma = Sigma_new, 
                   Lambda = Lambda_new, 
                   Loglikelihood = Loglikeli_new,
                   number = iter, 
                   index = index)
    return(result)
  }
  
  # 开始Simulation
  k <- 1
  
  result_fit <- list(Beta = array(dim = c(nrow(Beta), ncol(Beta), number)), 
                     Sigma = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                     Lambda = numeric(length = number), 
                     Loglikelihood = numeric(length = number), 
                     Iteration = numeric(length = number))
  
  while(k <= number) {
    Z1 <- rep(1, size)
    Z2 <- rnorm(size, 0, 1)
    Z3 <- rpois(size, 2)
    Z <- cbind(Z1, Z2, Z3)
    # Z <- matrix(rep(1, size), ncol = 1)
    X <- rRIGauN(size, Z, Beta, Sigma, Lambda)
    
    
    Res <- tryCatch({
      RIGauN.MLE(X, Z, Beta, Sigma, Lambda)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result_fit[[1]][, , k] <- Res$Beta
        result_fit[[2]][, , k] <- Res$Sigma
        result_fit[[3]][k] <- Res$Lambda
        result_fit[[4]][k] <- Res$Loglikelihood
        result_fit[[5]][k] <- Res$number
        k <- k + 1
      }
    }
  }
  
  result <- list(result_fit)
  
  return(result)
}

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 105
size <- 10000
Beta <- matrix(c(2, -1, 0, 0, 0.5, -3), ncol = 2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
Lambda <- 2

# 开始计时
start_time <- Sys.time()

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, Beta, Sigma, Lambda)
stopCluster(cl)

# 汇总拟合结果并保存
Result_Beta <- result_part[[1]][[1]]$Beta
Result_Sigma <- result_part[[1]][[1]]$Sigma
Result_Lambda <- result_part[[1]][[1]]$Lambda
Result_Loglikelihood <- result_part[[1]][[1]]$Loglikelihood
Result_Iteration <- result_part[[1]][[1]]$Iteration

for (i in 2:num_cores) {
  Result_Beta <- abind(Result_Beta, result_part[[i]][[1]]$Beta)
  Result_Sigma <- abind(Result_Sigma, result_part[[i]][[1]]$Sigma)
  Result_Lambda <- c(Result_Lambda, result_part[[i]][[1]]$Lambda)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]][[1]]$Loglikelihood)
  Result_Iteration <- c(Result_Iteration, result_part[[i]][[1]]$Iteration)
}

result <- list(MLE.Beta = Result_Beta, 
               MLE.Sigma = Result_Sigma, 
               MLE.Lambda = Result_Lambda, 
               MLE.Loglikelihood = Result_Loglikelihood, 
               MLE.Iteration = Result_Iteration)

MLE.Beta <- apply(result$MLE.Beta, c(1, 2), mean)
MLE.Sigma <- apply(result$MLE.Sigma, c(1, 2), mean)
MLE.Lambda <- mean(result$MLE.Lambda)
MLE.Loglikelihood <- mean(result$MLE.Loglikelihood)
MLE.Iteration <- mean(result$MLE.Iteration)

MLE.Beta.Bias <- MLE.Beta - Beta
MLE.Sigma.Bias <- MLE.Sigma - Sigma
MLE.Lambda.Bias <- MLE.Lambda - Lambda

b <- array(dim = dim(result$MLE.Beta))
for (i in 1:all_times) {
  b[, , i] <- result$MLE.Beta[, , i] - Beta
}
MLE.Beta.MSE <- apply(b^2, c(1, 2), mean)

b <- array(dim = dim(result$MLE.Sigma))
for (i in 1:all_times) {
  b[, , i] <- result$MLE.Sigma[, , i] - Sigma
}
MLE.Sigma.MSE <- apply(b^2, c(1, 2), mean)

MLE.Lambda.MSE <- sum((result$MLE.Lambda - Lambda)^2 / all_times)

MLE <- data.frame(Loglikelihood = MLE.Loglikelihood, 
                  Beta11 = MLE.Beta[1, 1], 
                  Beta12 = MLE.Beta[1, 2], 
                  Beta21 = MLE.Beta[2, 1], 
                  Beta22 = MLE.Beta[2, 2], 
                  Beta31 = MLE.Beta[3, 1], 
                  Beta32 = MLE.Beta[3, 2], 
                  Sigma11 = MLE.Sigma[1, 1], 
                  Sigma12 = MLE.Sigma[1, 2], 
                  Sigma22 = MLE.Sigma[2, 2], 
                  Lambda = MLE.Lambda, 
                  Beta11.MSE = MLE.Beta.MSE[1, 1], 
                  Beta12.MSE = MLE.Beta.MSE[1, 2], 
                  Beta21.MSE = MLE.Beta.MSE[2, 1], 
                  Beta22.MSE = MLE.Beta.MSE[2, 2], 
                  Beta31.MSE = MLE.Beta.MSE[3, 1], 
                  Beta32.MSE = MLE.Beta.MSE[3, 2], 
                  Sigma11.MSE = MLE.Sigma.MSE[1, 1], 
                  Sigma12.MSE = MLE.Sigma.MSE[1, 2], 
                  Sigma22.MSE = MLE.Sigma.MSE[2, 2], 
                  Lambda.MSE = MLE.Lambda.MSE, 
                  Iteration = MLE.Iteration)

MLE |> round(digits = 4)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)