library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

Simulation <- function(number, size, Mu, Sigma, Lambda) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressWarnings(suppressMessages(library(mvtnorm)))
  suppressWarnings(suppressMessages(library(e1071)))
  suppressWarnings(suppressMessages(library(Bessel)))
  
  rGaN <- function(n, Mu, Sigma, Lambda) {
    
    # 生成τ分布的随机数
    tau <- rgamma(n, Lambda, Lambda - 1)
    
    # 随机表示方法生成多元IGauN分布随机数
    Z <- rmvnorm(n, rep(0, nrow(Sigma)), Sigma)
    X <- matrix(rep(Mu,n), nrow = n, byrow = TRUE) + Z / sqrt(tau)
    return(X)
  }
  
  GaN.density <- function(x, Mu, Sigma, Lambda) {
    d <- length(Mu)
    
    # 计算密度函数
    term1 <- gamma(d / 2 + Lambda) * (Lambda - 1)^Lambda / (2 * pi)^(d / 2) / det(Sigma)^(0.5) / gamma(Lambda)
    term2 <- (as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) / 2 + Lambda - 1)^(-d / 2 - Lambda)
    
    return(term1 * term2)
  }
  
  GaN.loglikeli <- function(X, Mu, Sigma, Lambda) {
    # 返回密度函数的对数和
    return(sum(log(apply(X, 1, GaN.density, Mu = Mu, Sigma = Sigma, Lambda = Lambda))))
  }
  
  GaN.MLE <- function(X, Mu_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-8) {
    index <- FALSE
    iter <- 1
    
    n <- nrow(X)
    d <- ncol(X)
    
    # US算法求c+1/(s-1)+log(s)-ψ(s)=0的零点
    US <- function(s, c) {
      index <- 0
      k <- 0
      
      s_new <- s
      
      while (index == 0) {
        s <- s_new
        q <- c + log(s - 1) + 1 / (s - 1) - digamma(s) + pi^2 * s / 6 - 1 / s
        s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
        
        if (abs(s_new - s) < 1e-5) {
          index <- 1
          break
        }
        k <- k + 1
      }
      return(s_new)
    }
    
    Gamma.moment <- function(X, Mu, Sigma, Lambda) {
      # Gamma分布的参数α，β
      compute.beta <- function(x) {
        return(as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) / 2 + Lambda - 1)
      }
      alpha <- d / 2 + Lambda
      beta <- apply(X, 1, compute.beta)
      
      
      result1 <- alpha / beta
      result2 <- digamma(alpha) - log(beta)
      
      return(list(E1 = result1, Elog = result2))
    }
    
    Loglikeli_new <- GaN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
    
    while (iter <= max_iter) {
      Mu <- Mu_new
      Sigma <- Sigma_new
      Lambda <- Lambda_new
      Loglikeli <- Loglikeli_new
      
      Gamma.result <- Gamma.moment(X, Mu, Sigma, Lambda)
      
      a <- Gamma.result$E1
      Bd <- mean(Gamma.result$Elog)
      Ba <- mean(a)
      A <- diag(a)
      
      Mu_new <- colSums(a * X / sum(a))
      cat("Mu:", Mu_new, "\n")
      
      Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
      Sigma_new <- t(Q) %*% A %*% Q / n
      cat("Sigma:", Sigma_new, "\n")
      
      Lambda_new <- US(Lambda, 1 - Ba + Bd)
      cat("Lambda:", Lambda_new, "\n")
      
      Loglikeli_new <- GaN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
      cat("Loglikelihood:", Loglikeli_new, "\n")
      
      if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
        index <- TRUE
        break
      }
      
      cat("Iteration number:", iter, "\n")
      iter <- iter + 1
    }
    
    result <- list(Mu = Mu_new, 
                   Sigma = Sigma_new, 
                   Lambda = Lambda_new, 
                   Loglikelihood = Loglikeli_new,
                   number = iter, 
                   index = index)
    return(result)
  }
  
  # 开始Simulation
  k <- 1
  
  result_moment <- list(Mu = matrix(nrow = number, ncol = length(Mu)), 
                        Sigma = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                        Lambda = numeric(length = number), 
                        Loglikelihood = numeric(length = number))
  result_fit <- list(Mu = matrix(nrow = number, ncol = length(Mu)), 
                     Sigma = array(dim = c(nrow(Sigma), ncol(Sigma), number)), 
                     Lambda = numeric(length = number), 
                     Loglikelihood = numeric(length = number), 
                     Iteration = numeric(length = number))
  
  while(k <= number) {
    X <- rGaN(size, Mu, Sigma, Lambda)
    
    # 矩估计确定初值
    m <- colMeans(X)
    s <- var(X)
    kur <- apply(X, 2, kurtosis, type = 1) %>% mean()
    l <- (3 / kur) + 2
    
    
    Res <- tryCatch({
      GaN.MLE(X, m, s, l)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result_moment[[1]][k, ] <- m
        result_moment[[2]][, , k] <- s
        result_moment[[3]][k] <- l
        result_moment[[4]][k] <- GaN.loglikeli(X, m, s, l)
        result_fit[[1]][k, ] <- Res$Mu
        result_fit[[2]][, , k] <- Res$Sigma
        result_fit[[3]][k] <- Res$Lambda
        result_fit[[4]][k] <- Res$Loglikelihood
        result_fit[[5]][k] <- Res$number
        k <- k + 1
      }
    }
  }
  
  result <- list(result_fit, result_moment)
  
  return(result)
}

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 105
size <- 900
Mu <- c(-1, 2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
Lambda <- 2.5

# 开始计时
start_time <- Sys.time()

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, Mu, Sigma, Lambda)
stopCluster(cl)

# 汇总拟合结果并保存
Result_Mu <- result_part[[1]][[1]]$Mu
Result_Sigma <- result_part[[1]][[1]]$Sigma
Result_Lambda <- result_part[[1]][[1]]$Lambda
Result_Loglikelihood <- result_part[[1]][[1]]$Loglikelihood
Result_Iteration <- result_part[[1]][[1]]$Iteration
Result_Moment_Mu <- result_part[[1]][[2]]$Mu
Result_Moment_Sigma <- result_part[[1]][[2]]$Sigma
Result_Moment_Lambda <- result_part[[1]][[2]]$Lambda
Result_Moment_Loglikelihood <- result_part[[1]][[2]]$Loglikelihood

for (i in 2:num_cores) {
  Result_Mu <- rbind(Result_Mu, result_part[[i]][[1]]$Mu)
  Result_Sigma <- abind(Result_Sigma, result_part[[i]][[1]]$Sigma)
  Result_Lambda <- c(Result_Lambda, result_part[[i]][[1]]$Lambda)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]][[1]]$Loglikelihood)
  Result_Iteration <- c(Result_Iteration, result_part[[i]][[1]]$Iteration)
  Result_Moment_Mu <- rbind(Result_Moment_Mu, result_part[[i]][[2]]$Mu)
  Result_Moment_Sigma <- abind(Result_Moment_Sigma, result_part[[i]][[2]]$Sigma)
  Result_Moment_Lambda <- c(Result_Moment_Lambda, result_part[[i]][[2]]$Lambda)
  Result_Moment_Loglikelihood <- c(Result_Moment_Loglikelihood, result_part[[i]][[2]]$Loglikelihood)
}

result <- list(Moment.Mu = Result_Moment_Mu, 
               Moment.Sigma = Result_Moment_Sigma, 
               Moment.Lambda = Result_Moment_Lambda, 
               Moment.Loglikelihood = Result_Moment_Loglikelihood, 
               MLE.Mu = Result_Mu, 
               MLE.Sigma = Result_Sigma, 
               MLE.Lambda = Result_Lambda, 
               MLE.Loglikelihood = Result_Loglikelihood, 
               MLE.Iteration = Result_Iteration)

MLE.Mu <- colMeans(result$MLE.Mu)
MLE.Sigma <- apply(result$MLE.Sigma, c(1, 2), mean)
MLE.Lambda <- mean(result$MLE.Lambda)
MLE.Loglikelihood <- mean(result$MLE.Loglikelihood)
MLE.Iteration <- mean(result$MLE.Iteration)

Moment.Mu <- colMeans(result$Moment.Mu)
Moment.Sigma <- apply(result$Moment.Sigma, c(1, 2), mean)
Moment.Lambda <- mean(result$Moment.Lambda)
Moment.Loglikelihood <- mean(result$Moment.Loglikelihood)

MLE.Mu.MSE <- colSums((sweep(result$MLE.Mu, MARGIN = 2, STATS = Mu, FUN = "-"))^2 / all_times)
MLE.Lambda.MSE <- sum((result$MLE.Lambda - Lambda)^2 / all_times)

b <- array(dim = dim(result$MLE.Sigma))
for (i in 1:all_times) {
  b[, , i] <- result$MLE.Sigma[, , i] - Sigma
}
MLE.Sigma.MSE <- apply(b^2, c(1, 2), mean)

Moment.Mu.MSE <- colSums(((sweep(result$Moment.Mu, MARGIN = 2, STATS = Mu, FUN = "-"))^2)^2 / all_times)
Moment.Lambda.MSE <- sum((result$Moment.Lambda - Lambda)^2 / all_times)

b <- array(dim = dim(result$Moment.Sigma))
for (i in 1:all_times) {
  b[, , i] <- result$Moment.Sigma[, , i] - Sigma
}
Moment.Sigma.MSE <- apply(b^2, c(1, 2), mean)



MLE <- data.frame(Loglikelihood = MLE.Loglikelihood, 
                  Mu1 = MLE.Mu[1], 
                  Mu2 = MLE.Mu[2], 
                  Sigma11 = MLE.Sigma[1, 1], 
                  Sigma12 = MLE.Sigma[1, 2], 
                  Sigma22 = MLE.Sigma[2, 2], 
                  Lambda = MLE.Lambda, 
                  Mu1.MSE = MLE.Mu.MSE[1], 
                  Mu2.MSE = MLE.Mu.MSE[2], 
                  Sigma11.MSE = MLE.Sigma.MSE[1, 1], 
                  Sigma12.MSE = MLE.Sigma.MSE[1, 2], 
                  Sigma22.MSE = MLE.Sigma.MSE[2, 2], 
                  Lambda.MSE = MLE.Lambda.MSE, 
                  Iteration = MLE.Iteration)
Moment <- data.frame(Loglikelihood = Moment.Loglikelihood, 
                     Mu1 = Moment.Mu[1], 
                     Mu2 = Moment.Mu[2], 
                     Sigma11 = Moment.Sigma[1, 1], 
                     Sigma12 = Moment.Sigma[1, 2], 
                     Sigma22 = Moment.Sigma[2, 2], 
                     Lambda = Moment.Lambda, 
                     Mu1.MSE = Moment.Mu.MSE[1], 
                     Mu2.MSE = Moment.Mu.MSE[2], 
                     Sigma11.MSE = Moment.Sigma.MSE[1, 1], 
                     Sigma12.MSE = Moment.Sigma.MSE[1, 2], 
                     Sigma22.MSE = Moment.Sigma.MSE[2, 2], 
                     Lambda.MSE = Moment.Lambda.MSE)
Final <- list(MLE = MLE %>% round(digits = 4), 
              Moment = Moment %>% round(digits = 4))

Final

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)