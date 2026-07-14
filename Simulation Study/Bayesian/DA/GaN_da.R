GaN.DA <- function(X, Mu, Sigma, Lambda, B = 2500, G = 10000, 
                   Mu0 = NULL, Kappa0 = NULL, Lambda0 = NULL, Nu0 = NULL) {
  
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  # 存储结果的空数据集
  Mu.value = matrix(nrow = G + 1, ncol = d)
  Sigma.value = array(dim = c(d, d, G + 1))
  
  # 赋予初值
  Mu.value[1, ] <- Mu
  Sigma.value[, , 1] <- Sigma
  
  while (iter <= G) {
    # I-step: 生成潜在变量
    Beta <- apply(X, 1, function(x) as.numeric(t(x - Mu.value[iter, ]) %*% solve(Sigma.value[, , iter]) %*% (x - Mu.value[iter, ])) / 2 + Lambda - 1)
    tau <- sapply(Beta, function(beta) rgamma(1, Lambda + d / 2, beta))
    
    # P-step: 更新参数
    if (is.null(Mu0)) {
      # 生成μ随机数
      m <- colSums(tau * X / sum(tau))
      s <- Sigma.value[, , iter] / sum(tau)
      Mu.value[iter + 1, ] <- rmvnorm(1, m, s)
      
      # 生成Σ随机数
      Mu_new <- as.numeric(Mu.value[iter + 1, ])
      Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
      Psi <- t(Q) %*% diag(tau) %*% Q
      Sigma.value[, , iter + 1] <- riwish(n, Psi)
    } else {
      # 生成μ随机数
      m <- (colSums(tau * X) + Kappa0 * Mu0) / (sum(tau) + Kappa0)
      s <- Sigma.value[, , iter] / (sum(tau) + Kappa0)
      Mu.value[iter + 1, ] <- rmvnorm(1, m, s)
      
      # 生成Σ随机数
      Mu_new <- as.numeric(Mu.value[iter + 1, ])
      Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
      Psi <- t(Q) %*% diag(tau) %*% Q + Lambda0 + Kappa0 * (Mu_new - Mu0) %*% t(Mu_new - Mu0)
      Sigma.value[, , iter + 1] <- riwish(n + Nu0 + 1, Psi)
    }
    # cat("Iteration number:", iter, "\n")
    iter <- iter + 1
  }
  
  # 删去前B次迭代数据，保留后G-B次作为后验样本
  result <- list(Mu = Mu.value[-(1:(B+1)), ], Sigma = Sigma.value[, , -(1:(B+1))])
  return(result)
}
