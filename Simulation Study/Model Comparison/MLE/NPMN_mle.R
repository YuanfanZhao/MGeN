NPMN.MLE <- function(X, Z, Beta_new, Sigma_new, M, max_iter = 1000, ep = 1e-6) {
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  # 构造 τ 的 partition
  gl.quad <- function(M) {
    gl <- gauss.quad(M, kind = "laguerre")
    x <- gl$nodes      # 标准拉盖尔节点
    w <- gl$weights    # 对应权重
    return(list(node = x, weight = w))
  }
  gl_result <- gl.quad(M)
  tau <- gl_result$node
  
  # 定义混合变量的矩
  compute.G <- function(X, Z, Beta, Sigma, p) {
    # 位置参数Mu
    Mu <- Z %*% Beta
    
    compute.g <- function(data) {
      x <- data[1:d]
      mu <- data[-(1:d)]
      quad <- c(t(x - mu) %*% solve(Sigma) %*% (x - mu))
      
      h <- p * tau^(d / 2) * exp(-0.5 * tau * quad)
      
      return(h / sum(h))
    }
    
    G <- apply(cbind(X, Mu), 1, compute.g) |> t()
    return(G)
  }
  
  # 定义p初始化向量和似然函数
  p_new <- gl_result$weight
  Loglikeli_new <- NPMN.loglikeli(X, Z, Beta_new, Sigma_new, p_new)
  
  while (iter <= max_iter) {
    Beta <- Beta_new
    Sigma <- Sigma_new
    p <- p_new
    Loglikeli <- Loglikeli_new
    
    G <- compute.G(X, Z, Beta, Sigma, p)
    
    if (M == 1) {
      a <- as.numeric(G * tau)
    } else {
      a <- as.numeric(G %*% tau)
    }
    A <- diag(a)
    
    Beta_new <- solve(t(Z) %*% A %*% Z) %*% (t(Z) %*% A %*% X)
    # cat("Beta:", Beta_new, "\n")
    
    Q <- X - Z %*% Beta_new
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    p_new <- colMeans(G)
    # cat("p:", p, "\n")
    
    Loglikeli_new <- NPMN.loglikeli(X, Z, Beta_new, Sigma_new, p_new)
    # cat("Loglikelihood:", Loglikeli_new, "\n")
    
    if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
      index <- TRUE
      break
    }
    
    # cat("Iteration number:", iter, "\n")
    iter <- iter + 1
  }
  
  result <- list(Beta = Beta_new, 
                 Sigma = Sigma_new, 
                 p = p_new, 
                 Loglikelihood = Loglikeli_new,
                 number = iter, 
                 index = index)
}
