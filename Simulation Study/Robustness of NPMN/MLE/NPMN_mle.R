NPMN.MLE <- function(X, M, Mu_new = colMeans(X), Sigma_new = var(X), max_iter = 1000, ep = 1e-6) {
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
  compute.G <- function(X, Mu, Sigma, p) {
    compute.g <- function(x) {
      quad <- c(t(x - Mu) %*% solve(Sigma) %*% (x - Mu))
      
      h <- p * tau^(d / 2) * exp(-0.5 * tau * quad)
      
      return(h / sum(h))
    }
    
    G <- apply(X, 1, compute.g) |> t()
    return(G)
  }
  
  # 定义p初始化向量和似然函数
  p_new <- gl_result$weight
  Loglikeli_new <- NPMN.loglikeli(X, Mu_new, Sigma_new, p_new)
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    p <- p_new
    Loglikeli <- Loglikeli_new
    
    if (M == 1) {
      G <- compute.G(X, Mu, Sigma, p) %>% t()
    } else {
      G <- compute.G(X, Mu, Sigma, p)
    }
    
    a <- as.numeric(G %*% tau)
    A <- diag(a)
    
    Mu_new <- colSums(a * X / sum(a))
    # cat("Mu:", Mu_new, "\n")
    
    Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    p_new <- colMeans(G)
    # cat("p:", p_new, "\n")
    
    Loglikeli_new <- NPMN.loglikeli(X, Mu_new, Sigma_new, p_new)
    # cat("Loglikelihood:", Loglikeli_new, "\n")
    
    if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
      index <- TRUE
      break
    }
    
    # cat("Iteration number:", iter, "\n")
    iter <- iter + 1
  }
  
  result <- list(Mu = Mu_new, 
                 Sigma = Sigma_new, 
                 p = p_new, 
                 Loglikelihood = Loglikeli_new,
                 number = iter, 
                 index = index)
  
  return(result)
}
