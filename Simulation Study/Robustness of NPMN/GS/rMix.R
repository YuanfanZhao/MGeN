rMix <- function(n, Mu, Sigma, nv, Ratio) {
  
  # 生成多元Laplace数据的函数
  rmlaplace <- function(n, Mu, Sigma) {
    
    # 生成τ分布的随机数
    tau <- 1 / rgamma(n, 1, 0.5)
    
    # 随机表示方法生称多元GeN分布随机数
    z <- rmvnorm(n, rep(0, nrow(Sigma)), Sigma)
    x <- matrix(rep(Mu,n), nrow = n, byrow = TRUE) + z / sqrt(tau)
    return(x)
  }
  
  # 分别生成多元正态、多元t和多元Laplace并按照比例混合
  S <- nv / (nv - 2) * Sigma

  X1 <- LaplacesDemon::rmvt(n, Mu, Sigma, nv)
  X2 <- rmlaplace(n, Mu, S / 2)
  X3 <- mvrnorm(n, Mu, S)
  z <- rmultinom(n, 1, Ratio)
  X <- X1 * z[1, ] + X2 * z[2, ] + X3 * z[3, ]
  return(X)
}
