rRIGauN <- function(n, Mu, Sigma, Lambda) {
  
  # 生成τ分布的随机数
  tau <- 1 / rinvgauss(n, mean = 1, shape = Lambda)
  
  # 随机表示方法生成多元IGauN分布随机数
  Z <- rmvnorm(n, rep(0, nrow(Sigma)), Sigma)
  X <- matrix(rep(Mu,n), nrow = n, byrow = TRUE) + Z / sqrt(tau)
  return(X)
}