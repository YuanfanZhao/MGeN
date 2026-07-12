rNorm <- function(n, Z, Beta, Sigma) {
  # 位置参数μ（n×d）
  Mu <- Z %*% Beta
  
  Z <- rmvnorm(n, rep(0, nrow(Sigma)), Sigma)
  X <- Mu + Z
  return(X)
}
