NPMN.density <- function(x, z, Beta, Sigma, p) {
  M <- length(p)
  d <- ncol(Beta)
  
  # 位置参数μ和节点τ
  Mu <- t(Beta) %*% z
  tau_grid <- gauss.quad(M, kind = "laguerre")$nodes
  
  quad <- c(t(x - Mu) %*% solve(Sigma) %*% (x - Mu))
  result <- sum(p * tau_grid^(d / 2) * exp(-0.5 * tau_grid * quad)) / (2 * pi)^(d / 2) / det(Sigma)^(1/2)
  return(result)
}
