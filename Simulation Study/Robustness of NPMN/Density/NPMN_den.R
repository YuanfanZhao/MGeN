NPMN.density <- function(x, Mu, Sigma, p) {
  M <- length(p)
  d <- length(Mu)
  
  # 节点τ
  tau_grid <- gauss.quad(M, kind = "laguerre")$nodes
  
  quad <- c(t(x - Mu) %*% solve(Sigma) %*% (x - Mu))
  result <- sum(p * tau_grid^(d / 2) * exp(-0.5 * tau_grid * quad)) / (2 * pi)^(d / 2) / det(Sigma)^(1/2)
  return(result)
}
