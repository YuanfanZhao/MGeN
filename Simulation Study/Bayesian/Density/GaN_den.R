GaN.density <- function(x, Mu, Sigma, Lambda) {
  d <- length(Mu)
  
  # 计算密度函数
  term1 <- gamma(d / 2 + Lambda) * (Lambda - 1)^Lambda / (2 * pi)^(d / 2) / det(Sigma)^(0.5) / gamma(Lambda)
  term2 <- (as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) / 2 + Lambda - 1)^(-d / 2 - Lambda)
  
  return(term1 * term2)
}