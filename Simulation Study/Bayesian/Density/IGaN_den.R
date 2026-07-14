IGaN.density <- function(x, Mu, Sigma, Lambda) {
  d <- length(Mu)
  
  # GIGau分布的参数a，b，p
  a <- as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu))
  b <- 2 * Lambda
  p <- d / 2 - Lambda
  
  # 计算密度函数
  term1 <- Lambda^Lambda / (2 * pi)^(d / 2) / det(Sigma)^(0.5) / gamma(Lambda)
  term2 <- 2 * besselK(sqrt(a * b), p, expon.scaled = FALSE) / (a / b)^(p / 2)
  
  return(term1 * term2)
}