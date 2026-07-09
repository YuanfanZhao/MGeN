RIGauN.density <- function(x, z, Beta, Sigma, Lambda) {
  d <- ncol(Beta)
  
  # 位置参数Mu
  Mu <- t(Beta) %*% z
  
  # GIGau分布的参数a，b，p
  a <- as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) + Lambda
  b <- Lambda
  p <- (d + 1) / 2
  
  # 计算密度函数
  term1 <- exp(Lambda) * sqrt(Lambda) / (2 * pi)^((d + 1) / 2) / det(Sigma)^(0.5)
  term2 <- 2 * besselK(sqrt(a * b), p, expon.scaled = FALSE) / (a / b)^(p / 2)
  
  return(term1 * term2)
}