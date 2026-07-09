rMix <- function(size, Z, Beta, Sigma, Kurtosis, Ratio) {
  # 确定各分布的参数Lambda
  Lambda1 <- (3 / Kurtosis) + 2
  Lambda2 <- 3 / Kurtosis
  Lambda3 <- (3 + sqrt(12 * Kurtosis + 9)) / (2 * Kurtosis)
  Lambda4 <- 3 / Kurtosis
  
  X1 <- rGaN(size, Z, Beta, Sigma, Lambda1)
  X2 <- rIGaN(size, Z, Beta, Sigma, Lambda2)
  X3 <- rIGauN(size, Z, Beta, Sigma, Lambda3)
  X4 <- rRIGauN(size, Z, Beta, Sigma, Lambda4)
  
  mul <- rmultinom(size, 1, Ratio)
  X <- X1 * mul[1, ] + X2 * mul[2, ] + X3 * mul[3, ] + X4 * mul[4, ]
  return(X)
}
