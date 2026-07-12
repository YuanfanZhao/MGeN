Mix.density <- function(x, z, Beta, Sigma, Kurtosis, Ratio) {
  # 确定各分布的参数Lambda
  Lambda1 <- (3 / Kurtosis) + 2
  Lambda2 <- 3 / Kurtosis
  Lambda3 <- (3 + sqrt(12 * Kurtosis + 9)) / (2 * Kurtosis)
  Lambda4 <- 3 / Kurtosis
  
  Ratio[1] * GaN.density(x, z, Beta, Sigma, Lambda1) + Ratio[2] * IGaN.density(x, z, Beta, Sigma, Lambda2) + Ratio[3] * IGauN.density(x, z, Beta, Sigma, Lambda3) + 
    Ratio[4] * RIGauN.density(x, z, Beta, Sigma, Lambda4) + Ratio[5] * Norm.density(x, z, Beta, Sigma)
}
