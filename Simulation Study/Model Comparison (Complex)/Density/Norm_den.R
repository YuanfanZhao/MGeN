Norm.density <- function(x, z, Beta, Sigma) {
  d <- ncol(Beta)
  
  # 位置参数Mu
  Mu <- t(Beta) %*% z
  
  return(dmvnorm(x, Mu, Sigma, log = FALSE))
}