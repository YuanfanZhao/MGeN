Norm.MLE <- function(X) {
  n <- nrow(X)
  Mu <- colMeans(X)
  Sigma <- cov(X) * (n - 1) / n
  
  Loglikeli <- Norm.loglikeli(X, Mu, Sigma)
  
  list(Mu = Mu,
       Sigma = Sigma, 
       Loglikelihood = Loglikeli)
}
