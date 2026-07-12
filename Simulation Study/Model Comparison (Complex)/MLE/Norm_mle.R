Norm.MLE <- function(X, Z) {
  n <- nrow(X)
  Beta <- solve(t(Z) %*% Z) %*% t(Z) %*% X
  Sigma <- t(X - Z %*% Beta) %*% (X - Z %*% Beta) / n
  
  list(Beta = Beta, 
       Sigma = Sigma)
}