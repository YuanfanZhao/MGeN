Laplace.density <- function(x, Mu, Sigma) {
  d <- length(Mu)
  Q <- max(sqrt(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)), 1e-10)
  result <- Q^(1 - d / 2) * besselK(Q, d / 2 - 1, expo = FALSE) / (2 * pi)^(d/2) / det(Sigma)^0.5
  return(result)
}
