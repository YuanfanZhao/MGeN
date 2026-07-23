Mix.density <- function(x, Mu, Sigma, nv, Ratio) {
  S <- nv / (nv - 2) * Sigma
  
  Ratio[1] * t.density(x, Mu, Sigma, nv) + Ratio[2] * Laplace.density(x, Mu, S / 2) + Ratio[3] * Norm.density(x, Mu, S)
}
