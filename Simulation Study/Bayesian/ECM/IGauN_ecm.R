IGauN.ECM <- function(X, Mu_new, Sigma_new, Lambda, Mu0, Kappa0, Lambda0, Nu0, max_iter = 1000, ep = 1e-8) {
  
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  GIGau.moment <- function(X, Mu, Sigma) {
    
    # GIGau分布的参数a，b，p
    compute.a <- function(x) {
      return(as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) + (Lambda - 1)^2 / Lambda)
    }
    a <- apply(X, 1, compute.a)
    b <- Lambda
    p <- (d - 1) / 2
    
    return(BesselK(sqrt(a * b), p + 1) / BesselK(sqrt(a * b), p) * (b / a)^0.5)
  }
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    Loglikeli <- IGauN.loglikeli(X, Mu, Sigma, Lambda)
    
    a <- GIGau.moment(X, Mu, Sigma)
    A <- diag(a)
    
    if (is.null(Mu0)) {
      Mu_new <- colSums(a * X / sum(a))
      # cat("Mu:", Mu_new, "\n")
      
      Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
      Sigma_new <- t(Q) %*% A %*% Q / n
      # cat("Sigma:", Sigma_new, "\n")
    } else {
      Mu_new <- (colSums(a * X) + Kappa0 * Mu0) / (sum(a) + Kappa0)
      # cat("Mu:", Mu_new, "\n")
      
      Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
      Sigma_new <- (t(Q) %*% A %*% Q + Lambda0 + Kappa0 * (Mu_new - Mu0) %*% t(Mu_new - Mu0)) / (n + Nu0 + d + 2)
      # cat("Sigma:", Sigma_new, "\n")
    }
    
    Loglikeli_new <- IGauN.loglikeli(X, Mu_new, Sigma_new, Lambda)
    # cat("Loglikelihood:", Loglikeli_new, "\n")
    
    if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
      index <- TRUE
      break
    }
    
    # cat("Iteration number:", iter, "\n")
    iter <- iter + 1
  }
  
  result <- list(Mu = Mu_new, 
                 Sigma = Sigma_new, 
                 number = iter, 
                 index = index)
}
