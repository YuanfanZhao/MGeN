Laplace.MLE <- function(X, max_iter = 1000, ep = 1e-8) {
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  GIG.moment <- function(X, Mu, sigma) {
    
    compute.a <- function(x) {
      return(max(t(x - Mu) %*% solve(Sigma) %*% (x - Mu), 1e-10))
    }
    a <- apply(X, 1, compute.a)
    b <- 1
    p <- d / 2 - 1
    
    result1 <- BesselK(sqrt(a * b), p + 1, expo = TRUE) / BesselK(sqrt(a * b), p, expo = TRUE) * (b / a)^0.5
    
    return(result1)
  }
  
  Mu_new <- colMeans(X)
  Sigma_new <- var(X)
  Loglikeli_new <- Laplace.loglikeli(X, Mu_new, Sigma_new)
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    Loglikeli <- Loglikeli_new
    
    a <- GIG.moment(X, Mu, Sigma)
    A <- diag(a)
    
    Mu_new <- colSums(a * X / sum(a))
    # cat("Mu:", Mu_new, "\n")
    
    Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    Loglikeli_new <- Laplace.loglikeli(X, Mu_new, Sigma_new)
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
                 Loglikelihood = Loglikeli_new,
                 number = iter, 
                 index = index)
  return(result)
}
