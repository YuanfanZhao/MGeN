IGauN.MLE <- function(X, Mu_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-8) {
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  GIG.moment <- function(X, Mu, Sigma, Lambda) {
    # GIGau分布的参数a，b，p
    compute.a <- function(x) {
      return(as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) + (Lambda - 1)^2 / Lambda)
    }
    a <- apply(X, 1, compute.a)
    b <- Lambda
    p <- (d - 1) / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1) / BesselK(sqrt(a * b), p) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1) / BesselK(sqrt(a * b), p) * (b / a)^(-0.5)
    
    return(list(E1 = result1, Ei1 = result2))
  }
  
  Loglikeli_new <- IGauN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    Lambda <- Lambda_new
    Loglikeli <- Loglikeli_new
    
    GIG.result <- GIG.moment(X, Mu, Sigma, Lambda)
    
    a <- GIG.result$E1
    Bb <- mean(GIG.result$Ei1)
    Ba <- mean(a)
    A <- diag(a)
    
    Mu_new <- colSums(a * X / sum(a))
    # cat("Mu:", Mu_new, "\n")
    
    Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    Lambda_new <- (-1 - sqrt(1 - 4 * Ba * (2 - Ba - Bb))) / (2 * (2 - Ba - Bb))
    # cat("Lambda:", Lambda_new, "\n")
    
    Loglikeli_new <- IGauN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
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
                 Lambda = Lambda_new, 
                 Loglikelihood = Loglikeli_new,
                 number = iter, 
                 index = index)
  return(result)
}