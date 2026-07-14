GaN.ECM <- function(X, Mu_new, Sigma_new, Lambda, max_iter = 1000, ep = 1e-8, 
                    Mu0 = NULL, Kappa0 = NULL, Lambda0 = NULL, Nu0 = NULL) {
  
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  # US算法求c+1/(s-1)+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s - 1) + 1 / (s - 1) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  Gamma.moment <- function(X, Mu, Sigma) {
    
    # Gamma分布的参数α，β
    compute.beta <- function(x) {
      return(as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) / 2 + Lambda - 1)
    }
    alpha <- d / 2 + Lambda
    beta <- apply(X, 1, compute.beta)
    
    return(alpha / beta)
  }
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    Loglikeli <- GaN.loglikeli(X, Mu, Sigma, Lambda)
    
    a <- Gamma.moment(X, Mu, Sigma)
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
    
    Loglikeli_new <- GaN.loglikeli(X, Mu_new, Sigma_new, Lambda)
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
