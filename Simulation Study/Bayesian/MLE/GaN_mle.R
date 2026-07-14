GaN.MLE <- function(X, Mu_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-8) {
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
  
  Gamma.moment <- function(X, Mu, Sigma, Lambda) {
    # Gamma分布的参数α，β
    compute.beta <- function(x) {
      return(as.numeric(t(x - Mu) %*% solve(Sigma) %*% (x - Mu)) / 2 + Lambda - 1)
    }
    alpha <- d / 2 + Lambda
    beta <- apply(X, 1, compute.beta)
    
    
    result1 <- alpha / beta
    result2 <- digamma(alpha) - log(beta)
    
    return(list(E1 = result1, Elog = result2))
  }
  
  Loglikeli_new <- GaN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
  
  while (iter <= max_iter) {
    Mu <- Mu_new
    Sigma <- Sigma_new
    Lambda <- Lambda_new
    Loglikeli <- Loglikeli_new
    
    Gamma.result <- Gamma.moment(X, Mu, Sigma, Lambda)
    
    a <- Gamma.result$E1
    Bd <- mean(Gamma.result$Elog)
    Ba <- mean(a)
    A <- diag(a)
    
    Mu_new <- colSums(a * X / sum(a))
    # cat("Mu:", Mu_new, "\n")
    
    Q <- sweep(X, MARGIN = 2, STATS = Mu_new, FUN = "-")
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    Lambda_new <- US(Lambda, 1 - Ba + Bd)
    # cat("Lambda:", Lambda_new, "\n")
    
    Loglikeli_new <- GaN.loglikeli(X, Mu_new, Sigma_new, Lambda_new)
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