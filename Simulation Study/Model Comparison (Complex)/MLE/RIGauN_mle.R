RIGauN.MLE <- function(X, Z, Beta_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-6) {
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  GIG.moment <- function(X, Z, Beta, Sigma, Lambda) {
    
    # 位置参数Mu
    Mu <- Z %*% Beta
    
    # GIGau分布的参数a，b，p
    compute.a <- function(data) {
      x <- data[1:d]
      mu <- data[-(1:d)]
      return(as.numeric(t(x - mu) %*% solve(Sigma) %*% (x - mu)) + Lambda)
    }
    a <- apply(cbind(X, Mu), 1, compute.a)
    b <- Lambda
    p <- (d + 1) / 2
    
    result1 <- BesselK(sqrt(a * b), p + 1) / BesselK(sqrt(a * b), p) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1) / BesselK(sqrt(a * b), p) * (b / a)^(-0.5)
    
    return(list(E1 = result1, Ei1 = result2))
  }
  
  Loglikeli_new <- RIGauN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
  
  while (iter <= max_iter) {
    Beta <- Beta_new
    Sigma <- Sigma_new
    Lambda <- Lambda_new
    Loglikeli <- Loglikeli_new
    
    GIG.result <- GIG.moment(X, Z, Beta, Sigma, Lambda)
    
    a <- GIG.result$E1
    Bb <- mean(GIG.result$Ei1)
    Ba <- mean(a)
    A <- diag(a)
    
    Beta_new <- solve(t(Z) %*% A %*% Z) %*% (t(Z) %*% A %*% X)
    # cat("Beta:", Beta_new, "\n")
    
    Q <- X - Z %*% Beta_new
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    Lambda_new <- 1 / (Ba + Bb - 2)
    # cat("Lambda:", Lambda_new, "\n")
    
    Loglikeli_new <- RIGauN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
    # cat("Loglikelihood:", Loglikeli_new, "\n")
    
    if (abs((Loglikeli_new - Loglikeli) / Loglikeli) < ep) {
      index <- TRUE
      break
    }
    
    # cat("Iteration number:", iter, "\n")
    iter <- iter + 1
  }
  
  result <- list(Beta = Beta_new, 
                 Sigma = Sigma_new, 
                 Lambda = Lambda_new, 
                 Loglikelihood = Loglikeli_new,
                 number = iter, 
                 index = index)
  return(result)
}