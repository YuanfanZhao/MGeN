IGaN.MLE <- function(X, Z, Beta_new, Sigma_new, Lambda_new, max_iter = 1000, ep = 1e-6) {
  index <- FALSE
  iter <- 1
  
  n <- nrow(X)
  d <- ncol(X)
  
  # US算法求c+log(s)-ψ(s)=0的零点
  US <- function(s, c) {
    index <- 0
    k <- 0
    
    s_new <- s
    
    while (index == 0) {
      s <- s_new
      q <- c + log(s) - digamma(s) + pi^2 * s / 6 - 1 / s
      s_new <- (q + sqrt(q^2 + 2 * pi^2 / 3)) / (pi^2 / 3)
      
      if (abs(s_new - s) < 1e-5) {
        index <- 1
        break
      }
      k <- k + 1
    }
    return(s_new)
  }
  
  GIG.moment <- function(X, Z, Beta, Sigma, Lambda) {
    
    # 位置参数Mu
    Mu <- Z %*% Beta
    
    # GIGau分布的参数a，b，p
    compute.a <- function(data) {
      x <- data[1:d]
      mu <- data[-(1:d)]
      return(as.numeric(t(x - mu) %*% solve(Sigma) %*% (x - mu)))
    }
    a <- apply(cbind(X, Mu), 1, compute.a)
    b <- 2 * Lambda
    p <- d / 2 - Lambda
    
    result1 <- BesselK(sqrt(a * b), p + 1) / BesselK(sqrt(a * b), p) * (b / a)^0.5
    result2 <- BesselK(sqrt(a * b), p - 1) / BesselK(sqrt(a * b), p) * (b / a)^(-0.5)
    ep <- 1e-6
    result3 <- (log(BesselK(sqrt(a * b), p + ep, expo = FALSE)) - log(BesselK(sqrt(a * b), p, expo = FALSE))) / ep + 0.5 * log(b / a)
    
    return(list(E1 = result1, Ei1 = result2, Elog = result3))
  }
  
  Loglikeli_new <- IGaN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
  
  while (iter <= max_iter) {
    Beta <- Beta_new
    Sigma <- Sigma_new
    Lambda <- Lambda_new
    Loglikeli <- Loglikeli_new
    
    GIG.result <- GIG.moment(X, Z, Beta, Sigma, Lambda)
    
    a <- GIG.result$E1
    Bb <- mean(GIG.result$Ei1)
    Bd <- mean(GIG.result$Elog)
    A <- diag(a)
    
    Beta_new <- solve(t(Z) %*% A %*% Z) %*% (t(Z) %*% A %*% X)
    # cat("Beta:", Beta_new, "\n")
    
    Q <- X - Z %*% Beta_new
    Sigma_new <- t(Q) %*% A %*% Q / n
    # cat("Sigma:", Sigma_new, "\n")
    
    Lambda_new <- US(Lambda, 1 - Bb - Bd)
    # cat("Lambda:", Lambda_new, "\n")
    
    Loglikeli_new <- IGaN.loglikeli(X, Z, Beta_new, Sigma_new, Lambda_new)
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