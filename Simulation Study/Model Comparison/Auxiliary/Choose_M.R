# 基于交叉验证与MSE选择最优区间个数M
Choose.M <- function(X, Z, Beta, Sigma, Kurtosis, Ratio, K = 5, ep = 1e-6) {
  n <- nrow(X)
  M_grid <- seq(2, floor(sqrt(n)), by = 1)
  # cat("M:", M_grid, "\n")
  
  compute_loglik_over_M <- function(M_grid) {
    loglik_vec <- numeric(length(M_grid))
    
    # total_steps <- length(loglik_vec)
    # pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    # step_count <- 0 # 初始化步数计数器
    
    for (k in 1:length(M_grid)) {
      M <- M_grid[k]
      fit <- NPMN.MLE(X, Z, Beta, Sigma, M)
      loglik_vec[k] <- fit$Loglikelihood
      
      # cat("M:", M, " Loglikelihood:", loglik_vec[k], "\n")
      # step_count <- step_count + 1
      # setTxtProgressBar(pb, step_count)
    }
    
    # close(pb)
    
    data.frame(
      M = M_grid,
      loglik = loglik_vec
    )
  }
  
  select_M_by_stability <- function(loglik_df, threshold = ep) {
    M_vals <- loglik_df$M
    ll_vals <- loglik_df$loglik
    
    delta <- diff(ll_vals) / abs(ll_vals[-length(ll_vals)])
    
    stability_df <- data.frame(
      M = M_vals[-length(M_vals)],
      delta = delta
    )
    
    # 选择第一个进入稳定区间的 M
    idx <- which(delta < threshold)[1]
    
    if (is.na(idx)) {
      return(max(M_vals))
    } else {
      return(M_vals[idx])
    }
  }
  
  select_M_by_CV <- function(M_grid) {
    n <- nrow(X)
    d <- ncol(X)
    fold_id <- sample(rep(1:K, length.out = n))
    
    cv_MSE <- matrix(0, nrow = length(M_grid), ncol = K)
    
    # total_steps <- K * length(M_grid)
    # pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    # step_count <- 0 # 初始化步数计数器
    
    for (k in 1:K) {
      X_train <- X[fold_id != k, , drop = FALSE]
      X_test  <- X[fold_id == k, , drop = FALSE]
      Z_train <- Z[fold_id != k, drop = FALSE]
      Z_test <- Z[fold_id == k, drop = FALSE]
      # cat("Fold:", k, "\n")
      
      for (j in 1:length(M_grid)) {
        M <- M_grid[j]
        fit <- NPMN.MLE(X_train, Z_train, Beta, Sigma, M)
        
        result <- numeric(length = nrow(X_test))
        
        m <- function(data) {
          x <- data[1:d]
          z <- data[-(1:d)]
          (Mix.density(x, z, Beta, Sigma, Kurtosis, Ratio) - NPMN.density(x, z, fit$Beta, fit$Sigma, fit$p))^2
        }
        
        result <- apply(cbind(X_test, Z_test), 1, m)
        cv_MSE[j, k] <- mean(result)
        
        # step_count <- step_count + 1
        # setTxtProgressBar(pb, step_count)
      }
    }
    
    # close(pb)
    
    mean_MSE <- rowMeans(cv_MSE)
    
    data.frame(
      M = M_grid,
      cv_MSE = mean_MSE
    )
  }
  
  choose_best_M_CV <- function(cv_df) {
    cv_df$M[which.min(cv_df$cv_MSE)]
  }
  
  # Step1: All likelihood
  ll_df <- compute_loglik_over_M(M_grid)
  # print(ll_df)
  
  # Step 2: stability screening
  M_stable <- select_M_by_stability(ll_df, ep)
  # print(M_stable)
  
  M_candidates <- M_grid[M_grid >= M_stable]
  # print(M_candidates)
  
  cv_df <- select_M_by_CV(M_candidates)
  # print(cv_df %>% round(digits = 4))
  
  # Step 4: final choice
  M_final <- choose_best_M_CV(cv_df)
  
  return(M_final)
}
