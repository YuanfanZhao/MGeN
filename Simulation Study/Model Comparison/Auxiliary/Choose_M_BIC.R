# 基于BIC选择最优区间个数M
Choose.M.BIC <- function(X, Z, Beta, Sigma, Method = "AIC") {
  n <- nrow(X)
  d <- ncol(X)
  
  M_grid <- seq(2, floor(sqrt(n)), by = 1)
  cat("M:", M_grid, "\n")
  
  compute_BIC_over_M <- function(M_grid) {
    BIC_vec <- numeric(length(M_grid))
    
    # total_steps <- length(loglik_vec)
    # pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    # step_count <- 0 # 初始化步数计数器
    
    for (k in 1:length(M_grid)) {
      M <- M_grid[k]
      para <- d + d * (d + 1) / 2 + M
      fit <- NPMN.MLE(X, Z, Beta, Sigma, M)
      if (Method == "AIC") {
        BIC_vec[k] <- 2 * para - 2 * fit$Loglikelihood
      } else if (Method == "BIC") {
        BIC_vec[k] <- para * log(n) - 2 * fit$Loglikelihood
      }
      
      cat("M:", M, " BIC:", BIC_vec[k], "\n")
      # step_count <- step_count + 1
      # setTxtProgressBar(pb, step_count)
    }
    
    # close(pb)
    
    data.frame(
      M = M_grid,
      BIC = BIC_vec
    )
  }
  
  # Step1: All likelihood
  ll_df <- compute_BIC_over_M(M_grid)
  # print(ll_df)
  
  # Step 2: stability screening
  M_final <- ll_df$M[which.min(ll_df$BIC)]
  # print(M_stable)
  return(M_final)
}
