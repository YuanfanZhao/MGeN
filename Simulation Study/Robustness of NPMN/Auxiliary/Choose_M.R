# 基于交叉验证与(MSE, MAE和Likelihood)选择最优区间个数M
Choose.M <- function(X, Mu, Sigma, nv, Ratio, K = 5, ep = 1e-6) {
  n <- nrow(X)
  M_grid <- seq(1, floor(sqrt(n)), by = 1)
  # cat("M:", M_grid, "\n")
  
  compute_loglik_over_M <- function(M_grid) {
    loglik_vec <- numeric(length(M_grid))
    
    # total_steps <- length(loglik_vec)
    # pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    # step_count <- 0 # 初始化步数计数器
    
    for (k in 1:length(M_grid)) {
      M <- M_grid[k]
      fit <- NPMN.MLE(X, M)
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
    
    cv_criterion <- array(0, c(length(M_grid), K, 3))
    
    # total_steps <- K * length(M_grid)
    # pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
    # step_count <- 0 # 初始化步数计数器
    
    for (k in 1:K) {
      X_train <- X[fold_id != k, , drop = FALSE]
      X_test  <- X[fold_id == k, , drop = FALSE]
      # cat("Fold:", k, "\n")
      
      for (j in 1:length(M_grid)) {
        M <- M_grid[j]
        fit <- NPMN.MLE(X_train, M)
        
        m <- function(x) {
          MSE.value <- (Mix.density(x, Mu, Sigma, nv, Ratio) - NPMN.density(x, fit$Mu, fit$Sigma, fit$p))^2
          MAE.value <- abs(Mix.density(x, Mu, Sigma, nv, Ratio) - NPMN.density(x, fit$Mu, fit$Sigma, fit$p))
          Density.value <- NPMN.density(x, fit$Mu, fit$Sigma, fit$p)
          return(c(MSE.value, MAE.value, Density.value))
        }
        
        result <- apply(X_test, 1, m)
        cv_criterion[j, k, 1] <- mean(result[1, ])
        cv_criterion[j, k, 2] <- mean(result[2, ])
        cv_criterion[j, k, 3] <- -(result[3, ] |> log() |> sum())
        
        # step_count <- step_count + 1
        # setTxtProgressBar(pb, step_count)
      }
    }
    
    # close(pb)
    
    mean_criterion <- apply(cv_criterion, c(1, 3), mean)
    
    list(
      M = M_grid,
      cv_criterion = mean_criterion
    )
  }
  
  choose_best_M_CV <- function(cv_df) {
    M.MSE <- cv_df$M[which.min(cv_df$cv_criterion[, 1])]
    M.MAE <- cv_df$M[which.min(cv_df$cv_criterion[, 2])]
    M.Log <- cv_df$M[which.min(cv_df$cv_criterion[, 3])]
    data.frame(MSE = M.MSE, 
               MAE = M.MAE, 
               Log = M.Log)
  }
  
  # Step1: All likelihood
  ll_df <- compute_loglik_over_M(M_grid)
  # print(ll_df)
  
  # Step 2: stability screening
  M_stable <- select_M_by_stability(ll_df, ep)
  # print(M_stable)
  
  M_candidates <- M_grid[M_grid >= M_stable]
  # print(M_candidates)
  
  if (length(M_candidates) > 1) {
    cv_df <- select_M_by_CV(M_candidates)
    
    # Step 4: final choice
    M_final <- choose_best_M_CV(cv_df)
  } else {
    M_final <- M_candidates
  }
  
  return(M_final)
}
