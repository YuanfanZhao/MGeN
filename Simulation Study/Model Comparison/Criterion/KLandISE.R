KLandISE <- function(X_test, Z_test, para_true, para_est) {
  N_test <- nrow(X_test)
  
  # 初始化存储每个样本对数密度的向量/矩阵
  log_f_true <- numeric(N_test)
  log_f_hat <- matrix(0, nrow = N_test, ncol = 5)
  colnames(log_f_hat) <- c("GaN", "IGaN", "IGauN", "RIGauN", "NPMN")
  
  # 设置一个极小值防止 log(0) 导致数值溢出 (Numerical Underflow)
  epsilon <- 1e-300
  
  for (i in 1:N_test) {
    # 提取单行并保持为矩阵形式 (1 x d) 和 (1 x q+1)
    x_i <- X_test[i, ]
    z_i <- Z_test[i]
    
    # ================= 1. 计算真实的对数密度 =================
    f_true_val <- Mix.density(x_i, z_i, 
                              para_true$Beta, 
                              para_true$Sigma, 
                              para_true$Kurtosis, 
                              para_true$Ratio)
    log_f_true[i] <- log(max(f_true_val, epsilon))
    
    # ================= 2. 计算各个模型的估计对数密度 =================
    
    # 2.1 GaN
    f_GaN <- GaN.density(x_i, z_i, para_est$GaN$Beta, para_est$GaN$Sigma, para_est$GaN$Lambda)
    log_f_hat[i, "GaN"] <- log(max(f_GaN, epsilon))
    
    # 2.2 IGaN
    f_IGaN <- IGaN.density(x_i, z_i, para_est$IGaN$Beta, para_est$IGaN$Sigma, para_est$IGaN$Lambda)
    log_f_hat[i, "IGaN"] <- log(max(f_IGaN, epsilon))
    
    # 2.3 IGauN
    f_IGauN <- IGauN.density(x_i, z_i, para_est$IGauN$Beta, para_est$IGauN$Sigma, para_est$IGauN$Lambda)
    log_f_hat[i, "IGauN"] <- log(max(f_IGauN, epsilon))
    
    # 2.4 RIGauN
    f_RIGauN <- RIGauN.density(x_i, z_i, para_est$RIGauN$Beta, para_est$RIGauN$Sigma, para_est$RIGauN$Lambda)
    log_f_hat[i, "RIGauN"] <- log(max(f_RIGauN, epsilon))
    
    # 2.5 NPMN
    f_NPMN <- NPMN.density(x_i, z_i, para_est$NPMN$Beta, para_est$NPMN$Sigma, para_est$NPMN$p)
    log_f_hat[i, "NPMN"] <- log(max(f_NPMN, epsilon))
  }
  
  # ================= 3. 计算并返回 KL 散度 =================
  mean_log_f_true <- mean(log_f_true)
  mean_log_f_hat <- colMeans(log_f_hat)
  
  KLD_results <- mean_log_f_true - mean_log_f_hat
  
  # 2. 计算 积分平方误差 (ISE/MISE)
  # 公式: ISE = mean( (f_true - f_hat)^2 / f_true )
  f_true_vec <- exp(log_f_true) # 长度为 N_test 的向量
  ISE_results <- numeric(5)
  names(ISE_results) <- c("GaN", "IGaN", "IGauN", "RIGauN", "NPMN")
  
  # 用一个简单的循环对每个模型计算 ISE，保证 R 向量计算不串位
  for(k in 1:5) {
    f_hat_vec <- exp(log_f_hat[, k]) # 第 k 个模型在所有测试样本上的预测概率
    ISE_results[k] <- mean( (f_true_vec - f_hat_vec)^2 / f_true_vec )
  }
  
  # 返回 KL 散度与ISE
  return(list(
    KLD = KLD_results,
    ISE = ISE_results
  ))
}
