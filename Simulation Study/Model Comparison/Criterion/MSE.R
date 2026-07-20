MSE <- function(Z_test, Beta_true, Beta_est, X_test = NULL) {
  model_names <- c("GaN", "IGaN", "IGauN", "RIGauN")
  
  # 初始化存储结果的向量
  MSE_mu_results <- numeric(length(model_names))
  names(MSE_mu_results) <- model_names
  
  if (!is.null(X_test)) {
    MSE_pred_results <- numeric(length(model_names))
    names(MSE_pred_results) <- model_names
  }
  
  # 计算测试集上的真实理论均值矩阵
  Mu_true <- Z_test %*% Beta_true
  
  # 遍历每个模型，计算预测均值并对比
  for (mod in model_names) {
    Beta_hat <- Beta_est[[mod]]$Beta
    
    Mu_hat <- Z_test %*% Beta_hat
    
    # 指标 A：理论均值 MSE —— 评价 Beta 估计的纯粹精度
    MSE_mu_results[mod] <- mean((Mu_true - Mu_hat)^2)
    
    # 指标 B：经验预测 MSE —— 评价对新样本的预测能力
    if (!is.null(X_test)) {
      MSE_pred_results[mod] <- mean((X_test - Mu_hat)^2)
    }
  }
  
  # 汇总返回
  if (!is.null(X_test)) {
    return(list(MSE_mu = MSE_mu_results, MSE_pred = MSE_pred_results))
  } else {
    return(list(MSE_mu = MSE_mu_results))
  }
}
