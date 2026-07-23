GaN.loglikeli <- function(X, Mu, Sigma, Lambda) {
  # 返回密度函数的对数和
  return(sum(log(apply(X, 1, GaN.density, Mu = Mu, Sigma = Sigma, Lambda = Lambda))))
}