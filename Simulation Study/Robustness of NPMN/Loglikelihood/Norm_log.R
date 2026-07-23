Norm.loglikeli <- function(X, Mu, Sigma) {
  # 返回密度函数的对数和
  return(sum(log(apply(X, 1, Norm.density, Mu = Mu, Sigma = Sigma))))
}