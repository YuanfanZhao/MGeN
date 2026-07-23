NPMN.loglikeli <- function(X, Mu, Sigma, p) {
  # 返回密度函数的对数和
  return(sum(log(apply(X, 1, NPMN.density, Mu = Mu, Sigma = Sigma, p = p))))
}