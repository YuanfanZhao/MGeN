t.loglikeli <- function(X, Mu, Sigma, nv) {
  
  # 返回密度函数的对数和
  return(sum(log(apply(X, 1, t.density, Mu = Mu, Sigma = Sigma, nv = nv))))
}