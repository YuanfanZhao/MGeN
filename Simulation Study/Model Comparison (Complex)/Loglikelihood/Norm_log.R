Norm.loglikeli <- function(X, Z, Beta, Sigma) {
  d <- ncol(Beta)
  
  compute.density <- function(data) {
    x <- data[1:d]
    z <- data[-(1:d)]
    Norm.density(x, z, Beta, Sigma)
  }
  
  # 返回密度函数的对数和
  return(sum(log(apply(cbind(X, Z), 1, compute.density))))
}