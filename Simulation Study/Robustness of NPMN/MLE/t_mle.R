t.MLE <- function(X) {
  fit <- fit.mst(X)
  return(list(
    Mu = fit$mu,
    Sigma = fit$Sigma %>% as.matrix(),
    nv = fit$df, 
    loglik = fit$ll.max
  ))
}
