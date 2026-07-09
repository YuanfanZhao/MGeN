# 设定参数
size <- 1000
Beta <- matrix(c(1, -2.5, 0, 0, 1.5, -2), ncol = 2)
Sigma <- matrix(c(1, 0.5, 0.5, 2), nrow = 2)
Kurtosis <- 2
Ratio <- c(0.25, 0.25, 0.25, 0.25)
Z1 <- rep(1, size)
Z2 <- rnorm(size, 0, 1)
Z3 <- rpois(size, 2)
Z <- cbind(Z1, Z2, Z3)

# 生成混合数据
X <- rMix(size, Z, Beta, Sigma, Kurtosis, Ratio)


# 确定各分布的参数Lambda
Lambda1 <- (3 / Kurtosis) + 2
Lambda2 <- 3 / Kurtosis
Lambda3 <- (3 + sqrt(12 * Kurtosis + 9)) / (2 * Kurtosis)
Lambda4 <- 3 / Kurtosis
M <- 10
tau <- gauss.quad(M, kind = "laguerre")$nodes

GaN_fit <- GaN.MLE(X, Z, Beta, Sigma, Lambda1)
IGaN_fit <- IGaN.MLE(X, Z, Beta, Sigma, Lambda2)
IGauN_fit <- IGauN.MLE(X, Z, Beta, Sigma, Lambda3)
RIGauN_fit <- RIGauN.MLE(X, Z, Beta, Sigma, Lambda4)
NPMN_fit <- NPMN.MLE(X, Z, Beta, Sigma, M)

S <- sum(NPMN_fit$p / tau) * NPMN_fit$Sigma
