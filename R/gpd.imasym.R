#'GPD Asymptotic (Adjusted) Information Matrix Test
#'
#'Runs the IM Test using bootstrap estimated covariance matrix. Asymptotically (in sample size) follows an F(3, B-3) distribution (see reference for details).
#'@param data Data should be in vector form.
#'@param B Number of bootstrap replicates for the covariance estimate.
#'@param theta Estimate for theta in the vector form (scale, shape). If NULL, uses the MLE.
#'@references Dhaene, G., & Hoorelbeke, D. (2004). The information matrix test with bootstrap-based covariance matrix estimation. Economics Letters, 82(3), 341-347.
#'@examples
#'## Generate some data from GPD
#'dat <- rgpd(200, 0, 1, 0.2)
#'gpd.imasym(dat, 50)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Value of theta used in the test.
#'@return effective_bootnum Effective number of bootstrap replicates used for the covariance estimate (if some did not converge).
#'@export

gpd.imasym <- function(data, B, theta = NULL) {
  n <- length(data)
  if(is.null(theta)) {
    fit <- 9999
    try(fit <- gpd.fit(data, nextremes = n, method = "mle"), silent = TRUE)
    if (!is.list(fit))
      stop("Maximum likelihood failed to converge at initial step")
    scale <- fit$par.ests[1]
    shape <- fit$par.ests[2]
    theta <- c(scale, shape)
  }
  thresh <- findthresh(data, n)
  data <- data - thresh
  v <- gpd.imcov(data, B, theta)
  B <- v$boot_adj
  v <- v$cov
  u <- gpd.ind(data, theta)
  d <- colSums(u)
  stat <- (1/n) * t(d) %*% v %*% d
  stat <- as.vector(stat)
  stat <- stat*(B-3) / (3*B-3)
  p <- 1 - pf(stat, 3, (B-3))
  names(theta) <- c("scale", "shape")
  out <- list(stat, p, theta, B)
  names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  out
}
