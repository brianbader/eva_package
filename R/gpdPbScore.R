gpdPbGen <- function(n, theta, information) {
  data1 <- rgpd(n, loc = 0, scale = theta[1], shape = theta[2])
  fit1 <- 9999
  try(fit1 <- gpdFit(data1, nextremes = n, method = "mle"), silent = TRUE)
  if(!is.list(fit1)) {
    teststat <- NA
  } else {
    scale1 <- fit1$par.ests[1]
    shape1 <- fit1$par.ests[2]
    theta1 <- c(scale1, shape1)
    thresh1 <- findthresh(data1, n)
    data1 <- data1 - thresh1
    teststat <- gpdTestStat(data1, theta1, information)
  }
  teststat
}


#' GPD Parametric Bootstrap Score Test
#'
#' Parametric bootstrap score test procedure to assess goodness-of-fit to the Generalized Pareto distribution.
#' @param data Data should be in vector form.
#' @param B Number of bootstrap replicates.
#' @param information To use observed (default) or expected information in the test.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @examples
#' ## Generate some data from GPD
#' x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' gpdPbScore(x, 100)
#' @return statistic Test statistic.
#' @return p.value P-value for the test.
#' @return theta Estimated value of theta for the initial data.
#' @import parallel
#' @export

gpdPbScore <- function(data, B, information = c("observed", "expected"), allowParallel = FALSE, numCores = 1) {
  n <- length(data)
  information <-  match.arg(information)
  fit <- 9999
  try(fit <- gpdFit(data, nextremes = n, method = "mle"), silent = TRUE)
  if(!is.list(fit))
    stop("Maximum likelihood failed to converge at initial step")
  scale <- fit$par.ests[1]
  shape <- fit$par.ests[2]
  theta <- c(scale, shape)
  thresh <- findthresh(data, n)
  data <- data - thresh
  stat <- gpdTestStat(data, theta, information)
  if(allowParallel==TRUE) {
    cl <- makeCluster(numCores)
    fun <- function(cl) {
      parSapply(cl, 1:B, function(i,...) {gpdPbGen(n, theta, information)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  } else {
    teststat <- replicate(B, gpdPbGen(n, theta, information))
  }
  teststat <- teststat[!is.na(teststat)]
  B <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (B + 2)
  names(theta) <- c("scale", "shape")
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
