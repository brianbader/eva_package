gpd.adgen <- function(n, theta) {
  data1 <- rgpd(n, loc=0, scale=theta[1], shape=theta[2])
  fit1 <- 9999
  try(fit1 <- gpdfit(data1, nextremes=n, method="mle"), silent = TRUE)
  if(!is.list(fit1)){
    teststat <- NA
  }
  else{
    scale1 <- fit1$par.ests[1]
    shape1 <- fit1$par.ests[2]
    theta1 <- c(scale1, shape1)
    thresh1 <- min(data1)
    data1 <- data1 - thresh1 + 0.000001
    newdata1 <- pgpd(data1, loc = 0, scale = scale1, shape = shape1)
    newdata1 <- sort(newdata1)
    i <- seq(1, n, 1)
    teststat <- -n - (1/n)*sum((2*i - 1)*(log(newdata1) + log(1 - rev(newdata1))))
  }
  teststat
}


#'Generalized Pareto Distribution Anderson-Darling Test
#'
#'Anderson-Darling goodness-of-fit test for the Generalized Pareto distribution. Critical values are generated via parametric bootstrap.
#'@param data Data should be in vector form, assumed to be from the GP distribution.
#'@param B Number of bootstrap replicates.
#'@param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#'@param numCores If allowParallel is true, specify the number of cores to use.
#'@references Choulakian, V., & Stephens, M. A. (2001). Goodness-of-fit tests for the generalized Pareto distribution. Technometrics, 43(4), 478-484.
#'@examples
#'## Generate some data from GPD
#'dat <- rgpd(200, 0, 1, 0.2)
#'gpd.ad(dat, 100)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Estimated value of theta for the initial data.
#'@import parallel
#'@export
gpd.ad <- function(data, B, allowParallel=FALSE, numCores=1) {
  n <- length(data)
  fit <- 9999
  try(fit <- gpdfit(data, nextremes=n, method="mle"), silent = TRUE)
  if (!is.list(fit))
    stop("Maximum likelihood failed to converge at initial step")
  scale <- fit$par.ests[1]
  shape <- fit$par.ests[2]
  theta <- c(scale, shape)
  thresh <- min(data)
  data <- data - thresh + 0.000001
  newdata <- pgpd(data, loc = 0, scale = scale, shape = shape)
  newdata <- sort(newdata)
  i <- seq(1, n, 1)
  stat <- -n - (1/n)*sum((2*i - 1)*(log(newdata) + log(1 - rev(newdata))))
  if(allowParallel==TRUE){
    cl <- makeCluster(numCores)
    fun <- function(cl){
      parSapply(cl, 1:B, function(i,...) {gpd.adgen(n, theta)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  }
  else{
    teststat <- replicate(B, gpd.adgen(n, theta))
  }
  teststat <- teststat[!is.na(teststat)]
  B <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (B + 2)
  names(theta) <- c("scale", "shape")
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
