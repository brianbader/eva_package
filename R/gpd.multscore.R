#'GPD Multiplier Score Test
#'
#'Fast weighted bootstrap alternative to the parametric bootstrap procedure for the Generalized Pareto score test.
#'@param data Data should be in vector form. 
#'@param theta Estimate for theta in the vector form (scale, shape). If NULL, uses probability weighted moments as an initial estimate.
#'@param B Number of bootstrap replicates.
#'@param information To use observed or expected (default) information in the test.
#'@examples
#'data <- rgpd(100, 0, 1, 0.25)
#'gpd.multscore(data, 1000)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Value of theta used in the test.
#'@export
gpd.multscore <- function(data, B, theta = NULL, information = c("expected", "observed")) {
  n <- length(data)
  information <-  match.arg(information)
  if(is.null(theta)) {
    fit <- gpdfit(data, nextremes=n, method="pwm")
    scale <- fit$par.ests[1]
    shape <- fit$par.ests[2]
    theta <- c(scale, shape)
  }
  thresh <- min(data)
  data <- data - thresh + 0.000001
  u <- gpdscorectb(data, theta)
  w <- colSums(u)
  if(information == "observed"){
    info <- gpdfisherobs(data, theta)
  }
  else{
    info <- gpdfisher(data, theta)
  }
  stat <- t(w) %*% info %*% w
  stat <- as.vector(stat)
  h <- chol(info)
  u <- u %*% t(h)
  z <- matrix(rnorm(n*B, mean=0, sd=1), n, B)
  z <- scale(z, center = TRUE, scale = FALSE)
  v <- t(u) %*% z
  teststat <- diag(t(v) %*% v)
  p <- (sum(teststat > stat) + 1) / (B + 2)
  names(theta) <- c("scale", "shape")
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}