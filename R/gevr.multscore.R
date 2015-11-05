#'GEVr Multiplier Score Test
#'
#'Fast weighted bootstrap alternative to the parametric bootstrap procedure for the GEVr score test.
#'@param data Data should be contain n rows, each a GEVr observation.
#'@param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, an initial estimate is provided using PWM.
#'@param B Number of bootstrap replicates.
#'@param information To use observed (default) or expected information in the test.
#'@examples
#'data <- rgevr(500, 5, loc=0.5, scale=1, shape=0.3)
#'result <- gevr.multscore(data, 1000)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Value of theta used in the test.
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#'@importFrom ismev rlarg.fit
#'@export

gevr.multscore <- function(data, B, theta = NULL, information=c("observed", "expected"))
{
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  information <-  match.arg(information)
  if(is.null(theta)) {
    y <- 9999
    if(R == 1) {
      try(y <- gevr.fit(data, method = "pwm"), silent = TRUE)
      if (!is.list(y))
        stop("PWM failed to converge at initial step")
      theta <- y$par.ests
    }
    else{
      try(y <- rlarg.fit(as.matrix(data[, 1:(R-1)]), show = FALSE), silent = TRUE)
      if (!is.list(y))
        stop("Maximum likelihood failed to converge at initial step")
      theta <- y$mle
    }
  }
  u <- gevrscorectb(data, theta)
  w <- colSums(u)
  if(information == "observed"){
    info <- gevrfisherobs(data, theta)
  }
  else{
    info <- gevrfisher(data, theta)
  }
  stat <- (1/n) * t(w) %*% info %*% w
  stat <- as.vector(stat)
  h <- chol(info)
  u <- u %*% t(h)
  z <- matrix(rnorm(n*B, mean=0, sd=1), n, B)
  z <- scale(z, center = TRUE, scale = FALSE)
  v <- t(u) %*% z
  teststat <- (1/n) * diag(t(v) %*% v)
  p <- (sum(teststat > stat) + 1) / (B+2)
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
