#' GEVr Multiplier Score Test
#'
#' Fast weighted bootstrap alternative to the parametric bootstrap procedure for the GEVr score test.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, an initial estimate is provided.
#' @param B Number of bootstrap replicates.
#' @param information To use expected (default) or observed information in the test.
#' @examples
#' x <- rgevr(500, 5, loc = 0.5, scale = 1, shape = 0.3)
#' result <- gevrMultScore(x, 1000)
#' @return statistic Test statistic.
#' @return p.value P-value for the test.
#' @return theta Value of theta used in the test.
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation \eqn{i = 1, \ldots, n}.
#' @references Bader B., Jun Y., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#' @export

gevrMultScore <- function(data, B, theta = NULL, information = c("expected", "observed")) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  information <- match.arg(information)
  if(is.null(theta)) {
    if(R == 1) {
      y <- tryCatch(gevrFit(data, method = "pwm"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
      if(is.null(y))
        stop("PWM failed to converge at initial step")
      theta <- y$par.ests
    } else {
      y <- tryCatch(gevrFit(as.matrix(data[, 1:(R-1)]), method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
      if(is.null(y))
        stop("Maximum likelihood failed to converge at initial step")
      theta <- y$par.ests
    }
  }
  u <- gevrScore(data, theta)
  w <- colSums(u)
  if(information == "observed") {
    info <- gevrFisherObs(data, theta)
  } else {
    info <- gevrFisher(data, theta)
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
