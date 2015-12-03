#' GEVr Entropy Difference Test
#'
#' Goodness-of-fit test for GEVr using the difference in likelihood between GEVr and GEV(r-1).
#' This can be used sequentially to test for the choice of r.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, uses the MLE from the largest (r-1) order statistics.
#' @examples
#' ## This will test if the GEV2 distribution fits the data.
#' x <- rgevr(100, 2, loc = 0.5, scale = 1, shape = 0.5)
#'
#' ## Use the MLE from the GEV1 data by leaving theta input NULL.
#' result <- gevrEd(x)
#' @return statistic Test statistic.
#' @return p.value P-value for the test.
#' @return theta Value of theta used in the test.
#' @details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i = 1, ..., n.
#' @export
gevrEd <- function(data, theta = NULL) {
  data <- as.matrix(data)
  R <- ncol(data)
  if(R == 1) stop("R must be at least two")
  n <- nrow(data)
  if(is.null(theta)) {
    data1 <- as.matrix(data[, 1:(R-1)])
    y <- tryCatch(gevrFit(data1, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
    if(is.null(y))
      stop("Maximum likelihood failed to converge at initial step")
    theta <- y$par.ests
  }
  Diff <- dgevr(data[, 1:R], loc = theta[1], scale = theta[2], shape = theta[3], log.d = TRUE) -
          dgevr(data[, 1:(R-1)], loc = theta[1], scale = theta[2], shape = theta[3], log.d = TRUE)
  EstVar <- sum((Diff - mean(Diff))^2) / (n-1)
  FirstMom  <- - log(theta[2]) - 1 + (1 + theta[3])*digamma(R)
  Diff <- sum(Diff) / n
  Diff <- sqrt(n)*(Diff  - FirstMom) / sqrt(EstVar)
  p.value <- 2*(1-pnorm(abs(Diff)))
  out <- list(Diff, p.value, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
