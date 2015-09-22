#'Sequential Entropy Difference Tests for the GEVr Model
#'
#'Sequentially performs the entropy difference (ED) test for the dataset.
#'@param data Data should be contain n rows, each a GEVr observation.
#'@param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, uses the MLE from the largest (r-1) order statistics.
#'@examples
#'dat <- rgevr(500, 5, loc=0.5, scale=1, shape=0.5)
#'ed.seqtests(dat)
#'@return A matrix containing the p-value and test statistics of the sequential tests.
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#'@importFrom ismev rlarg.fit
#'@export

ed.seqtests <- function(data, theta = NULL) {
  R <- ncol(data)
  if(R==1) stop("R must be at least two")
  n <- nrow(data)
  result <- matrix(0, R-1, 3)
  for(i in 2:R){
    data1 <- as.matrix(data[, 1:(i-1)])
    data2 <- as.matrix(data[, 1:i])
    result[i-1, 1] <- i
    if(is.null(theta)){
      y <- 9999
      try(y <- rlarg.fit(data1, show = FALSE), silent = TRUE)
      if (!is.list(y))
        stop("Maximum likelihood failed to converge at one of the steps")
      theta <- y$mle
    }
    z <- ed.test(data2, theta)
    result[i-1, 2] <- z$p.value
    result[i-1, 3] <- z$statistic
  }
  colnames(result) <- c("r", "p.value", "statistic")
  result
}
