gev.pbgen <- function(n, R, theta, information) {
  data1 <- rgevr(n, R, theta[1], theta[2], theta[3])
  data1 <- as.matrix(data1)
  y1 <- 9999
  try(y1 <- rlarg.fit(data1, show = FALSE), silent = TRUE)
  if(!is.list(y1)){
    teststat <- NA
  }
  else{
    theta1 <- y1$mle
    teststat <- gevteststat(data1, theta1, information)
  }
  teststat
}


#'GEVr Parametric Bootstrap Score Test
#'
#'Parametric bootstrap score test procedure to assess goodness-of-fit to the GEVr distribution.
#'@param data Data should be contain n rows, each a GEVr observation.
#'@param B Number of bootstrap replicates.
#'@param information To use observed (default) or expected information in the test.
#'@param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#'@param numCores If allowParallel is true, specify the number of cores to use.
#'@examples
#'## Generate some data from GEVr
#'dat <- rgevr(200, 5, loc=0.5, scale=1, shape=0.5)
#'gev.pbscore(dat, 99)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Initial value of theta used in the test.
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#'@importFrom ismev rlarg.fit
#'@import parallel
#'@export
gev.pbscore <- function(data, B, information=c("observed", "expected"), allowParallel=FALSE, numCores=1) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  information <-  match.arg(information)
  y <- 9999
  try(y <- rlarg.fit(data, show = FALSE), silent = TRUE)
  if (!is.list(y))
    stop("Maximum likelihood failed to converge at initial step")
  theta <- y$mle
  stat <- gevteststat(data, theta, information)
  if(allowParallel==TRUE){
    cl <- makeCluster(numCores)
    fun <- function(cl){
      parSapply(cl, 1:B, function(i,...) {gev.pbgen(n, R, theta, information)})
    }
    teststat <- fun(cl)
    stopCluster(cl)    
  }
  else{
    teststat <- replicate(B, gev.pbgen(n, R, theta, information))
  }
  teststat <- teststat[!is.na(teststat)]
  B <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (B+2)
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
