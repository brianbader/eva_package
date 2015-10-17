#'Sequential Score Tests for the GEVr Model
#'
#'Sequentially performs either the multiplier or parametric bootstrap score test for the GEVr model.
#'@param data Data should be contain n rows, each a GEVr observation.
#'@param nsim Number of bootstrap simulations.
#'@param method Which test to run, multiplier or parametric (default) bootstrap.
#'@param information To use observed (default) or expected information in the test.
#'@param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, uses the MLE.
#'@param allowParallel Should the parametric boostrap procedure be run in parallel or not. Defaults to false.
#'@param numCores If allowParallel is true, specify the number of cores to use.
#'@examples
#'data <- rgevr(200, 3, loc=0.5, scale=1, shape=0.5)
#'gevscore.seqtests(data, 99)
#'@return A matrix containing the test statistics and p-value results of the sequential tests.
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#'@importFrom ismev rlarg.fit
#'@export

gevscore.seqtests <- function(data, nsim, method=c("pb", "mult"), information=c("observed", "expected"),
                              allowParallel=FALSE, numCores=1)
{
  data <- as.matrix(data)
  R <- ncol(data)
  n <- nrow(data)
  method <- match.arg(method)
  information <- match.arg(information)
  result <- matrix(0, R, 3)
  for(i in 1:R){
    result[i, 1] <- i
    x <- as.matrix(data[, 1:i])
    if(method == "mult") {
      fit <- gev.multscore(x, nsim, NULL, information)
    }
    else{
      fit <- gev.pbscore(x, nsim, information, allowParallel, numCores)
    }
    result[i, 2] <- fit$p.value
    result[i, 3] <- fit$statistic
  }
  colnames(result) <- c("r", "p.values", "statistic")
  as.data.frame(result)
}
