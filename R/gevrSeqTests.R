#' Sequential Tests for the GEVr Model
#'
#' Sequentially performs the entropy difference (ED) test or the multiplier or parametric bootstrap score tests for the GEVr model.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param method Which test to run: ED test (ed), multiplier (multscore) or parametric bootstrap (pbscore) score test.
#' @param nsim If method equals 'pbscore' or 'multscore', the number of bootstrap simulations to use.
#' @param information To use observed (default) or expected information in the score tests.
#' @param allowParallel If method equals 'pbscore', should the parametric boostrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @examples
#' x <- rgevr(200, 3, loc = 0.5, scale = 1, shape = 0.5)
#' gevrSeqTests(x, method = "ed")
#' @return A matrix containing the test statistics and p-value results of the sequential tests.
#' @details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#' @export

gevrSeqTests <- function(data, nsim = NULL, method = c("ed", "pbscore", "multscore"), information = c("observed", "expected"),
                              allowParallel = FALSE, numCores = 1) {
  data <- as.matrix(data)
  R <- ncol(data)
  method <- match.arg(method)
  if(method != "ed"){
    if(is.null(nsim))
      stop("Must enter the number of bootstrap replicates!")
    information <- match.arg(information)
    result <- matrix(0, R, 6)
    for(i in 1:R) {
      result[i, 1] <- i
      if(method == "multscore")
        fit <- gevrMultScore(data[, 1:i], nsim, NULL, information)
      if(method == "pbscore")
        fit <- gevrPbScore(data[, 1:i], nsim, information, allowParallel, numCores)
      result[i, 2] <- fit$p.value
      result[i, 3] <- fit$statistic
      result[i, 4:6] <- fit$theta
    }
  } else {
    if(R == 1)
      stop("R must be at least two")
    result <- matrix(0, R-1, 6)
    for(i in 2:R) {
      result[i-1, 1] <- i
      fit <- gevrEd(data[, 1:i])
      result[i-1, 2] <- fit$p.value
      result[i-1, 3] <- fit$statistic
      result[i-1, 4:6] <- fit$theta
    }
  }
  colnames(result) <- c("r", "p.values", "statistic", "est.loc", "est.scale", "est.shape")
  as.data.frame(result)
}
