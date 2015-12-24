#' GPD Multiple Threshold GoF Testing
#'
#' Wrapper functions to test multiple thresholds for goodness-of-fit to the Generalized Pareto model. Can choose which test to run.
#' @param data Orginal, full data in vector form.
#' @param thresholds A set of threshold values (either this or a set of the number of extremes must be given, but not both).
#' @param nextremes A set of the number of upper extremes to be used.
#' @param method Which test to run to sequentially test the thresholds. Must be one of 'cvm', 'ad', pbscore', 'multscore', 'imasym', or 'impb'.
#' @param nsim Number of boostrap replicates for the 'cvm', 'ad', 'pbscore', 'multscore', and 'imasym' tests.
#' @param inner Number of inner boostrap replicates if 'impb' test is chosen.
#' @param outer Number of outer boostrap replicates if 'impb' test is chosen.
#' @param information To use observed or expected (default) information for the 'pbscore' and 'multscore' tests.
#' @param allowParallel If selected, should the 'cvm', 'ad', 'pbscore', or 'impb' procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @details Function returns a matrix containing the thresholds used, the number of observations above each threshold, 
#' the corresponding test statistics, p-values (raw and transformed), and parameter estimates at each threshold.
#' @examples
#' set.seed(7)
#' x <- rgpd(10000, 0, 5, 0.2)
#' ## A vector of thresholds to test
#' threshes <- c(1, 2, 3, 4, 5)
#' gpdSeqTests(x, thresholds = threshes, method = "ad")
#' @return threshold The threshold used for the test.
#' @return num.above The number of observations above the given threshold.
#' @return p.values Raw p-values for the tresholds tested.
#' @return ForwardStop Transformed p-values according to the ForwardStop stopping rule.
#' @return StrongStop Transformed p-values according to the StrongStop stopping rule.
#' @return statistic Returned test statistics of each individual test.
#' @return est.scale Estimated scale parameter for the given r.
#' @return est.shape Estimated shape parameter for the given r.
#' @export

gpdSeqTests <- function(data, thresholds = NA, nextremes = NA, method = c("cvm", "ad", "pbscore", "multscore", "imasym", "impb"),
                        nsim = NULL, inner = NULL, outer = NULL, information = c("expected", "observed"),
                        allowParallel = FALSE, numCores = 1) {
  if(is.na(nextremes) && is.na(thresholds))
    stop("Enter either a set of thresholds or number of upper extremes")
  if(!is.na(nextremes) && !is.na(thresholds))
    stop("Enter EITHER a set of thresholds or number of upper extremes")
  information <- match.arg(information)
  method <- match.arg(method)
  if(((method == "pbscore" || method == "multscore" || method == "imasym") & is.null(nsim)) ||
     (method == "impb"  && (is.null(inner) || is.null(outer))))
    stop("Need to specify the number of bootstrap replicates")
  if(any(!is.na(nextremes)))
    thresholds <- findthresh(data, nextremes)
  num <- length(thresholds)
  result <- matrix(0, num, 9)
  for(i in 1:num) {
    x <- data[data>thresholds[i]]
    result[i, 1] <- i
    if(method == "cvm") {
      if(!is.null(nsim))
        fit <- gpdCvm(x, bootstrap = TRUE, B = nsim, allowParallel = allowParallel, numCores = numCores)
      else
        fit <- gpdCvm(x, allowParallel = allowParallel, numCores = numCores)
    }
    if(method == "ad") {
      if(!is.null(nsim))
        fit <- gpdAd(x, bootstrap = TRUE, B = nsim, allowParallel = allowParallel, numCores = numCores)
      else
        fit <- gpdAd(x, allowParallel = allowParallel, numCores = numCores)
    }
    if(method == "pbscore")
      fit <- gpdPbScore(x, nsim, information = information, allowParallel = allowParallel, numCores = numCores)
    if(method == "multscore")
      fit <- gpdMultScore(x, nsim, information = information)
    if(method == "imasym")
      fit <- gpdImAsym(x, nsim)
    if(method == "impb")
      fit <- gpdImPb(x, inner, outer, allowParallel = allowParallel, numCores = numCores)
    result[i, 2] <- thresholds[i]
    result[i, 3] <- length(x)
    result[i, 4] <- fit$p.value
    result[i, 7] <- fit$statistic
    result[i, 8:9] <- fit$theta
  }
  result[, 5] <- rev(pSeqStop(rev(result[, 4]))$ForwardStop)
  result[, 6] <- rev(pSeqStop(rev(result[, 4]))$StrongStop)
  colnames(result) <- c("testnum", "threshold", "num.above", "p.values", "ForwardStop", "StrongStop", "statistic", "est.scale", "est.shape")
  as.data.frame(result)
}
