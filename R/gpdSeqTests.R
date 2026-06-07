#' GPD Multiple Threshold Goodness-of-Fit Testing
#'
#' Wrapper function to test multiple thresholds for goodness-of-fit to the Generalized Pareto model. Can choose which test to
#' run from the available tests in this package.
#' @param data Original, full dataset in vector form.
#' @param thresholds A set of candidate threshold values (either this or a set of the number of extremes must be given, but not both). Must be provided as a vector.
#' @param nextremes A set of candidate numbers of upper extremes to be used, provided as a vector.
#' @param method Which test to run to sequentially test the thresholds. Must be one of `ad', `cvm', `pbscore', `multscore', `imasym', or `impb'.
#' @param nsim Number of boostrap replicates for the `ad', `cvm', `pbscore', `multscore', and `imasym' tests.
#' @param inner Number of inner boostrap replicates if `impb' test is chosen.
#' @param outer Number of outer boostrap replicates if `impb' test is chosen.
#' @param information To use observed or expected (default) information for the `pbscore' and `multscore' tests.
#' @param allowParallel If selected, should the `cvm', `ad', `pbscore', or `impb' procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @details Function returns a matrix containing the thresholds used, the number of observations above each threshold,
#' the corresponding test statistics, p-values (raw and transformed), and parameter estimates at each threshold. The user must provide
#' the data, a vector of thresholds or number of upper extremes to be used, and select the test.
#' Thresholds are tested in the order supplied. For sequential threshold selection, provide candidate thresholds
#' from lowest to highest. If using \code{nextremes}, provide the corresponding numbers of upper extremes in
#' descending order so that the implied thresholds are increasing.
#' At a chosen significance level, rejecting the first \eqn{k} sequential hypotheses means that the first
#' \eqn{k} candidate thresholds are rejected as too low. The next candidate threshold, if not rejected,
#' is the selected threshold. If no adjusted p-value rejects, the first candidate threshold is retained;
#' if all adjusted p-values reject, none of the candidate thresholds is high enough.
#' @examples
#' set.seed(7)
#' x <- rgpd(10000, loc = 0, scale = 5, shape = 0.2)
#' ## Candidate thresholds should be listed from lowest to highest.
#' candidate_thresholds <- c(1.5, 2.5, 3.5, 4.5, 5.5)
#' z <- gpdSeqTests(x, thresholds = candidate_thresholds, method = "ad")
#'
#' ## For a 5% ForwardStop rule, rejecting the first k tests rejects
#' ## candidate_thresholds[1:k] as too low.
#' k <- max(which(z$ForwardStop <= 0.05), 0)
#' rejected_thresholds <- candidate_thresholds[seq_len(k)]
#' selected_threshold <- if(k < length(candidate_thresholds)) candidate_thresholds[k + 1] else NA
#' @return
#' \item{threshold}{The threshold used for the test.}
#' \item{num.above}{The number of observations above the given threshold.}
#' \item{p.values}{Raw p-values for the thresholds tested.}
#' \item{ForwardStop}{Transformed p-values according to the ForwardStop stopping rule.}
#' \item{StrongStop}{Transformed p-values according to the StrongStop stopping rule.}
#' \item{statistic}{Returned test statistics of each individual test.}
#' \item{est.scale}{Estimated scale parameter for the given threshold.}
#' \item{est.shape}{Estimated shape parameter for the given threshold.}
#' @import parallel
#' @export

gpdSeqTests <- function(data, thresholds = NULL, nextremes = NULL, method = c("ad", "cvm", "pbscore", "multscore", "imasym", "impb"),
                        nsim = NULL, inner = NULL, outer = NULL, information = c("expected", "observed"),
                        allowParallel = FALSE, numCores = 1) {
  if(is.null(nextremes) && is.null(thresholds))
    stop("Enter either a set of thresholds or number of upper extremes")
  if(!is.null(nextremes) && !is.null(thresholds))
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
        fit <- gpdCvm(x, bootstrap = TRUE, bootnum = nsim, allowParallel = allowParallel, numCores = numCores)
      else
        fit <- gpdCvm(x, allowParallel = allowParallel, numCores = numCores)
    }
    if(method == "ad") {
      if(!is.null(nsim))
        fit <- gpdAd(x, bootstrap = TRUE, bootnum = nsim, allowParallel = allowParallel, numCores = numCores)
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
  seqStop <- pSeqStop(result[, 4])
  result[, 5] <- seqStop$ForwardStop
  result[, 6] <- seqStop$StrongStop
  colnames(result) <- c("testnum", "threshold", "num.above", "p.values", "ForwardStop", "StrongStop", "statistic", "est.scale", "est.shape")
  as.data.frame(result)
}
