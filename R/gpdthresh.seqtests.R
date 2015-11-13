#'GPD Multiple Threshold GoF Testing
#'
#'Wrapper functions to test multiple thresholds for goodness-of-fit to the Generalized Pareto model. Can choose which test to run.
#'@param data Orginal, full data in vector form.
#'@param thresholds A set of threshold values (either this or a set of the number of extremes must be given, but not both).
#'@param nextremes A set of the number of upper extremes to be used.
#'@param method Which test to run to sequentially test the thresholds. Must be one of 'cvm', 'ad', pbscore', 'multscore', 'imasym', or 'impb'.
#'@param nsim Number of boostrap replicates for the 'cvm', 'ad', 'pbscore', 'multscore', and 'imasym' tests.
#'@param inner Number of inner boostrap replicates if 'impb' test is chosen.
#'@param outer Number of outer boostrap replicates if 'impb' test is chosen.
#'@param information To use observed (default) or expected information for the 'pbscore' and 'multscore' tests.
#'@param allowParallel If selected, should the 'cvm', 'ad', 'pbscore', or 'impb' procedure be run in parallel or not. Defaults to false.
#'@param numCores If allowParallel is true, specify the number of cores to use.
#'@return A matrix containing the thresholds used, the number of observations above each threshold, and the corresponding test statistics and p-values.
#'@export
gpdthresh.seqtests <- function(data, thresholds = NA, nextremes = NA, method = c("cvm", "ad", "pbscore", "multscore", "imasym", "impb"),
                               nsim = NULL, inner = NULL, outer = NULL, information = c("observed", "expected"),
                               allowParallel = FALSE, numCores = 1) {
  if(is.na(nextremes) && is.na(thresholds))
    stop("Enter either a set of thresholds or number of upper extremes")
  if(!is.na(nextremes) && !is.na(thresholds))
    stop("Enter EITHER a set of thresholds or number of upper extremes")
  information <- match.arg(information)
  method <- match.arg(method)
  if(((method == "cvm" || method == "ad" || method == "pbscore" || method == "multscore" || method == "imasym") & is.null(nsim)) ||
       (method == "impb"  && (is.null(inner) || is.null(outer))))
    stop("Need to specify the number of bootstrap replicates")
  if(any(!is.na(nextremes)))
    thresholds <- findthresh(data, nextremes)
  num <- length(thresholds)
  result <- matrix(0, num, 5)
  for(i in 1:num) {
    x <- data[data>thresholds[i]]
    result[i, 1] <- i
    if(method == "cvm")
      fit <- gpd.cvm(x, nsim, allowParallel, numCores)
    if(method == "ad")
      fit <- gpd.ad(x, nsim, allowParallel, numCores)
    if(method == "pbscore")
      fit <- gpd.pbscore(x, nsim, information = information, allowParallel, numCores)
    if(method == "multscore")
      fit <- gpd.multscore(x, nsim, information = information)
    if(method == "imasym")
      fit <- gpd.imasym(x, nsim)
    if(method == "impb")
      fit <- gpd.impb(x, inner, outer, allowParallel, numCores)
    result[i, 2] <- thresholds[i]
    result[i, 3] <- length(x)
    result[i, 4] <- fit$p.value
    result[i, 5] <- fit$statistic
  }
  colnames(result) <- c("test", "threshold", "# above threshold", "p.values", "statistic")
  as.data.frame(result)
}
