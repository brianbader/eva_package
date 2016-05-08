gevCIfit <- function(z, locvec, scalevec, shapevec, umat, resampling, grouping, ordering) {
  if(resampling == "in") {
    x <- rgevr(z$n, 1, loc = locvec, scale = scalevec, shape = shapevec)
  }
  if(resampling == "ht") {
    nr <- nrow(umat)
    umat <- umat[sample(1:nr, nr, replace = TRUE), ]
    rkmat <- apply(umat, 2, rank, ties.method = "random")
    umat.new <- apply(rkmat, 2, function(x) sort(runif(nr))[x])
    x <- qgev(as.vector(umat.new)[order(order(grouping, ordering))], loc = locvec, scale = scalevec, shape = shapevec)
  }
  if(resampling == "sc") {
    nr <- nrow(umat)
    nc <- ncol(umat)
    gen.list <- rep(list(list(min = 0, max = 1)), nc)
    umat.new <- simulateMvMatrix(n = nr, cor.mat = cor(umat, method = "spearman"),
                                 distributions = rep("unif", nc), param.list = gen.list)
    x <- qgev(as.vector(umat.new)[order(order(grouping, ordering))], loc = locvec, scale = scalevec, shape = shapevec)
  }
  z1 <- tryCatch(gevrFit(x, method = z$method, information = z$information,
            locvars = as.data.frame(z$covars[[1]]), loclink = z$links[[1]], locform = z$forms[[1]],
            scalevars = as.data.frame(z$covars[[2]]), scalelink = z$links[[2]], scaleform = z$forms[[2]],
            shapevars = as.data.frame(z$covars[[3]]), shapelink = z$links[[3]], shapeform = z$forms[[3]],
            gumbel = z$gumbel, start = as.numeric(z$par.ests)),
            error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(z1)) rep(NA, (length(z$par.ests) + 1))
  else c(z1$par.ests, z1$converged)
}


#' Bootstrapped Standard Errors for Fitted GEV1 (Block Maxima) Models
#'
#' Computes bootstrap based standard errors for fitted block maxima models under independence and
#' two semi-parametric approaches to bootstrap from models with dependence.
#'
#' @param z A class object returned from `gevrFit'.
#' @param conf Confidence level. Defaults to 95 percent.
#' @param bootnum Number of bootstrap replicates.
#' @param resampling Whether to use independent, Heffernan-Tawn, or Spearman-Correlation resampling, respectively.
#' See details.
#' @param grouping Vector specifying grouping of the observations. Typically refers to location or site in environmental applications.
#' @param ordering Vector specifying the (temporal) ordering of the observations within each group.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#'
#' @details Under independence, resampling is entirely parametric -- that is, GEV1 data are generated
#' directly from the fitted distribution provided for each replicate. The Heffernan-Tawn approach is one
#' attempt to capture dependence within the resampling process, by preserving the temporal ordering across
#' groups when generating data. The Spearman-Correlation resampling approach, roughly speaking, generates
#' multivariate data according to the correlation structure of the original data residuals. Both are
#' considered semi-parametric approaches. See references and function `simulateMvMatrix' in the EnvStats
#' package for details.
#'
#' @references Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to inducing rank correlation among
#' input variables. Communications in Statistics-Simulation and Computation, 11(3), 311-334.
#' @references Heffernan, J. E., & Tawn, J. A. (2004). A conditional approach for multivariate extreme values (with
#' discussion). Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66(3), 497-546.
#' @examples
#' ## Not run
#' # n <- 100
#' # r <- 1
#' # x <- rgevr(n, r, loc = 100 + 1:n / 50,  scale = 1 + 1:n / 300, shape = 0)
#' # covs <- as.data.frame(seq(1, n, 1))
#' # names(covs) <- c("Trend1")
#' ## Create some unrelated covariates
#' # covs$Trend2 <- rnorm(n)
#' # covs$Trend3 <- 30 * runif(n)
#' # z <- gevrFit(data = x, method = "mle", locvars = covs, locform = ~ Trend1 + Trend2*Trend3,
#' #              scalevars = covs, scaleform = ~ Trend1)
#' # gevCIboot(z, bootnum = 100, resampling = "in")
#' ## Now create groupings for this data -- four sites
#' # grouping <- c(rep("A", n/4), rep("B", n/4), rep("C", n/4), rep("D", n/4))
#' # ordering <- rep(seq(1, n/4, 1), 4)
#' # gevCIboot(z, bootnum = 100, resampling = "sc", grouping = grouping, ordering = ordering)
#' ## See vignette for a detailed example with dependent data
#' @return
#' \item{CIs}{Confidence intervals for the GEV parameters.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates (only those that converged are used).}
#' \item{replicates}{Matrix of bootstrap replicates.}
#' @importFrom EnvStats simulateMvMatrix
#' @import parallel
#' @export
gevCIboot <- function(z, conf = .95, bootnum, resampling = c("in", "ht", "sc"), grouping = NULL, ordering = NULL,
                      allowParallel = FALSE, numCores = 1) {
  if(z$R > 1)
    stop("Bootstrap based parameter CIs are only available for R=1")
  resampling <- match.arg(resampling)
  dat.transformed <- NULL
  alpha <- (1 - conf) / 2
  locvec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  scalevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))
  if(!z$gumbel) {
    shapevec <- z$links[[3]](rowSums(t(z$par.ests[(z$parnum[1] + z$parnum[2] + 1):(z$parnum[1] + z$parnum[2] + z$parnum[3])] * t(z$covars[[3]]))))
  } else {
    shapevec <- rep(0, z$n)
  }
  if(resampling != "in") {
    if((length(grouping) != z$n) | (length(ordering) !=z$n))
      stop("Grouping vectors must be same length as data")
    if((z$n %% length(unique(grouping)) != 0) | (z$n %% length(unique(ordering)) != 0))
      stop("There must be the same number of observations in each group")
    if(length(unique(paste(grouping, ordering, sep=""))) != z$n)
      stop("Groupings are not unique")
    dat.transformed <- pgev(z$data, loc = locvec, scale = scalevec, shape = shapevec)
    dat.transformed <- matrix(dat.transformed[order(grouping, ordering)], ncol = length(unique(grouping)))
  }
  if(allowParallel == TRUE) {
    cl <- makeCluster(numCores)
    fun <- function(cl) {
      parSapply(cl, 1:bootnum, function(i,...) {gevCIfit(z, locvec, scalevec, shapevec, dat.transformed, resampling, grouping, ordering)})
    }
    bootsample <- fun(cl)
    stopCluster(cl)
  } else {
    bootsample <- replicate(bootnum, gevCIfit(z, locvec, scalevec, shapevec, dat.transformed, resampling, grouping, ordering))
  }
  bootsample <- bootsample[, (bootsample[nrow(bootsample), ] == 0 & !is.na(bootsample[nrow(bootsample), ]))]
  eff <- ncol(bootsample)
  bootsample <- bootsample[-nrow(bootsample), ]
  CIs <- t(apply(bootsample, 1, quantile, probs = c(alpha, 1 - alpha)))
  out <- list(CIs, eff, t(bootsample))
  names(out) <- c("CIs", "effective_bootnum", "replicates")
  out
}
