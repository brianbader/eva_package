gev.spatpbgen <- function(rkmat, n, m, loc_est, scale_est, shape_est){
  uu <- sort(runif(n*m))[rkmat]
  bootdat <- matrix(uu, nrow = n)
  teststat <- 0
  for(k in 1:m){
    bootdat[, k] <- (((-log(bootdat[, k]))^(-shape_est[k])) - 1) * (scale_est[k]/shape_est[k]) + loc_est[k]                   
    teststat <- teststat + gev.pbscore(bootdat[, k], 0, information=information)$statistic  
  }
  teststat
}

#'GEV1 (Spatially) Correlated Semi-Parametric Bootstrap Global Score Test
#'
#'Semi-parametric bootstrap score test procedure to assess global goodness-of-fit to the GEV1 distribution in 
#'the presence of (spatial) correlation. Each marginal distribution (column) is assumed to be from the GEV1 
#'distribution. The test statistic is a sum of marginal parametric bootstrap score test statistics (as in the 
#'function gev.pbscore). In the prescence of correlation among marginal distributions, a semi-parametric 
#'bootstrap procedure is done to preserve the dependence structure. See Heffernan and Tawn (2004) for more 
#'details on that procedure.
#'@param data Data should contain m columns, each data from a marginal GEV1 distribution. The number of rows 
#'corresponds to the number of observations.
#'@param B Number of bootstrap replicates.
#'@param information To use observed or expected (default) information in the test.
#'@param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#'@param numCores If allowParallel is true, specify the number of cores to use.
#'@examples
#'## Not run
#'## Generate data with max-stable dependence and GEV1 margins
#'## install.packages("SpatialExtremes")
#'## library(SpatialExtremes)
#'## n.site <- 20
#'## n.obs <- 50
#'## coord <- matrix(runif(2*n.site, 0, 10), ncol = 2)
#'## colnames(coord) <- c("lon", "lat")
#'## data <- rmaxstab(n.obs, coord, "gauss", cov11 = 100, cov12 = 25, cov22 = 220)
#'## param.loc <- -10 + 2 * coord[,2]
#'## param.scale <- 5 + 2 * coord[,1] + coord[,2]^2
#'## param.shape <- rep(0.2, n.site)
#'## for (i in 1:n.site) data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i], param.shape[i])
#'## gev.spatpbscore(data, 1000)
#'@return statistic Test statistic.
#'@return p.value P-value for the test.
#'@return theta Parameter estimates at each site.
#'@references Heffernan, J. E., & Tawn, J. A. (2004). A conditional approach for multivariate extreme values (with discussion). 
#'Journal of the Royal Statistical Society: Series B (Statistical Methodology), 66(3), 497-546.
#'@import parallel
#'@export
gev.spatpbscore <- function(data, B, information=c("expected", "observed"), allowParallel=FALSE, numCores=1) {
  m <- ncol(data)
  n <- nrow(data)
  information <- match.arg(information)
  ## Compute original global test statistic and get original parameter estimates
  stat <- 0
  loc_est <- rep(0, m)
  scale_est <- rep(0, m)
  shape_est <- rep(0, m)
  unifdata <- matrix(0, n, m)
  for(i in 1:m) {
    y <- gev.pbscore(data[, i], 0, information=information)
    stat <- stat + y$statistic
    loc_est[i] <- y$theta[1]
    scale_est[i] <- y$theta[2]
    shape_est[i] <- y$theta[3]
    ## Transform data to uniform scale
    unifdata[, i] <- exp( -(((data[, i] - loc_est[i]) * (shape_est[i]/scale_est[i]) + 1)^(-1/shape_est[i])) )
  }
  
  ## Rank the original dataset
  rkmat <- matrix(rank(unifdata, ties.method="random"), nrow = n)
  if(allowParallel==TRUE){
    cl <- makeCluster(numCores)
    fun <- function(cl){
      parSapply(cl, 1:B, function(i,...) {gev.spatpbgen(rkmat, n, m, loc_est, scale_est, shape_est)})
    }
    teststat <- fun(cl)
    stopCluster(cl)  
  }
  else{
    teststat <- replicate(B, gev.spatpbgen(rkmat, n, m, loc_est, scale_est, shape_est))
  }
  p <- (sum(teststat > stat) + 1) / (B + 2)
  theta <- cbind(loc_est, scale_est, shape_est)
  colnames(theta) <- c("location", "scale", "shape")
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "estimates")
  out
}
