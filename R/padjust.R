#'P-Value Sequential Adjustment
#'
#'Fast weighted bootstrap alternative to the parametric bootstrap procedure for the Generalized Pareto score test.
#'@param p Vector of ordered p-values.
#'@references G'Sell, M. G., Wager, S., Chouldechova, A., & Tibshirani, R. (2013). Sequential Selection Procedures and False Discovery Rate Control. arXiv preprint arXiv:1309.5352.
#'@references Bader B., Jun Y., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#'@examples
#'dat <- rgevr(500, 10, loc=0.5, scale=1, shape=0.5)
#'y <- gev.seqtests(dat, method = "ed")
#'padjust(rev(y$p.values))
#'@return StrongStop Vector of ordered p-values adjusted for the familywise error rate.
#'@return ForwardStop Vector of ordered p-values adjusted for the false discovery rate.
#'@return UnAdjusted Vector of non-transformed p-values.
#'@export
padjust <- function(p){
  m <- length(p)
  pFWER <- rep(0, m)
  pFDR <- rep(0, m)
  int <- seq(1, m, 1)
  for(k in 1:m){
    pFWER[k] <- (m/k)*exp(sum(log(p[k:m])/int[k:m]))
    pFDR[k] <- -(1/k)*sum(log(1-p[1:k]))
  }
  out <- list(pFWER, pFDR, p)
  names(out) <- c("StrongStop", "ForwardStop", "UnAdjusted")
  out
}
