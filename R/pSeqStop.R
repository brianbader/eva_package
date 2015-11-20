#' P-Value Sequential Adjustment
#'
#' Given a set of (ordered) p-values, returns p-values adjusted according to the ForwardStop and StrongStop stopping rules.
#' @param p Vector of ordered p-values.
#' @references G'Sell, M. G., Wager, S., Chouldechova, A., & Tibshirani, R. (2013). Sequential Selection Procedures and False Discovery Rate Control. arXiv preprint arXiv:1309.5352.
#' @references Bader B., Jun Y., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#' @examples
#' x <- rgevr(500, 10, loc = 0.5, scale = 1, shape = 0.5)
#' y <- gevrSeqTests(x, method = "ed")
#' pSeqStop(rev(y$p.values))
#' @return StrongStop Vector of ordered p-values adjusted for the familywise error rate.
#' @return ForwardStop Vector of ordered p-values adjusted for the false discovery rate.
#' @return UnAdjusted Vector of non-transformed p-values.
#' @export
pSeqStop <- function(p) {
  m <- length(p)
  int <- seq(1, m, 1)
  pFDR <- cumsum(-log(1-p[int])) / int
  pFWER <- rev(exp(cumsum(rev(log(p[int])) / rev(int)))) * (m / int)
  out <- cbind.data.frame(pFWER, pFDR, p)
  colnames(out) <- c("StrongStop", "ForwardStop", "UnAdjusted")
  out
}
