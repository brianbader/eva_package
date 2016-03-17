gpdRlPlot <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  p <- c(seq(0.001, 0.01, by = 0.005), seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
         0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, .999, seq(1, 100, by = 10))
  levels <- matrix(0, length(p), 3)
  for(i in 1:nrow(levels)) {
    y <- gpdRl(z, 1/p[i], conf = conf, method = method)
    levels[i, 1] <- y$Estimate
    levels[i, 2:3] <- y$CI
  }
  dat <- sort(z$data)
  plot((z$n + 1) / ((1:z$n) * z$npp), rev(dat), type = "n", log = "x",
       xlab = "Return Period", ylab = "Return Level", xlim = c(0.1, 1000), ylim = c(min(z$data, levels[, 1]), max(z$data, levels[, 1])))
  title("Return Level Plot")
  lines(1/p, levels[, 1])
  lines(1/p, levels[, 2], col = 4)
  lines(1/p, levels[, 3], col = 4)
  points((z$n + 1) / ((1:z$n) * z$npp), rev(dat))
}


gpdHist <- function(z) {
  excess <- z$data[z$data > z$threshold]
  h <- hist(excess, plot = FALSE)
  x <- seq(min(h$breaks), max(h$breaks), (max(h$breaks) - min(h$breaks))/1000)
  y <- dgpd(x, loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2])
  hist(excess, freq = FALSE, ylim = c(0, max(max(h$density), max(y))),
       xlab = "x", ylab = "Density", main = "Density Plot")
  points(excess, rep(0, length(excess)))
  lines(x, y, col = 4)
}


gpdPP <- function(z) {
  excess <- z$data[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  p <- pgpd(excess, loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2])
  p <- sort(p)
  plot(p, Series, xlab = "Empirical", ylab = "Model", xlim = c(0,1), ylim = c(0,1))
  title("Probability Plot")
  abline(0, 1, col = 4)
}


gpdQQ <- function(z) {
  excess <- z$data[z$data > z$threshold]
  n <- length(excess)
  Series <- seq(1, n, 1) / (n+1)
  emp <- qgpd(Series, loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2])
  plot(sort(excess), emp, xlab = "Empirical", ylab = "Model",
       xlim = c(min(excess, emp), max(excess, emp)), ylim = c(min(excess, emp), max(excess, emp)))
  title("Quantile Plot")
  abline(0, 1, col = 4)
}


#' Diagnostic plots for a fit to the Generalized Pareto distribution
#'
#' @param z A class object returned from 'gpdFit'.
#' @param conf Confidence level used in the return level plot.
#' @param method The method to compute the return level confidence interval - either delta method (default) or profile likelihood. Choosing profile likelihood may be quite slow.
#' @examples
#' ## Not run
#' # x <- rgpd(10000, loc = 0.5, scale = 1, shape = 0.1)
#' # z <- gpdFit(x, nextremes = 500)
#' # gpdDiag(z)
#' @return Provides return level, density, probability, and quantile plots for the GPD exceedances. The overlaid density is the 'true' density for the estimated parameters.
#' @details See the reference for details on how return levels are calculated.
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
#' @export
gpdDiag <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  par(ask = TRUE, mfcol = c(2, 2))
  try(gpdRlPlot(z, conf, method), silent = TRUE)
  try(gpdHist(z), silent = TRUE)
  try(gpdPP(z), silent = TRUE)
  try(gpdQQ(z), silent = TRUE)
  par(mfrow = c(1, 1))
}










