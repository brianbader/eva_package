gevrRlPlot <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  p <- c(seq(0.001, 0.01, by = 0.005), seq(0.01, 0.09, by = 0.01), 0.1, 0.2, 0.3, 0.4, 0.5,
         0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.995, 0.999)
  levels <- matrix(0, length(p), 3)
  for(i in 1:nrow(levels)) {
    y <- gevrRl(z, 1/p[i], conf = conf, method = method)
    levels[i, 1] <- y$Estimate
    levels[i, 2:3] <- y$CI
  }
  plot(-1/log((1:length(z$data[,1]))/(length(z$data[,1]) + 1)), sort(z$data[,1]), log = "x", type = "n",
       xlab = "Return Period", ylab = "Return Level",  xlim = c(0.1, 1000), ylim = c(min(z$data[, 1], levels[, 1]), max(z$data[, 1], levels[, 1])))
  title("Return Level Plot")
  lines(-1/log(1-p), levels[,1])
  lines(-1/log(1-p), levels[,2], col = 4)
  lines(-1/log(1-p), levels[,3], col = 4)
  points(-1/log((1:length(z$data[,1]))/(length(z$data[,1]) + 1)), sort(z$data[,1]))
}


pgevMarg1 <- function(z, j) {
  p <- rep(0, z$n)
  locvec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  scalevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))
  shapevec <- z$links[[3]](rowSums(t(z$par.ests[(z$parnum[1] + z$parnum[2] + 1):(z$parnum[1] + z$parnum[2] + z$parnum[3])] * t(z$covars[[3]]))))
  for(i in 1:z$n) {
    for(k in 0:(j-1)) {
      p[i] <- p[i] + nzsh((z$data[i, j] - locvec[i]) / scalevec[i], shapevec[i])^k / gamma(k+1)
    }
    p[i] <- p[i] * exp(-nzsh((z$data[i, j] - locvec[i]) / scalevec[i], shapevec[i]))
  }
  p
}


pgevMarg2 <- function(x, z, j, i) {
  p <- 0
  locvec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  scalevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))
  shapevec <- z$links[[3]](rowSums(t(z$par.ests[(z$parnum[1] + z$parnum[2] + 1):(z$parnum[1] + z$parnum[2] + z$parnum[3])] * t(z$covars[[3]]))))
  for(k in 0:(j-1)) {
    p <- p + nzsh((x - locvec[i]) / scalevec[i], shapevec[i])^k / gamma(k+1)
  }
  p * exp(-nzsh((x - locvec[i]) / scalevec[i], shapevec[i]))
}


gevrQQ <- function(z, j) {
  qgevMarg <- function(x, q, i) {
    q - pgevMarg2(x, z, j, i)
  }
  emp <- rep(0, z$n)
  Series <- seq(1, z$n, 1) / (z$n + 1)
  locvec <- z$links[[1]](rowSums(t(z$par.ests[1:z$parnum[1]] * t(z$covars[[1]]))))
  scalevec <- z$links[[2]](rowSums(t(z$par.ests[(z$parnum[1] + 1):(z$parnum[1] + z$parnum[2])] * t(z$covars[[2]]))))
  shapevec <- z$links[[3]](rowSums(t(z$par.ests[(z$parnum[1] + z$parnum[2] + 1):(z$parnum[1] + z$parnum[2] + z$parnum[3])] * t(z$covars[[3]]))))
  ztrans <- nzsh(-(z$data[, j] - locvec) / scalevec, - shapevec)
  Series <- Series[order(Series)[rank(ztrans)]]
  for(i in 1:z$n) {
    emp[i] <- uniroot(qgevMarg, interval = c(min(z$data[, j]) - 2, max(z$data[, j]) + 2), q = Series[i], i = i)$root
  }
  plot(z$data[, j], emp, xlab = "Empirical", ylab = "Model",
       xlim = c(min(z$data[, j], emp), max(z$data[, j], emp)), ylim = c(min(z$data[, j], emp), max(z$data[, j], emp)))
  if(z$stationary)
    title(paste("Quantile Plot, j=", j, sep = ""))
  if(!z$stationary)
    title(paste("Residual Quantile Plot, j=", j, sep = ""))
  abline(0, 1, col = 4)
}


## Provide it the fitted object and which marginal statistic to plot (j)
gevrPP <- function(z, j) {
  n <- z$n
  Series <- seq(1, z$n, 1) / (z$n + 1)
  p <- pgevMarg1(z, j)
  p <- sort(p)
  plot(p, Series, xlab = "Empirical", ylab = "Model", xlim = c(0,1), ylim = c(0,1))
  if(z$stationary)
    title(paste("Probability Plot, ", "j=", j, sep = ""))
  if(!z$stationary)
    title(paste("Residual Probability Plot, ", "j=", j, sep = ""))
  abline(0, 1, col = 4)
}


dgevMarg <- function(x, j, loc = loc, scale = scale, shape = shape) {
  if(length(shape) == 1)
    shape <- rep(shape, max(length(x), length(loc), length(scale)))
  w <- (x - loc) / scale
  ifelse(shape == 0,   exp(-exp(-w) - j*w) / (scale * factorial(j-1)),
         (nzsh(w, shape)^j / (scale * gamma(j))) * exp(-nzsh(w, shape)))
}


gevrHist <- function(z, j) {
  h <- hist(z$data[, j], plot = FALSE)
  x <- seq(min(h$breaks), max(h$breaks), (max(h$breaks) - min(h$breaks))/1000)
  if(j == 1)
    y <- dgevr(x, loc = z$par.ests[1], scale = z$par.ests[2], shape = z$par.ests[3])
  if(j > 1)
    y <- dgevMarg(x, j, loc = z$par.ests[1], scale = z$par.ests[2], shape = z$par.ests[3])
  hist(z$data[, j], freq = FALSE, ylim = c(0, max(max(h$density), max(y))),
       xlab = "x", ylab = "Density", main = paste("Density Plot, j=", j, sep = ""))
  points(z$data[, j], rep(0, length(z$data[, j])))
  lines(x, y, col = 4)
}


#' Diagnostic plots for a fit to the GEVr distribution.
#'
#' @param z A class object returned from 'gevrFit'.
#' @param conf Confidence level used in the return level plot.
#' @param method The method to compute the return level confidence interval - either delta method (default) or profile
#' likelihood. Choosing profile likelihood may be quite slow.
#' @examples
#' ## Not run
#' # x <- rgevr(500, 2, loc = 0.5, scale = 1, shape = 0.1)
#' # z <- gevrFit(x)
#' # gevrDiag(z)
#' @return For stationary models, provides return level plot and density, probability,
#' and quantile plots for each marginal order statistic. The overlaid density is the 'true' marginal
#' density for the estimated parameters. For nonstationary models, provides only residual probability and quantile plots.
#' @details In certain cases the quantile plot may fail, because it requires solving a root equation. See the references for details.
#' @references Tawn, J. A. (1988). An extreme-value theory model for dependent observations. Journal of Hydrology, 101(1), 227-250.
#' @references Smith, R. L. (1986). Extreme value theory based on the r largest annual events. Journal of Hydrology, 86(1), 27-43.
#' @export
gevrDiag <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  par(ask = TRUE, mfcol = c(2, 2))
  if(z$stationary)
    try(gevrRlPlot(z, conf, method), silent = TRUE)
  for(i in 1:z$R) {
    if(z$stationary)
      try(gevrHist(z, i), silent = TRUE)
    try(gevrPP(z, i), silent = TRUE)
    try(gevrQQ(z, i), silent = TRUE)
  }
  par(mfrow = c(1, 1))
}
