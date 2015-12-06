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


pgevMarg <- function(data, theta, j) {
  data <- as.matrix(data)
  n <- nrow(data)
  p <- rep(0, n)
  for(i in 1:n) {
    for(k in 0:(j-1)) {
      p[i] <- p[i] + nzsh((data[i, 1] - theta[1])/theta[2], theta[3])^k / gamma(k+1)
    }
    p[i] <- p[i] * exp(-nzsh((data[i, 1] - theta[1])/theta[2], theta[3]))
  }
  p
}


## Provide it the fitted object and which marginal statistic to plot (j)
gevrPP <- function(z, j) {
  n <- z$n
  Series <- seq(1, n, 1) / (n+1)
  p <- pgevMarg(z$data[,j], z$par.ests, j)
  p <- sort(p)
  plot(p, Series, xlab = "Empirical", ylab = "Model", xlim = c(0,1), ylim = c(0,1))
  title(paste("Probability Plot, ", "j=", j, sep = ""))
  abline(0, 1, col = 4)
}


gevrQQ <- function(z, j) {
  n <- z$n
  qgevMarg <- function(x, q) {
    q - pgevMarg(x, z$par.ests, j)
  }
  emp <- rep(0, n)
  Series <- seq(1, n, 1) / (n+1)
  for(i in 1:n) {
    emp[i] <- uniroot(qgevMarg, interval = c(min(z$data[, j]) - 2, max(z$data[, j]) + 2), q = Series[i])$root
  }
  plot(sort(z$data[, j]), emp, xlab = "Empirical", ylab = "Model",
       xlim = c(min(z$data[, j], emp), max(z$data[, j], emp)), ylim = c(min(z$data[, j], emp), max(z$data[, j], emp)))
  title(paste("Quantile Plot, j=", j, sep = ""))
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
#' @param z A class object returned from gevr.fit.
#' @param conf Confidence level used in the return level plot.
#' @param method The method to compute the return level confidence interval - either delta method (default) or profile likelihood. Choosing profile likelihood may be quite slow.
#' @examples
#' ## Not run
#' # x <- rgevr(500, 2, loc = 0.5, scale = 1, shape = 0.1)
#' # z <- gevrFit(x)
#' # gevrDiag(z)
#' @return Provides return level plot and density, probability, and quantile plots for each marginal order statistic. The overlaid density is the 'true' marginal density for the estimated parameters.
#' @details In certain cases the quantile plot may fail, because it requires solving a root equation. See the references for details.
#' @references Tawn, J. A. (1988). An extreme-value theory model for dependent observations. Journal of Hydrology, 101(1), 227-250.
#' @references Smith, R. L. (1986). Extreme value theory based on the r largest annual events. Journal of Hydrology, 86(1), 27-43.
#' @export
gevrDiag <- function(z, conf = 0.95, method = c("delta", "profile")) {
  method <- match.arg(method)
  opar <- par(ask = TRUE, mfcol = c(2, 2))
  try(gevrRlPlot(z, conf, method), silent = TRUE)
  for(i in 1:z$R) {
    try(gevrHist(z, i), silent = TRUE)
    try(gevrPP(z, i), silent = TRUE)
    try(gevrQQ(z, i), silent = TRUE)
  }
  par(opar)
}


