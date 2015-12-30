#' The GEVr Distribution
#'
#' Random number generation (rgevr) and density (dgevr) functions for the GEVr distribution with parameters loc, scale, and shape.
#' Also, quantile function (qgev) and cumulative distribution function (pgev) for the GEV1 distribution.
#' @name gevr
#' @param x Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations
#' @param r Number of order statistics for each observation.
#' @param log.d Logical: Whether or not to return the log density. (FALSE by default)
#' @param lower.tail Logical: If TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
#' @param log.p Logical: If TRUE, probabilities p are given as log(p). (FALSE by default)
#' @param loc Location parameter.
#' @param scale Scale parameter.
#' @param shape Shape Parameter.
#' @examples
#' ## Plot the densities of the heavy and bounded upper tail forms of GEVr
#' set.seed(7)
#' dat1 <- rgevr(1000, 1, loc = 0, scale = 1, shape = 0.25)
#' dat2 <- rgevr(1000, 1, loc = 0, scale = 1, shape = -0.25)
#' hist(dat1, col = rgb(1, 0, 0, 0.5), xlim = c(-5, 10), ylim = c(0, 0.4),
#' main = "Histogram of GEVr Densities", xlab = "Value", freq = FALSE)
#' hist(dat2, col = rgb(0, 0,1, 0.5), add = TRUE, freq = FALSE)
#' box()
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation
#' \eqn{i = 1, \ldots, n}. Note that currently the quantile and cdf functions are only for the GEV1 distribution. The GEVr
#' distribution is also known as the r-largest order statistics model and is a generalization of the block maxima model (GEV1).
#' The density function is given by \deqn{f_r (x_1, x_2, ..., x_r | \mu, \sigma, \xi) = \sigma^{-r} \exp\Big\{-(1+\xi z_r)^{-\frac{1}{\xi}}
#' - \left(\frac{1}{\xi}+1\right)\sum_{j=1}^{r}\log(1+\xi z_j)\Big\}} for some location parameter \eqn{\mu},
#' scale parameter \eqn{\sigma > 0}, and shape parameter \eqn{\xi}, where \eqn{x_1 > \cdots > x_r}, \eqn{z_j = (x_j - \mu) / \sigma},
#' and \eqn{1 + \xi z_j > 0} for \eqn{j=1, \ldots, r}. When \eqn{r = 1}, this distribution is exactly the GEV distribution.
#'
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
NULL

#' @rdname gevr
#' @export
dgevr <- function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE) {
  x <- as.matrix(x)
  r <- ncol(x)
  if (min(scale) <= 0)
    stop("invalid scale")
  if ((length(loc) != 1 & length(scale) != 1) | (length(loc) != 1 & length(shape) != 1) | (length(scale) != 1 & length(shape) != 1))
    stop("only one parameter argument can be a vector")
  if ((length(loc) != 1 | length(scale) != 1 | length(shape) != 1) & (r > 1))
    stop("parameter arguments cannot be a vector if r > 1")
  if (nrow(x) > 1 & (length(loc) != 1 | length(scale) != 1 | length(shape) != 1))
    stop("cannot have a vector of parameters AND observations")
  if(length(shape) == 1) shape <- rep(shape, max(nrow(x), length(loc), length(scale)))
  w <- matrix(((x - loc) / scale), ncol = r)
  z <- matrix(w*shape, ncol = r)
  z <- pmax(z, -1)
  log.density <- ifelse(shape == 0, rowSums(-log(scale) - w) - exp(-w[,r]),
                        rowSums(-log(scale) - ((1/shape) + 1) * log1p(z)) - exp((-1/shape) * log1p(z[,r])))
  log.density[is.nan(log.density) | is.infinite(log.density)] <- -Inf
  if (!log.d)
    log.density <- exp(log.density)
  log.density
}


#' @rdname gevr
#' @export
rgevr <- function(n, r, loc = 0, scale = 1, shape = 0) {
  if (min(scale) <= 0)
    stop("invalid scale")
  if ((length(loc) != 1 & length(scale) != 1) | (length(loc) != 1 & length(shape) != 1) | (length(scale) != 1 & length(shape) != 1))
    stop("only one parameter argument can be a vector")
  if ((length(loc) != 1 | length(scale) != 1 | length(shape) != 1) & (r > 1))
    stop("parameter arguments cannot be a vector if r > 1")
  if (n > 1 & (length(loc) != 1 | length(scale) != 1 | length(shape) != 1))
    stop("cannot have a vector of parameters AND observations")
  umat <- matrix(runif(n * r), n, r)
  if (r > 1) matrix(qgev(t(apply(umat, 1, cumprod)), loc, scale, shape), ncol = r)
  else qgev(umat, loc, scale, shape)
}


#' @rdname gevr
#' @export
qgev <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  p <- as.vector(p)
  if (log.p)
    p <- exp(p)
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) <= 0)
    stop("invalid scale")
  if ((length(loc) != 1 & length(scale) != 1) | (length(loc) != 1 & length(shape) != 1) | (length(scale) != 1 & length(shape) != 1))
    stop("only one parameter argument can be a vector")
  if (length(p) > 1 & (length(loc) != 1 | length(scale) != 1 | length(shape) != 1))
    stop("cannot have vectorized parameters and probabilities")
  if (!lower.tail)
    p <- 1 - p
  if(length(shape) == 1) shape <- rep(shape, max(length(p), length(loc), length(scale)))
  ifelse(shape == 0, loc - scale * log(-log(p)), loc + scale * expm1(log(-log(p)) * -shape) / shape)
}


#' @rdname gevr
#' @export
pgev <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  q <- as.vector(q)
  if (min(scale) <= 0)
    stop("invalid scale")
  if ((length(loc) != 1 & length(scale) != 1) | (length(loc) != 1 & length(shape) != 1) | (length(scale) != 1 & length(shape) != 1))
    stop("only one parameter argument can be a vector")
  if (length(q) > 1 & (length(loc) != 1 | length(scale) != 1 | length(shape) != 1))
    stop("cannot have vectorized parameters and quantiles")
  if(length(shape) == 1) shape <- rep(shape, max(length(q), length(loc), length(scale)))
  w <- (q - loc) / scale
  z <- pmax(shape*w, -1)
  p <- ifelse(shape == 0, exp(-exp(-w)), exp(-exp((-1/shape)*log1p(z))))
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  p
}
