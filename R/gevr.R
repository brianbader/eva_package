#'The GEVr Distribution
#'
#'Random number generation (rgevr) and density (dgevr) functions for the GEVr distribution with parameters loc, scale, and shape.
#'Also, quantile function (qgev) and cumulative distribution function (pgev) for the GEV1 distribution.
#'@param x Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#'@param q Vector of quantiles.
#'@param p Vector of probabilities.
#'@param n Number of observations
#'@param r Number of order statistics for each observation.
#'@param log.d Logical: Whether or not to return the log density. (FALSE by default)
#'@param lower.tail Logical: If TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
#'@param log.p Logical: If TRUE, probabilities p are given as log(p). (FALSE by default)
#'@param loc Location parameter.
#'@param scale Scale parameter.
#'@param shape Shape Parameter.
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i = 1, ..., n. Note
#'that currently the quantile and cdf functions are only for the GEV1 distribution.
#'
#'@rdname gevr
#'@export
dgevr <- function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE)
{
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


#'@rdname gevr
#'@export
rgevr <- function(n, r, loc = 0, scale = 1, shape = 0)
{
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


#'@rdname gevr
#'@export
qgev <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
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


#'@rdname gevr
#'@export
pgev <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
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
