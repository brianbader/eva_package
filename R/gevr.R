
## Handles the limit (1 + shape*w)^(-1/shape) as shape -> 0
xix <- function(w, shape, terms = 100){
  x <- pmax(w*shape, -1)
  n <- -1/shape
  if(shape == 0){
    out <- exp(-w)
  }
  else{
    out <- 1
    for(i in 1:terms){
      if(choose(n,i) != Inf & choose(n,i) != -Inf) out <- out + choose(n,i)*x^i
    }
    out <- ifelse(abs(x) < 1, out, (1 + x)^n)
  }
  out
}


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
  R <- ncol(x)
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(loc) != 1 | length(scale) != 1 | length(shape) != 1)
    stop("parameter arguments cannot be a vector")
  w <- (x - loc) / scale
  log.density <- rowSums( -log(scale) + log(xix(w, shape)) )
  log.density <- log.density - xix(w[,R], shape)
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
  if (length(loc) != 1 | length(scale) != 1 | length(shape) != 1)
    stop("parameter arguments cannot be a vector")
  umat <- matrix(runif(n * r), n, r)
  if (r > 1) umat <- t(apply(umat, 1, cumprod))
  qgev(umat, loc, scale, shape)[,,drop=TRUE]
}


#'@rdname gevr
#'@export
qgev <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p)
    p <- exp(p)
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  if (!lower.tail)
    p <- 1 - p
  if (shape == 0) loc - scale * log(-log(p))
  else {
    gev.stand <- expm1(-shape * log(-log(p))) / shape
    loc + scale * gev.stand
  }
}


#'@rdname gevr
#'@export
pgev <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  w <- (q - loc) / scale
  p <- exp(-xix(w, shape))
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  p
}
