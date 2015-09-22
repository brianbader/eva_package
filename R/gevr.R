## Helper function for the rgevr function
qgev <- function(p, loc, scale, shape) {
  if (shape == 0) -1 * scale * log(log(1 / p)) + loc
  else scale * ((log(1 / p))^(-1 * shape) - 1) / shape + loc
}


#'The GEVr Distribution
#'
#'Random number generation (rgevr) and density (dgevr) functions for the GEVr distribution with parameters loc, scale, and shape.
#'@param x Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#'@param n Number of observations
#'@param r Number of order statistics for each observation
#'@param log.d Logical: Whether or not to return the log density (FALSE by default)
#'@param loc Location parameter
#'@param scale Scale parameter
#'@param shape Shape Parameter
#'@details GEVr data (in matrix x) should be of the form x[i,1] > x[i, 2] > ... > x[i, r] for each observation i=1, ..., n.
#'
#'@rdname gevr
#'@export
dgevr <- function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(loc) != 1 | length(scale) != 1 | length(shape) != 1)
    stop("parameter arguments cannot be a vector")
  x <- as.matrix(x)
  R <- ncol(x)
  z <- pmax((shape / scale) * (x - loc), -1)
  w <- (1 / scale) * (x - loc)
  if(shape == 0) {
    log.density <- rowSums( -log(scale) - w )
    log.density <- log.density - exp(-w[,R])
  }
  else {
    log.density <- rowSums( -log(scale) - ((1/shape) + 1) * log1p(z) )
    log.density <- log.density - (1 + z[,R])^(-1 / shape)
    log.density[is.nan(log.density)] <- -Inf
  }
  if(!log.d) {
    exp(log.density)
  }
  else {
    log.density
  }
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
