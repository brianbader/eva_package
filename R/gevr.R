
## Accurately computes log(1 - exp(-a))
log1mexp <- function(a){
  if(a <= log(2)) out <- log(-expm1(-a))
  else out <-log1p(-exp(-a))
  out
}

## Accurately computes log(1 + exp(a))
log1pexp <- function(a){
  if(a <= -37) out <- exp(a)
  if(a > - 37 & a <= 18) out <- log1p(exp(a))
  if(a > 18 & a <= 33.3) out <- a + exp(-a)
  if(a > 33.3) out <- a
  out
}

## Handles the limit (1 + shape*w)^(-1/shape) as shape -> 0
xix <- function(w, shape){
  z <- pmax((1 + shape*w), 0)
  out <- z^(-1/shape)
  if(abs(shape) < 1e-3 & shape != 0){
    eps <- 1e-3
    shape1 <- shape + eps
    shape2 <- shape - eps
    check <- z^(-1/shape)
    check1 <- (1 + shape1*w)^(-1/shape1)
    check2 <- (1 + shape2*w)^(-1/shape2)
    if((check1 < check & check < check2) || (check2 < check & check < check1))
      out <- z^(-1/shape)
    else
      out <- exp(-w)
  }
  if(shape == 0) out <- exp(-w)
  out
}





#'@rdname gevr
#'@export
qgev <- function(p, loc, scale, shape)
{
  if (shape == 0) -1 * scale * log(-log(p)) + loc
  else {
    gev.stand <- expm1(-shape * log(-log(p))) / shape
    scale * gev.stand + loc
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
  x <- as.matrix(x)
  R <- ncol(x)
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(loc) != 1 | length(scale) != 1 | length(shape) != 1)
    stop("parameter arguments cannot be a vector")
  w <- (x - loc) / scale
  log.density <- rowSums( -log(scale) + log(xix(w, shape)) )
  log.density <- log.density -  xix(w[,R], shape)
  log.density[is.nan(log.density)] <- -Inf
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
