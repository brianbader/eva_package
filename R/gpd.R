#'The Generalized Pareto Distribution (GPD)
#'
#'Density, distribution function, quantile function and random number generation for the Generalized Pareto
#'distribution with location, scale, and shape parameters.
#'@param x, q Vector of quantiles.
#'@param p Vector of probabilities.
#'@param n Number of observations.
#'@param loc, scale, shape Location, scale and shape parameters, respectively; the shape argument cannot be a vector (must have length one).
#'@param log.d Logical; if TRUE, the log density is returned.
#'@param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#'@param log.p Logical; if TRUE, probabilities p are given as log(p).
#'@examples dgpd(2:4, 1, 0.5, 0.8)
#'@examples pgpd(2:4, 1, 0.5, 0.8)
#'@examples qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.8)
#'@examples rgpd(6, 1, 0.5, 0.8)
#'@examples p <- (1:9)/10
#'@examples pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
#'@examples ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#'
#'@details The Genralized Pareto distribution function is given (Pickands, 1975)
#'by \deqn{H(y) = 1 - \Big[1 + \frac{\xi (y - \mu)}{\sigma}\Big]^{-1/\xi}} defined
#'on \deqn{\{y : y > 0, (1 + \xi (y - \mu) / \sigma) > 0 \}}, with location \deqn{\mu},
#'scale \deqn{\sigma}, and shape parameter \deqn{\xi}. On an unrelated note, the core
#'basis for these functions come from R package 'evd' (Stephenson, 2010).
#'
#'@references Pickands III, J. (1975). Statistical inference using extreme order statistics. Annals of Statistics, 119-131.
#'@references Stephenson, A. (2010). Package 'evd', Functions for extreme value distributions.
#'
#'@rdname gpd
#'@export
rgpd <- function (n, loc = 0, scale = 1, shape = 0)
{
  if (min(scale) < 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  if (shape == 0)
    return(loc - scale * log(runif(n)))
  else return(loc + scale * (runif(n)^(-shape) - 1)/shape)
}


#'@rdname gpd
#'@export
pgpd <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  q <- pmax(q - loc, 0)/scale
  if (shape == 0)
    p <- 1 - exp(-q)
  else {
    p <- pmax(1 + shape * q, 0)
    p <- 1 - p^(-1/shape)
  }
  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- log(p)
  p
}


#'@rdname gpd
#'@export
dgpd <- function (x, loc = 0, scale = 1, shape = 0, log.d = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  d <- (x - loc)/scale
  nn <- length(d)
  scale <- rep(scale, length.out = nn)
  index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
  if (shape == 0) {
    d[index] <- log(1/scale[index]) - d[index]
    d[!index] <- -Inf
  }
  else {
    d[index] <- log(1/scale[index]) - (1/shape + 1) * log(1 + shape * d[index])
    d[!index] <- -Inf
  }
  if (!log.d)
    d <- exp(d)
  d
}


#'@rdname gpd
#'@export
qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p)
    p <- exp(p)
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) < 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  if (!lower.tail)
    p <- 1 - p
  if (shape == 0)
    return(loc - scale * log(p))
  else return(loc + scale * (p^(-shape) - 1)/shape)
}
