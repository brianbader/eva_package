#'The Generalized Pareto Distribution (GPD)
#'
#'Density, distribution function, quantile function and random number generation for the Generalized Pareto
#'distribution with location, scale, and shape parameters.
#'@param x Vector of observations.
#'@param q Vector of quantiles.
#'@param p Vector of probabilities.
#'@param n Number of observations.
#'@param loc Location parameter.
#'@param scale Scale parameter.
#'@param shape Shape parameter.
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
#'scale \deqn{\sigma}, and shape parameter \deqn{\xi}.
#'
#'@references Pickands III, J. (1975). Statistical inference using extreme order statistics. Annals of Statistics, 119-131.
#'
#'@rdname gpd
#'@export
rgpd <- function (n, loc = 0, scale = 1, shape = 0)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(loc) != 1 | length(scale) != 1 | length(shape) != 1)
    stop("invalid shape")
  qgpd(runif(n), loc, scale, shape)
}


#'@rdname gpd
#'@export
pgpd <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  w <- (q - loc) / scale
  p <- 1 - xix(w, shape)
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
  w <- (x - loc) / scale
  log.density <- -log(scale) + log(xix(w, shape))
  log.density[is.nan(log.density) | is.infinite(log.density)] <- -Inf
  if (!log.d)
    log.density <- exp(log.density)
  log.density
}


#'@rdname gpd
#'@export
qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE)
{
  if (log.p)
    p <- exp(p)
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 1)
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) <= 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  if (lower.tail)
    p <- 1 - p
  if (shape == 0) loc - scale * log(p)
  else {
    gpd.stand <- expm1(-shape * log(p)) / shape
    loc + scale * gpd.stand
  }
}
