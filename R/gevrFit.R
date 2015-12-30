#' Parameter estimation for the GEVr distribution model
#'
#' This function provides maximum likelihood estimation for the GEVr model, with the option of probability weighted moment and maximum product
#' spacing estimation for block maxima (GEV1) data.
#' @param data Data should be a matrix from the GEVr distribution.
#' @param method Method of estimation - maximum likelihood (mle), probability weighted moments (pwm), and maximum product spacings (mps). Uses mle by default. For \eqn{r > 1}, only mle can be used.
#' @param information Whether standard errors should be calculated via observed or expected (default) information. For probability weighted moments, only expected information will be used if possible.
#' @param start Option to provide a set of starting parameters to optim; a vector of location, scale, and shape, in that order. Otherwise, the routine attempts to find good starting parameters.
#' @examples
#' x <- rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
#' result <- gevrFit(x, method = "mps")
#' @return A list describing the fit, including parameter estimates and standard errors for the mle and mps methods. Returns as a class object 'gevrFit' to be used with diagnostic plots.
#' @import stats graphics
#' @export
gevrFit <- function (data, method = c("mle", "mps", "pwm"),
                     information = c("expected", "observed"), start = NULL) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if(R > 1 & (method == "mps" | method == "pwm"))
     stop("If R > 1, MLE must be used")
  ## Probability Weighted Moments
  ## Also use this as the intial estimates for other methods
  y <- function(x, w0, w1, w2) {
    (3^x - 1)/(2^x - 1) - (3 * w2 - w0)/(2 * w1 - w0)
  }
  nmom <- 3
  x <- rev(sort(as.vector(data[,1])))
  moments <- rep(0, nmom)
  moments[1] <- mean(x)
  for (i in 1:n) {
    weight <- 1/n
    for (j in 2:nmom) {
      weight <- weight * (n - i - j + 2)/(n - j + 1)
      moments[j] <- moments[j] + weight * x[i]
    }
  }
  w0 <- moments[1]
  w1 <- moments[2]
  w2 <- moments[3]
  shape0 <- uniroot(f = y, interval = c(-5, +5), w0 = w0, w1 = w1,
                w2 = w2)$root
  scale0 <- (2 * w1 - w0) * shape0/gamma(1 - shape0)/(2^shape0 - 1)
  loc0 <- w0 + scale0 * (1 - gamma(1 - shape0))/shape0
  theta0 <- c(loc0, scale0, shape0)

  if(is.null(start))
    start <- theta0

  if(method == "pwm") {
    out <- list(n = n, data = data, type = "pwm",
                par.ests = theta0, par.ses = NA, varcov = NA,
                converged = NA, nllh.final = NA, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
  }

  if(method == "mle") {
    negloglik <- function(theta, x) {
      loc <- theta[1]
      scale <- theta[2]
      shape <- theta[3]
      z <- (shape / scale) * (x - loc)
      if ((scale < 0) || (min(1+z) < 0))
        out <- .Machine$double.xmax
      else {
        out <- - sum(dgevr(x, loc = loc, scale = scale, shape = shape, log.d = TRUE))
      }
      out
    }
    fit <- optim(start, negloglik, hessian = TRUE, x = data)
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed") {
      varcov <- solve(fit$hessian)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = data, type = "mle",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, nllh.final = fit$value, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
    names(out$par.ses) <- c("Location", "Scale", "Shape")
  }

  if(method == "mps") {
    data <- sort(as.vector(data))
    negloglik <- function(theta, x) {
      loc <- theta[1]
      scale <- theta[2]
      shape <- theta[3]
      z <- (shape / scale) * (x - loc)
      if ((scale < 0) || (min(1+z) < 0))
        out <- .Machine$double.xmax
      else {
        cdf <- pgev(x, loc = loc, scale = scale, shape = shape)
        cdf <- c(0, cdf, 1)
        D <- diff(cdf)
        ## Check if any differences are zero due to rounding and adjust
        D <- ifelse(D == 0, .Machine$double.eps, D)
        out <- - sum(log(D))
      }
      out
    }
    fit <- optim(start, negloglik, hessian = TRUE, x = data)
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed") {
      varcov <- solve(fit$hessian)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = as.matrix(data), type = "mps",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, moran = fit$value, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
    names(out$par.ses) <- c("Location", "Scale", "Shape")
  }
  class(out) <- "gevrFit"
  out
}
