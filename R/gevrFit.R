#' Parameter estimation for the GEVr distribution model
#'
#' This function provides maximum likelihood estimation for the GEVr model, with the option of probability weighted moment and maximum product
#' spacing estimation for block maxima (GEV1) data. It also allows generalized linear modeling of the parameters.
#' @param data Data should be a matrix from the GEVr distribution.
#' @param method Method of estimation - maximum likelihood (mle), probability weighted moments (pwm), and maximum product spacings (mps). Uses mle by default.
#' For \eqn{r > 1}, only mle can be used.
#' @param information Whether standard errors should be calculated via observed or expected (default) information. For probability weighted moments,
#' only expected information will be used if possible.
#' @param covars A matrix or dataframe of covariates to use. Parameter intercepts are automatically handled by the function. It is highly recommended to center
#' and scale the covariates to ease the optimization process.
#' @param locvars A vector specifying the columns of the covariates to use for modeling of the location parameter. Defaults to NULL.
#' @param scalevars A vector specifying the columns of the covariates to use for modeling of the scale parameter. Defaults to NULL.
#' @param shapevars A vector specifying the columns of the covariates to use for modeling of the shape parameter. Defaults to NULL.
#' @param loclink A link function specifying the relationship between the covariates and location parameter. Defaults to the identity function.
#' @param scalelink A link function specifying the relationship between the covariates and scale parameter. Defaults to the identity function.
#' @param shapelink A link function specifying the relationship between the covariates and shape parameter. Defaults to the identity function.
#' @param start Option to provide a set of starting parameters to optim; a vector of location, scale, and shape, in that order. Otherwise, the routine attempts
#' to find good starting parameters. See details.
#' @param opt Optimization method to use with optim.
#' @param maxit Number of iterations to use in optimization, passed to optim. Defaults to 10,000.
#' @param ... Additional arguments to pass to optim.
#' @examples
#' set.seed(7)
#' x1 <- rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
#' result1 <- gevrFit(x1, method = "mps")
#'
#' ## A linear trend in the location parameter
#' n <- 80
#' r <- 5
#' x2 <- matrix(0, ncol = r, nrow = n)
#' for(i in 1:n) {
#'   x2[i, ] <- rgevr(1, r, loc = 100 + i/40, scale = 1, shape = 0)
#' }
#'
#' covs <- scale(as.matrix(seq(1, n, 1)))
#' result2 <- gevrFit(data = x2, covars = covs, method = "mle", locvars = c(1))
#'
#' @return A list describing the fit, including parameter estimates and standard errors for the mle and mps methods. Returns as a class
#' object 'gevrFit' to be used with diagnostic plots.
#' @details If covariates are used, it is highly recommended that they are centered and scaled, for example, using R function `scale'.
#' In the stationary case (no covariates), starting parameters for mle and mps estimation are the probability weighted moment estimates.
#' In the case where covariates are used, the starting intercept parameters are the probability weighted moment estimates from the stationary case
#' and the parameters based on covariates are initially set to zero. For non-stationary parameters, the first reported estimate refers to the
#' intercept term. \cr
#' Intercept terms are automatically handled by the function. By default, the link functions are the identity function and the covariate dependent
#' scale parameter estimates are forced to be positive. For some link function \eqn{f(\cdot)} and for example, location parameter \eqn{\mu}, the
#' link is written as \eqn{\mu = f(\mu_1 x_1 + \mu_2 x_2 + \ldots + \mu_k x_k)}. \cr
#' Maximum likelihood estimation can be used in all cases. Probability weighted moment estimation can only be used if \eqn{r = 1} and data is
#' assumed to be stationary. Maximum product spacings estimation can be used in the non-stationary case, but only if \eqn{r = 1}.
#'
#' @import stats graphics
#' @export
gevrFit <- function(data, method = c("mle", "mps", "pwm"), information = c("expected", "observed"), covars = NULL, locvars = NULL,
                    scalevars = NULL, shapevars = NULL, loclink = identity, scalelink = identity, shapelink = identity, start = NULL,
                    opt = "Nelder-Mead", maxit = 10000, ...) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if(R > 1 & (method == "mps" | method == "pwm"))
    stop("If R > 1, MLE must be used")
  if(!is.null(covars) & method == "pwm")
    stop("Probability weighted moments can only be fitted for stationary data")

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

  if(is.null(start)) {
    locinit <- c(loc0, rep(0, length(locvars)))
    scaleinit <- c(scale0, rep(0, length(scalevars)))
    shapeinit <- c(shape0, rep(0, length(shapevars)))
    init <- c(locinit, scaleinit, shapeinit)
  }

  if(method == "pwm") {
    out <- list(n = n, data = data, type = "pwm",
                par.ests = theta0, par.ses = NA, varcov = NA,
                converged = NA, nllh.final = NA, R = R,
                stationary = TRUE)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
  }

  if(method == "mle") {
    negloglik <- function(a) {
      loc <- a[1:length(locinit)]
      scale <- a[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
      shape <- a[(length(locinit) + length(scaleinit) + 1):length(a)]

      locmat <- t(loc * t(cbind(rep(1, nrow(data)), covars[, c(locvars)])))
      scalemat <- t(scale * t(cbind(rep(1, nrow(data)), covars[, c(scalevars)])))
      shapemat <- t(shape * t(cbind(rep(1, nrow(data)), covars[, c(shapevars)])))

      locvec <- loclink(rowSums(locmat))
      scalevec <- scalelink(rowSums(scalemat))
      shapevec <- shapelink(rowSums(shapemat))

      w <- matrix(((data - locvec) / scalevec), ncol = R)
      z <- matrix(w * shapevec, ncol = R)
      z <- pmax(z, -1)

      log.density <- ifelse(shapevec == 0, rowSums(-log1p(scalevec -1) - w) - exp(-w[,R]),
                            rowSums(-log1p(scalevec - 1) - ((1/shapevec) + 1) * log1p(z)) - exp((-1/shapevec) * log1p(z[,R])))
      log.density[is.nan(log.density) | is.infinite(log.density)] <- -Inf

      if(any(scalevec < 0)) {
        out <- .Machine$double.xmax
      } else {
        out <- - sum(log.density)
      }
      out
    }

    fit <- optim(init, negloglik, hessian = TRUE, method = opt, control = list(maxit = maxit, ...))
    if(fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed" || !is.null(covars)) {
      varcov <- solve(fit$hessian)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = data, type = "mle",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, nllh.final = fit$value, R = R,
                stationary = is.null(covars))

    if(is.null(locvars))   lab1 <- "Location" else lab1 <- paste("Location", seq(1, length(locvars) + 1, 1), sep = "")
    if(is.null(scalevars)) lab2 <- "Scale"    else lab2 <- paste("Scale", seq(1, length(scalevars) + 1, 1), sep = "")
    if(is.null(shapevars)) lab3 <- "Shape"    else lab3 <- paste("Shape", seq(1, length(shapevars) + 1, 1), sep = "")

    names(out$par.ests) <- c(lab1, lab2, lab3)
    names(out$par.ses) <- names(out$par.ests)
  }

  if(method == "mps") {
    mpsobj <- function(a) {
      loc <- a[1:length(locinit)]
      scale <- a[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
      shape <- a[(length(locinit) + length(scaleinit) + 1):length(a)]

      locmat <- t(loc * t(cbind(rep(1, nrow(data)), covars[, c(locvars)])))
      scalemat <- t(scale * t(cbind(rep(1, nrow(data)), covars[, c(scalevars)])))
      shapemat <- t(shape * t(cbind(rep(1, nrow(data)), covars[, c(shapevars)])))

      locvec <- loclink(rowSums(locmat))
      scalevec <- scalelink(rowSums(scalemat))
      shapevec <- shapelink(rowSums(shapemat))

      w <- (as.vector(data) - locvec) / scalevec
      z <- pmax(w * shapevec, -1)
      cdf <- ifelse(shapevec == 0, exp(-exp(-w)), exp(-exp((-1/shapevec)*log1p(z))))
      cdf <- sort(cdf)

      if(any(scalevec < 0)) {
        out <- .Machine$double.xmax
      } else {
        cdf <- c(0, cdf, 1)
        D <- diff(cdf)
        ## Check if any differences are zero due to rounding and adjust
        D <- ifelse(D == 0, .Machine$double.eps, D)
        out <- - sum(log1p(D - 1))
      }
      out
    }

    fit <- optim(init, mpsobj, hessian = TRUE, method = opt, control = list(maxit = maxit, ...))
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed" || !is.null(covars)) {
      varcov <- solve(fit$hessian)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = as.matrix(data), type = "mps",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, moran = fit$value, R = R,
                stationary = is.null(covars))

    if(is.null(locvars))   lab1 <- "Location" else lab1 <- paste("Location", seq(1, length(locvars) + 1, 1), sep = "")
    if(is.null(scalevars)) lab2 <- "Scale"    else lab2 <- paste("Scale", seq(1, length(scalevars) + 1, 1), sep = "")
    if(is.null(shapevars)) lab3 <- "Shape"    else lab3 <- paste("Shape", seq(1, length(shapevars) + 1, 1), sep = "")

    names(out$par.ests) <- c(lab1, lab2, lab3)
    names(out$par.ses) <- names(out$par.ests)
  }
  class(out) <- "gevrFit"
  out
}
