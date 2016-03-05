#' Parameter estimation for the GEVr distribution model
#'
#' This function provides maximum likelihood estimation for the GEVr model, with the option of probability weighted moment and maximum product
#' spacing estimation for block maxima (GEV1) data. It also allows generalized linear modeling of the parameters.
#' @param data Data should be a matrix from the GEVr distribution.
#' @param method Method of estimation - maximum likelihood (mle), maximum product spacings (mps), and probability weighted moments (pwm). Uses mle by default.
#' For \eqn{r > 1}, only mle can be used.
#' @param information Whether standard errors should be calculated via observed or expected (default) information. For probability weighted moments,
#' only expected information will be used if possible. In the case with covariates, only observed information is available.
#' @param locvars A dataframe of covariates to use for modeling of the location parameter. Parameter intercepts are automatically handled by the function. Defaults to NULL.
#' @param scalevars A dataframe of covariates to use for modeling of the scale parameter. Parameter intercepts are automatically handled by the function. Defaults to NULL.
#' @param shapevars A dataframe of covariates to use for modeling of the shape parameter. Parameter intercepts are automatically handled by the function. Defaults to NULL.
#' @param locform An object of class `formula' (or one that can be coerced into that class), specifying the model of the location
#' parameter. If NULL, assumes stationary (intercept only) model. See details.
#' @param scaleform An object of class `formula' (or one that can be coerced into that class), specifying the model of the scale
#' parameter. If NULL, assumes stationary (intercept only) model. See details.
#' @param shapeform An object of class `formula' (or one that can be coerced into that class), specifying the model of the shape
#' parameter. If NULL, assumes stationary (intercept only) model. See details.
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
#' covs <- as.data.frame(seq(1, n, 1))
#' names(covs) <- c("Trend")
#' result2 <- gevrFit(data = x2, method = "mle", locvars = covs, locform = ~ Trend)
#'
#' ## Calculate p-values for the parameters
#' 2 * pnorm(abs(result2$par.ests / result2$par.ses), lower.tail = FALSE)
#'
#' @return A list describing the fit, including parameter estimates and standard errors for the mle and mps methods. Returns as a class
#' object 'gevrFit' to be used with diagnostic plots.
#' @details In the stationary case (no covariates), starting parameters for mle and mps estimation are the probability weighted moment estimates.
#' In the case where covariates are used, the starting intercept parameters are the probability weighted moment estimates from the stationary case
#' and the parameters based on covariates are initially set to zero. For non-stationary parameters, the first reported estimate refers to the
#' intercept term. \cr
#' Formulas for generalized linear modeling of the parameters should be given in the form `~ var1 + var2 + $\eqn{\ldots}'. Essentially, specification
#' here is the same as would be if using function `lm' for only the right hand side of the equation. Interactions, polynomials, etc. can be
#' handled as in the `formula' class. \cr
#' Intercept terms are automatically handled by the function. By default, the link functions are the identity function and the covariate dependent
#' scale parameter estimates are forced to be positive. For some link function \eqn{f(\cdot)} and for example, location parameter \eqn{\mu}, the
#' link is written as \eqn{\mu = f(\mu_1 x_1 + \mu_2 x_2 + \ldots + \mu_k x_k)}. \cr
#' Maximum likelihood estimation can be used in all cases. Probability weighted moment estimation can only be used if \eqn{r = 1} and data is
#' assumed to be stationary. Maximum product spacings estimation can be used in the non-stationary case, but only if \eqn{r = 1}.
#'
#' @import stats graphics
#' @export
gevrFit <- function(data, method = c("mle", "mps", "pwm"), information = c("expected", "observed"), locvars = NULL, scalevars = NULL,
                    shapevars = NULL, locform = NULL, scaleform = NULL, shapeform = NULL, loclink = identity, scalelink = identity,
                    shapelink = identity, start = NULL, opt = "Nelder-Mead", maxit = 10000, ...) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  information <- match.arg(information)
  if(!is.null(locvars))
    if(nrow(locvars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(scalevars))
    if(nrow(scalevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if(!is.null(shapevars))
    if(nrow(shapevars) != n)
      stop("Dimension of covariates does not match dimension of responses!")
  if((!is.null(locform) & is.null(locvars)) | (!is.null(scaleform) & is.null(scalevars)) | (!is.null(shapeform) & is.null(shapevars)))
    stop("Need to specify covariates!")
  if(R > 1 & (method == "mps" | method == "pwm"))
    stop("If R > 1, MLE must be used")
  if((!is.null(locvars) | !is.null(scalevars) | !is.null(shapevars)) & method == "pwm")
    stop("Probability weighted moments can only be fitted for stationary data")

  locvars.model <- as.matrix(rep(1, n))
  scalevars.model <- as.matrix(rep(1, n))
  shapevars.model <- as.matrix(rep(1, n))

  locnames <- colnames(locvars.model)
  scalenames <- colnames(scalevars.model)
  shapenames <- colnames(shapevars.model)

  loctrans1 <- 0
  scaletrans1 <- 0
  shapetrans1 <- 0

  loctrans2 <- 1
  scaletrans2 <- 1
  shapetrans2 <- 1

  if(!is.null(locform)) {
    locvars.model <- model.matrix(locform, data = locvars)
    locnames <- colnames(locvars.model)[2:ncol(locvars.model)]
    locvars.model <- scale(model.matrix(locform, data = locvars)[, -1])
    loctrans1 <- c(loctrans1, attr(locvars.model, "scaled:center"))
    loctrans2 <- c(loctrans2, attr(locvars.model, "scaled:scale"))
    locvars.model <- cbind(rep(1, n), locvars.model)
  }

  if(!is.null(scaleform)) {
    scalevars.model <- model.matrix(scaleform, data = scalevars)
    scalenames <- colnames(scalevars.model)[2:ncol(scalevars.model)]
    scalevars.model <- scale(model.matrix(scaleform, data = scalevars)[, -1])
    scaletrans1 <- c(scaletrans1, attr(scalevars.model, "scaled:center"))
    scaletrans2 <- c(scaletrans2, attr(scalevars.model, "scaled:scale"))
    scalevars.model <- cbind(rep(1, n), scalevars.model)
  }

  if(!is.null(shapeform)) {
    shapevars.model <- model.matrix(shapeform, data = shapevars)
    shapenames <- colnames(shapevars.model)[2:ncol(shapevars.model)]
    shapevars.model <- scale(model.matrix(shapeform, data = shapevars)[, -1])
    shapetrans1 <- c(shapetrans1, attr(shapevars.model, "scaled:center"))
    shapetrans2 <- c(shapetrans2, attr(shapevars.model, "scaled:scale"))
    shapevars.model <- cbind(rep(1, n), shapevars.model)
  }

  trans1 <- c(loctrans1, scaletrans1, shapetrans1)
  trans2 <- c(loctrans2, scaletrans2, shapetrans2)

  locvars.model1 <- t((t(locvars.model) * loctrans2) + loctrans1)
  scalevars.model1 <- t((t(scalevars.model) * scaletrans2) + scaletrans1)
  shapevars.model1 <- t((t(shapevars.model) * shapetrans2) + shapetrans1)

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
    locinit <- c(loc0, rep(0, ncol(locvars.model) - 1))
    scaleinit <- c(scale0, rep(0, ncol(scalevars.model) - 1))
    shapeinit <- c(shape0, rep(0, ncol(shapevars.model) - 1))
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
    negloglik <- function(vars, locvars1, scalevars1, shapevars1) {
      loc <- vars[1:length(locinit)]
      scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
      shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)]

      locmat <- t(loc * t(locvars1))
      scalemat <- t(scale * t(scalevars1))
      shapemat <- t(shape * t(shapevars1))

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

    init.fit <- optim(init, negloglik, hessian = TRUE, method = opt, control = list(maxit = maxit, ...),
                      locvars1 = locvars.model, scalevars1 = scalevars.model, shapevars1 = shapevars.model)

    fit <- optim(init.fit$par / trans2, negloglik, hessian = TRUE, method = opt, control = list(maxit = maxit, ...),
                 locvars1 = locvars.model1, scalevars1 = scalevars.model1, shapevars1 = shapevars.model1)

    if(fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed" | (!is.null(locvars) | !is.null(scalevars) | !is.null(shapevars))) {
      varcov <- unname(solve(fit$hessian), force = TRUE)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = data, type = "mle",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, nllh.final = fit$value, R = R,
                stationary = (is.null(locvars) & is.null(scalevars) & is.null(shapevars)))

    names(out$par.ests) <- c("Location Intercept", locnames, "Scale Intercept", scalenames, "Shape Intercept", shapenames)
    names(out$par.ses) <- names(out$par.ests)
  }

  if(method == "mps") {
    mpsobj <- function(vars, locvars1, scalevars1, shapevars1) {
      loc <- vars[1:length(locinit)]
      scale <- vars[(length(locinit) + 1):(length(locinit) + length(scaleinit))]
      shape <- vars[(length(locinit) + length(scaleinit) + 1):length(vars)]

      locmat <- t(loc * t(locvars1))
      scalemat <- t(scale * t(scalevars1))
      shapemat <- t(shape * t(shapevars1))

      locvec <- loclink(rowSums(locmat))
      scalevec <- scalelink(rowSums(scalemat))
      shapevec <- shapelink(rowSums(shapemat))

      w <- as.vector((data - locvec) / scalevec)
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

    init.fit <- optim(init, mpsobj, hessian = TRUE, method = opt, control = list(maxit = maxit, ...),
                      locvars1 = locvars.model, scalevars1 = scalevars.model, shapevars1 = shapevars.model)

    fit <- optim(init.fit$par / trans2, mpsobj, hessian = TRUE, method = opt, control = list(maxit = maxit, ...),
                 locvars1 = locvars.model1, scalevars1 = scalevars.model1, shapevars1 = shapevars.model1)

    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    if(information == "observed" | (!is.null(locvars) | !is.null(scalevars) | !is.null(shapevars))) {
      varcov <- unname(solve(fit$hessian), force = TRUE)
    } else {
      varcov <- gevrFisher(data, par.ests) / n
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = as.matrix(data), type = "mps",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, moran = fit$value, R = R,
                stationary = (is.null(locvars) & is.null(scalevars) & is.null(shapevars)))

    names(out$par.ests) <- c("Location Intercept", locnames, "Scale Intercept", scalenames, "Shape Intercept", shapenames)
    names(out$par.ses) <- names(out$par.ests)
  }
  class(out) <- "gevrFit"
  out
}
