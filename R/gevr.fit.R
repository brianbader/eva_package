num.decimals.max <- function(x) {
  n <- length(x)
  out <- rep(0, n)
  for(i in 1:n) {
    if((x[i] %% 1) != 0)
      out[i] <- nchar(strsplit(sub('0+$', '', as.character(x[i])), ".", fixed=TRUE)[[1]][[2]])
  }
  max(out)
}

#'Fits generalized extreme value distribution (GEV) to block maxima data.
#'
#'@param data Data should be a numeric vector from the GEV distribution.
#'@param method Method of estimation - maximum likelihood (mle), probability weighted moments (pwm), and maximum product spacings (mps). Uses mle by default.
#'@examples
#'x <- rgevr(500, 1, loc = 0.5, scale = 1, shape = 0.3)
#'result <- gevr.fit(x, "mps")
#'@return A list describing the fit, including parameter estimates and standard errors for the mle and mps methods. Returns as a class object 'gevr.fit' to be used with diagnostic plots.
#'@export

gevr.fit <- function (data, method = c("mle", "mps", "pwm")) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  method <- match.arg(method)
  if(R > 1 & (method == "mps" | method == "pwm"))
     stop("If R > 1, MLE must be used")
  ## Probability Weighted Moments.
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
  shape <- uniroot(f = y, interval = c(-5, +5), w0 = w0, w1 = w1,
                w2 = w2)$root
  scale <- (2 * w1 - w0) * shape/gamma(1 - shape)/(2^shape - 1)
  loc <- w0 + scale * (1 - gamma(1 - shape))/shape
  theta <- c(loc, scale, shape)

  if(method == "pwm"){
    out <- list(n = n, data = data, type = "pwm",
                par.ests = theta, par.ses = NA, varcov = NA,
                converged = NA, nllh.final = NA, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
  }

  if(method == "mle"){
    negloglik <- function(theta, x) {
      loc <- theta[1]
      scale <- theta[2]
      shape <- theta[3]
      z <- (shape / scale) * (x - loc)
      if ((scale < 0) || (min(1+z) < 0))
        out <- 1e+06
      else {
        out <- - sum(dgevr(x, loc = loc, scale = scale, shape = shape, log.d = TRUE))
      }
      out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, x = data)
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = data, type = "mle",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, nllh.final = fit$value, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
    names(out$par.ses) <- c("Location", "Scale", "Shape")
  }

  if(method == "mps"){
    data <- sort(as.vector(data))
    negloglik <- function(theta, x) {
      loc <- theta[1]
      scale <- theta[2]
      shape <- theta[3]
      z <- (shape / scale) * (x - loc)
      if ((scale < 0) || (min(1+z) < 0))
        out <- 1e+06
      else {
        cdf <- pgev(x, loc = loc, scale = scale, shape = shape)
        cdf <- c(0, cdf, 1)
        D <- diff(cdf)
        ## Check if any values are zero due to rounding and adjust
        len <- num.decimals.max(cdf)
        D <- ifelse(D == 0, 1/(2*(10^len)), D)
        out <- - sum(log(D))
      }
      out
    }
    fit <- optim(theta, negloglik, hessian = TRUE, x = data)
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    varcov <- solve(fit$hessian)
    par.ses <- sqrt(diag(varcov))
    out <- list(n = n, data = as.matrix(data), type = "mps",
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                converged = fit$convergence, moran = fit$value, R = R)
    names(out$par.ests) <- c("Location", "Scale", "Shape")
    names(out$par.ses) <- c("Location", "Scale", "Shape")
  }
  class(out) <- "gevr.fit"
  out
}
