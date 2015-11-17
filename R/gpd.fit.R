findthresh <- function(data, ne) {
  data <- rev(sort(as.numeric(data)))
  thresholds <- unique(data)
  indices <- match(data[ne], thresholds)
  indices <- pmin(indices + 1, length(thresholds))
  thresholds[indices] - min(c(1e-10, abs(diff(data))))
}

#'Fits the generalized pareto distribution to data
#'
#'The base code for this function is taken from the R package evir. See citation below. The addition here includes estimation via Maximum Product Spacings (MPS).
#'@param data Data should be a numeric vector from the GPD distribution.
#'@param threshold A threshold value (either this or the number of extremes must be given, but not both).
#'@param nextremes Number of upper extremes to be used.
#'@param npp Length of each period (typically year). Is used in return level estimation. Defaults to 365.
#'@param information Whether standard errors should be calculated via observed or expected information. For probability weighted moments, only expected information will be used if possible.
#'@param method Method of estimation - maximum likelihood (mle), probability weighted moments (pwm), and maximum product spacings (mps). Uses mle by default.
#'@return A list describing the fit, including parameter estimates and standard errors.
#'@references Pfaff, Bernhard, Alexander McNeil, and A. Stephenson. "evir: Extreme Values in R." R package version (2012): 1-7.
#'@export

gpd.fit <- function(data, threshold = NA, nextremes = NA, npp = 365, method = c("mle", "mps", "pwm"),
                    information = c("observed", "expected")) {
  data <- as.numeric(data)
  n <- length(data)
  if(is.na(nextremes) && is.na(threshold))
    stop("Enter either a threshold or the number of upper extremes")
  if(!is.na(nextremes) && !is.na(threshold))
    stop("Enter EITHER a threshold or the number of upper extremes")
  if(!is.na(nextremes))
    threshold <- findthresh(data, nextremes)
  exceedances <- data[data > threshold]
  excess <- exceedances - threshold
  Nu <- length(excess)
  p.less.thresh <- 1 - Nu/n
  method <- match.arg(method)
  xbar <- mean(excess)

  a0 <- xbar
  gamma <- -0.35
  delta <- 0
  pvec <- ((1:Nu) + gamma)/(Nu + delta)
  a1 <- mean(sort(excess) * (1 - pvec))
  shape0 <- 2 - a0/(a0 - 2 * a1)
  scale0 <- (2 * a0 * a1)/(a0 - 2 * a1)
  start <- c(scale0, shape0)

  if(method == "pwm") {
    denom <- Nu * (1 - 2 * shape0) * (3 - 2 * shape0)
    if(shape0 > 0.5) {
      denom <- NA
      warning("Asymptotic standard errors not available for",
              "PWM Method when shape > 0.5")
    }
    one <- (7 - 18 * shape0 + 11 * shape0^2 - 2 * shape0^3) * scale0^2
    two <- (1 - shape0) * (1 - shape0 + 2 * shape0^2) * (2 - shape0)^2
    cov <- scale0 * (2 - shape0) * (2 - 6 * shape0 + 7 * shape0^2 - 2 *
                                shape0^3)
    varcov <- matrix(c(one, cov, cov, two), 2)/denom
    par.ses <- sqrt(diag(varcov))
    information <- "expected"
    out <- list(n = length(data), data = data, threshold = threshold,
                p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
                par.ests = start, par.ses = par.ses, varcov = varcov,
                information = information, npp = npp, rate = 1 - p.less.thresh)
  }

  if (method == "mle") {
    mle_est <- function(theta, dat) {
      scale <- theta[1]
      shape <- theta[2]
      cond1 <- scale <= 0
      cond2 <- (shape <= 0) && (max(dat) > (-scale/shape))
      if (cond1 || cond2) {
        out <- 1e+06
      } else {
        out <- - sum(dgpd(dat, loc = 0, scale = scale, shape = shape, log.d = TRUE))
      }
      out
    }
    fit <- optim(start, mle_est, hessian = TRUE, dat = excess)
    if(fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    converged <- fit$convergence
    nllh.final <- fit$value
    information <- match.arg(information)
    if(information == "observed")
      varcov <- solve(fit$hessian)
    if(information == "expected") {
      one <- (2 * (1 + par.ests[2]) * par.ests[1]^2)/Nu
      two <- (1 + par.ests[2])^2/Nu
      cov <- -((1 + par.ests[2]) * par.ests[1])/Nu
      varcov <- matrix(c(one, cov, cov, two), 2)
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = length(data), data = data, threshold = threshold,
                p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                information = information, converged = converged, nllh.final = nllh.final,
                npp = npp,  rate = 1 - p.less.thresh)
  }

  if (method == "mps") {
    mps_est <- function(theta, dat) {
      scale <- theta[1]
      shape <- theta[2]
      z <- 1 + shape*(dat / scale)
      cond1 <- scale <= 0
      cond2 <- (shape <= 0) && (max(dat) > (-scale/shape))
      if(cond1 || cond2) {
        out <- 1e+06
      } else {
        cdf <- pgpd(x, loc = 0, scale = scale, shape = shape)
        cdf <- sort(cdf)
        cdf <- c(0, cdf, 1)
        D <- diff(cdf)
        ## Check if any values are zero due to rounding and adjust
        len <- num.decimals.max(cdf)
        D <- ifelse(D==0, 1/(2*(10^len)), D)
        out <- - sum(log(D))
      }
      out
    }
    fit <- optim(start, mps_est, hessian = TRUE, dat = excess)
    if (fit$convergence)
      warning("optimization may not have succeeded")
    par.ests <- fit$par
    converged <- fit$convergence
    nllh.final <- fit$value
    information <- match.arg(information)
    if(information == "observed") {
      varcov <- solve(fit$hessian)
    } else {
      one <- (2 * (1 + par.ests[2]) * par.ests[1]^2)/Nu
      two <- (1 + par.ests[2])^2/Nu
      cov <- -((1 + par.ests[2]) * par.ests[1])/Nu
      varcov <- matrix(c(one, cov, cov, two), 2)
    }
    par.ses <- sqrt(diag(varcov))
    out <- list(n = length(data), data = data, threshold = threshold,
                p.less.thresh = p.less.thresh, n.exceed = Nu, method = method,
                par.ests = par.ests, par.ses = par.ses, varcov = varcov,
                information = information, converged = converged, moran = fit$value,
                npp = npp, rate = 1 - p.less.thresh)
  }

  names(out$par.ests) <- c("Scale", "Shape")
  names(out$par.ses) <- c("Scale", "Shape")
  class(out) <- "gpd.fit"
  out
}

