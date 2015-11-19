#' GPD Return Level Estimate and Confidence Interval
#'
#' Computes m-period return level estimates for the generalized pareto distribution, using either the delta method or profile likelihood.
#'
#' @param z An object of class gpd.fit.
#' @param period The number of periods to use for the return level.
#' @param conf Confidence level. Defaults to 95 percent.
#' @param method The method to compute the confidence interval - either delta method (default) or profile likelihood.
#' @param opt Optimization method to maximize the profile likelihood if that is selected. The default method is Nelder-Mead.
#'
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
#' @examples
#' x <- rgpd(5000, loc = 0, scale = 1, shape = 0.1)
#' ## Compute 50-period return level.
#' z <- gpd.fit(x, nextremes = 200)
#' gpd.rl(z, 50, method = "delta")
#' @return Estimate Estimated m-period return level.
#' @return CI Confidence interval for the m-period return level.
#' @return Period The period length used.
#' @details Caution: The profile likelihood optimization may be slow for large datasets.
#' @export
gpd.rl <- function(z, period, conf = .95, method = c("delta", "profile"),
                             opt = c("Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B", "Brent")) {
  method <- match.arg(method)
  m <- period * z$npp
  est <- z$threshold + (z$par.ests[1] / z$par.ests[2]) * ((m * z$rate)^z$par.ests[2] - 1)
  est <- as.numeric(est)
  if(method == "delta") {
    cov <- matrix(0, 3, 3)
    cov[2:3, 2:3] <- z$varcov
    cov[1, 1] <- z$rate * (1 - z$rate) / z$n
    del <- matrix(0, 3, 1)
    del[1, 1] <- z$par.ests[1] * (m^z$par.ests[2]) * (z$rate^(z$par.ests[2] - 1))
    del[2, 1] <- (1 / z$par.ests[2]) * ((m * z$rate)^z$par.ests[2] - 1)
    del[3, 1] <- (-z$par.ests[1] / z$par.ests[2]^2) * ((m * z$rate)^z$par.ests[2] - 1) +
      (z$par.ests[1] / z$par.ests[2]) * ((m * z$rate)^z$par.ests[2]) * log(m * z$rate)
    se <- sqrt(t(del) %*% cov %*% del)
    se <- as.vector(se)
    alpha <- (1 - conf) / 2
    lower <- est - qnorm(1-alpha)*se
    upper <- est + qnorm(1-alpha)*se
    CI <- as.numeric(c(lower, upper))
  } else {
    opt <- match.arg(opt)
    sol <- z$par.ests[2]
    gpd.lik <- function(shape, xp) {
      if(shape == 0) {
        scale <- (xp - z$threshold) / log(m * z$rate)
      } else {
        scale <- ((xp - z$threshold) * shape) / ((m * z$rate)^shape - 1)
      }
      if(scale <= 0) {
        out <- 1e6
      } else {
        out <- dgpd(z$data[z$data > z$threshold], loc = z$threshold, scale = scale, shape = shape, log.d = TRUE)
        out <- - sum(out)
        if(out == Inf)
          out <- 1e6
      }
      out
    }
    cutoff <- qchisq(conf, 1)
    prof <- function(xp) {
      lmax <- dgpd(z$data[z$data > z$threshold], loc = z$threshold, scale = z$par.ests[1], shape = z$par.ests[2], log.d = TRUE)
      lmax <- sum(lmax)
      yes <- optim(sol, gpd.lik, method = opt, xp = xp)
      sol <- yes$par
      lci <- -yes$value
      2*(lmax-lci) - cutoff
    }
    prof <- Vectorize(prof)
    suppressWarnings(out1 <- uniroot(prof, c(est - 1e-6, est), extendInt="downX"))
    suppressWarnings(out2 <- uniroot(prof, c(est, est + 1e-6), extendInt="upX"))
    CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  }
  out <- list(est, CI, period)
  names(out) <- c("Estimate", "CI", "Period")
  out
}
