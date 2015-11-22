#' GPD Shape Parameter Profile Likelihood
#'
#' Computes the profile likelihood based confidence interval for the shape parameter of the generalized pareto model.
#'
#' @param z A class object returned from gpdFit.
#' @param conf Confidence level to use. Defaults to 95 percent.
#' @param opt Optimization method to maximize the profile likelihood if that is selected. The default method is Nelder-Mead.
#'
#' @examples
#' x <- rgpd(500, loc = 0, scale = 1, shape = 0.4)
#' z <- gpdFit(x, threshold = 0)
#' gpdProfShape(z)
#' @return Estimate Estimated shape parameter.
#' @return CI Profile likelihood based confidence interval for the shape parameter.
#' @return ConfLevel The confidence level used.
#' @export
gpdProfShape <- function(z, conf = .95, opt = c("Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B", "Brent")) {
  data <- z$data
  threshold <- z$threshold
  theta <- as.numeric(z$par.ests)
  opt <- match.arg(opt)
  sol <- theta[1]
  gpdLikShape <- function(a, sh) {
    if(a[1] <= 0) {
      out <- .Machine$double.xmax
    } else {
      out <- dgpd(data[data > threshold], loc = threshold, scale = a[1], shape = sh, log.d = TRUE)
      out <- - sum(out)
      if(out == Inf)
        out <- .Machine$double.xmax
    }
    out
  }
  cutoff <- qchisq(conf, 1)
  prof <- function(sh) {
    lmax <- dgpd(data[data > threshold], loc = threshold, scale = theta[1], shape = theta[2], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gpdLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    2*(lmax-lci) - cutoff
  }
  prof <- Vectorize(prof)
  suppressWarnings(out1 <- uniroot(prof, c(theta[2] - 1e-6, theta[2]), extendInt="downX"))
  suppressWarnings(out2 <- uniroot(prof, c(theta[2], theta[2] + 1e-6), extendInt="upX"))
  CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  out <- list(theta[2], CI, conf)
  names(out) <- c("Estimate", "CI", "ConfLevel")
  out
}
