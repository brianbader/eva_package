#'GEVr Return Level Estimate and Confidence Interval
#'
#'Computes m-period return level estimates, using either the delta method or profile likelihood. Can be for block maxima (GEV1) or GEVr.
#'
#'@param data Data should be a numeric vector (GEV1) or matrix from the GEVr distribution.
#'@param period The number of periods to use for the return level.
#'@param conf Confidence level. Defaults to 95 percent.
#'@param method The method to compute the confidence interval - either delta method (default) or profile likelihood.
#'@param opt Optimization method to maximize the profile likelihood if that is selected. The default method is Nelder-Mead.
#'
#'@references http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8391/background/chapter2.pdf
#'@examples
#'data <- rgevr(100, 2, loc=0.5, scale=1, shape=0.3)
#'## Compute 50-period return level.
#'gevr.returnlevel(data, 50, method = "profile")
#'@return Estimate Estimated m-period return level.
#'@return CI Confidence interval for the m-period return level.
#'@return Period The period length used.
#'@details Caution: The profile likelihood optimization may be slow (on the order of minutes).
#'@importFrom ismev rlarg.fit
#'@export

gevr.returnlevel <- function(data, period, conf=.95, method=c("delta", "profile"),
                             opt=c("Nelder-Mead", "SANN", "BFGS", "CG", "L-BFGS-B", "Brent"))
{
  method <- match.arg(method)
  data <- as.matrix(data)
  y <- NA
  try(y <- rlarg.fit(data, show = FALSE), silent = TRUE)
  if (!is.list(y))
    stop("Maximum likelihood failed to converge at initial step")
  cov <- y$cov
  theta <- y$mle
  z <- -log1p(-(1/period))
  est <- theta[1] - (theta[2]/theta[3])*(-expm1(-theta[3]*log(z)))
  if(method == "delta") {
    del <- matrix(ncol=1, nrow=3)
    del[1,1] <- 1
    del[2,1] <- -((theta[3])^(-1))*(-expm1(-theta[3]*log(z)))
    del[3,1] <- ((theta[2])*(theta[3]^(-2))*(-expm1(-theta[3]*log(z)))) - ((theta[2])*((theta[3])^(-1))*((1+expm1(-theta[3]*log(z)))*log(z)))
    se <- sqrt(t(del) %*% cov %*% del)
    se <- as.vector(se)
    alpha <- (1-conf)/2
    lower <- est - qnorm(1-alpha)*se
    upper <- est + qnorm(1-alpha)*se
    CI <- c(lower, upper)
  }
  else{
    opt <- match.arg(opt)
    sol <- c(theta[2], theta[3])
    gevr.lik <- function(a, xp) {
      loc <- xp + (a[1]/a[2])*(-expm1(-a[2]*log(z)))
      if(a[1] <= 0){
       out <- 1e6
      }
      else{
        out <- dgevr(data, loc = loc, scale = a[1], shape = a[2], log.d = TRUE)
        out <- - sum(out)
        if(out == Inf) out <- 1e6
      }
      out
    }
    cutoff <- qchisq(conf, 1)
    prof <- function(xp){
      lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
      lmax <- sum(lmax)
      yes <- optim(sol, gevr.lik, method = opt, xp = xp)
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
