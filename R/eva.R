#' eva: Extreme Value Analysis with Goodness-of-Fit Testing
#'
#' The focus of this package is to provide much needed automated diagnostic
#' tools (in the form of statistical hypothesis testing) to extreme value
#' models. Other useful functionality is efficient and user-friendly non-stationary
#' model fitting, profile likelihood confidence intervals, data generation
#' in the r-largest order statistics model (GEVr), and ordered p-value multiplicity
#' adjustments. Also, all routines are implemented to efficiently handle
#' the near-zero shape parameter, which may cause numerical issues in other
#' packages. Functions can be roughly assigned to the following topics:
#'
#' @section Formal (Automated) Goodness-of-Fit Testing:
#'
#' \code{\link{gevrSeqTests}} is a wrapper function that performs sequential testing
#' for r in the GEVr distribution, with adjusted p-values. It can implement three
#' tests:
#'
#' \code{\link{gevrEd}} An entropy difference test, which uses an asymptotic normal central limit theorem result.
#'
#' \code{\link{gevrPbScore}} A score test, implemented using parametric bootstrap and can be run in parallel.
#'
#' \code{\link{gevrMultScore}} An asymptotic approximation to the score test (computationally effciient).
#'
#' \code{\link{gpdSeqTests}} is a wrapper function that performs sequential testing
#' for thresholds in the Generalized Pareto distribution (GPD), with adjusted
#' p-values. It can implement the following tests:
#'
#' \code{\link{gpdAd}} The Anderson-Darling test, with log-linear interpolated p-values. Can also be bootstrapped
#' (with a parallel option).
#'
#' \code{\link{gpdCvm}} The Cramer-Von Mises test, with log-linear interpolated p-values. Can also be bootstrapped
#' (with a parallel option).
#'
#' \code{\link{gpdImAsym}} An asymptotic information matrix test, with bootstrapped covariance estimates.
#'
#' \code{\link{gpdImPb}} A full bootstrap version of information matrix test, with bootstrapped covariance estimates and critical values.
#'
#' \code{\link{gpdPbScore}} A score test, implemented using parametric bootstrap and can be run in parallel.
#'
#' \code{\link{gpdMultScore}} An asymptotic approximation to the score test (computationally effciient).
#'
#' \code{\link{pSeqStop}} A simple function that reads in raw, ordered p-values and returns two sets that adjust
#' for the familywise error rate and false discovery rate.
#'
#' @section Data generation and model fitting:
#'
#' All the functions in this section (and package) efficiently handle a near-zero
#' value of the shape parameter, which can cause numerical instability in similar
#' functions from other packages. See the vignette for an example.
#'
#' Data generation, density, quantile, and distribution functions can handle
#' non-stationarity and vectorized inputs.
#'
#' \code{\link{gevr}} Data generation and density function for the GEVr distribution,
#' with distribution function and quantile functions available for GEV1 (block maxima).
#'
#' \code{\link{gpd}} Data generation, distribution, quantile, and density functions
#' for the GPD distribution.
#'
#' \code{\link{gevrFit}} Non-stationary fitting of the GEVr distribution, with the option
#' of maximum product spacings estimation when r=1. Uses formula statements for
#' user friendliness and automatically centers/scales covariates when appropriate
#' to speed up optimization.
#'
#' \code{\link{gpdFit}} Non-stationary fitting of the GP distribution, with same
#' options and implementation as \code{\link{gevrFit}}. Allows non-stationary
#' threshold to be used.
#'
#' \code{\link{gevrProfShape}} Profile likelihood estimation for the shape
#' parameter of the stationary GEVr distribution.
#'
#' \code{\link{gpdProfShape}} Profile likelihood estimation for the shape
#' parameter of the stationary GP distribution.
#'
#' \code{\link{gevrRl}} Profile likelihood estimation for return levels
#' of the stationary GEVr distribution.
#'
#' \code{\link{gpdRl}} Profile likelihood estimation for return levels
#' of the stationary GP distribution.
#'
#' \code{\link{gevCIboot}} Bootstrapped standard errors for fitted GEV1
#' (block maxima) models, with options to resample under assumed dependence.
#'
#' @section Visual Diagnostics:
#'
#' \code{\link{gevrDiag}}, \code{\link{gpdDiag}} Diagnostic plots for a fit to
#' the GEVr (GP) distribution. For stationary models, return level, density, quantile,
#' and probability plots are returned. For non-stationary models, residual quantile,
#' residual probability, and residuals versus covariate plots are returned.
#'
#' \code{\link{mrlPlot}} Plots the empirical mean residual life, with
#' confidence intervals. Visual diagnostic tool to choose a threshold
#' for exceedances.
#'
#' @section Data:
#'
#' \code{\link{fortmax}} Top ten annual precipitation events (inches) for
#' one rain gauge in Fort Collins, Colorado from 1900 through 1999.
#'
#' \code{\link{lowestoft}} Top ten annual sea levels at the LoweStoft Station
#' tide gauge from 1964 - 2014.
#'
#' @import stats graphics
#' @docType package
#' @name eva
NULL
