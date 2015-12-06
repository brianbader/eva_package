
## Functions to generate system data
# setwd("E:/MyPapers/GPD Tests/R Scripts")
# ADQuantiles <- read.csv("AD_Quantiles.csv", header = TRUE)
# CVMQuantiles <- read.csv("CVM_Quantiles.csv", header = TRUE)
# ADQuantiles <- ADQuantiles[, -1]
# CVMQuantiles <- CVMQuantiles[, -1]
# colnames(ADQuantiles) <- round(rev(seq(.001, .999, .001)), 3)
# colnames(CVMQuantiles) <- round(rev(seq(.001, .999, .001)), 3)
# rownames(ADQuantiles) <- round(seq( -0.5, 1, 0.01), 2)
# rownames(CVMQuantiles) <- round(seq( -0.5, 1, 0.01), 2)
# setwd("E:/MyPapers/EVA Package")
# devtools::use_data(ADQuantiles, CVMQuantiles, internal = TRUE)


## S3 functions for class gpdFit
plot.gpdFit <- function(x, ...) {
  gpdDiag(x, ...)
}


summary.gpdFit <- function(object, ...) {
  object$data <- NULL
  object
}


## Returns matrix of indicators needed for the GPD IM Test
gpdInd <- function(data, theta) {
  scale <- theta[1]
  shape <- theta[2]
  n <- length(data)
  w <- 1 + (1/shape)
  z <- 1 + (shape*data/scale)
  p12 <- - data/(shape*(scale^2)*z) + (w*data)/((scale^2)*z) - (w*(data^2)*shape)/((scale^3)*(z^2))
  p22 <- - (2*log(z))/(shape^3) + (2*data)/((shape^2)*scale*z) + (w*(data^2))/((scale^2)*(z^2))
  p11 <- (1/(scale^2)) - (2*w*shape*data)/((scale^3)*z) + (w*(shape^2)*(data^2))/((scale^4)*(z^2))
  p2 <- log(z)/(shape^2) - (w*data)/(scale*z)
  p1 <- - (1/scale) + (w*shape*data)/((scale^2)*z)
  D11 <- p1*p1 + p11
  D12 <- p1*p2 + p12
  D22 <- p2*p2 + p22
  D <- matrix(0, n, 3)
  D[ , 1] <- p1*p1 + p11
  D[ , 2] <- p1*p2 + p12
  D[ , 3] <- p2*p2 + p22
  D
}


## Feed a 3 by 3 matrix in, outputs the inverse (done analytically)
solve3by3 <- function(V) {
  a <- V[1,1]
  b <- V[1,2]
  c <- V[1,3]
  d <- V[2,1]
  e <- V[2,2]
  f <- V[2,3]
  g <- V[3,1]
  h <- V[3,2]
  i <- V[3,3]
  W <- matrix(0, 3, 3)
  W[1,1] <- e*i - f*h
  W[1,2] <- c*h - b*i
  W[1,3] <- b*f - c*e
  W[2,1] <- f*g - d*i
  W[2,2] <- a*i - c*g
  W[2,3] <- c*d - a*f
  W[3,1] <- d*h - e*g
  W[3,2] <- b*g - a*h
  W[3,3] <- a*e - b*d
  m <- a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g)
  W <- W / m
  W
}

## Helper function for gpd.imcov
gpdImCovGen <- function(n, theta) {
  scale <- theta[1]
  shape <- theta[2]
  y <- rgpd(n, loc = 0, scale = scale, shape = shape)
  fit1 <- tryCatch(gpdFit(y, nextremes = n, method = "mle", information = "expected"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit1)) {
    temp <- rep(NA, 3)
  } else {
    scale1 <- fit1$par.ests[1]
    shape1 <- fit1$par.ests[2]
    theta1 <- c(scale1, shape1)
    y <- y - findthresh(y, n)
    D1 <- gpdInd(y, theta1)
    D1 <- colSums(D1) / sqrt(n)
    temp <- D1
  }
temp
}


## Function returns GPD bootstrapped indicator covariance matrix
gpdImCov<- function(data, B, theta) {
  n <- length(data)
  temp <- t(replicate(B, gpdImCovGen(n, theta)))
  temp <- temp[complete.cases(temp), ]
  B <- nrow(temp)
  Dbar <- colMeans(temp)
  temp[,1] <- temp[,1] - Dbar[1]
  temp[,2] <- temp[,2] - Dbar[2]
  temp[,3] <- temp[,3] - Dbar[3]
  V <- (1/(B-1)) * t(temp) %*% temp
  V <- solve3by3(V)
  out <- list(V, B)
  names(out) <- c("cov", "boot_adj")
  out
}


## Returns expected inverse fisher information matrix
gpdFisher <- function(data, theta) {
  scale <- theta[1]
  shape <- theta[2]
  n <- length(data)
  one <- (2 * (1 + shape) * scale^2)/n
  two <- (1 + shape)^2/n
  cov <- -((1 + shape) * scale)/n
  varcov <- matrix(c(one, cov, cov, two), 2)
  varcov
}


## Returns observed inverse fisher information matrix
gpdFisherObs <- function(data, theta) {
  scale <- theta[1]
  shape <- theta[2]
  w <- 1 + (1/shape)
  z <- 1 + (shape*data/scale)
  p12 <- - data/(shape*(scale^2)*z) + (w*data)/((scale^2)*z) - (w*(data^2)*shape)/((scale^3)*(z^2))
  p22 <- - (2*log(z))/(shape^3) + (2*data)/((shape^2)*scale*z) + (w*(data^2))/((scale^2)*(z^2))
  p11 <- (1/(scale^2)) - (2*w*shape*data)/((scale^3)*z) + (w*(shape^2)*(data^2))/((scale^4)*(z^2))
  p11 <- sum(p11)
  p12 <- sum(p12)
  p22 <- sum(p22)
  varcov <- matrix(c(p11, p12, p12, p22), 2)
  varcov <- - varcov
  varcov <- solve(varcov)
  varcov
}


## Outputs matrix with row contributions to score. Need to sum over the columns to get full score.
gpdScore <- function(data, theta) {
  scale <- theta[1]
  shape <- theta[2]
  n <- length(data)
  w <- 1 + (1/shape)
  z <- 1 + (shape*data/scale)
  p1 <- - (1/scale) + (w*shape*data)/((scale^2)*z)
  p2 <- log(z)/(shape^2) - (w*data)/(scale*z)
  w <- matrix(0, n, 2)
  w[, 1] <- p1
  w[, 2] <- p2
  w
}


## Returns test statistic for the score test
gpdTestStat <- function(data, theta, information) {
  n <- length(data)
  w <- gpdScore(data, theta)
  w <- colSums(w)
  if(information == "observed") {
    info <- gpdFisherObs(data, theta)
  } else {
    info <- gpdFisher(data, theta)
  }
  score <- t(w) %*% info %*% w
  score <- as.vector(score)
  score
}
