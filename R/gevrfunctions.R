## Helper function to handle (1 + x*shape)^(-1/shape) as shape -> 0
nzsh <- function(x, shape) {
  exp((-1/shape) * log1p(x * shape))
}


## S3 functions for class gevrFit
plot.gevrFit <- function(x, ...) {
  gevrDiag(x, ...)
}


summary.gevrFit <- function(object, ...) {
  object$data <- NULL
  object
}


## Returns expected inverse fisher information matrix for GEV distribution.
gevrFisher <- function(data, theta) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  gr1 <- gamma(R + shape + 1) / gamma(R)
  gr2 <- gamma(R + 2*shape + 1) / gamma(R)
  A <- (((1+shape)^2)*gr2)/((scale^2)*(1+2*shape))
  B <- (-1/((scale^2)*shape*(1+2*shape)))*(((1+shape)^2)*gr2 - (1+2*shape)*gr1)
  C <- (1/(scale*(shape^2)*(1+2*shape)))*(((1+2*shape)*gr1)*(shape*digamma(R+shape+1) + (shape^2+shape+1)/(1+shape)) - ((1+shape)^2)*gr2)
  D <- (1/((scale^2)*(shape^2)*(1+2*shape)))*(R*(1+2*shape) - 2*(1+2*shape)*gr1 + ((1+shape)^2)*gr2)
  E <- (1/(scale*((-shape)^3)*(1+2*shape)))*(R*(-shape)*(1+2*shape)*digamma(R+1) + (1+2*shape)*gr1*(shape*digamma(R+shape+1) + (1+(1+shape)^2)/(1+shape)) - ((1+shape)^2)*gr2 - R*(1+2*shape))
  F <- (1/(((-shape)^4)*(1+2*shape)))*((2*(1+2*shape)*gr1)*((-shape)*digamma(R+shape+1) - (shape^2+shape+1)/(1+shape)) + ((1+shape)^2)*gr2 + (R*(1+2*shape))*(1 + 2*shape*digamma(R+1) + (shape^2)*(1 + trigamma(R+1) + ((digamma(R+1))^2))))
  info <- matrix(c(A, B, C, B, D, E, C, E, F), nrow=3, ncol=3)
  info <- solve3by3(info)
  info
}


## Returns observed inverse fisher information matrix for GEV distribution.
gevrFisherObs <- function(data, theta) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  z <- (data - loc) / scale
  A <- rowSums( -((1+shape)/scale)*((1+shape*z)^(-1)) - ((shape*(1+shape))/scale^2)*((1+shape*z)^(-2)) )
  A <- A + (1/scale) * nzsh(z[,R], shape)^(-1) + (1/scale^2) * (1+shape) * nzsh(z[,R], shape)^(-2)
  A <- sum(A)
  B <- rowSums( ((1+shape)/(scale^2))*((1+shape*z)^(-1)) - ((shape*(1+shape))/scale^2)*z*((1+shape*z)^(-2)) )
  B <- B - (1/(scale^2)) * nzsh(z[,R], shape)^(-1) + (1/scale^2) * (1+shape) * z[,R] * nzsh(z[,R], shape)^(-2)
  B <- sum(B)
  C <- rowSums( (scale^(-1))*((1+shape*z)^(-1)) - ((1+shape)/scale)*z*((1+shape*z)^(-2)) )
  C <- C - (1/scale) * (1/shape^2) * nzsh(z[,R], shape)^(-1) * log(1+shape*z[,R]) - ((1+shape)/(-1*scale*shape)) * z[,R] * nzsh(z[,R], shape)^(-2)
  C <- sum(C)
  D <- rowSums( (2*(1+shape)/scale^2)*(z/(1+shape*z)) - (shape*(1+shape)/scale^2)*((z/(1+shape*z))^2) )
  D <- D - (R/scale^2) - (2/scale^2) * z[,R] * nzsh(z[,R], shape)^(-1) +  ((1+shape)/scale^2) * z[,R]^2 * nzsh(z[,R], shape)^(-2)
  D <- sum(D)
  E <- rowSums( (1/scale)*(z/(1+shape*z)) - ((1+shape)/scale)*((z/(1+shape*z))^2) )
  E <- E - (1/scale) * (1/shape^2) * nzsh(z[,R], shape)^(-1) * z[,R] * log(1+shape*z[,R]) - ((1+shape)/(-1*scale*shape)) * z[,R]^2 * nzsh(z[,R], shape)^(-2)
  E <- sum(E)
  F <- rowSums( (2/shape^3)*log(1+shape*z) - (2/shape^2)*(z/(1+shape*z)) - ((1+shape)/shape)*((z/(1+shape*z))^2) )
  F <- F - (2/shape^3) * nzsh(z[,R], shape) * log(1+shape*z[,R]) + (2/shape^2) * z[,R] * nzsh(z[,R], shape)^(-1) + ((1+shape)/shape^2) * z[,R]^2 * nzsh(z[,R], shape)^(-2) + (1/shape^4) * nzsh(z[,R], shape) * log(1+shape*z[,R])^2  - (2/shape^3) * z[,R] * nzsh(z[,R], shape)^(-1) * log(1+shape*z[,R])
  F <- sum(F)
  info <- matrix(c(A, B, C, B, D, E, C, E, F), nrow = 3, ncol = 3)
  info <- (1/N) * info
  info <- solve3by3(info)
  info
}


## Outputs matrix with row contributions to score. Need to sum over the columns to get full score.
gevrScore <- function(data, theta) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  z <- (data - loc) / scale
  dLoc <- rowSums((((1/shape)+1)*shape) / (scale*(1+shape*z)))
  dLoc <- dLoc - nzsh(z[,R], shape) / (scale*(1+(shape*z[,R])))
  dScale <- rowSums((((1/shape)+1)*shape*z) / (scale*(1+shape*z)))
  dScale <- dScale - (R/scale) - z[,R] * nzsh(z[,R], shape) / (scale*(1+shape*z[,R]))
  dShape <-  rowSums(log(1+shape*z) / shape^2 - (((1/shape)+1)*z) / (1+shape*z))
  dShape <- dShape - nzsh(z[,R], shape) * ((log(1+shape*z[,R]) / shape^2) - z[,R] / (shape*(1+shape*z[,R])))
  score <- matrix(c(dLoc, dScale, dShape), ncol = 3)
  score
}


## Calculates the score test statistic
gevrTestStat <- function(data, theta, information) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  u <- gevrScore(data, theta)
  u <- colSums(u)
  if(information == "observed") {
    info <- gevrFisherObs(data, theta)
  } else {
    info <- gevrFisher(data, theta)
  }
  stat <- (1/N) * t(u) %*% info %*% u
  stat <- as.vector(stat)
  stat
}
