## Returns expected inverse fisher information matrix for GEV distribution.
gevrfisher <- function(dat, theta) {
  R <- ncol(dat)
  N <- nrow(dat)
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
gevrfisherobs <- function(dat, theta) {
  R <- ncol(dat)
  N <- nrow(dat)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  z <- (dat - loc) / scale
  A <- rowSums( -((1+shape)/scale)*((1+shape*z)^(-1)) - ((shape*(1+shape))/scale^2)*((1+shape*z)^(-2)) )
  A <- A + (1/scale)*((1+shape*z[,R])^((-1/shape)-1)) + (1/scale^2)*(1+shape)*((1+shape*z[,R])^((-1/shape)-2))
  A <- sum(A)
  B <- rowSums( ((1+shape)/(scale^2))*((1+shape*z)^(-1)) - ((shape*(1+shape))/scale^2)*z*((1+shape*z)^(-2)) )
  B <- B - (1/(scale^2))*((1+shape*z[,R])^((-1/shape)-1)) + (1/scale^2)*(1+shape)*z[,R]*((1+shape*z[,R])^((-1/shape)-2))
  B <- sum(B)
  C <- rowSums( (scale^(-1))*((1+shape*z)^(-1)) - ((1+shape)/scale)*z*((1+shape*z)^(-2)) )
  C <- C - (1/scale)*(1/shape^2)*((1+shape*z[,R])^((-1/shape)-1))*log(1+shape*z[,R]) - ((1+shape)/(-1*scale*shape))*z[,R]*((1+shape*z[,R])^((-1/shape)-2))
  C <- sum(C)
  D <- rowSums( (2*(1+shape)/scale^2)*(z/(1+shape*z)) - (shape*(1+shape)/scale^2)*((z/(1+shape*z))^2) )
  D <- D - (R/scale^2) - (2/scale^2)*z[,R]*((1+shape*z[,R])^((-1/shape)-1)) +  ((1+shape)/scale^2)*(z[,R]^2)*((1+shape*z[,R])^((-1/shape)-2))
  D <- sum(D)
  E <- rowSums( (1/scale)*(z/(1+shape*z)) - ((1+shape)/scale)*((z/(1+shape*z))^2) )
  E <- E - (1/scale)*(1/shape^2)*((1+shape*z[,R])^((-1/shape)-1))*z[,R]*log(1+shape*z[,R]) - ((1+shape)/(-1*scale*shape))*(z[,R]^2)*((1+shape*z[,R])^((-1/shape)-2))
  E <- sum(E)
  F <- rowSums(  (2/shape^3)*log(1+shape*z) - (2/shape^2)*(z/(1+shape*z)) - ((1+shape)/shape)*((z/(1+shape*z))^2)  )
  F <- F - (2/shape^3)*((1+shape*z[,R])^(-1/shape))*log(1+shape*z[,R]) + (2/shape^2)*z[,R]*((1+shape*z[,R])^((-1/shape)-1)) + ((1+shape)/shape^2)*(z[,R]^2)*((1+shape*z[,R])^((-1/shape)-2)) + (1/shape^4)*((1+shape*z[,R])^(-1/shape))*((log(1+shape*z[,R]))^2)  - (2/shape^3)*z[,R]*((1+shape*z[,R])^((-1/shape)-1))*log(1+shape*z[,R])
  F <- sum(F)
  info <- matrix(c(A, B, C, B, D, E, C, E, F), nrow=3, ncol=3)
  info <- (1/N)*info
  info <- solve3by3(info)
  info
}

## Outputs matrix with row contributions to score. Need to sum over the columns to get full score.
gevrscorectb <- function(dat, theta) {
  R <- ncol(dat)
  N <- nrow(dat)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  z <- (dat - loc) / scale
  dLoc <- rowSums((((1/shape)+1)*shape) / (scale*(1+shape*z)))
  dLoc <- dLoc - ((1+(shape*z[,R]))^(-1/shape)) / (scale*(1+(shape*z[,R])))
  dScale <- rowSums((((1/shape)+1)*shape*z) / (scale*(1+shape*z)))
  dScale <- dScale - (R/scale) - z[,R]*((1+shape*z[,R])^(-1/shape)) / (scale*(1+shape*z[,R]))
  dShape <-  rowSums(log(1+shape*z) / shape^2 - (((1/shape)+1)*z) / (1+shape*z))
  dShape <- dShape - ((1+shape*z[,R])^(-1/shape))*((log(1+shape*z[,R]) / shape^2) - z[,R] / (shape*(1+shape*z[,R])))
  score <- matrix(c(dLoc, dScale, dShape), ncol=3)
  score
}

## Calculates the score test statistic
gevrteststat <- function(dat, theta, information) {
  R <- ncol(dat)
  N <- nrow(dat)
  u <- gevrscorectb(dat, theta)
  u <- colSums(u)
  if(information == "observed"){
    info <- gevrfisherobs(dat, theta)
  }
  else{
    info <- gevrfisher(dat, theta)
  }
  stat <- (1/N) * t(u) %*% info %*% u
  stat <- as.vector(stat)
  stat
}
