# Functions for building and evaluating forecasts

## dispersal
# D = matrix of distances between grid cells (should already be squared)
# tau = dispersal constant
dispersal = function(D, tau) { # D should be distance matrix already squared
  M = exp(-(D/tau))# output is dispersal matrix
  return(M)
}

## growth

# growth base
growthsimple = function(a, b, nt) {
  G = exp(a+b*nt) # output is growth matrix
  return(G)
}

# growth with climate covariates
growthClim = function(a0,a1,a2,b0,X,nt) {
  G = exp((a0+X[,1]*a1+X[,2]*a2)+b0*nt)
  return(G)
}

# growth with topographic covariates
growthTopo = function(a0,a1,a2,b0,b1,b2,X,nt) {
  G = exp((a0+X[,3]*a1+X[,4]*a2)+(b0+X[,3]*b1+X[,4]*b2)*nt)
  return(G)
}

# growth with topo and climate covariates
growthTopoClim = function(a0,a1,a2,a3,a4,b0,b1,b2,X,nt) {
  G = exp((a0+X[,1]*a1+X[,2]*a2+X[,3]*a3+X[,4]*a4)+(b0+X[,3]*b1+X[,4]*b2)*nt)
  return(G)
}

# rmse
rmsefunc <- function(pred,obs) {
  n <- length(pred)
  RMSE <- sqrt((sum((pred-obs)^2))/n)
  return(RMSE)
}

# bias
biasfunc <- function(pred,obs) {
  n <- length(pred)
  bias <- (sum(pred-obs))/n
  return(bias)
}

# density
densefunc <- function(pred,obs,sig_o) {
  dense <- exp(sum(dnorm(obs,pred,sig_o,log=T)))
  return(dense)
}

