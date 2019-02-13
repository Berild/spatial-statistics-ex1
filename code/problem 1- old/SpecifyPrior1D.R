#Prior distribution 1D
######################################################################
library(fields)
#Grid
n = 50
L = 1:50

#Specify expected value
mu = 0
E = matrix(rep(mu,50), ncol = 1)

#Possible model parameters
sigma2 = c(1,5)
nu_m=c(1,3)
nu_e = c(1,1.9)

##########################################
#Covariance matrix
tau = rdist(L,L)/10

covMatPrior = function(tau, sigmasq, nu, func="m"){
  if (tolower(func)== "e"){
    corr = exp(-tau^nu)
  }else if(tolower(func) == "m"){
    corr = Matern(tau, range = 1, nu = nu) #What should range be?
  }else{
    warning("Unknown correlation function")
    return(NULL)
  }
  cov = sigmasq*corr
  return(cov)
}
