#Posterior Distribution 1D

library(fields) #rdist, Matern is here
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS) #mvrnorm is here
library(akima)

######################################################
#Make mean and covariance matrix for prior
source('SpecifyPrior1D.R')
sigmar = sigma2[2]
nu = nu_e[2]
func = "e"

meanPrior = E
covPrior = covMatPrior(tau, sigmar, nu, "e")

######################################################
#Likelihood

source('specifyLikelihood1D.R')

H = obsMatrix(x, n)
covLik1 = diag(sigmad[1], ncol = length(x), nrow=length(x))
covLik2 = diag(sigmad[2], ncol = length(x), nrow=length(x))

#Draw measurements
d1 = simulateMeasurements(savedRealization, H, sigmad[1])
d2 = simulateMeasurements(savedRealization, H, sigmad[2])

######################################################
#Mean and variance for posterior

postModel = function(priorMean, priorCov, obsMat, covLik, d){
  C= solve(H %*%(priorCov)%*%t(H)+covLik)
  postMean=priorMean +priorCov%*% t(H) %*% C%*%(d-H%*%priorMean)
  postCov=priorCov-priorCov %*% t(H)%*% C%*%H%*%priorCov
  
  model = list(postMean = postMean, postCov = postCov)
  return(model)
}

posterior1 = postModel(meanPrior, covPrior, H, covLik1, d1)
posterior2 = postModel(meanPrior, covPrior, H, covLik2, d2)



