#Problem 1 - Functions
#######################################################
#For prior
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

#######################################################
#For likelihood
#Observation matrix
obsMatrix = function(x,n){
  H = matrix(0,ncol = n, nrow = length(x))
  for (i in 1:length(x)){
    H[i,x[i]] = 1
  }
  return(H)
}

#Draw measurements 
simulateMeasurements = function(r_true, obsMat, sigmaLik){
  m = nrow(obsMat)
  measurements = obsMat%*%r_true + sigmaLik*rnorm(m) 
  return(measurements)
}

#################################################################
#For posterior 

postModel = function(priorMean, priorCov, obsMat, covLik, d){
  C= solve(H %*%(priorCov)%*%t(H)+covLik)
  postMean=priorMean +priorCov%*% t(H) %*% C%*%(d-H%*%priorMean)
  postCov=priorCov-priorCov %*% t(H)%*% C%*%H%*%priorCov
  
  model = list(postMean = postMean, postCov = postCov)
  return(model)
}
