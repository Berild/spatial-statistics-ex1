#Likelihood 1D

#Realization that is used to draw realizations from
load('savedRealization.RData')


sigmad = c(0,0.25)
x = c(10,25,30)

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
