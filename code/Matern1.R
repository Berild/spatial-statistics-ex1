library(fields)
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS)
library(akima)
library(ggplot2)

distance=0:0.01:50/10
nuMat=seq(1,3, length = 3)
CorrelationFunctionMat=Matern(distance , range = 1,nu=nuMat[1]) #What should range be?

nuExp = seq(1,1.9, length = 3)
CorrelationFunctionExp = exp(-distance^nuExp[1])

#plot(distance,CorrelationFunction,type="l")
#lines(distance,CorrelationFunction,type="l",col="blue")

ggplot(data.frame(distance, CorrelationFunctionMat, CorrelationFunctionExp))+
  geom_path(aes(x = distance, y = CorrelationFunctionMat), col = 'blue') + 
  geom_path(aes(x = distance, y = CorrelationFunctionExp), col = 'red')


