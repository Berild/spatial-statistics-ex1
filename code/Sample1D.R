#Sample 1D
library(fields) #rdist, Matern is here
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS) #mvrnorm is here
library(akima)
library(ggplot2)
library(reshape2) #melt

#Grid
n = 50
L = 1:50


#Specify expected value
mu = 0
E = matrix(rep(mu,50), ncol = 1)

#Covariance matrix
tau = rdist(L,L)/10
sigma2 = c(1,5)

nuMat=c(1,3)
CorrelationFunctionMat=Matern(tau, range = 1, nu=nuMat[1]) #What should range be?

nuExp = c(1,1.9)
CorrelationFunctionExp = exp(-tau^nuExp[1])

Cov = sigma2[1]*CorrelationFunctionMat

#Make samples
sample = mvrnorm(10, E, Cov)


sampleDF = as.data.frame(t(sample))
sampleDF$L = L
#Changing format to be able to plot all realisations
long_sample = melt(data, id = "L")

#Display
ggplot(long_data,
       aes(x=L, y=value, colour=variable)) +
  geom_line()+
  ggtitle("Realizations of 1D Gaussian RF") +
  xlab("x")+
  ylab("sampled value") +
  #annotate("text", x = 5, y = max(long_data$value), label = "some text")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
