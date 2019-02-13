#Problem 1 - code
library(fields)#Matern is here
library(geoR) #Need to install XQuartz. cov.spatial() belongs here
library(MASS)
library(akima)
library(ggplot2)
library(reshape2)#melt
library(gridExtra)#grid.arrange

source('problem1_functions.R')
set.seed(42)
###################################################################################
###################################################################################
#a) Correlation function and variograms
###################################################################################
#Plot the Matern correlation function, the Powered Exponential function
#And the correponding variograms
#This code is inspired by Henning Omre's "Matern1.R"

#variance, for the variogram
sigma2 = c(1,5)

distance=0:0.01:49
nuMat=c(1,3)
CorrMatern1=Matern(distance, range = 10, nu=nuMat[1]) #What should range be?
CorrMatern2=Matern(distance, range = 10, nu=nuMat[2])

nuExp = c(1,1.9)
CorrExp1 = exp(-(distance/10)^nuExp[1])
CorrExp2 = exp(-(distance/10)^nuExp[2])


#plot(distance,CorrelationFunction,type="l")
#lines(distance,CorrelationFunction,type="l",col="blue")
CorrFuncDF = data.frame(distance, CorrMatern1, CorrMatern2,
                        CorrExp1, CorrExp2)
longDF = melt(CorrFuncDF, id = 'distance')

ggplot(longDF) + 
  geom_line(aes(x = distance, y = value, colour = variable)) + 
  scale_colour_discrete(labels = c("Matern\n nu = 1\n", "Matern\n nu = 3\n", 
                                   "Powered exp.\n nu = 1\n", "Powered exp.\n nu = 1.9\n")) +
  ggtitle("Correlation function") +
  theme(plot.title = element_text(hjust = 0.5))

#Save figure:
#ggsave("../figures/corrfunc.pdf", plot = last_plot(), device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 4, units = "in",
#       dpi = 300, limitsize = TRUE)

#####################################
#Plotting the variogram functions
plots = list()
spec = paste("Sigma^2 =", sigma2)

#####################################
plots[[1]] = ggplot(longDF) + 
  geom_line(aes(x = distance, y = sigma2[1]*(1-value), colour = variable)) + 
  scale_colour_discrete(name = spec[1], labels = 
                          c("Matern\n nu = 1\n", "Matern\n nu = 3\n", 
                            "Powered exp.\n nu = 1\n", "Powered exp.\n nu = 1.9\n")) +
  ggtitle("Variogram functions") +
  ylab('value') +
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]] = ggplot(longDF) + 
  geom_line(aes(x = distance, y = sigma2[2]*(1-value), colour = variable)) + 
  scale_colour_discrete(name = spec[2], labels = 
                          c("Matern\n nu = 1\n", "Matern\n nu = 3\n", 
                            "Powered exp.\n nu = 1\n", "Powered exp.\n nu = 1.9\n")) +
  #ggtitle("Variogram functions") +
  ylab('value') +
  theme(plot.title = element_text(hjust = 0.5))
#######################################

#Save plot in one figure
plotGrid = grid.arrange(grobs = plots, ncol = 1)

#ggsave("../figures/variograms.pdf", plot = plotGrid, device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 2*4, units = "in",
#       dpi = 300, limitsize = TRUE)

###################################################################################
###################################################################################
#b) Prior distribution
###################################################################################
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

#Covariance matrix
tau = rdist(L,L)/10

############################
#Make vectors of parameter values to iterate through
sigmaLong = rep(sigma2, times = 4)
nuLong = c(rep(nu_e, each = 2), rep(nu_m, each = 2))
funcLong = c(rep("e", 4), rep("m", 4))
f = ifelse(funcLong=="e", "Pow Exp: ", "Matern: ")
specifications = paste(f,"nu=",nuLong,", sigma=",sigmaLong, sep = "")

#############################
#Make samples and display them
nsamps = 10

##############################
#Choosing parameters for the "True system"
#Unsure whether this is the right way to do it
sigma_true = 5
func_true = "e"
nu_true = nu_e[2]
which_save = nsamps

###############################
#plots = list()
for(i in 1:8){
  Cov = covMatPrior(tau,sigmaLong[i], nuLong[i], funcLong[i])
  realizations = mvrnorm(nsamps, E, Cov)
  ############################
  #Save one realization for sigma^2 = 5 and the others as specified
  if(sigmaLong[i]==sigma_true & nuLong[i] == nu_true & funcLong[i]==func_true){
    savedRealization = realizations[which_save, ]
    save(savedRealization, file = "savedRealization.RData", ascii=TRUE)
  }
  ############################
  #Displaying:
  sampleDF = as.data.frame(t(realizations))
  sampleDF$L = L
  #Changing format to be able to plot all realisations
  long_realis = melt(sampleDF, id = "L")
  
  plot= ggplot(long_realis,
                aes(x=L, y=value, colour=variable)) +
    geom_line()+
    ggtitle("Realizations of 1D Gaussian RF") +
    xlab("x")+
    ylab("Realizations") +
    annotate("text", x = 10, y = max(long_realis$value)*1.1, label = specifications[i])+ 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  ########################
  #Save plots
  name = paste("../figures/sample1conf",i,".pdf", sep = "")
  #ggsave(name, plot = plot, device = NULL, path = NULL,
  #     scale = 1, width = 5.5, height = 4, units = "in",
  #   dpi = 300, limitsize = TRUE)
  
  print(plot)
  #plots[[i]] = plot
}

###################################################################################
###################################################################################
#c) Likelihood
###################################################################################

#Realization that is used to draw realizations from
load('savedRealization.RData')


sigmad = c(0,0.25)
x = c(10,25,30)
###################################################################################
###################################################################################
#d) Posterior distribution and prediction
###################################################################################

#Prior
sigmar = sigma2[2]
nu = nu_e[2]
func = "e"
meanPrior = E
covPrior = covMatPrior(tau, sigmar, nu, "e")

#Likelihood
H = obsMatrix(x, n)
covLik1 = diag(sigmad[1], ncol = length(x), nrow=length(x))
covLik2 = diag(sigmad[2], ncol = length(x), nrow=length(x))

#Draw measurements
set.seed(1)
d1 = simulateMeasurements(savedRealization, H, sigmad[1])
d2 = simulateMeasurements(savedRealization, H, sigmad[2])

#Compute posterior
posterior1 = postModel(meanPrior, covPrior, H, covLik1, d1)
posterior2 = postModel(meanPrior, covPrior, H, covLik2, d2)

################################
#Predicted values
predicted1 = posterior1$postMean
predicted2 = posterior2$postMean

#Prediction interval
alpha = 0.1
z = qnorm(alpha/2, lower.tail = FALSE)

#Got some rounding error giving negative variance...?
stdDev1 = sqrt(diag(abs(posterior1$postCov))) 

stdDev2 = sqrt(diag(posterior2$postCov))

intLower1 = predicted1-z*stdDev1
intUpper1 = predicted1+z*stdDev1

intLower2 = predicted2-z*stdDev2
intUpper2 = predicted2+z*stdDev2

prediction1 = data.frame(predicted1, intLower1, intUpper1)
prediction2 = data.frame(predicted2, intLower2, intUpper2)

#Display
specifications = paste("Measurement error = ", sigmad, sep = "")

plots = list()
plots[[1]] = ggplot(prediction1, aes(x = L)) + 
  geom_errorbar(aes(ymin = prediction1$intLower1, ymax = prediction1$intUpper1), col = "red")+
  geom_line(aes(y = prediction1$predicted1)) + 
  geom_point(data = as.data.frame(x), aes(x = x, y = d1)) + 
  ggtitle("Posterior predictions with intervals") +
  xlab("x")+
  ylab("Predicted values") +
  annotate("text", x = 9, y = max(prediction1$intUpper1)*1.1, label = specifications[1])+ 
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]] = ggplot(prediction2, aes(x = L)) + 
  geom_errorbar(aes(ymin = prediction2$intLower2, ymax = prediction2$intUpper2), col = "red")+
  geom_line(aes(y = prediction2$predicted2)) + 
  geom_point(data = as.data.frame(x), aes(x = x, y = d2)) + 
  xlab("x")+
  ylab("Predicted values") +
  annotate("text", x = 10, y = max(prediction2$intUpper2)*1.1, label = specifications[2])

plotGrid = grid.arrange(grobs = plots, nrow = 2)

#ggsave("../figures/predictions.pdf", plot = plotGrid, device = NULL, path = NULL,
#       scale = 1, width = 5.5, height = 2*4, units = "in",
#       dpi = 300, limitsize = TRUE)

###################################################################################
###################################################################################
#e) Simulations
###################################################################################
#Numerical approximation

nsamps = 100

postSamps1 = mvrnorm(nsamps, posterior1$postMean, posterior1$postCov)
postSamps2 = mvrnorm(nsamps, posterior2$postMean, posterior2$postCov)

predictedNum1 = colMeans(postSamps1)
predictedNum2 = colMeans(postSamps2)

stdDevNum1 = sqrt((1/(nsamps -1))*colSums((sweep(postSamps1,2,predictedNum1))^2))
stdDevNum2 = sqrt((1/(nsamps -1))*colSums((sweep(postSamps2,2,predictedNum2))^2))

#How is the prediction interval now?
#t-distributed?
t = qt(alpha/2, df = nsamps-1, lower.tail = FALSE)

predictionNum1 = data.frame(pred = predictedNum1, 
                            lower = predictedNum1 - t*stdDevNum1*sqrt(1+1/nsamps),
                            upper = predictedNum1 + t*stdDevNum1*sqrt(1+1/nsamps))

predictionNum2 = data.frame(pred = predictedNum2, 
                            lower = predictedNum2 - t*stdDevNum2*sqrt(1+1/nsamps),
                            upper = predictedNum2 + t*stdDevNum2*sqrt(1+1/nsamps))

#Displaying:

sampleDF1 = as.data.frame(t(postSamps1))
sampleDF1$L = L
#Changing format to be able to plot all realisations
long_realis1 = melt(sampleDF1, id = "L")

sampleDF2 = as.data.frame(t(postSamps2))
sampleDF2$L = L
#Changing format to be able to plot all realisations
long_realis2 = melt(sampleDF2, id = "L")

plots= list()

plots[[1]] = 
  ggplot(long_realis1,aes(x=L, y=value, colour = variable)) +
  geom_line()+
  geom_errorbar(data = predictionNum1,aes(x = L, ymin=predictionNum1$lower, 
                                          ymax = predictionNum1$upper),inherit.aes = FALSE, col = "grey30")+
  geom_line(data = predictionNum1, aes(x = L, y = predictionNum1$pred), inherit.aes = FALSE)+
  ggtitle("Realizations and empirical estimated predictions")+
  xlab("x")+
  ylab("Posterior values") +
  annotate("text", x = 10, y = max(long_realis1$value)*1.1, label = specifications[1])+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

plots[[2]]= 
  ggplot(long_realis2,aes(x=L, y=value, colour = variable)) +
  geom_line()+
  geom_errorbar(data = predictionNum2,aes(x = L, ymin=predictionNum2$lower, 
                                          ymax = predictionNum2$upper),inherit.aes = FALSE, col = "grey30")+
  geom_line(data = predictionNum2, aes(x = L, y = predictionNum2$pred), inherit.aes = FALSE)+
  xlab("x")+
  ylab("Posterior values") +
  annotate("text", x = 10, y = max(long_realis2$value)*1.1, label = specifications[2])+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")


plotGrid = grid.arrange(grobs = plots, nrow = 2)

#ggsave("../figures/posteriorSamps.pdf", plot = plotGrid, device = NULL, path = NULL,
#      scale = 1, width = 5.5, height = 2*4, units = "in",
#      dpi = 300, limitsize = TRUE)


###################################################################################
###################################################################################
#f) Function predictors
###################################################################################
#Using realizations without measurement errors.
realizations = postSamps1
N = nrow(realizations)
Ahat = mean(rowSums(realizations>2))

#Need prediction variance
varA = (1/(N-1))*sum((rowSums(realizations>2)-Ahat)^2)
#var(rowSums(realizations>2) gives the same :) 
varAhat = varA/100

#Using posterior Mean
Atilde = sum(posterior1$postMean>2)
