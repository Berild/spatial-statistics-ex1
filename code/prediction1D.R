#Prediction and prediction interval 1D
library(gridExtra)

#Get posterior model
source('ComputePosterior1D.R')

########################################################
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
  annotate("text", x = 10, y = max(prediction1$intUpper1)*1.1, label = specifications[1])+ 
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]] = ggplot(prediction2, aes(x = L)) + 
  geom_errorbar(aes(ymin = prediction2$intLower2, ymax = prediction2$intUpper2), col = "red")+
  geom_line(aes(y = prediction2$predicted2)) + 
  geom_point(data = as.data.frame(x), aes(x = x, y = d2)) + 
  xlab("x")+
  ylab("Predicted values") +
  annotate("text", x = 10, y = max(prediction2$intUpper2)*1.1, label = specifications[2])

plotGrid = grid.arrange(grobs = plots, nrow = 2)

ggsave("../figures/predictions.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

############################################################
#Numerical approximation
library(MASS)
set.seed(42)
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
                            lower = predictedNum1 - t*stdDevNum1,
                            upper = predictedNum1 + t*stdDevNum1)

predictionNum2 = data.frame(pred = predictedNum2, 
                            lower = predictedNum2 - t*stdDevNum2,
                            upper = predictedNum2 + t*stdDevNum2)

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

ggsave("../figures/posteriorSamps.pdf", plot = plotGrid, device = NULL, path = NULL,
      scale = 1, width = 5.5, height = 2*4, units = "in",
      dpi = 300, limitsize = TRUE)