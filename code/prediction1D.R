#Prediction and prediction interval 1D
library(gridExtra)

#Get posterior model
source('ComputePosterior1D.R')


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
  ylab("predictions") +
  annotate("text", x = 10, y = max(prediction2$intUpper2)*1.1, label = specifications[2])

plotGrid = grid.arrange(grobs = plots, nrow = 2)

ggsave("../figures/predictions.pdf", plot = plotGrid, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 2*4, units = "in",
       dpi = 300, limitsize = TRUE)

