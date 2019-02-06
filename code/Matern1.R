library(fields)
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS)
library(akima)
library(ggplot2)
library(reshape2)

distance=0:0.01:49/10
nuMat=c(1,3)
CorrMatern1=Matern(distance , range = 1, nu=nuMat[1]) #What should range be?
CorrMatern2=Matern(distance , range = 1, nu=nuMat[2])

nuExp = c(1,1.9)
CorrExp1 = exp(-distance^nuExp[1])
CorrExp2 = exp(-distance^nuExp[2])

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

ggsave("../figures/corrfunc.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

