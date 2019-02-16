#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(fields) #rdist, Matern is here
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS) #mvrnorm is here
library(akima)
library(ggplot2)
library(reshape2) #melt
library(gridExtra) #grid.arrange
source("multiplot.R")

# data
topo.data <- read.table(
  "topo.dat", 
  header=TRUE
  )

# using akima library
topo.interp <- interp(
  x = topo.data$x, 
  y = topo.data$y, 
  z = topo.data$z
  )
topo.zmin <- min(topo.interp$z,na.rm=TRUE)
topo.zmax <- max(topo.interp$z,na.rm=TRUE)
breaks <- pretty(c(topo.zmin,topo.zmax),20)
colors <- terrain.colors(length(breaks)-1)
pdf("../figures/terrain2a.pdf",width = 5.5, height = 5.5)
with(topo.interp, image  (x,y,z, breaks=breaks, col=colors))
with(topo.interp,contour(x,y,z, levels=breaks, add=TRUE))
dev.off()

# Creating a geodata variable to use in geoR calculations
topogeo <- as.geodata(topo.data)

# Creating the grid of predicted values
grd <-expand.grid(
  x = seq(from = 0, to = 315, by = 1), 
  y = seq(from = 0, to = 315, by = 1)
  )

# Running krige prediction with 2nd order expectation function
kc2 <- krige.conv(
  topogeo, 
  locations = grd, 
  krige = krige.control(
    type.krige = "ok",
    cov.model = "sph",
    cov.pars = c(2500,100),
    nugget = 0, 
    trend.l = "2nd", 
    trend.d = "2nd"
    )
  )

# Creating dataframe
kc2.df = data.frame(
  x = grd[,1],
  y = grd[,2], 
  pred = kc2$predict, 
  var = kc2$krige.var
  )

# Plotting predicted values
kc2.plot.pred <- ggplot(kc2.df, 
                        aes(x=x,y=y,fill = pred)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "green3",
                       mid = "sandybrown",
                       high = "white",
                       midpoint = (max(kc2.df$pred) + min(kc2.df$pred))/2,
                       limits = range(kc2.df$pred),
                       name = "prediction")+ 
  geom_contour(aes(x = x, y =y, z = kc2.df$pred), 
               stat = "contour", 
               position = "identity", 
               lineend = "butt", 
               linejoin = "round",
               color = "gray52") + 
  xlab("x") +
  ylab("y") + 
  ggtitle("UK 2nd order") +
  theme(plot.title = element_text(hjust = 0.5))
kc2.plot.pred
ggsave("../figures/uk2.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


# Plotting variance of predicted values
kc2.plot.se <- ggplot(kc2.df, aes(x=x,y=y,fill = var)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "Spectral")+
  xlab("x") +
  ylab("y") + 
  ggtitle("UK Variance 2nd order") +
  theme(plot.title = element_text(hjust = 0.5))
kc2.plot.se
ggsave("../figures/uk2se.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
#scale_fill_gradient(high = "white",low = "blue",limits = range(kc2.df$var), name = "Var")+ 
# --------------------------------------
# 1st degree expectation function
kc1 <- krige.conv(
  topogeo, 
  locations = grd, 
  krige = krige.control(
    type.krige = "ok",
    cov.model = "sph",
    cov.pars = c(2500,100),
    nugget = 0, 
    trend.l = "1st", 
    trend.d = "1st"
    )
  )

# Creating dataframe
kc1.df = data.frame(
  x = grd[,1],
  y = grd[,2], 
  pred = kc1$predict, 
  var = kc1$krige.var
  )

# Creating a dataframe for the contour plots
kc.contour <- data.frame(
  x = rep(grd[,1],2),
  y = rep(grd[,2],2),
  pred_contour = c(kc1.df$pred,kc2.df$pred),
  group = c(rep("1st",length(kc1.df$pred)),rep("2nd",length(kc2.df$pred)))
  )
                                                                                                     
# looking at the difference between 1st and 2nd degree expectation function
# plotting the contours of the predicted value
ggplot(kc.contour)+
  geom_contour(aes(x=x,y=y,z = pred_contour,color=group),size=1)+
  scale_color_manual(
    name = "g(x) order",
    labels = c("1st", "2nd"),
    values = c("royalblue3", "tan3")
  )+
  xlab("x") +
  ylab("y") + 
  ggtitle("1st and 2nd order expectation function") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("../figures/uk1n2.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
# Adding observation errors to the krigeing
# => adding nugget tot the kige.conv function

nugget5 <- krige.conv(
  topogeo, 
  locations = grd, 
  krige = krige.control(
    type.krige = "ok",
    cov.model = "sph",
    cov.pars = c(2500,100),
    nugget = 5, # Here we add the observation error
    trend.l = "1st", 
    trend.d = "1st"
    )
  )

nugget25 <- krige.conv(
  topogeo, 
  locations = grd, 
  krige = krige.control(
    type.krige = "ok",
    cov.model = "sph",
    cov.pars = c(2500,100),
    nugget = 25, 
    trend.l = "1st", 
    trend.d = "1st"
    )
  )

# Creating data.frame of the two predictions with nuggets
nugget5.df = data.frame(
  x = grd[,1],
  y = grd[,2], 
  pred = nugget5$predict, 
  se = sqrt(nugget5$krige.var)
  )

nugget25.df = data.frame(
  x = grd[,1],
  y = grd[,2], 
  pred = nugget25$predict, 
  se = sqrt(nugget25$krige.var)
  )         

# Creating a collected data.frame for compairing in ggplot
nugget.df <- data.frame(
  x= rep(grd[,1],2),
  y=rep(grd[,2],2),
  pred=c(nugget5.df$pred,nugget25.df$pred),
  se=c(nugget5.df$se,nugget25.df$se),
  group=c(rep("5",length(nugget5.df$pred)),
          rep("25",length(nugget25.df$pred)))
  )

# Plotting the contours of the standard error
# of both predictions with nuggets
ggplot(nugget.df,aes(x=x,y=y))+
  geom_contour(aes(z = se,color=group),size=0.75)+
  scale_color_manual(
    name = "Nugget",
    labels = c("5", "25"),
    values = c("royalblue3", "tan3")
  )+
  xlab("x") +
  ylab("y") + 
  ggtitle("Observation error") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("../figures/nugget.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)
