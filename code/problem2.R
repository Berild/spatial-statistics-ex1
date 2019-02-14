setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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
topo.data <- read.table("topo.dat", header=TRUE)

# using akima library
topo.interp <- interp(x = topo.data$x, y = topo.data$y, z = topo.data$z)
topo.zmin <- min(topo.interp$z,na.rm=TRUE)
topo.zmax <- max(topo.interp$z,na.rm=TRUE)
breaks <- pretty(c(topo.zmin,topo.zmax),20)
colors <- terrain.colors(length(breaks)-1)
pdf("../figures/terrain2a.pdf",width = 5.5, height = 4)
with(topo.interp, image  (x,y,z, breaks=breaks, col=colors))
with(topo.interp,contour(x,y,z, levels=breaks, add=TRUE))
points (topo.data, pch = 3)
dev.off()

# Creating a geodata variable to use in geoR calculations
topogeo <- as.geodata(topo.data)

# Creating the grid of predicted values
grd <-expand.grid(x = seq(from = 0, to = 315, by = 1), y = seq(from = 0, to = 315, by = 1))


kc2 <- krige.conv(topogeo, locations = grd, krige = krige.control(type.krige = "ok",
        cov.model = "sph",cov.pars = c(2500,100),nugget = 0, trend.l = "2nd", trend.d = "2nd"))

# Creating dataframe
kc2.df = data.frame(x = grd[,1],y = grd[,2], pred = kc2$predict, se = sqrt(kc2$krige.var))

# Plotting predicted values
kc2.plot.pred <- ggplot(kc2.df, aes(x=x,y=y,fill = pred)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "navy",
                       mid = "yellow",
                       high = "red",
                       midpoint = (max(kc2.df$pred) + min(kc2.df$pred))/2,
                       limits = range(kc2.df$pred),
                       name = "prediction")+ 
  xlab("x") +
  ylab("y") + 
  ggtitle("UK 2nd order") +
  theme(plot.title = element_text(hjust = 0.5))
kc2.plot.pred
ggsave("../figures/uk2.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


# Plotting variance of predicted values
kc2.plot.se <- ggplot(kc2.df, aes(x=x,y=y,fill = se)) + 
  geom_tile() + 
  scale_fill_gradient2(high = "lightpink2",
                       mid = "steelblue2",
                       low = "navyblue",
                       midpoint = (max(kc2.df$se)+min(kc2.df$se))/2,
                       limits = range(kc2.df$se),
                       name = "SE")+ 
  xlab("x") +
  ylab("y") + 
  ggtitle("UK Variance 2nd order") +
  theme(plot.title = element_text(hjust = 0.5))
kc2.plot.se
ggsave("../figures/uk2se.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

# --------------------------------------
# 1st degree expectation function
kc1 <- krige.conv(topogeo, locations = grd, krige = krige.control(type.krige = "ok",
        cov.model = "sph",cov.pars = c(2500,100),nugget = 0, trend.l = "1st", trend.d = "1st"))

# Creating dataframe
kc1.df = data.frame(x = grd[,1],y = grd[,2], pred = kc2$predict, se = sqrt(kc2$krige.var))

# Plotting predicted values
kc1.plot.pred <- ggplot(kc1.df, aes(x=x,y=y,fill = pred)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "navyblue",
                       mid = "yellow",
                       high = "red",
                       midpoint = (max(kc1.df$pred) + min(kc1.df$pred))/2,
                       limits = range(kc1.df$pred),
                       name = "prediction")+ 
  xlab("x") +
  ylab("y") + 
  ggtitle("UK 1st order") +
  theme(plot.title = element_text(hjust = 0.5))
kc1.plot.pred
ggsave("../figures/uk1.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


# Plotting variance of predicted values
kc1.plot.se <- ggplot(kc1.df, aes(x=x,y=y,fill = se)) + 
  geom_tile() + 
  scale_fill_gradient2(high = "lightpink2",
                       mid = "steelblue2",
                       low = "navyblue",
                       midpoint = (max(kc1.df$se)+min(kc1.df$se))/2,
                       limits = range(kc1.df$se),
                       name = "SE")+ 
  xlab("x") +
  ylab("y") + 
  ggtitle("UK Variance 1st order") +
  theme(plot.title = element_text(hjust = 0.5))
kc1.plot.se
ggsave("../figures/uk1se.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

# multiplot of everything
multiplot(kc1.plot.pred, kc1.plot.se, kc2.plot.pred, kc2.plot.se, cols = 2)
ggsave("../figures/ukmulti.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


# ----------------------------------------------
# calculation of probablities d-f

pnorm(700, mean = kc2.df$pred[(kc2.df$x==100)&(kc2.df$y==100)], sd = (kc2.df$se[(kc2.df$x==100)&(kc2.df$y==100)]),lower.tail = FALSE)
