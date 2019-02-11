setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())
library(fields) #rdist, Matern is here
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS) #mvrnorm is here
library(akima)
library(ggplot2)
library(reshape2) #melt
library(gridExtra) #grid.arrange

# data
topo.data <- read.table("topo.dat", header=TRUE)

# using akima library
topo.interp <- interp(x = topo.data$x, y = topo_data$y, z = topo_data$z)
topo.zmin <- min(topo.interp$z,na.rm=TRUE)
topo.zmax <- max(topo.interp$z,na.rm=TRUE)
breaks <- pretty(c(topo.zmin,topo.zmax),20)
colors <- heat.colors(length(breaks)-1)
with(topo.interp, image  (x,y,z, breaks=breaks, col=colors))
with(topo.interp,contour(x,y,z, levels=breaks, add=TRUE))
points (topo.data, pch = 3)

topogeo <- as.geodata(topo.data)
grd <-expand.grid(x = seq(from = 0, to = 315, by = 1), y = seq(from = 0, to = 315, by = 1))
kc <- krige.conv(topogeo, locations = grd, krige = krige.control(type.krige = "ok",
        cov.model = "sph",cov.pars = c(30000,60),nugget = 10000, trend.l = "1st", trend.d = "1st"))
#kc1 <- krige.conv(topogeo, locations = data.frame(30,40), krige = krige.control(type.krige ="ok",
#      cov.model = "exp",cov.pars =c(10,3.33),nugget = 0,trend.l = "1st", trend.d = "1st"))
kc_melt = data.frame(x = grd[,1],y = grd[,2], z = kc$predict)
ggplot(kc_melt, aes(x=x,y=y,fill = z)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "navy",
                       mid = "yellow",
                       high = "red",
                       midpoint = (max(kc_melt$z) + min(kc_melt$z))/2,
                       limits = range(kc_melt$z),
                       name = "prediction")+ 
  xlab("x") +
  ylab("y") + 
  ggtitle("Universal Kriging") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("../figures/uk1.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)