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



loci <- expand.grid(seq(0,315), seq(0,315))
krige.conv(topo.interp, locations = c(145,120), krige)
krige.control(type.krige = "uk", trend.d = "cte", trend.l = "cte",
             obj.model = NULL, beta, cov.model, cov.pars, kappa,
             nugget, micro.scale = 0, dist.epsilon = 1e-10,
             aniso.pars, lambda)
