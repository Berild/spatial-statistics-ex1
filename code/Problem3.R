rm(list=ls())
library(fields)
library(MASS)
library(ggplot2)
library(geoR)

#Make prior, draw realisation and locations to measure.

###########################################################################
#Set grid L
buildGrid = function(N1, N2){
  N = N1*N2
  x1_val = 1:N1
  x2_val = 1:N2
  x1_pos = rep(x1_val,N2)
  x2_pos = rep(x2_val, each= N1)
  vec_pos = 1:N
  grid = list(N1 = N1, N2 = N2, x1 = x1_pos, x2 = x2_pos, vec = vec_pos)
  return(grid)
}

N1 = 30
N2 = 30
L = buildGrid(N1, N2)

###########################################################################
#Build prior model

buildPrior = function(L, mu_r, sigma_r, tau, xi){
  N = L$N1*L$N2
  tau = rdist(cbind(L$x1, L$x2))

  Sigma = sigma_r*exp(-(1/xi)*tau)
  E = rep(mu_r, N)
  prior = list(Sigma = Sigma, E = E)
  return(prior)
}

#Specifications:
mu_r = 0
sigma_r = 2
xi = 15

prior = buildPrior(L, mu_r, sigma_r, tau, xi)

###########################################################################
# Draw realization
set.seed(12)
realization = mvrnorm(n=1, mu = prior$E, Sigma = prior$Sigma)


# Display realization
# Remove stat_contour if you don't want to display the contour lines
plot_realization = ggplot(data.frame(x1 = L$x1, x2 = L$x2, realization = realization),
       aes(x=x1, y=x2)) + geom_raster(aes(fill = realization)) +
  scale_fill_continuous(high='blue 4', low='white') +
  stat_contour(aes(z = realization), col= "grey 50") +
  ggtitle("Sample from prior model") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_realization)

# Save plot
name = paste("../figures/3a_realization.pdf", sep = "")
ggsave(name, plot = plot_realization, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

###########################################################################
# Compute variogram and empirical variogram
variogram = function(sigma, xi, tau) {
  return(sigma_r - exp(-tau/xi))
}

coords = cbind(L$x1, L$x2)
tau = unique(c(rdist(coords))) # unique distances
gamma = variogram(sigma_r, xi, tau)
gamma_empirical = variog(coords=coords, data=realization)

# Plot variogram and empirical variogram
df_var = data.frame(
  tau=c(tau, gamma_empirical$u),
  gamma=c(gamma, gamma_empirical$v),
  var=c(
    rep("Theoretical", length(tau)),
    rep("Empirical", length(gamma_empirical$v))
  )
)

plot_variogram = ggplot(df_var) +
  geom_line(aes(x=tau, y=gamma, col=var)) +
  ggtitle("Variogram function") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram)

# Save plot
name = paste("../figures/3b_variogram.pdf", sep = "")
ggsave(name, plot = plot_variogram, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)


###########################################################################
#Draw positions to measure

randomPositions= function(num, L){
  N=L$N1*L$N2
  positions = ceiling(N*runif(num))
  return(positions)
}

obsMatrix = function(positions, L){
  N=L$N1*L$N2
  H = matrix(0,ncol = N, nrow = length(positions))
  for (i in 1:length(positions)){
    H[i,positions[i]] = 1
  }
  return(H)
}

num = 36
pos = randomPositions(num, L)
H = obsMatrix(pos, L)

#Display positions
#Same plot as above, only with the drawn positions
ggplot(data.frame(x1 = L$x1, x2 = L$x2, realization = realization),
       aes(x=x1, y=x2)) + geom_raster(aes(fill = realization)) +
  scale_fill_continuous(high='blue 4', low='white') +
  #stat_contour(aes(z = realization), col= "grey 50") +
  geom_point(data = data.frame(x1=L$x1[pos], x2 = L$x2[pos]),
             aes(x = x1, y=x2), shape = 21, col = 'red') +
  ggtitle("Sample from prior model") +
  theme(plot.title = element_text(hjust = 0.5))

#################################################
#Draw exact measurements in these positions
simulateMeasurements = function(r_true, obsMat, measure_error){
  m = nrow(obsMat)
  measurements = obsMat%*%r_true + measure_error*rnorm(m)
  return(measurements)
}

measurements = simulateMeasurements(realization, H, 0)

