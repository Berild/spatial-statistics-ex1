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

buildPrior = function(L, mu_r, sigma_sqrd, tau, xi){
  N = L$N1*L$N2
  tau = rdist(cbind(L$x1, L$x2))

  Sigma = sigma_sqrd*exp(-tau/xi)
  E = rep(mu_r, N)
  prior = list(Sigma = Sigma, E = E)
  return(prior)
}

#Specifications:
mu_r = 0
sigma_r_sqrd = 2
xi = 15

prior = buildPrior(L, mu_r, sigma_r_sqrd, tau, xi)

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
  #ggtitle("Sample from prior model") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_realization)

# Save plot
ggsave("../figures/3a_realization.pdf",
       plot = plot_realization, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 4, units = "in",
       dpi = 300, limitsize = TRUE)

###########################################################################
# Compute variogram and empirical variogram
variogram = function(sigma_sqrd, xi, tau) {
  return(sigma_sqrd * (1 - exp(-tau/xi)))
}

coords = cbind(L$x1, L$x2)
tau = unique(c(rdist(coords))) # unique distances
gamma = variogram(sigma_r_sqrd, xi, tau)
gamma_empirical = variog(coords=coords, data=realization)

#gamma_empirical_mean = gamma_empirical$v
#k = 100
#for(i in 1:(k-1)) {
#  realization_i = mvrnorm(n=1, mu = prior$E, Sigma = prior$Sigma)
#  gamma_empirical_mean = gamma_empirical_mean + variog(coords=coords,
#                                                       data=realization_i)$v
#}
#gamma_empirical_mean = gamma_empirical_mean / k

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
  #ggtitle("Variogram") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram)

# Save plot
ggsave("../figures/3b_variogram.pdf",
       plot = plot_variogram, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)


###########################################################################
#Draw positions to measure

randomPositions= function(num, L){
  N=L$N1*L$N2
  positions = sample(1:N, num)
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
  #ggtitle("Random observations from sample") +
  theme(plot.title = element_text(hjust = 0.5))

###########################################################################
#Draw exact measurements in these positions
simulateMeasurements = function(r_true, obsMat, measure_error){
  m = nrow(obsMat)
  measurements = obsMat%*%r_true + measure_error*rnorm(m)
  return(measurements)
}

measurements = simulateMeasurements(realization, H, 0)

# Compute empirical variogram from observations
coords_obs = coords[pos,]
gamma_empirical_obs = variog(coords=coords_obs, data=measurements)

# Plot variogram and empirical variogram
df_var_obs = data.frame(
  tau=c(tau, gamma_empirical_obs$u),
  gamma=c(gamma, gamma_empirical_obs$v),
  var=c(
    rep("Theoretical", length(tau)),
    rep("Empirical", length(gamma_empirical_obs$v))
  )
)

plot_variogram_obs = ggplot(df_var_obs) +
  geom_line(aes(x=tau, y=gamma, col=var)) +
  #ggtitle("Variogram") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram_obs)

# Save plot
ggsave("../figures/3c_variogram_obs.pdf",
       plot = plot_variogram_obs, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)

###########################################################################
# Estimate sigma_r^2 and xi_r from the full realization and the 36 observations
initial_values = c(1, 1)
fit_full = likfit(coords=coords, data=realization,
                  ini.cov.pars=initial_values, cov.model="exponential",
                  fix.nugget=)#, trend=rep(0, length(realization)))
fit_obs = likfit(coords=coords_obs, data=measurements,
                 ini.cov.pars=initial_values, cov.model="exponential",
                 fix.nugget=TRUE)#, trend=rep(0, length(measurements)))
print(fit_full$cov.pars)
print(fit_obs$cov.pars)

# Plot theoretical varigram and the two estimates
gamma_est_full = variogram(fit_full$sigmasq, fit_full$phi, tau)
gamma_est_obs = variogram(fit_obs$sigmasq, fit_obs$phi, tau)

df_var_est = data.frame(
  tau=c(tau, tau, tau),
  gamma=c(gamma, gamma_est_full, gamma_est_obs),
  var=c(
    rep("Theoretical", length(tau)),
    rep("Full estimate", length(tau)),
    rep("Observation estimate", length(tau))
  )
)

plot_variogram_est = ggplot(df_var_est) +
  geom_line(aes(x=tau, y=gamma, col=var)) +
  #ggtitle("Variogram") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram_est)

# Save plot
ggsave("../figures/3c_variogram_estimates.pdf",
       plot = plot_variogram_est, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)

###########################################################################
# Procedure above repeated with different number of observations
nums = c(9, 64, 100)
poss = lapply(nums, function(x)(randomPositions(x, L)))
Hs = lapply(poss, function(x)(obsMatrix(x, L)))

measurementss = lapply(Hs, function(x)(simulateMeasurements(realization, x, 0)))

# Compute empirical variogram from observations
coords_obss = lapply(poss, function(x)(coords[x,]))
gamma_empirical_obss = vector("list", length(nums))
for(i in 1:length(nums)) {
  gamma_empirical_obss[[i]] = variog(coords=coords_obss[[i]],
                                     data=measurementss[[i]])
}

# Plot variogram and empirical variogram
df_var_obss = data.frame(
  tau=c(
    tau,
    gamma_empirical_obss[[1]]$u,
    gamma_empirical_obss[[2]]$u,
    gamma_empirical_obss[[3]]$u
    ),
  gamma=c(
    gamma,
    gamma_empirical_obss[[1]]$v,
    gamma_empirical_obss[[2]]$v,
    gamma_empirical_obss[[3]]$v
    ),
  var=c(
    rep("Theoretical", length(tau)),
    rep(paste(c("Empirical, ", nums[1], " observations"), collapse=""),
        length(gamma_empirical_obss[[1]]$v)),
    rep(paste(c("Empirical, ", nums[2], " observations"), collapse=""),
        length(gamma_empirical_obss[[2]]$v)),
    rep(paste(c("Empirical, ", nums[3], " observations"), collapse=""),
        length(gamma_empirical_obss[[3]]$v))
    )
)

plot_variogram_obss = ggplot(df_var_obss) +
  geom_line(aes(x=tau, y=gamma, col=var)) +
  #ggtitle("Variogram") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram_obss)

# Save plot
ggsave("../figures/3d_variogram_obss.pdf",
       plot = plot_variogram_obss, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)

# Estimate sigma_r^2 and xi_r
initial_values = c(1, 1)
fit1 = likfit(coords=coords_obss[[1]], data=measurementss[[1]],
              ini.cov.pars=initial_values, cov.model="exponential",
              fix.nugget=)#, trend=rep(0, length(realization)))
fit2 = likfit(coords=coords_obss[[2]], data=measurementss[[2]],
              ini.cov.pars=initial_values, cov.model="exponential",
              fix.nugget=TRUE)#, trend=rep(0, length(measurements)))
fit3 = likfit(coords=coords_obss[[3]], data=measurementss[[3]],
              ini.cov.pars=initial_values, cov.model="exponential",
              fix.nugget=TRUE)#, trend=rep(0, length(measurements)))
print(fit1$cov.pars)
print(fit2$cov.pars)
print(fit3$cov.pars)

# Plot theoretical varigram and the two estimates
gamma_est1 = variogram(fit1$sigmasq, fit1$phi, tau)
gamma_est2 = variogram(fit2$sigmasq, fit2$phi, tau)
gamma_est3 = variogram(fit3$sigmasq, fit3$phi, tau)

df_var_ests = data.frame(
  tau=c(tau, tau, tau, tau),
  gamma=c(gamma, gamma_est1, gamma_est2, gamma_est3),
  var=c(
    rep("Theoretical", length(tau)),
    rep(paste(c("Estimate, ", nums[1], " observations"), collapse=""),
        length(tau)),
    rep(paste(c("Estimate, ", nums[2], " observations"), collapse=""),
        length(tau)),
    rep(paste(c("Estimate, ", nums[3], " observations"), collapse=""),
        length(tau))
  )
)

plot_variogram_ests = ggplot(df_var_ests) +
  geom_line(aes(x=tau, y=gamma, col=var)) +
  #ggtitle("Variogram") +
  labs(y="value", col="") +
  theme(plot.title = element_text(hjust = 0.5))

# Show plot
print(plot_variogram_ests)

# Save plot
ggsave("../figures/3d_variogram_estimatess.pdf",
       plot = plot_variogram_ests, device = NULL, path = NULL,
       scale = 1, width = 5.5, height = 3, units = "in",
       dpi = 300, limitsize = TRUE)
