#Sample 1D
library(fields) #rdist, Matern is here
library(geoR) #Need to install XQuartz? cov.spatial() belongs here
library(MASS) #mvrnorm is here
library(akima)
library(ggplot2)
library(reshape2) #melt
library(gridExtra) #grid.arrange

######################################################################
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

##########################################
#Covariance matrix
tau = rdist(L,L)/10

CovMatPrior = function(tau, sigmasq, nu, func="m"){
  if (tolower(func)== "e"){
    corr = exp(-tau^nu)
  }else if(tolower(func) == "m"){
    corr = Matern(tau, range = 1, nu = nu) #What should range be?
  }else{
    warning("Unknown correlation function")
    return(NULL)
  }
  Cov = sigmasq*corr
  return(Cov)
}

############################
#Make vectors of parameter values to iterate through
sigmaLong = rep(sigma2, times = 4)
nuLong = c(rep(nu_e, each = 2), rep(nu_m, each = 2))
funcLong = c(rep("e", 4), rep("m", 4))
f = ifelse(funcLong=="e", "Pow Exp: ", "Matern: ")
specifications = paste(f,"nu=",nuLong,", sigma=",sigmaLong, sep = "")

####################################################################################
#Make samples and display them
set.seed(42)
nsamps = 10

##########################################
#Choosing parameters for the "True system"
#Unsure whether this is the right way to do it
sigma_true = 5
func_true = "e"
nu_true = nu_e[2]
which_save = nsamps

###########################################
plots = list()
for(i in 1:8){
  Cov = CovMatPrior(tau,sigmaLong[i], nuLong[i], funcLong[i])
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
  
  plots= ggplot(long_realis,
         aes(x=L, y=value, colour=variable)) +
    geom_line()+
    ggtitle("Realizations of 1D Gaussian RF") +
    xlab("x")+
    ylab("Realizations") +
    annotate("text", x = 10, y = max(long_realis$value)*1.1, label = specifications[i])+ 
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  ########################
  #Save plot
  #name = paste("../figures/sample1conf",i,".pdf", sep = "")
  #ggsave(name, plot = plot, device = NULL, path = NULL,
    #     scale = 1, width = 5.5, height = 4, units = "in",
     #   dpi = 300, limitsize = TRUE)
  
  #print(plot)
  plots[[i]] = plot
}

#########################################################################
#Display more than one plot at once:
#plotGrid1 = grid.arrange(grobs = plots[1:4], ncol = 2)
#plotGrid2 = grid.arrange(grobs = plots[5:8], ncol = 2)

