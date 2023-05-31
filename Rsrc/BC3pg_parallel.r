library(r3PG)
library(dplyr)
library(ggplot2)
library(multidplyr)
library(tidyr)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(data.table)
library(parallel)

source("Rsrc/functions.r")
load('myData/input3pgLists.rdata')

sites <- c(1:26,28:57)#length(allInputs$site)
pErr <- c(0.1,0.001)
nCores = 1
climate <- climate[1:1000,]

## calibration settings
iterations=3e3
thin = 1 #100
nChains <- 3

# ####runModel for multiple sites and report likelihood
# ll <- numeric(nSites)
# ll <-mclapply(sites, function(i,ll){
#   print(i)
#   ll[i] <- multi_r3pg(allInputs$site[[i]],allInputs$species[[i]],
#                      allInputs$thinning[[i]],allInputs$climate[[i]],
#                      allInputs$obsData[[i]],parameters)
# },ll=ll,mc.cores = nCores)
# 
# 
# ####runModel for multiple sites and return a list of databases for each site with the simulated ond observed data
# datX <- list()
# datX <-mclapply(sites, function(i,datX){
#   print(i)
#   datX[[i]] <- multi_r3pg(allInputs$site[[i]],allInputs$species[[i]],
#                       allInputs$thinning[[i]],allInputs$climate[[i]],
#                       allInputs$obsData[[i]],parameters,outType="datX")
# },datX=datX,mc.cores = nCores)
# 


# calibration
### Define the prior
pIds <- 1:20
par <- list()
par$best <- parameters$`Pinus sylvestris`[pIds]
par$min <- parameters$`Pinus sylvestris`[pIds] * 0.8
par$max <- parameters$`Pinus sylvestris`[pIds] * 1.2
par$names <- parameters$parameter[pIds]

pY <- which(par$best<0)
par$min[pY] <- par$best[pY] * 1.2
par$max[pY] <- par$best[pY] * 0.8
pY <- which(par$best==0)
par$min[pY] <- par$best[pY] -0.1
par$max[pY] <- par$best[pY] + 0.1 

prior <- createUniformPrior(lower = par$min, upper = par$max)
### Create Bayesian Setup
BCmod <- createBayesianSetup(logLike, prior, best = par$best, 
                                names = par$names)

settings <- list(iterations = iterations, nrChains = nChains,thin=thin)

calibration <- runMCMC(BCmod, sampler="DREAMzs", settings = settings)

save(calibration, file="calOut/calibration.rdata")
