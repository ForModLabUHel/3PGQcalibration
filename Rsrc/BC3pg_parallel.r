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
library(readxl)
library(tictoc)

source("Rsrc/functions.r")
load('myData/dataInputs.rdata')
parameters <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parameters')
data_sizeDist <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_sizeDist')
climateData <- data.table(read.csv('myData/monthly_weather.csv'))

sites <- 1:length(site_list)

load("NAsites.rdata")
sites <- sites[!sites %in% NAs]

nSites <- length(sites)
pErr <- rep(c(0.1,0.001),7)
nCores = 1
# climate <- climate[1:1000,]

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
par$best <- c(parameters$`Pinus sylvestris`[pIds],parameters$`Picea abies`[pIds],
              parameters$`Pinus contorta`[pIds],parameters$`Betula alba`[pIds],
              parameters$`other deciduous`[pIds],pErr)
par$min <- par$best * 0.8
par$max <- par$best * 1.2
par$names <- c(paste0(parameters$parameter[pIds],"_pine"),
               paste0(parameters$parameter[pIds],"_spruce"),
               paste0(parameters$parameter[pIds],"_contorta"),
               paste0(parameters$parameter[pIds],"_birch"),
               paste0(parameters$parameter[pIds],"_decid"),
               "aN","bN","aB","bB","aV","bV","aD","bD",
               "aWs","bWs","aWr","bWr","aWf","bWf")

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




# nCores = 5 #30
# cl <- parallel::makeCluster(nCores)
# pLikelihood <- function(param) parallel::parApply(cl = cl, X = param, MARGIN = 1, FUN = logLike)
# parallel::clusterExport(cl, varlist = ls())#c('parameters', 'par',
#                                         # 'pIds', 'multi_r3pg', 'site_list',
#                                         # "species_list","climateData","obsAll",
#                                         # 'run_3PG',# 'settings_3pg', #'n3PGpars', 
#                                         # "calc_ll", "Sivia_log", 'nSites', 'sites'))
# 
# prior <- createUniformPrior(lower = par$min, upper = par$max)
# ### Create Bayesian Setup
# BCmod <- createBayesianSetup(pLikelihood, prior, best = par$best, 
#                              names = par$names, 
#                              parallel = "external")
# iterations=100
# settings <- list(iterations = iterations, nrChains = nChains, thin=thin)
# 
# # https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html
# # rm(i)
# # calibration <- runMCMC(BCmod, sampler="DREAMzs", settings = settings)
# tictoc::tic(); set.seed(1234); calibration <- runMCMC(BCmod, sampler="DREAMzs", settings = settings); t = tictoc::toc()
# stopCluster(cl)
# 
# x = getSample(calibration, coda = F, parametersOnly = F) %>% data.frame()
# dim(x)
# 
# x %>% 
#   pivot_longer(everything()) %>% 
#   ggplot(aes(x = value)) +
#   geom_density() +
#   facet_wrap(~name, scales = "free") +
#   theme_bw()
# 
# # readr::write_rds(calibration, "calibration_50000_iter_10_plots.rds")
# # readr::write_rds(calibration, "calibration_200000_iter_10_plots.rds")
# best_param = MAP(calibration)
