library(r3PG)
library(dplyr)
library(ggplot2)
# library(multidplyr)
library(tidyr)
library(purrr)
library(BayesianTools)
library(sensitivity)
library(data.table)
library(parallel)
library(readxl)
library(tictoc)

###read data
load('myData/dataInputs.rdata')
parameters <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parameters')
pars_soilQlitter <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parsQlitter')
pars_SoilInit <- read_excel('myData/INPUT_cal.xlsx', sheet = 'pars_SoilInit')
data_sizeDist <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_sizeDist')
climateData <- data.table(read.csv('myData/monthly_weather.csv'))
###remove 0 radiation
climateData$srad[which(climateData$srad==0)]  <- 0.5
par <- fread('myData/parameters.csv')

source(url('https://raw.githubusercontent.com/ForModLabUHel/3PGQcalibration/master/Rsrc/functions.r'))

if(exists("initFile")){
  load(paste0("calOut/",initFile))
  startValue <- rbind(MAP(calibration[[1]])$parametersMAP,
                      MAP(calibration[[2]])$parametersMAP,
                      MAP(calibration[[3]])$parametersMAP)
}

####select sites forcalibration
nCalSites <-  floor(length(siteData)*2/3)
set.seed(1234)
calSites <- sort(sample(siteData,nCalSites))
valSites <- siteData[!siteData %in% calSites]
sites <- calSites
# 
# load("NAsites.rdata")
# sites <- sites[!sites %in% NAs]

nSites <- length(sites)
pErr <- rep(c(0.1,0.001),7)
# climate <- climate[1:1000,]

## calibration settings
if(!exists("calN")) calN <- 0
if(!exists("iterations")) iterations=3e3
if(!exists("nCores")) nCores = 1
if(!exists("thin")) thin = 100
if(!exists("nChains")) nChains <- 3

print(thin)
print(iterations)

sets <-split(sites, ceiling(seq_along(sites)/(nSites/21)))

###soilParamters
pSoil <- c(7,0.36,0.25,0.5,1.1,0.3) #"beta", "eta", "e0", "fc", "q0", "z
pErrSoil <- c(0,10) ##bias error, standard deviation

# calibration
### Define the prior
pIds <- c(1:7,10,13:15,31:32,39,41,44,47,49)
noDecid_par <- 7
# par <- list()
# par$best2 <- c(parameters$`Pinus sylvestris`[pIds],parameters$`Picea abies`[pIds],
#               parameters$`Pinus contorta`[pIds],parameters$`Betula alba`[pIds[-noDecid_par]],
#               parameters$`other deciduous`[pIds[-noDecid_par]],pErr,pSoil,pErrSoil)
# par$min <- par$best * 0.8
# par$max <- par$best * 1.2
# par$names <- c(paste0(parameters$parameter[pIds],"_pine"),
#                paste0(parameters$parameter[pIds],"_spruce"),
#                paste0(parameters$parameter[pIds],"_contorta"),
#                paste0(parameters$parameter[pIds[-noDecid_par]],"_birch"),
#                paste0(parameters$parameter[pIds[-noDecid_par]],"_decid"),
#                "aN","bN","aB","bB","aV","bV","aD","bD",
#                "aWs","bWs","aWr","bWr","aWf","bWf",
#                "beta", "eta", "e0", "fc", "q0", "z",
#                "biasSoil","sdSoil")

# pY <- which(par$best<0)
# par$min[pY] <- par$best[pY] * 1.2
# par$max[pY] <- par$best[pY] * 0.8
# pY <- which(par$best==0)
# par$min[pY] <- par$best[pY] -0.1
# par$max[pY] <- par$best[pY] + 0.1 
# par$min[121] <- -10
# par$max[121] <- 10

if(!exists("startValue")){
  prior <- createUniformPrior(lower = par$min, upper = par$max)
  ### Create Bayesian Setup
  BCmod <- createBayesianSetup(likelihood, prior, best = par$best,
                               names = par$names)
  
  startValue <- rbind(runif(length(par$best),par$min,par$max),
                      runif(length(par$best),par$min,par$max),
                      pMAP[1:110])
  
  settings <- list(iterations = iterations,thin=thin,startValue=startValue,
                   message=FALSE,consoleUpdates=1000)
  
  tic(paste0("calibration time."," iteratios: ",iterations))
  calibration <- runMCMC(BCmod, sampler="DEzs", settings = settings)
  toc()
}else{
  if(exists("fx")){
    calibration$settings$f <- fx
  } 
  tic(paste0("calibration time."," iteratios: ",iterations))
  calibration = runMCMC(calibration)
  toc()
}

save(calibration, file=paste0("calOut/calibration_",calN,".rdata"))
