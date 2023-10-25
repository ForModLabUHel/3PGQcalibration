library(Metrics)
library(gridExtra)
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
setwd("C:/Users/minunno/Documents/github/3PGQcalibration")
load('myData/dataInputs.rdata')
load("myData/pMAP.rdata")

parameters <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parameters')
pars_soilQlitter <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parsQlitter')
pars_SoilInit <- read_excel('myData/INPUT_cal.xlsx', sheet = 'pars_SoilInit')
data_sizeDist <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_sizeDist')
climateData <- data.table(read.csv('myData/monthly_weather.csv'))
###remove 0 radiation
climateData$srad[which(climateData$srad==0)]  <- 0.5
par <- fread('myData/parameters.csv')
npars <- length(par$best)
pMAP <- pMAP[1:npars]
source(url('https://raw.githubusercontent.com/ForModLabUHel/3PGQcalibration/master/Rsrc/functions.r'))
#source("Rsrc/functions.r")
# if(exists("initFile")){
#   load(paste0("calOut/",initFile))
#   startValue <- rbind(MAP(calibration[[1]])$parametersMAP,
#                       MAP(calibration[[2]])$parametersMAP,
#                       MAP(calibration[[3]])$parametersMAP)
# }

####select sites forcalibration
nCalSites <-  floor(length(siteData)*2/3)
set.seed(1234)
calSites <- sort(sample(siteData,nCalSites))
valSites <- siteData[!siteData %in% calSites]
sites <- valSites #calSites

# 
# load("NAsites.rdata")
# sites <- sites[!sites %in% NAs]

nSites <- length(sites)
# pErr <- rep(c(0.1,0.001),7)
# climate <- climate[1:1000,]

# ## calibration settings
# if(!exists("calN")) calN <- 0
# if(!exists("iterations")) iterations=3e3
# if(!exists("nCores")) nCores = 1
# if(!exists("thin")) thin = 100
# if(!exists("nChains")) nChains <- 3
# 
# print(thin)
# print(iterations)
# 
# sets <-split(siteData, ceiling(seq_along(siteData)/100))
# 
# ###soilParamters
# pSoil <- c(7,0.36,0.25,0.5,1.1,0.3) #"beta", "eta", "e0", "fc", "q0", "z
# pErrSoil <- c(0,10) ##bias error, standard deviation

# calibration
### Define the prior
pIds <- c(1:7,10,13:15,31:32,39,41,44,47,49)
noDecid_par <- 7

dataOr <- runModelOut(par$best,sites)
dataMAP <- runModelOut(pMAP,sites)

for(i in 1:length(dataMAP)) dataMAP[[i]]$siteID <- sites[i]
for(i in 1:length(dataOr)) dataOr[[i]]$siteID <- sites[i]

dataMAP <- do.call(rbind,dataMAP)
dataOr <- do.call(rbind,dataOr)

dataMAP$sims <- as.numeric(dataMAP$sims)
dataMAP$obs <- as.numeric(dataMAP$obs)


###add species
dataMAP$species <- dataOr$species <- "some"
for(i in sites){
  species <- inputs[[i]]$species[,.(species, planted, fertility, stems_n, biom_stem, biom_root, biom_foliage)]
  layerX <- as.numeric(which(species$biom_stem>0))
  speciesX <- species[layerX,]
  nLayers <- nrow(speciesX)
  
  for(ijx in 1:nLayers){
    dataMAP[siteID==i & layerID==ijx]$species <- speciesX[ijx]$species
    dataOr[siteID==i & layerID==ijx]$species <- speciesX[ijx]$species
  }
}


vars <- unique(dataMAP$var_name)
vars <- vars[-4]
plotListMAP <- plotListOr <- list()
for(i in vars){
  plotListMAP[[i]] <-  ggplot(data=dataMAP[var_name==i], aes(x=sims, y=obs,col=species)) +
    geom_point() + geom_abline(slope=1,intercept=0)+
    ggtitle(i)
  
  plotListOr[[i]] <- ggplot(data=dataOr[var_name==i], aes(x=sims, y=obs,col=species)) +
    geom_point() + geom_abline(slope=1,intercept=0) +
    ggtitle(i)
}

pdf("resCal/resultsCal.pdf")
  grid.arrange(grobs = c(plotListOr[1:4], plotListMAP[1:4]), ncol = 2, as.table = FALSE)
  grid.arrange(grobs = c(plotListOr[5:7], plotListMAP[5:7]), ncol = 2, as.table = FALSE)
dev.off()


rmseTab <- matrix(NA, length(vars),2)
for(i in 1:length(vars)){
  rmseTab[i,1] <- rmse(dataOr[var_name==vars[i]]$obs, dataOr[var_name==vars[i]]$sims)
  rmseTab[i,2] <- rmse(dataMAP[var_name==vars[i]]$obs, dataMAP[var_name==vars[i]]$sims)
}
colnames(rmseTab) <- c("original","calibrated")
rownames(rmseTab) <- vars

write.csv(rmseTab,file="resCal/rmse.csv")


for(i in 1:7){
  png(file=paste0("resCal/plots/",vars[i],"_pDef.png"))
    print(plotListOr[i])
  dev.off()

  png(file=paste0("resCal/plots/",vars[i],"_pCal.png"))
    print(plotListMAP[i])
  dev.off()
}
