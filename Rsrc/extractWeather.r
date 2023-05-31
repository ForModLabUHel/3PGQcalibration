# this is needed to run python inside R
library(reticulate)
library(readxl)
library(data.table)
library(plyr)
library(tidyverse)
library(sf)
library(raster)
##########
# coordPlots = coordinates of plots; must have columns long and lat
# sY = start year, sM = start month, sD = start day, eY = end year, eM = end month, eD = end day
getWDpoints <- function(coordPlots, sY=1951, sM=1, sD=1, eY=2100, eM=12, eD=31) {
  # here we store the weather data
  wDs <- list()
  
  # this is the number of coordinate points
  nDFs <- nrow(coordPlots)
  
  for (i in 1:nDFs) {
    lon <- coordPlots$long[i]
    lat <- coordPlots$lat[i]
    
    # getting the weather data from clipick
    result <- getWeatherData(lon, lat, sY, sM, sD, eY, eM, eD)
    
    # weather datas as data.table
    wD <- data.table(t(sapply(result,c)))
    
    # fix the column names
    names(wD) <- as.character(wD[1,])
    wD <- wD[-1,]
    
    # this is the number of the day
    wD$rday <- 1:nrow(wD)
    
    wDs[[i]] <- wD
  }
  
  return(wDs)
}

#####################

setwd("C:/Users/minunno/Documents/github/3PGQcalibration/")
source("Rsrc/extractWeatherFun.r")
dataX <- data.table(read_excel("myData/U22128_fm.xlsx", sheet = "U22128_MI"))

coords <- dataX[,.(`Unique plot-ID`,`Coordinate SWEREF 99 TM East. Distorted +-200-800 m`,`Coordinate SWEREF 99 TM North. Distorted +-200-800 m`)]
setnames(coords,c("plotID","long","lat"))
coords <- unique(coords)
temp <- coords[,2:3]
temp <- SpatialPoints (temp, proj4string = CRS ('+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
coords[,2:3] <- as.data.frame (spTransform(temp, CRS('+proj=longlat +datum=WGS84 +no_defs')))

# specify here the path of python in your computer
# the code works with python 2 but not with python 3
use_python("C:/Python27")

source_python("Rsrc/clipick.py")
##if error:
# devtools::install_github("rstudio/reticulate")
# reticulate::install_miniconda()

##extract weather data
WD <- getWDpoints(coords[,2:3])

# coord <- data.table(lat=60.2833333,long=14.969444444444445)
# clim <- getWDpoints(coord)
# x = as.data.frame(clim)

# write.csv(clim,file = "C:/Users/39348/OneDrive/Desktop/LAVORO_CHECCO/CALIBRAZIONE/Allometry/DATI_CLIMATICI/Dati_giornalieri_estratti/climate_data_Nyhammar.csv",)

###convert weather data
for(i in length(WD)){
  WD[[i]]$PAR<-as.numeric(WD[[i]]$rss)*0.44*4.56
  # MJ / sq meter TO mol/m2 # 4.56 (?mol/J* or ?mol s???1W*???1)
  # 0.44 Percents of the energy of solar radiation is actually in the 400 - 700 nm range  
  WD[[i]]$TAir<-(as.numeric(WD[[i]]$tasmax)+
                   as.numeric(WD[[i]]$tasmin))/2
  SVP <-  610.7*10^(7.5*as.numeric(WD[[i]]$TAir)/
                      (237.3+as.numeric(WD[[i]]$TAir)))
  WD[[i]]$VPD <- SVP*(1-(as.numeric(WD[[i]]$hursmax)+
                as.numeric(WD[[i]]$hursmin))/2/100)/1000
  WD[[i]]$Precip<-as.numeric(WD[[i]]$pr)
}
save(WD,file="weatherInputs.rdata")
