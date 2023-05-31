library(rgdal)
library(data.table)
library(readxl)

extractWeather=F
# setwd("C:/Users/minunno/Documents/research/3pgq/data")
dataX <- data.table(read_excel("myData/U22128_fm.xlsx", sheet = "U22128_MI"))
setnames(dataX,c("Unique plot-ID","Year of inventory","Dry weight kg/ha stump incl roots > 2mm all species"),c("Plot_ID","year","Wroot"))
####convert coordinates
coords <- dataX[,.(Plot_ID,year,`Coordinate SWEREF 99 TM East. Distorted +-200-800 m`,`Coordinate SWEREF 99 TM North. Distorted +-200-800 m`)]
setnames(coords,c("Plot_ID","year","long","lat"))
# coords <- unique(coords)
temp <- coords[,3:4]
temp <- SpatialPoints (temp, proj4string = CRS ('+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
coords[,3:4] <- as.data.frame (spTransform(temp, CRS('+proj=longlat +datum=WGS84 +no_defs')))
setkey(dataX,Plot_ID,year)
setkey(coords,Plot_ID,year)
dataX <- merge(dataX,coords,by=c("Plot_ID","year"))

selData <- dataX[,.SD[which(length(year)>1)],by=Plot_ID]$Plot_ID
selData <- dataX[Plot_ID %in% selData]

###this is a bit stupid
initVals <- selData[, .SD[which.min(year)], by = Plot_ID]
obs <- selData[, .SD[which.max(year)], by = Plot_ID]
initVals$dataUse="init"
obs$dataUse="cal"
selData <- rbind(initVals,obs)
setkey(selData,Plot_ID,year)

d_site  <- data.table(siteName=selData[dataUse=="init"]$`Plot ID within cluster`, 
                      Site_ID=selData[dataUse=="init"]$`Plot ID within cluster`,
                      Plot_Name= selData[dataUse=="init"]$Plot_ID,
                      Plot_ID=selData[dataUse=="init"]$Plot_ID,
                      latitude= selData[dataUse=="init"]$lat,
                      altitude=100,
                      soil_class = 3,
                      asw_i=999,
                      asw_min=0,
                      asw_max=160,
                      from=paste0(selData[dataUse=="init"]$year,"-01"),
                      to=paste0(selData[dataUse=="cal"]$year+1,"-01")
                      
)

selData[,planted:=paste0(year-`Stand age`,"-01")]

###avoid NaN
# selData[`Dry weight kg/ha branches incl small branches and needles all species`==0]$`Dry weight kg/ha branches incl small branches and needles all species`=0.0001
# selData[`Dry weight kg/ha stem incl bark all species`==0]$`Dry weight kg/ha stem incl bark all species` = 0.001
# selData[Wroot==0]$Wroot=0.0001

selData[,WstemPine:=`Vol/ha m3sk/ha Pine`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WstemSpruce:=`Vol/ha m3sk/ha Spruce`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WstemBirch:=`Vol/ha m3sk/ha Birch`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WstemContorta:=`Vol/ha m3sk/ha Contorta`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WstemDecid:=`Vol/ha m3sk/ha Other dec.`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]

selData[,WfoliagePine:=`BA/ha m2/ha Pine`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfoliageSpruce:=`BA/ha m2/ha Spruce`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfoliageBirch:=`BA/ha m2/ha aBirch`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfoliageContorta:=`BA/ha m2/ha Contorta`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfoliageDecid:=`BA/ha m2/ha Other dec.`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]

selData[,WrootPine:=`BA/ha m2/ha Pine`/`BA/ha m2/ha all species`*Wroot/1000]
selData[,WrootSpruce:=`BA/ha m2/ha Spruce`/`BA/ha m2/ha all species`*Wroot/1000]
selData[,WrootBirch:=`BA/ha m2/ha aBirch`/`BA/ha m2/ha all species`*Wroot/1000]
selData[,WrootContorta:=`BA/ha m2/ha Contorta`/`BA/ha m2/ha all species`*Wroot/1000]
selData[,WrootDecid:=`BA/ha m2/ha Other dec.`/`BA/ha m2/ha all species`*Wroot/1000]

setnames(selData,"Plot ID within cluster","siteName")
Wstem <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                        WstemPine,
                                        WstemSpruce,
                                        WstemBirch,
                                        WstemContorta,
                                        WstemDecid)],id.vars=1:3)
Wstem[variable=="WstemPine",layerID:=1]
Wstem[variable=="WstemPine",species:="Pinus sylvestris"]
Wstem[variable=="WstemSpruce",layerID:=2]
Wstem[variable=="WstemSpruce",species:="Picea abies"]
Wstem[variable=="WstemContorta",layerID:=3]
Wstem[variable=="WstemContorta",species:="Pinus contorta"]
Wstem[variable=="WstemBirch",layerID:=4]
Wstem[variable=="WstemBirch",species:="Betula alba"]
Wstem[variable=="WstemDecid",layerID:=5]
Wstem[variable=="WstemDecid",species:="other deciduous"]

Wroot <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                        WrootPine,
                                        WrootSpruce,
                                        WrootBirch,
                                        WrootContorta,
                                        WrootDecid)],id.vars=1:3)
Wroot[variable=="WrootPine",layerID:=1]
Wroot[variable=="WrootPine",species:="Pinus sylvestris"]
Wroot[variable=="WrootSpruce",layerID:=2]
Wroot[variable=="WrootSpruce",species:="Picea abies"]
Wroot[variable=="WrootContorta",layerID:=3]
Wroot[variable=="WrootContorta",species:="Pinus contorta"]
Wroot[variable=="WrootBirch",layerID:=4]
Wroot[variable=="WrootBirch",species:="Betula alba"]
Wroot[variable=="WrootDecid",layerID:=5]
Wroot[variable=="WrootDecid",species:="other deciduous"]

Wfoliage <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                           WfoliagePine,
                                           WfoliageSpruce,
                                           WfoliageBirch,
                                           WfoliageContorta,
                                           WfoliageDecid)],id.vars=1:3)
Wfoliage[variable=="WfoliagePine",layerID:=1]
Wfoliage[variable=="WfoliagePine",species:="Pinus sylvestris"]
Wfoliage[variable=="WfoliageSpruce",layerID:=2]
Wfoliage[variable=="WfoliageSpruce",species:="Picea abies"]
Wfoliage[variable=="WfoliageContorta",layerID:=3]
Wfoliage[variable=="WfoliageContorta",species:="Pinus contorta"]
Wfoliage[variable=="WfoliageBirch",layerID:=4]
Wfoliage[variable=="WfoliageBirch",species:="Betula alba"]
Wfoliage[variable=="WfoliageDecid",layerID:=5]
Wfoliage[variable=="WfoliageDecid",species:="other deciduous"]

stems <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                        `no of stems/ha  Pine`,
                                        `no of stems/ha  Spruce`,
                                        `no of stems/ha Birch`,
                                        `no of stems/ha  Contorta`,
                                        `no of stems/ha Other dec.`)],id.vars=1:3)
stems[variable=="no of stems/ha  Pine",layerID:=1]
stems[variable=="no of stems/ha  Pine",species:="Pinus sylvestris"]
stems[variable=="no of stems/ha  Spruce",layerID:=2]
stems[variable=="no of stems/ha  Spruce",species:="Picea abies"]
stems[variable=="no of stems/ha  Contorta",layerID:=3]
stems[variable=="no of stems/ha  Contorta",species:="Pinus contorta"]
stems[variable=="no of stems/ha Birch",layerID:=4]
stems[variable=="no of stems/ha Birch",species:="Betula alba"]
stems[variable=="no of stems/ha Other dec.",layerID:=5]
stems[variable=="no of stems/ha Other dec.",species:="other deciduous"]

setnames(Wstem,"value","biom_stem");Wstem$variable <- NULL
setnames(Wroot,"value","biom_root");Wroot$variable <- NULL
setnames(Wfoliage,"value","biom_foliage");Wfoliage$variable <- NULL
setnames(stems,"value","stems_n");stems$variable <- NULL

setkey(Wstem,"siteName","Plot_ID","planted","layerID","species")
setkey(Wroot,"siteName","Plot_ID","planted","layerID","species")
setkey(Wfoliage,"siteName","Plot_ID","planted","layerID","species")
setkey(stems,"siteName","Plot_ID","planted","layerID","species")

d_species <- merge(Wstem,Wroot,by=c("siteName","Plot_ID","planted","layerID","species"))
d_species <- merge(d_species,Wfoliage,by=c("siteName","Plot_ID","planted","layerID","species"))
d_species <- merge(d_species,stems,by=c("siteName","Plot_ID","planted","layerID","species"))

###process weather inputs
data_climate2 <- data.table(read.csv('myData/prebas_sweden_may23_monthly_weather.csv'))
coords2 <- data.table(read.csv('myData/coords_climid.csv'))
# setkey(data_climate2, climID)
# setkey(coords2, climID)
coords2$siteID <- as.numeric(as.factor(coords2$climID))
climateData <- merge(data_climate2,unique(coords2[, .(climID,siteID)]),by="climID")
dataX2 <- merge(dataX,coords2[, .(Plot_ID,siteID)],by="Plot_ID")
length(unique(dataX2$siteID))

library(r3PG)
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(ggpubr)


# setwd(pathX)
# setwd("C:/Users/minunno/Documents/github/3PGQcalibration")
data_site <- d_site#read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_site')
data_species <- d_species#read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_species')
data_climate1 <- read_excel('myData/INPUT_cal.xlsx', sheet= 'd_climate')
data_thinning <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_thinnings')
data_parameters <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parameters')
data_sizeDist <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_sizeDist')

###select randomly 3 sites for testing simulations
# data_siteSel <- c(122,123,1234)
#!!! filter NA
data_species <- data_species[-unique(which(is.na(data_species),arr.ind=T)[,1])]
plotIDx <- unique(data_species$Plot_ID)
data_site <- data_site[Plot_ID %in% plotIDx]

# plotIDx <- data_site[data_siteSel,]$Plot_ID
data_site$plotID0 <- as.numeric(as.factor(data_site$Plot_ID))
data_species <- merge(data_species,data_site[,.(Plot_ID,plotID0)],by="Plot_ID")
data_species$Plot_ID <- data_species$plotID0 
data_site$Plot_ID <- data_site$plotID0 
data_site$plotID0 <- NULL
data_species$plotID0 <- NULL

##filter NAs


plotIDsel <- 1:10
climIDs <- 1:5
data_species <- data_species[Plot_ID %in% plotIDsel,] 
data_site <- data_site[Plot_ID %in% plotIDsel]
data_site$Site_ID <- sample(climIDs,length(plotIDsel),replace = T)
# data_site$Plot_ID <-  match(plotIDx,data_site$Plot_ID)
# data_species$Plot_ID <- match(data_species$Plot_ID,plotIDx)#,data_site$Plot_ID)
climateData <- climateData[climateData$siteID %in% unique(data_site$Site_ID)]

###obsData 
obsData <- data.table(read_excel("myData/TABELLA_OBS.xlsx"))

Plot_ID <- data_site$Plot_ID
data_species$planted = as.character(data_species$planted)
data_species$fertility = 0.5#as.double(data_species$fertility)


##### set data #####

all_site = data_site[,5:12]
indX <- which(names(data_species) %in% c("species", "planted", "fertility", "stems_n", "biom_stem", "biom_root", "biom_foliage"))
all_species = data_species [,..indX]
all_climate = data_climate1 [,c(3:6,8:10)]
#all_thinning = data_thinning [,3:8]
all_parameters = data_parameters
all_sizeDist = data_sizeDist

##
my_out <- list()
for(i in Plot_ID){
  my_site <- data_site[data_site$Plot_ID==i,5:12]
  my_species <- data_species[data_species$Plot_ID == i,..indX]
  my_climate = climateData[,.(year,month,tmp_min,tmp_max,prcp,srad,frost_days)]#all_climate
  #my_thinning = data_thinning[data_thinning$Plot_ID==i,3:8]
  my_parameters = all_parameters
  my_sizeDist = all_sizeDist
  my_out[[i]] = run_3PG(
    site        = my_site,
    species     = my_species[which(my_species$biom_stem>0),c(2,1,7,6,3:5)],
    climate     = my_climate,
    # thinning    = my_thinning,
    parameters  = my_parameters,
    size_dist   = my_sizeDist,
    settings    = list(light_model = 2, transp_model = 2, phys_model = 2,
                       height_model = 1, correct_bias = 0, calculate_d13c = 0),
    check_input = TRUE, df_out = F)
  
  print(paste("siteID",i))
}

i_output[i_output$group_id==2,]
plot(my_out[[1]][,1,2,10],col=1)
points(my_out[[1]][,2,2,10],col=2)
points(my_out[[1]][,3,2,10],col=3)

plot(my_out[[1]][,1,6,1])
plot(my_out[[2]][,1,6,1])
plot(my_out[[3]][,1,6,1])




#' Title extractSims
#'extract data from a list of simulations and add to the data table
#' @param plotID plot ID
#' @param modOut list of model output
#' @param obsData tables with observed data and the "coordinates of data points
#' @param dataType data type: total trees, remaining, thinned, dead trees
#' @param colSobsData columns of obsData where the coordinates are stored
#'
#' @return
#' @export
#'
#' @examples
extractSims <- function(plotID, modOut,obsData,
                        dataType="total",colSobsData = 3:6){
  # select a plot
  plotX <- plotID
  # select 3pg out of that plot
  modOutX <- modOut[[plotX]]
  #select the index argument
  dataX <- obsData[Plot_ID == plotX & data_type==dataType]
  dataIndX <- as.matrix(dataX[,..colSobsData])
  ###extractvector of Data
  simsX <- modOutX[dataIndX]
  ###add data to table
  dataX$sims <- simsX
  dataX$obs <- as.numeric(dataX$obs)
  return(dataX)
}

#example:
##plot all sites
dataForPlot <- data.table()
for(plotID in 1:57){
  dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
  dataForPlot <- rbind(dataForPlot,dataX)
}
varXs <- unique(dataTest$var_name)
for(varX in varXs){
  print(ggplot(dataForPlot[var_name==varX],aes(x=obs, y=sims)) + ggtitle(varX)+
          geom_point() + geom_abline(slope = 1,intercept = 0))
}


# extract data and plot for a single site
plotID = 1
dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)

coordX <- unique(dataX[,5:6])
for(i in 1:nrow(coordX)){
  groupX <- unlist(coordX[i,1])
  variableX <- unlist(coordX[i,2])
  plot(my_out[[plotID]][,,groupX,variableX],main=varXs[i])
  dataX[group==groupX & variable==variableX,points(n_month,obs,col=2,pch=20)]
}


#' Title extractSims
#'extract data from a list of simulations and add to the data table
#' @param plotID plot ID
#' @param modOut list of model output
#' @param obsData tables with observed data and the "coordinates of data points
#' @param dataType data type: total trees, remaining, thinned, dead trees
#' @param colSobsData columns of obsData where the coordinates are stored
#'
#' @return
#' @export
#'
#' @examples
extractSims <- function(plotID, modOut,obsData,
                        dataType="total",colSobsData = 3:6){
  # select a plot
  plotX <- plotID
  # select 3pg out of that plot
  modOutX <- modOut[[plotX]]
  #select the index argument
  dataX <- obsData[Plot_ID == plotX & data_type==dataType]
  dataIndX <- as.matrix(dataX[,..colSobsData])
  ###extractvector of Data
  simsX <- modOutX[dataIndX]
  ###add data to table
  dataX$sims <- simsX
  dataX$obs <- as.numeric(dataX$obs)
  return(dataX)
}

#example:
##plot all sites
dataForPlot <- data.table()
for(plotID in 1:57){
  dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
  dataForPlot <- rbind(dataForPlot,dataX)
}
varXs <- unique(dataForPlot$var_name)
for(varX in varXs){
  print(ggplot(dataForPlot[var_name==varX],aes(x=obs, y=sims)) + ggtitle(varX)+
          geom_point() + geom_abline(slope = 1,intercept = 0))
}


# extract data and plot for a single site
plotID = 1
dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)

coordX <- unique(dataX[,5:6])
for(i in 1:nrow(coordX)){
  groupX <- unlist(coordX[i,1])
  variableX <- unlist(coordX[i,2])
  plot(my_out[[plotID]][,,groupX,variableX],main=varXs[i])
  dataX[group==groupX & variable==variableX,points(n_month,obs,col=2,pch=20)]
}

#' extractData3PG
#'
#' @param out array with 3PG output
#' @param siteX site ID for which to extract the data
#' @param varX character string indicating the variable to extract
#'
#' @return vector of model output
#' @export
#'
#' @examples
extractData3PG <- function(out,plotX,varX){
  indX <- which(i_output==varX,arr.ind=T)
  groupID <- i_output[indX[1],]$group_id
  varID <- i_output[indX[1],]$variable_id
  outX <- out[[plotX]][,,groupID,varID]
  if(all(is.na(outX))){
    stop("check output variable. ",
         "choose variable betwee these: ",
         i_output[,3])
  }
  return(outX)
}

######### DATA'S SIMULATION #############
plotX = c(1:length(my_out))
varXs <- unique(obsData$var_name)
n_months =read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'n_month')

datax_tab = data.table()
for (i in varXs) {   
  g = which(i_output[,3]==i) #variable row
  j = i_output[g,c(1,2)] #####group and variable ID####
  for (n in plotX) {
    p = which(n_months$Plot_ID==n)
    
    datax_tab = rbind(datax_tab, data.table(Plot_ID = n,
                                            n_month = n_months$n_mese[which(n_months$Plot_ID==n)],
                                            group = j[1],
                                            variable = j[2],
                                            value = my_out[[n]][n_months$n_mese[which(n_months$Plot_ID==n)],,as.numeric(j[1]),as.numeric(j[2])],
                                            var_name = i_output[g,3]))
    
    print(paste("Plot_ID",i))
  }
}

setnames(datax_tab,c('Plot_ID', 'n_month','group','variable','sim','var_name'))
datax_tab$data_type = "simulated"
obsXtab = obsData[,c(2:3,5:9)]
setnames(obsXtab,c('Plot_ID', 'n_month','group','variable','obs','var_name','data_type'))
setkey(datax_tab, Plot_ID, n_month,group,variable,var_name)
setkey(obsXtab, Plot_ID, n_month,group,variable,var_name)
all_data_tab = merge(datax_tab,obsXtab[data_type=='total'])
plot(all_data_tab$obs, all_data_tab$sim)
all_data_tab[var_name=='stems_n',plot(obs,sim)]
i = 1
all_data_tab[var_name ==varXs[i],plot(obs,sim,main=varXs[i])]


#per ogni dato metti in tab mese, layer, gruppo ,variabile
#calcola il residual= estrai tutti i dati simulati che corrispondono agli osservati (usa il n_mese)
#### my_out[[Plot_n]][month,layer,group_id,var_id] #####

################### RESIDUALS ##################################
####### Select a variable to calculate residuals ########
print(varXs)
variablex = 'height'

#Residuals_tab =data.table()

g = which(i_output[,3]==variablex) 
c1 =  as.numeric(i_output[g,1])
c2 = as.numeric(i_output[g,2])
sim_data=data.table()
for (p in Plot_ID) {
  value_month = obsData[obsData$Plot_ID==p & obsData$var_name == variablex  &obsData$data_type == 'total',3]
  sim_tab = my_out[[p]]
  
  sim_data = rbind(sim_data, data.table(Plot_ID = p,
                                        n_month = value_month,
                                        value = my_out[[p]][unlist(value_month),1,c1,c2])
  )
}


obs_variablex = as.numeric(unlist(obsData[obsData$var_name == variablex & obsData$data_type == 'total',7]))
variablex_res = data.table(Plot_ID = sim_data$Plot_ID,
                           n_months = sim_data$n_month.n_month,
                           var_name = variablex,
                           sim_value= sim_data$value,
                           obs_value= obs_variablex)

variablex_res$residual_value =  variablex_res$obs_value -variablex_res$sim_value
Residuals_tab = rbind(Residuals_tab,variablex_res)
write.csv(Residuals_tab, 'myData/Residuals_tab.csv')
