library(r3PG)
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(ggpubr)
library(rgdal)

###read the data
dataX <- data.table(read_excel("myData/U22128_fm.xlsx", sheet = "U22128_MI"))
setnames(dataX,c("Unique plot-ID","Year of inventory","Dry weight kg/ha stump incl roots > 2mm all species"),c("Plot_ID","year","Wr"))

####convert coordinates
coords <- dataX[,.(Plot_ID,year,`Coordinate SWEREF 99 TM East. Distorted +-200-800 m`,`Coordinate SWEREF 99 TM North. Distorted +-200-800 m`)] #original coordinates
setnames(coords,c("Plot_ID","year","long","lat"))
temp <- coords[,3:4]
temp <- SpatialPoints (temp, proj4string = CRS ('+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

coords[,3:4] <- as.data.frame (spTransform(temp, CRS('+proj=longlat +datum=WGS84 +no_defs'))) ####convert
setkey(dataX,Plot_ID,year)
setkey(coords,Plot_ID,year)
dataX <- merge(dataX,coords,by=c("Plot_ID","year")) #add new coordinates to the data

###read and process weather inputs
climateData <- data.table(read.csv('myData/monthly_weather.csv'))
coords2 <- data.table(read.csv('myData/coords_climid.csv'))
# setkey(data_climate2, climID)
# setkey(coords2, climID)
# coords2$siteID <- as.numeric(as.factor(coords2$climID))
# coords2$siteID <- as.factor(coords2$climID)
# climateData <- merge(data_climate2,unique(coords2[, .(climID,siteID)]),by="climID")
dataX <- merge(dataX,coords2[, .(Plot_ID,climID)],by="Plot_ID")
length(unique(dataX$climID))

###select only the plots that have a remeasurement
selPlots <- dataX[,.SD[which(length(year)>1)],by=Plot_ID]$Plot_ID
selData <- dataX[Plot_ID %in% selPlots]

###this is a bit stupid and could be more elegant but works
initVals <- selData[, .SD[which.min(year)], by = Plot_ID] ####it's the data of the 1st measurement for each plot that will be used to initialize the model
obs <- selData[, .SD[which.max(year)], by = Plot_ID] ####it's the data of the 2nd measurement for each plot that will be used to test/calibrate the model
initVals$dataUse="init"
obs$dataUse="cal"
selData <- rbind(initVals,obs)
setkey(selData,Plot_ID,year)

selData[, initYear:=.SD[which.min(year)], by = Plot_ID] ####it's the data of the 1st measurement for each plot that will be used to initialize the model
selData[, simMonth:=(year-initYear)*12, by = Plot_ID] ####it's the data of the 1st measurement for each plot that will be used to initialize the model


##sites with 0 data
site0 <- unique(selData$Plot_ID[which(selData$`BA/ha m2/ha all species`==0)])
length(site0)
selData <- selData[! Plot_ID %in% site0]
#####prepare the data in the model format input
d_site  <- data.table(siteName=selData[dataUse=="init"]$`Plot ID within cluster`, 
                      Site_ID=selData[dataUse=="init"]$`Plot ID within cluster`,
                      Plot_Name= selData[dataUse=="init"]$Plot_ID,
                      Plot_ID=selData[dataUse=="init"]$Plot_ID,
                      latitude= selData[dataUse=="init"]$lat,
                      climID= selData[dataUse=="init"]$climID,
                      altitude=100,
                      soil_class = 3,
                      asw_i=999,
                      asw_min=0,
                      asw_max=160,
                      from=paste0(selData[dataUse=="init"]$year,"-01"),
                      to=paste0(selData[dataUse=="cal"]$year+1,"-01")
                      
)

selData[,planted:=paste0(year-`Stand age`,"-01")] ###year of planting

###avoid NaN
# selData[`Dry weight kg/ha branches incl small branches and needles all species`==0]$`Dry weight kg/ha branches incl small branches and needles all species`=0.0001
# selData[`Dry weight kg/ha stem incl bark all species`==0]$`Dry weight kg/ha stem incl bark all species` = 0.001
# selData[Wr==0]$Wr=0.0001

####partitioning biomasses by species (start)
selData[,WsPine:=`Vol/ha m3sk/ha Pine`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WsSpruce:=`Vol/ha m3sk/ha Spruce`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WsBirch:=`Vol/ha m3sk/ha Birch`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WsContorta:=`Vol/ha m3sk/ha Contorta`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]
selData[,WsDecid:=`Vol/ha m3sk/ha Other dec.`/`Vol/ha m3sk/ha all species`*`Dry weight kg/ha stem incl bark all species`/1000]

selData[,WfPine:=`BA/ha m2/ha Pine`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfSpruce:=`BA/ha m2/ha Spruce`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfBirch:=`BA/ha m2/ha aBirch`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfContorta:=`BA/ha m2/ha Contorta`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]
selData[,WfDecid:=`BA/ha m2/ha Other dec.`/`BA/ha m2/ha all species`*`Dry weight kg/ha branches incl small branches and needles all species`/1000]

selData[,WrPine:=`BA/ha m2/ha Pine`/`BA/ha m2/ha all species`*Wr/1000]
selData[,WrSpruce:=`BA/ha m2/ha Spruce`/`BA/ha m2/ha all species`*Wr/1000]
selData[,WrBirch:=`BA/ha m2/ha aBirch`/`BA/ha m2/ha all species`*Wr/1000]
selData[,WrContorta:=`BA/ha m2/ha Contorta`/`BA/ha m2/ha all species`*Wr/1000]
selData[,WrDecid:=`BA/ha m2/ha Other dec.`/`BA/ha m2/ha all species`*Wr/1000]

setnames(selData,"Plot ID within cluster","siteName")
Ws <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                        WsPine,
                                        WsSpruce,
                                        WsBirch,
                                        WsContorta,
                                        WsDecid)],id.vars=1:3)
Ws[variable=="WsPine",layerID:=1]
Ws[variable=="WsPine",species:="Pinus sylvestris"]
Ws[variable=="WsSpruce",layerID:=2]
Ws[variable=="WsSpruce",species:="Picea abies"]
Ws[variable=="WsContorta",layerID:=3]
Ws[variable=="WsContorta",species:="Pinus contorta"]
Ws[variable=="WsBirch",layerID:=4]
Ws[variable=="WsBirch",species:="Betula alba"]
Ws[variable=="WsDecid",layerID:=5]
Ws[variable=="WsDecid",species:="other deciduous"]

Wr <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                        WrPine,
                                        WrSpruce,
                                        WrBirch,
                                        WrContorta,
                                        WrDecid)],id.vars=1:3)
Wr[variable=="WrPine",layerID:=1]
Wr[variable=="WrPine",species:="Pinus sylvestris"]
Wr[variable=="WrSpruce",layerID:=2]
Wr[variable=="WrSpruce",species:="Picea abies"]
Wr[variable=="WrContorta",layerID:=3]
Wr[variable=="WrContorta",species:="Pinus contorta"]
Wr[variable=="WrBirch",layerID:=4]
Wr[variable=="WrBirch",species:="Betula alba"]
Wr[variable=="WrDecid",layerID:=5]
Wr[variable=="WrDecid",species:="other deciduous"]

Wf <- melt(selData[dataUse=="init",.(siteName,Plot_ID,planted,
                                           WfPine,
                                           WfSpruce,
                                           WfBirch,
                                           WfContorta,
                                           WfDecid)],id.vars=1:3)
Wf[variable=="WfPine",layerID:=1]
Wf[variable=="WfPine",species:="Pinus sylvestris"]
Wf[variable=="WfSpruce",layerID:=2]
Wf[variable=="WfSpruce",species:="Picea abies"]
Wf[variable=="WfContorta",layerID:=3]
Wf[variable=="WfContorta",species:="Pinus contorta"]
Wf[variable=="WfBirch",layerID:=4]
Wf[variable=="WfBirch",species:="Betula alba"]
Wf[variable=="WfDecid",layerID:=5]
Wf[variable=="WfDecid",species:="other deciduous"]
####partitioning biomasses by species (end)

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

setnames(Ws,"value","biom_stem");Ws$variable <- NULL
setnames(Wr,"value","biom_root");Wr$variable <- NULL
setnames(Wf,"value","biom_foliage");Wf$variable <- NULL
setnames(stems,"value","stems_n");stems$variable <- NULL

setkey(Ws,"siteName","Plot_ID","planted","layerID","species")
setkey(Wr,"siteName","Plot_ID","planted","layerID","species")
setkey(Wf,"siteName","Plot_ID","planted","layerID","species")
setkey(stems,"siteName","Plot_ID","planted","layerID","species")

d_species <- merge(Ws,Wr,by=c("siteName","Plot_ID","planted","layerID","species"))
d_species <- merge(d_species,Wf,by=c("siteName","Plot_ID","planted","layerID","species"))
d_species <- merge(d_species,stems,by=c("siteName","Plot_ID","planted","layerID","species"))


###still preparing the data
data_site <- d_site
data_species <- d_species
# data_thinning <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_thinnings')
data_parameters <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_parameters')
data_sizeDist <- read_excel('myData/INPUT_cal.xlsx', sheet = 'd_sizeDist')

#!!! filter NA
# data_species <- data_species[-unique(which(is.na(data_species),arr.ind=T)[,1])]
# plotIDx <- unique(data_species$Plot_ID)
# data_site <- data_site[Plot_ID %in% plotIDx]

data_site$plotID0 <- as.numeric(as.factor(data_site$Plot_ID))
data_species <- merge(data_species,data_site[,.(Plot_ID,plotID0)],by="Plot_ID")
data_species$Plot_ID <- data_species$plotID0 
data_site$Plot_ID <- data_site$plotID0 
data_site$plotID0 <- NULL
data_species$plotID0 <- NULL

###process observed data
selData <- selData[Plot_ID %in% data_site$Plot_Name]
setnames(selData,"Plot_ID","Plot_Name")
selData <- merge(selData,data_site[,.(Plot_Name,Plot_ID)],by="Plot_Name")
setnames(selData,c("Vol/ha m3sk/ha Pine", "Vol/ha m3sk/ha Contorta","Vol/ha m3sk/ha Spruce", "Vol/ha m3sk/ha Birch","Vol/ha m3sk/ha Other dec."),
         c("Vpine","Vcont", "Vspruce", "Vbirch","Vdecid"))# 
baNames <- c( "BA/ha m2/ha all species", "BA/ha m2/ha Pine", "BA/ha m2/ha Contorta", "BA/ha m2/ha Spruce", "BA/ha m2/ha aBirch","BA/ha m2/ha Other dec.")
dNames <- c("Mean DBH (mm) all species", "Mean DBH (mm) Pine", "Mean DBH (mm) Contorta",
            "Mean DBH (mm) Spruce", "Mean DBH (mm) Birch", "Mean DBH (mm) Other dec.") 
nNames <- c("no of stems/ha all species", "no of stems/ha  Pine", "no of stems/ha  Contorta", 
          "no of stems/ha  Spruce", "no of stems/ha Birch", "no of stems/ha Other dec." )
setnames(selData,baNames, c("BAtot","BApine","BAcont","BAspruce","BAbirch","BAdecid"))
setnames(selData,dNames, c("Dtot","Dpine","Dcont","Dspruce","Dbirch","Ddecid"))
setnames(selData,nNames, c("Ntot","Npine","Ncont","Nspruce","Nbirch","Ndecid"))

obsBA <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                         BApine,BAcont,BAspruce,BAbirch,BAdecid)],id.vars = c("Plot_ID","simMonth"))
obsV <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                       Vpine,Vcont,Vspruce,Vbirch,Vdecid)],id.vars = c("Plot_ID","simMonth"))
obsN <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                       Npine,Ncont,Nspruce,Nbirch,Ndecid)],id.vars = c("Plot_ID","simMonth"))
obsD <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                      Dpine,Dcont,Dspruce,Dbirch,Ddecid)],id.vars = c("Plot_ID","simMonth"))

obsWs <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                      WsPine,WsContorta,WsSpruce,WsBirch,WsDecid)],id.vars = c("Plot_ID","simMonth"))
obsWf <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                       WfPine,WfContorta,WfSpruce,WfBirch,WfDecid)],id.vars = c("Plot_ID","simMonth"))
obsWr <- melt(selData[dataUse=="cal",.(Plot_ID,simMonth,
                                       WrPine,WrContorta,WrSpruce,WrBirch,WrDecid)],id.vars = c("Plot_ID","simMonth"))

obsD$value <- obsD$value * 0.1 #converts mm to cm 

###BA
obsBA[variable=="BApine",layer:=1]
obsBA[variable=="BAspruce",layer:=2]
obsBA[variable=="BAcont",layer:=3]
obsBA[variable=="BAbirch",layer:=4]
obsBA[variable=="BAdecid",layer:=5]

###V
obsV[variable=="Vpine",layer:=1]
obsV[variable=="Vspruce",layer:=2]
obsV[variable=="Vcont",layer:=3]
obsV[variable=="Vbirch",layer:=4]
obsV[variable=="Vdecid",layer:=5]
###N
obsN[variable=="Npine",layer:=1]
obsN[variable=="Nspruce",layer:=2]
obsN[variable=="Ncont",layer:=3]
obsN[variable=="Nbirch",layer:=4]
obsN[variable=="Ndecid",layer:=5]
###D
obsD[variable=="Dpine",layer:=1]
obsD[variable=="Dspruce",layer:=2]
obsD[variable=="Dcont",layer:=3]
obsD[variable=="Dbirch",layer:=4]
obsD[variable=="Ddecid",layer:=5]
###Ws
obsWs[variable=="WsPine",layer:=1]
obsWs[variable=="WsSpruce",layer:=2]
obsWs[variable=="WsContorta",layer:=3]
obsWs[variable=="WsBirch",layer:=4]
obsWs[variable=="WsDecid",layer:=5]
###Wf
obsWf[variable=="WfPine",layer:=1]
obsWf[variable=="WfSpruce",layer:=2]
obsWf[variable=="WfContorta",layer:=3]
obsWf[variable=="WfBirch",layer:=4]
obsWf[variable=="WfDecid",layer:=5]
###Wr
obsWr[variable=="WrPine",layer:=1]
obsWr[variable=="WrSpruce",layer:=2]
obsWr[variable=="WrContorta",layer:=3]
obsWr[variable=="WrBirch",layer:=4]
obsWr[variable=="WrDecid",layer:=5]


obsBA[,groupID:=2]; obsBA[,variableID:=3]
obsBA <- obsBA[value!=0]
obsBA[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]

obsN[,groupID:=2]; obsN[,variableID:=2]
obsN <- obsN[value!=0]
obsN[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]
obsV[,groupID:=2]; obsV[,variableID:=10]
obsV <- obsV[value!=0]
obsV[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]
obsD[,groupID:=2]; obsD[,variableID:=5]
obsD <- obsD[value!=0]
obsD[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]
obsWs[,groupID:=4]; obsWs[,variableID:=1]
obsWs <- obsWs[value!=0]
obsWs[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]
obsWf[,groupID:=4]; obsWf[,variableID:=3]
obsWf <- obsWf[value!=0]
obsWf[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]
obsWr[,groupID:=4]; obsWr[,variableID:=2]
obsWr <- obsWr[value!=0]
obsWr[,layerID:=as.numeric(as.factor(layer)),by=Plot_ID]

obsAll <- rbind(obsBA,obsD);obsAll <- rbind(obsAll,obsD);obsAll <- rbind(obsAll,obsN)
obsAll <- rbind(obsAll,obsV);obsAll <- rbind(obsAll,obsWs);obsAll <- rbind(obsAll,obsWf);obsAll <- rbind(obsAll,obsWr)
# setnames(obsAll,"value","obs")
# plotIDsel <- 1:10
# climIDs <- 1:5
# data_species <- data_species[Plot_ID %in% plotIDsel,] 
# data_site <- data_site[Plot_ID %in% plotIDsel]
# data_site$Site_ID <- sample(climIDs,length(plotIDsel),replace = T)
# data_site$Plot_ID <-  match(plotIDx,data_site$Plot_ID)
# data_species$Plot_ID <- match(data_species$Plot_ID,plotIDx)#,data_site$Plot_ID)
# climateData <- climateData[climateData$siteID %in% unique(data_site$Site_ID)]

###obsData 
# obsData <- data.table(read_excel("myData/TABELLA_OBS.xlsx"))

Plot_ID <- data_site$Plot_ID
data_species$planted = as.character(data_species$planted)
data_species$fertility = 0.5#as.double(data_species$fertility)


##### set data #####

# all_site = data_site[,5:12]
# all_species = data_species [,..indX]
# all_climate = data_climate2[,c(3:6,8:10)]
# #all_thinning = data_thinning [,3:8]
all_parameters = data_parameters
all_sizeDist = data_sizeDist

site_list <- split(data_site,data_site$Plot_ID)
species_list <- split(data_species,data_species$Plot_ID)
indXsite <- which(names(data_site) %in% c("latitude","altitude","soil_class","asw_i","asw_min","asw_max","from","to"))
indXspecies <- which(names(data_species) %in% c("species", "planted", "fertility", "stems_n", "biom_stem", "biom_root", "biom_foliage"))

save(obsAll,site_list,species_list,file="myData/dataInputs.rdata")

##
nLayersInit <- d_species[,length(which(stems_n>0)),by=Plot_ID]
setnames(nLayersInit,"Plot_ID","Plot_Name")
nLayersInit <- merge(nLayersInit,selData[dataUse=="init",.(Plot_Name,Plot_ID)])
nLayersObs <- obsAll[variableID==2 & groupID==2,length(value),by=Plot_ID]
setnames(nLayersInit,"V1","nLayersInit")
setnames(nLayersObs,"V1","nLayersObs")
nLayers <- merge(nLayersInit,nLayersObs,by="Plot_ID")
nLayers[,dNlayers:=nLayersInit-nLayersObs]
hist(nLayers$dNlayers)

plotToCheck <- nLayers[dNlayers!=0]$Plot_ID
Plot_IDsel <- which(!Plot_ID %in% plotToCheck)
my_out <- my_obs <- list()
obsAll$sims <- -99999
for(i in Plot_IDsel){
  climIDi <- site_list[[i]]$climID
  my_site <- site_list[[i]][,..indXsite]
  my_species <- species_list[[i]][,..indXspecies]
  my_climate = climateData[climID==climIDi,.(year,month,tmp_min,tmp_max,prcp,srad,frost_days)]#all_climate
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
  
  # my_obs[[i]]
  simsX <- my_out[[i]][as.matrix(obsAll[Plot_ID==i,.(simMonth,layerID,groupID,variableID)])]
  obsAll[Plot_ID==i]$sims=unlist(simsX)
  
  print(paste("siteID",i))
}

obsAll[sims == -99999,sims:=NA]
obsAll[variable=="BAspruce",plot(sims,value)]
obsAll[variable=="BApine",hist(sims-value)]

varX=3
groupX=2
plotX=1
nameX <- i_output[i_output$group_id==groupX & i_output$variable_id==varX,]$variable_name
ylimX <- range(my_out[[plotX]][,,groupX,varX])
# i_output[i_output$group_id==groupX,]
plot(my_out[[plotX]][,1,groupX,varX],col=1,ylim=ylimX,main=nameX)
points(my_out[[plotX]][,2,groupX,varX],col=2)
points(my_out[[plotX]][,3,groupX,varX],col=3)

#' 
#' plot(my_out[[1]][,1,6,1])
#' plot(my_out[[2]][,1,6,1])
#' plot(my_out[[3]][,1,6,1])
#' 
#' 
#' 
#' 
#' #' Title extractSims
#' #'extract data from a list of simulations and add to the data table
#' #' @param plotID plot ID
#' #' @param modOut list of model output
#' #' @param obsData tables with observed data and the "coordinates of data points
#' #' @param dataType data type: total trees, remaining, thinned, dead trees
#' #' @param colSobsData columns of obsData where the coordinates are stored
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' extractSims <- function(plotID, modOut,obsData,
#'                         dataType="total",colSobsData = 3:6){
#'   # select a plot
#'   plotX <- plotID
#'   # select 3pg out of that plot
#'   modOutX <- modOut[[plotX]]
#'   #select the index argument
#'   dataX <- obsData[Plot_ID == plotX & data_type==dataType]
#'   dataIndX <- as.matrix(dataX[,..colSobsData])
#'   ###extractvector of Data
#'   simsX <- modOutX[dataIndX]
#'   ###add data to table
#'   dataX$sims <- simsX
#'   dataX$obs <- as.numeric(dataX$obs)
#'   return(dataX)
#' }
#' 
#' #example:
#' ##plot all sites
#' dataForPlot <- data.table()
#' for(plotID in 1:57){
#'   dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
#'   dataForPlot <- rbind(dataForPlot,dataX)
#' }
#' varXs <- unique(dataTest$var_name)
#' for(varX in varXs){
#'   print(ggplot(dataForPlot[var_name==varX],aes(x=obs, y=sims)) + ggtitle(varX)+
#'           geom_point() + geom_abline(slope = 1,intercept = 0))
#' }
#' 
#' 
#' # extract data and plot for a single site
#' plotID = 1
#' dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
#' 
#' coordX <- unique(dataX[,5:6])
#' for(i in 1:nrow(coordX)){
#'   groupX <- unlist(coordX[i,1])
#'   variableX <- unlist(coordX[i,2])
#'   plot(my_out[[plotID]][,,groupX,variableX],main=varXs[i])
#'   dataX[group==groupX & variable==variableX,points(n_month,obs,col=2,pch=20)]
#' }
#' 
#' 
#' #' Title extractSims
#' #'extract data from a list of simulations and add to the data table
#' #' @param plotID plot ID
#' #' @param modOut list of model output
#' #' @param obsData tables with observed data and the "coordinates of data points
#' #' @param dataType data type: total trees, remaining, thinned, dead trees
#' #' @param colSobsData columns of obsData where the coordinates are stored
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' extractSims <- function(plotID, modOut,obsData,
#'                         dataType="total",colSobsData = 3:6){
#'   # select a plot
#'   plotX <- plotID
#'   # select 3pg out of that plot
#'   modOutX <- modOut[[plotX]]
#'   #select the index argument
#'   dataX <- obsData[Plot_ID == plotX & data_type==dataType]
#'   dataIndX <- as.matrix(dataX[,..colSobsData])
#'   ###extractvector of Data
#'   simsX <- modOutX[dataIndX]
#'   ###add data to table
#'   dataX$sims <- simsX
#'   dataX$obs <- as.numeric(dataX$obs)
#'   return(dataX)
#' }
#' 
#' #example:
#' ##plot all sites
#' dataForPlot <- data.table()
#' for(plotID in 1:57){
#'   dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
#'   dataForPlot <- rbind(dataForPlot,dataX)
#' }
#' varXs <- unique(dataForPlot$var_name)
#' for(varX in varXs){
#'   print(ggplot(dataForPlot[var_name==varX],aes(x=obs, y=sims)) + ggtitle(varX)+
#'           geom_point() + geom_abline(slope = 1,intercept = 0))
#' }
#' 
#' 
#' # extract data and plot for a single site
#' plotID = 1
#' dataX <- extractSims(plotID = plotID,obsData = obsData,modOut = my_out)
#' 
#' coordX <- unique(dataX[,5:6])
#' for(i in 1:nrow(coordX)){
#'   groupX <- unlist(coordX[i,1])
#'   variableX <- unlist(coordX[i,2])
#'   plot(my_out[[plotID]][,,groupX,variableX],main=varXs[i])
#'   dataX[group==groupX & variable==variableX,points(n_month,obs,col=2,pch=20)]
#' }
#' 
#' #' extractData3PG
#' #'
#' #' @param out array with 3PG output
#' #' @param siteX site ID for which to extract the data
#' #' @param varX character string indicating the variable to extract
#' #'
#' #' @return vector of model output
#' #' @export
#' #'
#' #' @examples
#' extractData3PG <- function(out,plotX,varX){
#'   indX <- which(i_output==varX,arr.ind=T)
#'   groupID <- i_output[indX[1],]$group_id
#'   varID <- i_output[indX[1],]$variable_id
#'   outX <- out[[plotX]][,,groupID,varID]
#'   if(all(is.na(outX))){
#'     stop("check output variable. ",
#'          "choose variable betwee these: ",
#'          i_output[,3])
#'   }
#'   return(outX)
#' }
#' 
#' ######### DATA'S SIMULATION #############
#' plotX = c(1:length(my_out))
#' varXs <- unique(obsData$var_name)
#' n_months =read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'n_month')
#' 
#' datax_tab = data.table()
#' for (i in varXs) {   
#'   g = which(i_output[,3]==i) #variable row
#'   j = i_output[g,c(1,2)] #####group and variable ID####
#'   for (n in plotX) {
#'     p = which(n_months$Plot_ID==n)
#'     
#'     datax_tab = rbind(datax_tab, data.table(Plot_ID = n,
#'                                             n_month = n_months$n_mese[which(n_months$Plot_ID==n)],
#'                                             group = j[1],
#'                                             variable = j[2],
#'                                             value = my_out[[n]][n_months$n_mese[which(n_months$Plot_ID==n)],,as.numeric(j[1]),as.numeric(j[2])],
#'                                             var_name = i_output[g,3]))
#'     
#'     print(paste("Plot_ID",i))
#'   }
#' }
#' 
#' setnames(datax_tab,c('Plot_ID', 'n_month','group','variable','sim','var_name'))
#' datax_tab$data_type = "simulated"
#' obsXtab = obsData[,c(2:3,5:9)]
#' setnames(obsXtab,c('Plot_ID', 'n_month','group','variable','obs','var_name','data_type'))
#' setkey(datax_tab, Plot_ID, n_month,group,variable,var_name)
#' setkey(obsXtab, Plot_ID, n_month,group,variable,var_name)
#' all_data_tab = merge(datax_tab,obsXtab[data_type=='total'])
#' plot(all_data_tab$obs, all_data_tab$sim)
#' all_data_tab[var_name=='stems_n',plot(obs,sim)]
#' i = 1
#' all_data_tab[var_name ==varXs[i],plot(obs,sim,main=varXs[i])]
#' 
#' 
#' #per ogni dato metti in tab mese, layer, gruppo ,variabile
#' #calcola il residual= estrai tutti i dati simulati che corrispondono agli osservati (usa il n_mese)
#' #### my_out[[Plot_n]][month,layer,group_id,var_id] #####
#' 
#' ################### RESIDUALS ##################################
#' ####### Select a variable to calculate residuals ########
#' print(varXs)
#' variablex = 'height'
#' 
#' #Residuals_tab =data.table()
#' 
#' g = which(i_output[,3]==variablex) 
#' c1 =  as.numeric(i_output[g,1])
#' c2 = as.numeric(i_output[g,2])
#' sim_data=data.table()
#' for (p in Plot_ID) {
#'   value_month = obsData[obsData$Plot_ID==p & obsData$var_name == variablex  &obsData$data_type == 'total',3]
#'   sim_tab = my_out[[p]]
#'   
#'   sim_data = rbind(sim_data, data.table(Plot_ID = p,
#'                                         n_month = value_month,
#'                                         value = my_out[[p]][unlist(value_month),1,c1,c2])
#'   )
#' }
#' 
#' 
#' obs_variablex = as.numeric(unlist(obsData[obsData$var_name == variablex & obsData$data_type == 'total',7]))
#' variablex_res = data.table(Plot_ID = sim_data$Plot_ID,
#'                            n_months = sim_data$n_month.n_month,
#'                            var_name = variablex,
#'                            sim_value= sim_data$value,
#'                            obs_value= obs_variablex)
#' 
#' variablex_res$residual_value =  variablex_res$obs_value -variablex_res$sim_value
#' Residuals_tab = rbind(Residuals_tab,variablex_res)
#' write.csv(Residuals_tab, 'myData/Residuals_tab.csv')
