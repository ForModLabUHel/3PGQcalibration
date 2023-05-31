library(r3PG)
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(ggpubr)

###select path for data and set to working directory
pathX <- "C:/Users/39348/OneDrive/Documents/Github/3PGQslu"
# pathX <- "C:/Users/minunno/Documents/Github/3PGQslu"
setwd(pathX)

data_site <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_site')
data_species <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_species')
data_climate <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet= 'd_climate')
data_thinning <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_thinnings')
data_parameters <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_parameters')
data_sizeDist <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_sizeDist')

###obsData 
obsData <- data.table(read_excel("myData/TABELLA_OBS.xlsx"))

Plot_ID <- data_site$Plot_ID
data_species$planted = as.character(data_species$planted)
data_species$fertility = as.double(data_species$fertility)


##### set data #####

all_site = data_site[,5:12]
all_species = data_species [,4:9]
all_climate = data_climate [,2:9]
all_thinning = data_thinning [,3:8]
all_parameters = data_parameters
all_sizeDist = data_sizeDist

##
my_out <- list()
for(i in Plot_ID){
  my_site <- data_site[data_site$Plot_ID==i,5:12]
  my_species <- data_species[data_species$Plot_ID == i,4:10]
  my_climate = all_climate
  my_thinning = data_thinning[data_thinning$Plot_ID==i,3:8]
  my_parameters = all_parameters
  my_sizeDist = all_sizeDist
  my_out[[i]] = run_3PG(
    site        = my_site,
    species     = my_species,
    climate     = my_climate,
    thinning    = my_thinning,
    parameters  = my_parameters,
    size_dist   = my_sizeDist,
    settings    = list(light_model = 2, transp_model = 2, phys_model = 2,
                       height_model = 1, correct_bias = 0, calculate_d13c = 0),
    check_input = TRUE, df_out = F)
  
  print(paste("siteID",i))
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
