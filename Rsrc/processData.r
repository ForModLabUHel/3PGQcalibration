library(r3PG)
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(ggpubr)
library(multidplyr)
library(tidyr)
library(purrr)

###select path for data and set to working directory
# pathX <- "C:/Users/39348/OneDrive/Documents/Github/3PGQslu"
# pathX <- "C:/Users/minunno/Documents/Github/3PGQslu"
# setwd(pathX)

data_site <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_site')
data_species <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_species')
data_climate <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet= 'd_climate')
data_thinning <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_thinnings')
data_thinning$age <- data_thinning$age + 0.3
data_parameters <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_parameters')
data_sizeDist <- read_excel('myData/INPUT_R_ALL_SITES.xlsx', sheet = 'd_sizeDist')

###obsData 
obsData <- data.table(read_excel("myData/TABELLA_OBS.xlsx"))
nas <- which(is.na(as.numeric(obsData$obs)))
obsData <- obsData[-nas]
varNames <- unique(obsData$var_name)

Plot_ID <- data_site$Plot_ID
data_species$planted = as.character(data_species$planted)
data_species$fertility = as.double(data_species$fertility)


##### set data #####
all_site = data_site[,5:12]
all_species = data_species [,4:9]
all_climate = data_climate [,c(2:6,8:10)]
all_thinning = data_thinning [,3:8]
all_parameters = data_parameters
all_sizeDist = data_sizeDist

climate.list <- site.list <- species.list <- thinning.list <- obsData.list <- list()
parameters <- data_parameters
sizeDist <- data_sizeDist
climate <- all_climate
settings_3pg  = list(light_model = 2, transp_model = 2, phys_model = 2,
                        height_model = 1, correct_bias = 0, calculate_d13c = 0)

##
for(i in Plot_ID){
  site.list[[i]] <- data_site[data_site$Plot_ID==i,5:12]
  species.list[[i]] <- data_species[data_species$Plot_ID == i,4:10]
  thinning.list[[i]] = data_thinning[data_thinning$Plot_ID==i,3:8]
  obsData.list[[i]] <- obsData[data_type=="total" & Plot_ID==i]
  climate.list[[i]] <- all_climate[all_climate$Plot_ID==i,]
}

nSites <- length(site.list)
site.gridX <- tibble(
  grid_id = numeric(nSites),
  site = numeric(nSites),
)
site.gridX$grid_id <- 1:nSites
site.gridX$site <- site.list

thinning.gridX <- tibble(
  grid_id = numeric(nSites),
  thinning = numeric(nSites),
)
thinning.gridX$grid_id <- 1:nSites
thinning.gridX$thinning <- thinning.list

species.gridX <- tibble(
  grid_id = numeric(nSites),
  species = numeric(nSites),
)
species.gridX$grid_id <- 1:nSites
species.gridX$species <- species.list

climate.gridX <- tibble(
  grid_id = numeric(nSites),
  climate = numeric(nSites),
)
climate.gridX$grid_id <- 1:nSites
climate.gridX$climate <- climate.list

obsData.gridX <- tibble(
  grid_id = numeric(nSites),
  obsData = numeric(nSites),
)
obsData.gridX$grid_id <- 1:nSites
obsData.gridX$obsData <- obsData.list

allInputs <- inner_join(site.gridX,species.gridX,by="grid_id")
allInputs <- inner_join(allInputs,thinning.gridX,by="grid_id")
allInputs <- inner_join(allInputs,obsData.gridX,by="grid_id")
allInputs <- inner_join(allInputs,climate.gridX,by="grid_id")

save(allInputs,parameters,settings_3pg,sizeDist,varNames,
     file="myData/input3pgLists.rdata")


