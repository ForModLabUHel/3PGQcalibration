###functions

# Heavy tailed noraml distribution
Sivia_log<-function(diff,sd){
  # sd[which(sd<=0)]<-11e-6e-6
  diff[which(abs(diff)<=1e-6)]<-1e-4
  R2<-(diff/sd)^2
  prob<-1/(sd*(pi*2)^0.5)*(1-exp(-R2/2))/R2
  log(prob)
}
calc_ll <- function(dataX,varNameX,pErr){
  llX <- 0
  if(nrow(dataX[var_name==varNameX])>0){
    subData <- dataX[var_name==varNameX]
    resX <- as.numeric(subData$sims) - as.numeric(subData$obs)
    sdX <- pErr[1] + pErr[2] * as.numeric(subData$sims)
    llX <- sum(Sivia_log(resX,sdX))
  }
  return(llX)
}




multi_r3pg <- function(inputs, climate, obsData,parameters,pErr,
                         pars_soilQlitter,pars_SoilInit,pErrSoil,outType="ll"){
  #' @description simulate the n runs for a given site with the drawn parameter combination
  
  #' r3pg_int <- function( par_df){
  #'   #' @description function to run one site and return required output on standing biomass
  #'   #' @param par_df a data.frame of parameters

  # species <- species[which(stems_n>0)]
  climIDi <- inputs$site$climID
  Plot_IDi <- inputs$site$Plot_ID
  climate = climate[climID==climIDi,.(year,month,tmp_min,tmp_max,prcp,srad,frost_days)]
  obsData = obsAll[Plot_ID==Plot_IDi]
  
  site <- inputs$site[,.(latitude,altitude,soil_class, asw_i,asw_min, asw_max, from,to)]
  species <- inputs$species[,.(species, planted, fertility, stems_n, biom_stem, biom_root, biom_foliage)]
  
  layerX <- as.numeric(which(species$biom_stem>0))
  speciesX <- species[layerX,]
  
  out <- run_3PG(site = site, species = speciesX, climate = climate,
                 parameters = parameters,
                 settings= list(light_model = 2, transp_model = 2, phys_model = 2,
                                height_model = 1, correct_bias = 0, calculate_d13c = 0),
                 parsQlitter = pars_soilQlitter,
                 parsQsoc = pars_SoilInit,
                 # check_input = T,
                 df_out = F)
  out[which(is.na(out))] <- 0
  dataX <- obsData[,.(simMonth,layerID,groupID,variableID,obs,var_name)]
  #remove NAs
  naObs <- which(is.na(as.numeric(dataX$obs)))
  if(length(naObs>1)) dataX <- dataX[-naObs]
  dataX$sims <- out[as.matrix(dataX[,1:4])]
  
  ####Soil C (initialization and processing)
  initSoilCx <- as.numeric(initSoil[Plot_IDi]$value)
  obsSoilCx <- as.numeric(obsSoil[Plot_IDi]$value)
  
  if(!is.na(initSoilCx) & !is.na(obsSoilCx)){
    nMonths <- dim(out)[1]
    soilCsoil <- Q_soc_m(SOC_inputs=rep(0,nMonths),
                         beta=as.numeric(pars_SoilInit[1,2]),
                         eta_11=as.numeric(pars_SoilInit[2,2]),
                         e0=as.numeric(pars_SoilInit[3,2]),
                         fc=as.numeric(pars_SoilInit[4,2]),
                         u0=0.113,
                         q0=as.numeric(pars_SoilInit[6,2]),
                         initSoil=initSoilCx)$SOCt
    soilCtrees <- apply(out [,,9,6:9],1,sum)
    soilC <- soilCtrees + soilCsoil
    soilCsim <- soilCsoil[120]
    resSoil <- soilCsim - obsSoilCx
    ll_soil <- dnorm(resSoil,mean = pErrSoil[1],sd=pErrSoil[2],log = T)
    if(outType=="datX"){
      dataSoil <- data.table(
        simMonth=120,
        layerID=NA,
        groupID=NA,
        variableID=90,
        obs = obsSoilCx,
        var_name ="soilC",
        sims = soilCsim
      )
      dataX <- rbind(dataX,dataSoil)    
    }
  }else{
    ll_soil <- 0
  }
  ll_stem <- calc_ll(dataX,"N",pErr[1:2])
  ll_ba <- calc_ll(dataX,"BA",pErr[3:4])
  ll_v <- calc_ll(dataX,"V",pErr[5:6])
  ll_d <- calc_ll(dataX,"D",pErr[7:8])
  ll_Wstem <- calc_ll(dataX,"Ws",pErr[9:10])
  ll_Wroot <- calc_ll(dataX,"Wr",pErr[11:12])
  ll_Wfol <- calc_ll(dataX,"Wf",pErr[13:14])
  
  ll <- ll_soil + ll_stem + ll_ba + ll_v + ll_d + ll_Wstem + ll_Wroot + ll_Wfol
  if(outType=="ll") return(ll)
  if(outType=="modOut") return(out)
  if(outType=="datX") return(dataX)
}

###likelihood
logLike1 <- function(siteXs,pValues){
  # print("here")

  parameters$`Pinus sylvestris`[pIds] <- pValues[1:length(pIds)]
  parameters$`Picea abies`[pIds] <- pValues[(length(pIds)+1):(length(pIds)*2)]
  parameters$`Pinus contorta`[pIds] <- pValues[(length(pIds)*2+1):(length(pIds)*3)]
  parameters$`Betula alba`[pIds[-noDecid_par]] <- pValues[(length(pIds)*3+1):(length(pIds)*4-1)]
  parameters$`other deciduous`[pIds[-noDecid_par]] <- pValues[(length(pIds)*4+1-1):(length(pIds)*5-2)]
  pErr <- pValues[(length(pIds)*5 +1-2):(length(pIds)*5 + 14-2)]
  
  ###soilParameters
  ind_pSoil <- (length(pIds)*5-2 + length(pErr)+1):(length(pIds)*5-2 + length(pErr)+8)
  soil_pars <- pValues[ind_pSoil]
  pSoilX <- c(soil_pars[1:4],0.,soil_pars[5:6])
  pars_soilQlitter[,2:6] <- matrix(pSoilX,nrow=28,ncol=5)
  pars_SoilInit[,2] <- pSoilX[1:6]
  pErrSoil <- soil_pars[7:8]
  
  ll <- lapply(inputs[siteXs],FUN=multi_r3pg,
               climate=climateData,
               obsData=obsAll,
               parameters=parameters,
               pErr=pErr,
               pars_soilQlitter=pars_soilQlitter,
               pars_SoilInit = pars_SoilInit,
               pErrSoil = pErrSoil
         )

  loglike <- sum(unlist(ll))
  return(loglike)
}

likelihood <- function(pValues){
  logLike <- mclapply(sets, function(jx) {
    logLike1(jx,pValues)  
  }, mc.cores = nCores)  
  llP <- sum(unlist(logLike))
  return(llP)
}


runModelOut <- function(pValues,siteXs){
  ####runModel for multiple sites and return a list of databases for each site with the simulated ond observed data
  parameters$`Pinus sylvestris`[pIds] <- pValues[1:length(pIds)]
  parameters$`Picea abies`[pIds] <- pValues[(length(pIds)+1):(length(pIds)*2)]
  parameters$`Pinus contorta`[pIds] <- pValues[(length(pIds)*2+1):(length(pIds)*3)]
  parameters$`Betula alba`[pIds[-noDecid_par]] <- pValues[(length(pIds)*3+1):(length(pIds)*4-1)]
  parameters$`other deciduous`[pIds[-noDecid_par]] <- pValues[(length(pIds)*4+1-1):(length(pIds)*5-2)]
  pErr <- pValues[(length(pIds)*5 +1-2):(length(pIds)*5 + 14-2)]
  
  ###soilParameters
  ind_pSoil <- (length(pIds)*5-2 + length(pErr)+1):(length(pIds)*5-2 + length(pErr)+8)
  soil_pars <- pValues[ind_pSoil]
  pSoilX <- c(soil_pars[1:4],0.,soil_pars[5:6])
  pars_soilQlitter[,2:6] <- matrix(pSoilX,nrow=28,ncol=5)
  pars_SoilInit[,2] <- pSoilX[1:6]
  pErrSoil <- soil_pars[7:8]
  
  ll <- lapply(inputs[siteXs],FUN=multi_r3pg,
               climate=climateData,
               obsData=obsAll,
               parameters=parameters,
               pErr=pErr,
               pars_soilQlitter=pars_soilQlitter,
               pars_SoilInit = pars_SoilInit,
               pErrSoil = pErrSoil,
               outType = "datX"
  )
}


