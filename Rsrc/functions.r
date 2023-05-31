###functions

# Heavy tailed noraml distribution
Sivia_log<-function(diff,sd){
  # sd[which(sd<=0)]<-11e-6e-6
  diff[which(abs(diff)<=1e-6)]<-1e-6
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




#' Title
#'
#' @param site 
#' @param species 
#' @param thinning 
#' @param climate #weather inputs
#' @param obsData 
#' @param parameters 
#' 
#' @param outType 
#' modOut returns model output;
#' ll returns the log likelihood
#' datX returns the database of observed and simulated 
#'
#' @return
#' @export
#'
#' @examples
multi_r3pg <- function(site, species, thinning,climate, obsData,parameters,outType="ll"){
  #' @description simulate the n runs for a given site with the drawn parameter combination
  
  
  #' r3pg_int <- function( par_df){
  #'   #' @description function to run one site and return required output on standing biomass
  #'   #' @param par_df a data.frame of parameters
  
  out <- run_3PG(site = site, species = species, climate = climate,
                 thinning = thinning, parameters = parameters,
                 settings=settings_3pg, check_input = TRUE,
                 df_out = F)
  
  dataX <- obsData
  #remove NAs
  naObs <- which(is.na(as.numeric(dataX$obs)))
  if(length(naObs>1)) dataX <- dataX[-naObs]
  dataX$sims <- out[as.matrix(dataX[,3:6])]
  
  ll_stem <- calc_ll(dataX,"stems_n",pErr)
  ll_ba <- calc_ll(dataX,"basal_area",pErr)
  ll_v <- calc_ll(dataX,"volume",pErr)
  ll_d <- calc_ll(dataX,"dbh",pErr)
  ll_h <- calc_ll(dataX,"height",pErr)
  ll_Wstem <- calc_ll(dataX,"biom_stem",pErr)
  ll_Wroot <- calc_ll(dataX,"biom_root",pErr)
  ll_Wfol <- calc_ll(dataX,"biom_foliage",pErr)
  
  ll <- ll_stem + ll_ba + ll_v + ll_d + ll_h + ll_Wstem + ll_Wroot + ll_Wfol
  if(outType=="ll") return(ll)
  if(outType=="modOut") return(out)
  if(outType=="datX") return(dataX)
}


###likelihood 
logLike <- function(pX){
  
  # print("here")
  parameters$`Pinus sylvestris`[pIds] <- pX
  ll <- numeric(nSites)
  ll <-mclapply(sites, function(i,ll){
    # print(i)
    ll[i] <- multi_r3pg(allInputs$site[[i]],allInputs$species[[i]],
                        allInputs$thinning[[i]],allInputs$climate[[i]],
                        allInputs$obsData[[i]],parameters)
  },ll=ll,mc.cores = nCores)
  loglike <- sum(unlist(ll))
  return(loglike)
}

