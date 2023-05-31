library(BayesianTools)
library(coda)
library(r3PG)


source("Rsrc/functions.r")
load("calOut/calibration.rdata")
load('myData/input3pgLists.rdata')


###make traceplot
pdf("outCal/tracePlots.pdf")
tracePlot(calibration$chain,start = 5000)
dev.off()

pdf("outCal/marginalPlots.pdf")
marginalPlot(calibration,start = 5000,prior=T)
dev.off()

###first select the parameters that are highly correlated
corrX <- cor(calibration$Z[1e5:dim(calibration$Z)[1],])
diag(corrX) <- NA
indPars <- unique(which(abs(corrX)>0.5,arr.ind=T)[,1])
###create correlation plot
pdf("outCal/correlationPlots.pdf")
correlationPlot(calibration,start=5000,whichParameters=indPars)
dev.off()


####runModel for multiple sites and report likelihood
ll <- numeric(nSites)
ll <-mclapply(sites, function(i,ll){
  print(i)
  ll[i] <- multi_r3pg(allInputs$site[[i]],allInputs$species[[i]],
                      allInputs$thinning[[i]],allInputs$climate[[i]],
                      allInputs$obsData[[i]],parameters)
},ll=ll,mc.cores = nCores)


####runModel for multiple sites and return a list of databases for each site with the simulated ond observed data
datX <- list()
datX <-mclapply(sites, function(i,datX){
  print(i)
  datX[[i]] <- multi_r3pg(allInputs$site[[i]],allInputs$species[[i]],
                          allInputs$thinning[[i]],allInputs$climate[[i]],
                          allInputs$obsData[[i]],parameters,outType="datX")
},datX=datX,mc.cores = nCores)

