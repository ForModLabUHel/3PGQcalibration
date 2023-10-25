.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library(data.table)
library("runjags")
library("coda")

setwd("/scratch/project_2000994/calibrations/3PGQcalibration")
#setwd("C:/Users/minunno/Documents/yucatrote/SLU")

calSets <- 9

load(paste0("calOut/calibration_",calSets[1],".1.rdata"))
namesX <- c(calibration$setup$names,"lp","ll","pr")

###settings
npar <- calibration$setup$numPars
indRun <-9 #number of independent calibration runs
thin=1
pChain <- mcmc.list()
pMAP <- NULL

for(i in 1:indRun){
  ij=calSets# for(ij in calSets){
    load(paste0("calOut/calibration_",ij,".",i,".rdata"))
    lChain <- dim(calibration$chain[[1]])[1]
    seqX <- (lChain*0.7):lChain
    # seqX <- seq(round(lChain/2),lChain,by=1)
  #   if(ij==calSets[1]){
       chain1 <- calibration$chain[[1]][seqX,]
       colnames(chain1) <- namesX
       chain2 <- calibration$chain[[2]][seqX,]
       colnames(chain2) <- namesX
       chain3 <- calibration$chain[[3]][seqX,]
       colnames(chain3) <- namesX
  #   }else{
  #     chain1 <- rbind(chain1,calibration$chain[[1]])
  #     chain2 <- rbind(chain2,calibration$chain[[2]])
  #     chain3 <- rbind(chain3,calibration$chain[[3]])
  #   }
  # # }
  # #create a list of 3 chains. Note those 3 chains are the 3 parallel chains that are used in the DE (differential evolution algorithms)
  # mcmcList <- mcmc.list(mcmc(chain1),mcmc(chain2),mcmc(chain3))
  # 
  #combine the 3 parallel chain in one unique chain
  # pChain[[i]] <- mcmc(getSample(calibration,start=600,numSamples = 30,parametersOnly = F))
  pChain[[i]] <- mcmc(rbind(chain1,chain2,chain3))
  rm(calibration);gc()
  LPind <- dim(pChain[[i]])[2]-2 #index to identify the prior in the chains
  pMAPx <- pChain[[i]][which.max(pChain[[i]][,LPind]),] ###check which one is the maximum a posteriori parameter vector
  if(i ==1){
    pMAP <- pMAPx
  }else{
    if(pMAPx[LPind] > pMAP[LPind]) pMAP <- pMAPx  
  }
  print(i)
}

pChainSample <- data.table()
dim(pChain[[1]])
# tracePlot(sampler = pChain[[1]], whichParameters = 2,thin=1)
# mcmcList <- as.mcmc.list(pChain)
save(pChain, file="calOut/allChain.rdata")
pChainSample <- getSample(pChain,numSamples = 1000)
save(pChainSample,pMAP,file="calOut/pMAP.rdata")


# pdf(paste0("calOut/tracePlots_",min(calSets),"to",max(calSets),".pdf"))
pdf("calOut/tracePlots.pdf")
tracePlot(pChain[c(2:4,6:9)])
dev.off()

####remove LL,LP, LPrior from the chains
pChain2 <- list()
for(i in 1:indRun) pChain2[i] <- pChain[i][,1:npar]

### combine the independent calibrations
# Note that you first you need to identify and remove those chains that got stuck in a local maxima or clearly did not converged.
# to do that look at the marginal trace plots "outCal/tracePlots.pdf" starting from the loglikelihhod

# combined chains
pChainComb <- list()
xx <- mcmc.list(pChain[1:3])
pChainComb[[1]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain[4:6])
pChainComb[[2]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain[6:7])
pChainComb[[3]] <- combine.mcmc(xx)


gelman.diag(pChain2,multivariate = T,autoburnin = F)
# gelman.diag(pChain,multivariate = T,autoburnin = F)

###filter the chains
###find out the chains that got stuck in local maxima looking at the loglikelihood
MAPx <- rep(0,indRun)
for(i in 2:indRun) MAPx[i] <- max(pChain[[i]][,(npar+1)])

indX <- sort.int(MAPx, decreasing = T, index.return = TRUE)
filtChain <- indX$ix[(2:7)]

# combined chains
pChainCombFilt <- list()
set1 <- filtChain[c(1,3,6)]#seq(2,length(filtChain),by=2)
set2 <- filtChain[c(2,4,5)]#seq(3,length(filtChain),by=2)


dim(pChain[[2]])
set1
dim(pChain[[2]])
chain1 <- matrix(0,(516*3),110)
chain1[seq(from=1,to=516*3,by=3),] <- pChain[[2]][,1:110]
chain1[seq(from=2,to=516*3,by=3),] <- pChain[[6]][,1:110]
chain1[seq(from=3,to=516*3,by=3),] <- pChain[[8]][,1:110]

set2

chain2 <- matrix(0,(516*3),110)
chain2[seq(from=1,to=516*3,by=3),] <- pChain[[3]][,1:110]
chain2[seq(from=2,to=516*3,by=3),] <- pChain[[4]][,1:110]
chain2[seq(from=3,to=516*3,by=3),] <- pChain[[7]][,1:110]

pChainCombFilt <- mcmc.list(mcmc(chain1),mcmc(chain2))

# xx <- mcmc.list(pChain[filtChain[set1]])
# pChainCombFilt[[1]] <- as.mcmc(combine.mcmc(xx))
# xx <- mcmc.list(pChain[filtChain[set2]][,1:110])
# pChainCombFilt[[2]] <- as.mcmc(combine.mcmc(xx))

# ciao <- mcmc.list(pChainCombFilt[[1]],pChainCombFilt[[2]])
gelman.diag(pChain2[[filtChain]],multivariate = T,autoburnin = F)

pchainX <- list()
pchainX[[i]] <- pChain2[[filtChain[i]]]

# pdf(paste0("calOut/tracePlots_",min(calSets),"to",max(calSets),".pdf"))
pdf("calOut/tracePlotsFilt.pdf")
tracePlot(pChainCombFilt)
dev.off()

