.libPaths(c("/scratch/project_2000994/project_rpackages", .libPaths()))
library("BayesianTools")
library("runjags")
library("coda")

setwd("/scratch/project_2000994/calibrations/3PGQcalibration")

load("calOut/calibration_0.1.rdata")
namesX <- colnames(calibration[[1]]$chain[[1]])

calSets <- 5:6
###settings
npar <- calibration[[1]]$setup$numPars
indRun <-9 #number of independent calibration runs
thin=1
pChain <- mcmc.list()
pMAP <- NULL

for(i in 1:indRun){
  for(ij in calSets){
    load(paste0("calOut/calibration_",ij,".",i,".rdata"))
    lChain <- dim(calibration[[1]]$chain[[1]])[1]
    seqX <- seq(thin,lChain,by=thin)
    if(ij==calSets[1]){
      chain1 <- rbind(calibration[[1]]$chain[[1]][seqX,],
                      calibration[[1]]$chain[[2]][seqX,],
                      calibration[[1]]$chain[[3]][seqX,])
      colnames(chain1) <- namesX
      chain2 <- rbind(calibration[[2]]$chain[[1]][seqX,],
                      calibration[[2]]$chain[[2]][seqX,],
                      calibration[[2]]$chain[[3]][seqX,])
      colnames(chain2) <- namesX
      chain3 <- rbind(calibration[[3]]$chain[[1]][seqX,],
                      calibration[[3]]$chain[[2]][seqX,],
                      calibration[[3]]$chain[[3]][seqX,])
      colnames(chain3) <- namesX
    }else{
      chain1 <- rbind(chain1,
                      calibration[[1]]$chain[[1]][seqX,],
                      calibration[[1]]$chain[[2]][seqX,],
                      calibration[[1]]$chain[[3]][seqX,])
      chain2 <- rbind(chain2,
                      calibration[[2]]$chain[[1]][seqX,],
                      calibration[[2]]$chain[[2]][seqX,],
                      calibration[[2]]$chain[[3]][seqX,])
      chain3 <- rbind(chain3,
                      calibration[[3]]$chain[[1]][seqX,],
                      calibration[[3]]$chain[[2]][seqX,],
                      calibration[[3]]$chain[[3]][seqX,])
    }
  }
  #create a list of 3 chains. Note those 3 chains are the 3 parallel chains that are used in the DE (differential evolution algorithms)
  mcmcList <- mcmc.list(mcmc(chain1),mcmc(chain2),mcmc(chain3))
  
  #combine the 3 parallel chain in one unique chain
  pChain[[i]] <- combine.mcmc(mcmc.objects=mcmcList)
  
  LPind <- dim(pChain[[i]])[2]-2 #index to identify the prior in the chains
  pMAPx <- pChain[[i]][which.max(pChain[[i]][,LPind]),] ###check which one is the maximum a posteriori parameter vector
  if(i ==1){
    pMAP <- pMAPx
  }else{
    if(pMAPx[LPind] > pMAP[LPind]) pMAP <- pMAPx  
  }
}

save(pChain, file="calOut/allChain.rdata")
pChainSample <- getSample(pChain,numSamples = 1000)
save(pChainSample,pMAP,file="calOut/pMAP.rdata")


pdf("calOut/tracePlots_0to3.pdf")
tracePlot(pChain)
dev.off()

####remove LL,LP, LPrior from the chains
pChain2 <- list()
for(i in 1:indRun) pChain2[i] <- pChain[i][,1:npar]

### combine the independent calibrations
# Note that you first you need to identify and remove those chains that got stuck in a local maxima or clearly did not converged.
# to do that look at the marginal trace plots "outCal/tracePlots.pdf" starting from the loglikelihhod

# combined chains
pChainComb <- list()
xx <- mcmc.list(pChain2[1:3])
pChainComb[[1]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[4:6])
pChainComb[[2]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[7:9])
pChainComb[[3]] <- combine.mcmc(xx)


gelman.diag(pChainComb,multivariate = T,autoburnin = F)

###filter the chains
###find out the chains that got stuck in local maxima looking at the loglikelihood
MAPx <- rep(0,indRun)
for(i in 1:indRun) MAPx[i] <- max(pChain[[i]][,(npar+1)])

indX <- sort.int(MAPx, decreasing = T, index.return = TRUE)
filtChain <- indX$ix[(1:6)]

# combined chains
pChainCombFilt <- list()
set1 <- seq(1,length(filtChain),by=2)
set2 <- seq(2,length(filtChain),by=2)

xx <- mcmc.list(pChain2[filtChain[set1]])
pChainCombFilt[[1]] <- combine.mcmc(xx)
xx <- mcmc.list(pChain2[filtChain[set2]])
pChainCombFilt[[2]] <- combine.mcmc(xx)

gelman.diag(pChainCombFilt,multivariate = T,autoburnin = F)
