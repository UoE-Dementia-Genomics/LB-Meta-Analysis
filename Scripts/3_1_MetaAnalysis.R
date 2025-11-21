
setwd("/lustre/projects/Research_Project-T112069/")
# Script to run meta analyses across cohorts


.libPaths("/lustre/projects/Research_Project-T112069/packages")
#Setting up parallel processors
library(doParallel)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(limma)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(Haplin)
library(sva)
library(reshape)
library(meta)

######################
#### Full Cohort #####
######################


load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex.Rdata")
res_UKBBN <- res_UKBBN[,-3]
colnames(res_UKBBN) <- c("Values","SE","T","P")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/BDR/BDR_Res_SV1.Rdata")

resUK <- res_UKBBN
resNBB <- res_NBB_SV3
resBDR <- BDR_res_LB_SV1

overlap <-Reduce(intersect, list(rownames(resUK),rownames(resNBB),rownames(resBDR)))
resUK <- resUK[overlap,]
resNBB <- resNBB[overlap,]

resBDR <- resBDR[overlap,]

colnames(resUK) <- c("Beta_UK","SE_UK","T_UK","P_UK")
colnames(resNBB) <- c("Beta_NBB","SE_NBB","T_NBB","P_NBB")
colnames(resBDR) <- c("Beta_BDR","SE_BDR","T_BDR","P_BDR")

resAll <- cbind(resUK, resNBB,resBDR)
resAll <- as.data.frame(resAll)



metaPar <- function(x){
  .libPaths("/lustre/projects/Research_Project-T112069/packages")
  library(meta)
  out<-try(meta::metagen(c(x[1], x[5],x[9]), c(x[2], x[6],x[10])))
  print(c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q)))
} 

cl<- makeCluster(14)
allSumstat <- t(parApply(cl,resAll,1,metaPar))
stopCluster(cl)


allSumstat <- as.data.frame(allSumstat)
allSumstat$FDR <- p.adjust(allSumstat[,3],method = "fdr")
allSumstat <- allSumstat[order(allSumstat[,3]),]

colnames(allSumstat)[1:8] <- c("Fixed_Effect","Fixed_SE","Fixed_P","Random_Effect","Random_SE","Random_P","I2","Het_P")
save(allSumstat, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")


######################
#### Pure LB #########
######################
#Note, clear environment before running
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex_Filtered.Rdata")
res_UKBBN_Filtered <- res_UKBBN_Filtered[,-3]
colnames(res_UKBBN_Filtered) <- c("Values","SE","T","P")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3_filtered.Rdata")
resUK <- res_UKBBN_Filtered
resNBB <- NBB_res_SV3_filtered


overlap <-Reduce(intersect, list(rownames(resUK),rownames(resNBB)))
resUK <- resUK[overlap,]
resNBB <- resNBB[overlap,]


colnames(resUK) <- c("Beta_UK","SE_UK","T_UK","P_UK")
colnames(resNBB) <- c("Beta_NBB","SE_NBB","T_NBB","P_NBB")

resAll <- cbind(resUK, resNBB)
resAll <- as.data.frame(resAll)

metaPar <- function(x){
  .libPaths("/lustre/projects/Research_Project-T112069/packages")
  library(meta)
  out<-try(meta::metagen(c(x[1], x[5]), c(x[2], x[6])))
  print(c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q)))
} 

cl<- makeCluster(14)
lowSumstat <- t(parApply(cl,resAll,1,metaPar))
stopCluster(cl)


lowSumstat <- as.data.frame(lowSumstat)
lowSumstat$FDR <- p.adjust(lowSumstat[,3],method = "fdr")
lowSumstat <- lowSumstat[order(lowSumstat[,3]),]

colnames(lowSumstat)[1:8] <- c("Fixed_Effect","Fixed_SE","Fixed_P","Random_Effect","Random_SE","Random_P","I2","Het_P")
save(lowSumstat, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_PureLB.Rdata")



##############################
#### Thal Controlled #########
##############################
#! NOTE: Clear environment before running

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex_Thal_SV3.Rdata")
UKBBN_res_Thal_SV3 <- UKBBN_res_Thal_SV3[,-3]
colnames(UKBBN_res_Thal_SV3) <- c("Values","SE","T","P")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3_Thal.Rdata")
resUK <- UKBBN_res_Thal_SV3
resNBB <- NBB_res_SV3_Thal
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/BDR/BDR_Res_Thal_SV1.Rdata")
resBDR <- BDR_res_LB_Thal_SV1 

overlap <-Reduce(intersect, list(rownames(resUK),rownames(resNBB),rownames(resBDR)))
resUK <- resUK[overlap,]
resNBB <- resNBB[overlap,]

resBDR <- resBDR[overlap,]

colnames(resUK) <- c("Beta_UK","SE_UK","T_UK","P_UK")
colnames(resNBB) <- c("Beta_NBB","SE_NBB","T_NBB","P_NBB")
colnames(resBDR) <- c("Beta_BDR","SE_BDR","T_BDR","P_BDR")

resAll <- cbind(resUK, resNBB,resBDR)
resAll <- as.data.frame(resAll)


metaPar <- function(x){
  .libPaths("/lustre/projects/Research_Project-T112069/packages")
  library(meta)
  out<-try(meta::metagen(c(x[1], x[5],x[9]), c(x[2], x[6],x[10])))
  print(c(out$TE.fixed, out$seTE.fixed, out$pval.fixed, out$TE.random, out$seTE.random, out$pval.random, out$I2, 1-pchisq(out$Q, out$df.Q)))
} 

cl<- makeCluster(14)
thalSumstat <- t(parApply(cl,resAll,1,metaPar))
stopCluster(cl)

thalSumstat <- as.data.frame(thalSumstat)
thalSumstat$FDR <- p.adjust(thalSumstat[,3],method = "fdr")
thalSumstat <- thalSumstat[order(thalSumstat[,3]),]
colnames(thalSumstat)[1:8] <- c("Fixed_Effect","Fixed_SE","Fixed_P","Random_Effect","Random_SE","Random_P","I2","Het_P")

save(thalSumstat, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_ThalControlled.Rdata")

