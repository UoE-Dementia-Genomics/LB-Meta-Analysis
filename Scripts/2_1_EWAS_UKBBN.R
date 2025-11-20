setwd("/lustre/projects/Research_Project-T112069/")
.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(doParallel)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(limma)
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(Haplin)
library(devtools)
library(qqman)
library(fastDummies)
library(CETYGO)
library(sva)

#########################
### Full Cohort EWAS ####
#########################

##### Load data #####
# Load in data
load("Meth/Methylation/QC/FinalData/betasPFC.Rdata")
load("Meth/Methylation/QC/ewasPheno.Rdata")
load("Meth/Methylation/QC/FinalData/betasCNG.Rdata")

# Merge data and index
betas <- cbind(betas.cng,betas.pfc)
betas <- betas[,rownames(phenoP)]
phenoP$count <- 1
betasP <-  betas[,phenoP$Basename]

# Defined EWAS Variables
BraakLB <- as.numeric(phenoP$Braak.LB.Numeric)
Sex <- as.numeric(as.factor(phenoP$Gender))
Age <- as.numeric(phenoP$Age)
PMI <- as.numeric(phenoP$Post.Mortem.Delay)
NeuNNeg_Sox10Neg_IRF8Pos <- phenoP$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- phenoP$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- phenoP$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- phenoP$NeuNNeg_Sox10Neg_IRF8Neg
Plate <- as.factor(phenoP$plateNumber)
BB <- as.factor(phenoP$Brain.Bank)
BR <- as.factor(phenoP$Brain_Region)
ID <- as.factor(phenoP$ID)
BraakNFT <- as.numeric(phenoP$Braak.Tangle)

EWAS <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR , random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 


LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

##### Run EWAS in parallel & save results #####
cl<- makeCluster(16)
res_UKBBN <-t(parApply(cl, betasP, 1, EWAS,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR))
stopCluster(cl)

colnames(res_UKBBN) <- c("Effect","SE","DF","T","P")
res_UKBBN <- as.data.frame(res_UKBBN)

LambdaInf(na.omit(res_UKBBN$P))
# [1] 1.094549

save(res_UKBBN, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex.Rdata")


#######################
### "Pure LB" EWAS ####
#######################
##### Load data #####
# Load in data
load("Meth/Methylation/QC/FinalData/betasPFC.Rdata")
load("Meth/Methylation/QC/ewasPheno.Rdata")
load("Meth/Methylation/QC/FinalData/betasCNG.Rdata")

# Merge data and index
rownames(phenoP) <- phenoP$Basename
betas <- cbind(betas.cng,betas.pfc)
betas <- betas[,rownames(phenoP)]
phenoP$count <- 1
betasP <-  betas[,phenoP$Basename]

phenoP <- phenoP[-which(phenoP$Braak.Tangle > 2),]
phenoP$count <- 1
betasP <-  betas[,phenoP$Basename]

# Defined EWAS Variables
BraakLB <- as.numeric(phenoP$Braak.LB.Numeric)
Sex <- as.numeric(as.factor(phenoP$Gender))
Age <- as.numeric(phenoP$Age)
PMI <- as.numeric(phenoP$Post.Mortem.Delay)
NeuNNeg_Sox10Neg_IRF8Pos <- phenoP$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- phenoP$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- phenoP$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- phenoP$NeuNNeg_Sox10Neg_IRF8Neg
Plate <- as.factor(phenoP$plateNumber)
BB <- as.factor(phenoP$Brain.Bank)
BR <- as.factor(phenoP$Brain_Region)
ID <- as.factor(phenoP$ID)
BraakNFT <- as.numeric(phenoP$Braak.Tangle)

##### Define EWAS functions #####

EWAS <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR , random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

##### Run EWAS in parallel & save results #####
cl<- makeCluster(16)
res_UKBBN_Filtered <-t(parApply(cl, betasP, 1, EWAS,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR))
stopCluster(cl)
save(res_UKBBN_Filtered, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex_Filtered.Rdata")

colnames(res_UKBBN_Filtered) <- c("Effect","SE","DF","T","P")
res_UKBBN_Filtered <- as.data.frame(res_UKBBN_Filtered)

# LambdaInf(na.omit(res_UKBBN_Filtered$P))
# [1] 1.038432


#############################
### Thal controlled EWAS ####
#############################

##### Load data #####
# Load in data
load("Meth/Methylation/QC/FinalData/betasPFC.Rdata")
load("Meth/Methylation/QC/ewasPheno.Rdata")
load("Meth/Methylation/QC/FinalData/betasCNG.Rdata")

# Merge data and index
betas <- cbind(betas.cng,betas.pfc)
betas <- betas[,rownames(phenoP)]
phenoP$count <- 1
betasP <-  betas[,phenoP$Basename]

phenoP <- phenoP[-which(is.na(phenoP$Thal)),]
betasP <- betasP[,phenoP$Basename]

##### Define variables for EWAS #####
sampleSet <- rownames(phenoP)

BraakLB <- as.numeric(phenoP$Braak.LB.Numeric)
Sex <- as.numeric(as.factor(phenoP$Gender))
Age <- as.numeric(phenoP$Age)
PMI <- as.numeric(phenoP$Post.Mortem.Delay)
NeuNNeg_Sox10Neg_IRF8Pos <- phenoP$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <-phenoP$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <-phenoP$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- phenoP$NeuNNeg_Sox10Neg_IRF8Neg
Plate <- as.factor(phenoP$plateNumber)
BB <- as.factor(phenoP$Brain.Bank)
BR <- as.factor(phenoP$Brain_Region)
ID <- as.factor(phenoP$ID)
BraakNFT <- as.numeric(phenoP$Braak.Tangle)
Thal <- phenoP$Thal
svFrame <- data.frame(BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal)
svModel <- model.matrix(~ 0 + BraakLB + Sex + Age + PMI + BraakNFT + NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR + Thal, data = svFrame)
nullModel <- model.matrix(~ 0 +  Sex + Age + PMI + BraakNFT + NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR + Thal, data = svFrame)
svSet <- sva(betasP[complete.cases(betasP),],svModel,nullModel, n.sv = 10)  

SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2]
SV3 <- svSet$sv[,3]

##### Define EWAS functions #####

EWAS <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR+Thal , random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

EWASsv1 <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR+Thal+SV1 , random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

EWASsv2 <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1,SV2){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR+Thal +SV1+SV2, random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

EWASsv3 <- function(x, BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1,SV2,SV3){
  library(nlme)
  res<-try(nlme::lme(as.numeric(x) ~ BraakLB +Sex+Age+PMI+BraakNFT+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+BB+BR+Thal +SV1+SV2+SV3, random=~1|ID), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,5))
  return(coef(summary(res))[2,])
} 

LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}



##### Run EWAS in parallel & save results #####
#cl<- makeCluster(16)
#UKBBN_res_Thal <-t(parApply(cl, betasP, 1, EWAS,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal))
#stopCluster(cl)
# LambdaInf(na.omit(UKBBN_res_Thal[,5]))
# [1] 1.276167

#save(resFiveCell_Thal, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/Multi/FiveCell/highPath_CrossCortex_FiveCell_ThalSub.Rdata")


#cl<- makeCluster(16)
#UKBBN_res_Thal_SV1 <-t(parApply(cl, betasP, 1, EWASsv1,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1))
#stopCluster(cl)
# LambdaInf(na.omit(UKBBN_res_Thal_SV1[,5]))
# [1] 1.186803
#save(resFiveCell_Thal_SV1, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/Multi/FiveCell/highPath_CrossCortex_FiveCell_ThalSub_SV1.Rdata")


# cl<- makeCluster(16)
# UKBBN_res_Thal_SV2 <-t(parApply(cl, betasP, 1, EWASsv2,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1,SV2))
# stopCluster(cl)
# save(resFiveCell_Thal_SV2, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/Multi/FiveCell/highPath_CrossCortex_FiveCell_ThalSub_SV2.Rdata")

# LambdaInf(na.omit(UKBBN_res_Thal_SV2[,5]))
# [1] 1.187224

# cl<- makeCluster(16)
# UKBBN_res_Thal_SV3 <-t(parApply(cl, betasP, 1, EWASsv3,BraakLB,ID, Sex, Age, PMI,BraakNFT,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,BB,BR,Thal,SV1,SV2,SV3))
# stopCluster(cl)
# 
# LambdaInf(na.omit(UKBBN_res_Thal_SV3[,5]))
# [1] 1.105494

save(UKBBN_res_Thal_SV3, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex_Thal_SV3.Rdata")

