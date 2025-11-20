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

save(res_UKBBN, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex")
