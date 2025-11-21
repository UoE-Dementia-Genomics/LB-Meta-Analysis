#Setting up parallel processors
.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(dplyr)
library(ggplot2)
library(cowplot)
library(nlme)
library(Haplin)
library(multcomp)
library(car)
library(parallel)
library(Rcpp)
library(fastcluster)
library(devtools)
library(WGCNA)
library(splitstackshape)
library(sva)
library(data.table)
library(CETYGO)


######################
##### Full Cohort ####
######################

load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Pheno.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Betas.Rdata")



# Calculate top 10 surrogate variables for Braak aSyn model:
design <- model.matrix(~ 0 + LB_stage+Age+SexNum+PMI+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+as.factor(PlateNum)+BraakTangle_numeric, data = lbPheno)
null <- model.matrix(~ 0 + Age+SexNum+PMI+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+as.factor(PlateNum)+BraakTangle_numeric, data = lbPheno)

svSet<- sva(betasLB[complete.cases(betasLB),],design,null, n.sv = 10)

Braak_aSyn_stage<- lbPheno$LB_stage
Sex <- as.factor(lbPheno$SexNum)
Age_death<- lbPheno$Age
PMD_min<- lbPheno$PMI
NeuNNeg_Sox10Neg_IRF8Pos <- lbPheno$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- lbPheno$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- lbPheno$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- lbPheno$NeuNNeg_Sox10Neg_IRF8Neg
Plate <- as.factor(lbPheno$Plate)
braak_nft_stage<- lbPheno$BraakTangle_numeric
SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2]
SV3 <- svSet$sv[,3]
SV4 <- svSet$sv[,4]
SV5 <- svSet$sv[,5]


EWAS <- function(x, Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage){
  res<-try(lm(as.numeric(x) ~   Braak_aSyn_stage+Sex+Age_death+PMD_min+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+braak_nft_stage), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv1 <- function(x, Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,SV1){
  res<-try(lm(as.numeric(x) ~   Braak_aSyn_stage+Sex+Age_death+PMD_min+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+braak_nft_stage+SV1), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}


LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}


# cl<- makeCluster(12)
# BDR_res_LB <- t(parApply(cl,betasLB,1,EWAS,Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage))
# stopCluster(cl)
# LambdaInf(na.omit(BDR_res_LB[,4]))
# [1] 1.738183
# save(BDR_res_LB, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/PFC/BDR/fiveCell/BDR_HighNFT_FiveCell.Rdata")

cl<- makeCluster(12)
BDR_res_LB_SV1 <- t(parApply(cl,betasLB,1,EWAS_sv1,Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,SV1))
stopCluster(cl)
LambdaInf(na.omit(BDR_res_LB_SV1[,4]))
# [1] 1.060798
save(BDR_res_LB_SV1, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/BDR/BDR_Res_SV1.Rdata")

##########################
##### Thal Controlled ####
##########################
load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Pheno.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Betas.Rdata")


lbPheno <- lbPheno[-which(is.na(lbPheno$Thal_amyloid)),]
betasLB <- betasLB[,rownames(lbPheno)]

# Calculate top 10 surrogate variables for Braak aSyn model:
design <- model.matrix(~ 0 + LB_stage+Age+SexNum+PMI+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+as.factor(PlateNum)+BraakTangle_numeric+Thal_amyloid, data = lbPheno)
null <- model.matrix(~ 0 + Age+SexNum+PMI+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+as.factor(PlateNum)+BraakTangle_numeric+Thal_amyloid, data = lbPheno)

svSet<- sva(betasLB[complete.cases(betasLB),],design,null, n.sv = 10)

Braak_aSyn_stage<- lbPheno$LB_stage
Sex <- as.factor(lbPheno$SexNum)
Age_death<- lbPheno$Age
PMD_min<- lbPheno$PMI
NeuNNeg_Sox10Neg_IRF8Pos <- lbPheno$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- lbPheno$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- lbPheno$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- lbPheno$NeuNNeg_Sox10Neg_IRF8Neg
Plate <- as.factor(lbPheno$Plate)
braak_nft_stage<- lbPheno$BraakTangle_numeric
Thal <- lbPheno$Thal_amyloid
SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2]
SV3 <- svSet$sv[,3]
SV4 <- svSet$sv[,4]
SV5 <- svSet$sv[,5]


EWAS <- function(x, Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,Thal){
  res<-try(lm(as.numeric(x) ~   Braak_aSyn_stage+Sex+Age_death+PMD_min+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+braak_nft_stage + Thal), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv1 <- function(x, Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,Thal,SV1){
  res<-try(lm(as.numeric(x) ~   Braak_aSyn_stage+Sex+Age_death+PMD_min+NeuNNeg_Sox10Neg_IRF8Pos+NeuNPos_SOX6Neg+NeuNPos_SOX6Pos+NeuNNeg_Sox10Neg_IRF8Neg+Plate+braak_nft_stage+Thal+SV1), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}


LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

# 
# cl<- makeCluster(12)
# BDR_res_LB_Thal <- t(parApply(cl,betasLB,1,EWAS,Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,Thal))
# stopCluster(cl)
# BDR_res_LB_Thal(na.omit(BDR_res_LB_Thal[,4]))
# # [1] 1.440808
# 
# save(BDR_res_LB_Thal, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/PFC/BDR/fiveCell/BDR_HighNFT_FiveCell_Thal.Rdata")

cl<- makeCluster(12)
BDR_res_LB_Thal_SV1 <- t(parApply(cl,betasLB,1,EWAS_sv1,Braak_aSyn_stage,Sex,Age_death,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos,NeuNPos_SOX6Neg,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,Plate,braak_nft_stage,Thal,SV1))
stopCluster(cl)
LambdaInf(na.omit(BDR_res_LB_Thal_SV1[,4]))
# [1] 1.034281

save(BDR_res_LB_Thal_SV1, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/BDR/BDR_Res_Thal_SV1.Rdata")
