##### Load packages and set working directory #####
setwd("/lustre/projects/Research_Project-T112069/Meth/NBB/")
.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(minfi)
library(parallel)
library(wateRmelon)
library(corrplot)
library(sva)
library(GEOquery)
library(ggplot2)
library(dplyr)
library(data.table)
library(FlowSorted.DLPFC.450k)
library(CETYGO)

############################
##### Full Cohort EWAS #####
############################

##### Load data and filter samples #####
# Load in data
NBB<-getGEO(filename='GSE203332_series_matrix.txt',
            GSEMatrix=TRUE,getGPL = FALSE)
NBBpheno <- pData(object = NBB)

# Process the pheno filename
NBBpheno <- NBBpheno[,c(49:65)]
colnames(NBBpheno) <-  gsub(pattern = ":ch1", replacement = "",colnames(NBBpheno)) # sort out the colnames
NBBpheno <- NBBpheno[-which(NBBpheno$neuropathological_diagnosis %in% c("AD","mixed_AD_LBD")),]
rownames(NBBpheno) <- paste(NBBpheno$`slide`,NBBpheno$`slide_position`,sep = "_")
NBBgeo <- NBBpheno
NBBpheno <- read.table("NBB_pheno.txt", header = T)
rownames(NBBpheno) <- NBBpheno$PlatePos_ID
NBBpheno <- cbind(NBBpheno, NBBgeo[rownames(NBBpheno),c("braak_nft_stage","cerad_stage")])
NBBpheno$braak_nft_stage <- as.numeric(NBBpheno$braak_nft_stage)

NBBpheno$PlateNum <- as.numeric(as.factor(NBBpheno$Plate))
NBBpheno$SexNum <- as.numeric(as.factor(NBBpheno$Sex))


#Filter beta matrix
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
filter <- read.table("/lustre/projects/Research_Project-T112069/Meth/reference/SNPProbes_McCartney.txt")
row.names(epicMani) <- epicMani[,1]
NBBbetas <- fread("GSE203332_Matrix_processed_dasen_betas.txt")
NBBbetas <- as.data.frame(NBBbetas)
rownames(NBBbetas) <- NBBbetas$ID_REF
NBBbetas$ID_REF <- NULL
NBBbetas <- NBBbetas[,NBBpheno$PlatePos_ID]
xyProbes <- rownames(epicMani[which(epicMani$CHR %in% c("X","Y")),])
NBBbetas <- as.data.frame(NBBbetas)
NBBbetas <- NBBbetas[-which(rownames(NBBbetas) %in% xyProbes),]

##### Run new celltype deconvolution #####
cellSet<-projectCellTypeWithError(NBBbetas, modelBrainCoef[["IDOL"]][[8]])
cellSet <- as.data.frame(cellSet)

##### Compute SVs #####

betas <- NBBbetas
betas <- as.matrix(betas)
pheno <- cbind(NBBpheno[,c("Braak_aSyn_stage" , "Age_death" , "SexNum" , "PMD_min" , "PlateNum" , "braak_nft_stage")],
               cellSet[,c("NeuNNeg_Sox10Neg_IRF8Pos", "NeuNPos_SOX6Neg", "NeuNPos_SOX6Pos","NeuNNeg_Sox10Neg_IRF8Neg")])


# Calculate top 10 surrogate variables for Braak aSyn model:
design <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                         NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage, data = pheno)
null_model <- model.matrix(~ 0 + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                             NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage, data = pheno)
svSet <- sva(betas[complete.cases(betas),],design,null_model, n.sv = 10)

##### Define variables for EWAS #####
Braak_aSyn_stage <- pheno$Braak_aSyn_stage
Age_death <- pheno$Age_death
SexNum <- pheno$SexNum
PMD_min <- pheno$PMD_min
PlateNum <- as.factor(pheno$PlateNum)
braak_nft_stage <- pheno$braak_nft_stage
NeuNNeg_Sox10Neg_IRF8Pos <- pheno$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- pheno$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- pheno$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- pheno$NeuNNeg_Sox10Neg_IRF8Neg
SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2] 
SV3 <- svSet$sv[,3] 
SV4 <- svSet$sv[,4] 
SV5 <- svSet$sv[,5] 


##### Define EWAS functions #####



EWAS <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv1 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv2 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1+SV2), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv3 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2,SV3){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1 + SV2 + SV3), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}


LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda)
}

##### Run EWAS in parallel & save results #####

#cl<- makeCluster(32)
#res_NBB <- t(parApply(cl,betas,1,EWAS, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage))
#stopCluster(cl)
#save(res_NBB, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/PFC/NBB/FiveCell/NBB_knownCovar.Rdata")
# LambdaInf(na.omit(res_NBB[,4]))
# [1] 1.724735


#cl<- makeCluster(32)
#res_NBB_SV1 <- t(parApply(cl,betas,1,EWAS_sv1, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1))
#stopCluster(cl)
#save(res_NBB_SV1, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/PFC/NBB/FiveCell/NBB_knownCovar_SV1.Rdata")
# LambdaInf(na.omit(res_NBB_SV1[,4]))
# 1.659109

#cl<- makeCluster(32)
#res_NBB_SV2 <- t(parApply(cl,betas,1,EWAS_sv2, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2))
#stopCluster(cl)
#save(res_NBB_SV2, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Pathology/PFC/NBB/FiveCell/NBB_knownCovar_SV2.Rdata")
# LambdaInf(na.omit(res_NBB_SV2[,4]))
# [1] 2.104532

cl<- makeCluster(32)
res_NBB_SV3 <- t(parApply(cl,betas,1,EWAS_sv3, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2,SV3))
stopCluster(cl)
save(res_NBB_SV3, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3.Rdata")
# LambdaInf(na.omit(res_NBB_SV3[,4]))
# [1] 1.088889


##########################
##### "Pure LB" EWAS #####
##########################

# NOTE! Clear environment before running this section 
##### Load data and filter samples #####
# Load in data
NBB<-getGEO(filename='GSE203332_series_matrix.txt',
            GSEMatrix=TRUE,getGPL = FALSE)
NBBpheno <- pData(object = NBB)

# Process the pheno filename
NBBpheno <- NBBpheno[,c(49:65)]
colnames(NBBpheno) <-  gsub(pattern = ":ch1", replacement = "",colnames(NBBpheno)) # sort out the colnames
NBBpheno <- NBBpheno[-which(NBBpheno$neuropathological_diagnosis %in% c("AD","mixed_AD_LBD")),]
rownames(NBBpheno) <- paste(NBBpheno$`slide`,NBBpheno$`slide_position`,sep = "_")
NBBgeo <- NBBpheno
NBBpheno <- read.table("NBB_pheno.txt", header = T)
rownames(NBBpheno) <- NBBpheno$PlatePos_ID
NBBpheno <- cbind(NBBpheno, NBBgeo[rownames(NBBpheno),c("braak_nft_stage","cerad_stage")])
NBBpheno$braak_nft_stage <- as.numeric(NBBpheno$braak_nft_stage)
NBBpheno <- NBBpheno[-which(NBBpheno$braak_nft_stage > 2),]


NBBpheno$PlateNum <- as.numeric(as.factor(NBBpheno$Plate))
NBBpheno$SexNum <- as.numeric(as.factor(NBBpheno$Sex))


#Filter beta matrix
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
filter <- read.table("/lustre/projects/Research_Project-T112069/Meth/reference/SNPProbes_McCartney.txt")
row.names(epicMani) <- epicMani[,1]
NBBbetas <- fread("GSE203332_Matrix_processed_dasen_betas.txt")
NBBbetas <- as.data.frame(NBBbetas)
rownames(NBBbetas) <- NBBbetas$ID_REF
NBBbetas$ID_REF <- NULL
NBBbetas <- NBBbetas[,NBBpheno$PlatePos_ID]
xyProbes <- rownames(epicMani[which(epicMani$CHR %in% c("X","Y")),])
NBBbetas <- as.data.frame(NBBbetas)
NBBbetas <- NBBbetas[-which(rownames(NBBbetas) %in% xyProbes),]

##### Run new celltype deconvolution #####
cellSet<-projectCellTypeWithError(NBBbetas, modelBrainCoef[["IDOL"]][[8]])
cellSet <- as.data.frame(cellSet)

##### Compute SVs #####

betas <- NBBbetas
betas <- as.matrix(betas)
pheno <- cbind(NBBpheno[,c("Braak_aSyn_stage" , "Age_death" , "SexNum" , "PMD_min" , "PlateNum" , "braak_nft_stage")],
               cellSet[,c("NeuNNeg_Sox10Neg_IRF8Pos", "NeuNPos_SOX6Neg", "NeuNPos_SOX6Pos","NeuNNeg_Sox10Neg_IRF8Neg")])


# Calculate top 10 surrogate variables for Braak aSyn model:
design <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                         NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage, data = pheno)
null_model <- model.matrix(~ 0 + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                             NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage, data = pheno)
svSet <- sva(betas[complete.cases(betas),],design,null_model, n.sv = 10)

##### Define variables for EWAS #####
Braak_aSyn_stage <- pheno$Braak_aSyn_stage
Age_death <- pheno$Age_death
SexNum <- pheno$SexNum
PMD_min <- pheno$PMD_min
PlateNum <- as.factor(pheno$PlateNum)
braak_nft_stage <- pheno$braak_nft_stage
NeuNNeg_Sox10Neg_IRF8Pos <- pheno$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- pheno$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- pheno$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- pheno$NeuNNeg_Sox10Neg_IRF8Neg
SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2] 
SV3 <- svSet$sv[,3] 
SV4 <- svSet$sv[,4] 
SV5 <- svSet$sv[,5] 


##### Define EWAS functions #####



EWAS <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv1 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv2 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1+SV2), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv3 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2,SV3){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + SV1 + SV2 + SV3), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}




##### Run EWAS in parallel & save results #####



cl<- makeCluster(32)
NBB_res_SV3_filtered <- t(parApply(cl,betas,1,EWAS_sv3, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,SV1,SV2,SV3))
stopCluster(cl)
save(NBB_res_SV3_filtered, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3_filtered.Rdata")
# LambdaInf(na.omit(NBBresSV3[,4]))
# [1] 1.041707



################################
##### Thal controlled EWAS #####
################################
#! NOTE, clear environment before running this step 

##### Load data and filter samples #####
# Load in data
NBB<-getGEO(filename='GSE203332_series_matrix.txt',
            GSEMatrix=TRUE,getGPL = FALSE)
NBBpheno <- pData(object = NBB)

# Process the pheno filename
NBBpheno <- NBBpheno[,c(49:65)]
colnames(NBBpheno) <-  gsub(pattern = ":ch1", replacement = "",colnames(NBBpheno)) # sort out the colnames
NBBpheno <- NBBpheno[-which(NBBpheno$neuropathological_diagnosis %in% c("AD","mixed_AD_LBD")),]
rownames(NBBpheno) <- paste(NBBpheno$`slide`,NBBpheno$`slide_position`,sep = "_")
NBBgeo <- NBBpheno
NBBpheno <- read.table("NBB_pheno.txt", header = T)
rownames(NBBpheno) <- NBBpheno$PlatePos_ID
NBBpheno <- cbind(NBBpheno, NBBgeo[rownames(NBBpheno),c("braak_nft_stage","cerad_stage","thal_amyloid_phase")])
NBBpheno$braak_nft_stage <- as.numeric(NBBpheno$braak_nft_stage)
NBBpheno$thal_amyloid_phase <- as.numeric(NBBpheno$thal_amyloid_phase)

NBBpheno$PlateNum <- as.numeric(as.factor(NBBpheno$Plate))
NBBpheno$SexNum <- as.numeric(as.factor(NBBpheno$Sex))


#Filter beta matrix
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
filter <- read.table("/lustre/projects/Research_Project-T112069/Meth/reference/SNPProbes_McCartney.txt")
row.names(epicMani) <- epicMani[,1]
NBBbetas <- fread("GSE203332_Matrix_processed_dasen_betas.txt")
NBBbetas <- as.data.frame(NBBbetas)
rownames(NBBbetas) <- NBBbetas$ID_REF
NBBbetas$ID_REF <- NULL
NBBbetas <- NBBbetas[,NBBpheno$PlatePos_ID]
xyProbes <- rownames(epicMani[which(epicMani$CHR %in% c("X","Y")),])
NBBbetas <- as.data.frame(NBBbetas)
NBBbetas <- NBBbetas[-which(rownames(NBBbetas) %in% xyProbes),]

##### Run new celltype deconvolution #####
cellSet<-projectCellTypeWithError(NBBbetas, modelBrainCoef[["IDOL"]][[8]])
cellSet <- as.data.frame(cellSet)

##### Subset to those with Thal
NBBpheno <- NBBpheno[-which(is.na(NBBpheno$thal_amyloid_phase)),]
# 262  24

##### Compute SVs #####



betas <- NBBbetas
betas <- as.matrix(betas)
betas <- betas[,rownames(NBBpheno)]
pheno <- cbind(NBBpheno[,c("Braak_aSyn_stage" , "Age_death" , "SexNum" , "PMD_min" , "PlateNum" , "braak_nft_stage","thal_amyloid_phase")],
               cellSet[rownames(NBBpheno),c("NeuNNeg_Sox10Neg_IRF8Pos", "NeuNPos_SOX6Neg", "NeuNPos_SOX6Pos","NeuNNeg_Sox10Neg_IRF8Neg")])


# Calculate top 10 surrogate variables for Braak aSyn model:
design <- model.matrix(~ 0 + Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                         NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + thal_amyloid_phase, data = pheno)
null_model <- model.matrix(~ 0 + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos +  
                             NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage+ thal_amyloid_phase, data = pheno)
svSet <- sva(betas[complete.cases(betas),],design,null_model, n.sv = 10)

##### Define variables for EWAS #####
Braak_aSyn_stage <- pheno$Braak_aSyn_stage
Age_death <- pheno$Age_death
SexNum <- pheno$SexNum
PMD_min <- pheno$PMD_min
PlateNum <- as.factor(pheno$PlateNum)
braak_nft_stage <- pheno$braak_nft_stage
NeuNNeg_Sox10Neg_IRF8Pos <- pheno$NeuNNeg_Sox10Neg_IRF8Pos
NeuNPos_SOX6Neg <- pheno$NeuNPos_SOX6Neg
NeuNPos_SOX6Pos <- pheno$NeuNPos_SOX6Pos
NeuNNeg_Sox10Neg_IRF8Neg <- pheno$NeuNNeg_Sox10Neg_IRF8Neg
Thal <- pheno$thal_amyloid_phase
SV1 <- svSet$sv[,1]
SV2 <- svSet$sv[,2] 
SV3 <- svSet$sv[,3] 
SV4 <- svSet$sv[,4] 
SV5 <- svSet$sv[,5] 


##### Define EWAS functions #####



EWAS <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,Thal){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + Thal), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv1 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,Thal,SV1){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage +Thal + SV1), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv2 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,Thal,SV1,SV2){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage +Thal+ SV1+SV2), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}

EWAS_sv3 <- function(x, Braak_aSyn_stage, Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,Thal,SV1,SV2,SV3){
  res<-try(lm(as.numeric(x) ~  Braak_aSyn_stage + Age_death + SexNum + PMD_min + NeuNNeg_Sox10Neg_IRF8Pos + NeuNPos_SOX6Neg + NeuNPos_SOX6Pos + NeuNNeg_Sox10Neg_IRF8Neg + PlateNum + braak_nft_stage + Thal + SV1 + SV2 + SV3), silent=T)
  if(inherits(res,'try-error')) return(rep(NA,4))
  return(coef(summary(res))[2,])
}



##### Run EWAS in parallel & save results #####
cl<- makeCluster(32)
NBB_res_SV3_Thal <- t(parApply(cl,betas,1,EWAS_sv3, Braak_aSyn_stage,Age_death,SexNum,PMD_min,NeuNNeg_Sox10Neg_IRF8Pos, NeuNPos_SOX6Neg ,NeuNPos_SOX6Pos,NeuNNeg_Sox10Neg_IRF8Neg,PlateNum,braak_nft_stage,Thal,SV1,SV2,SV3))
stopCluster(cl)
# > LambdaInf(na.omit(NBB_res_SV3_Thal[,4]))
# [1] 1.023153

save(NBB_res_SV3_Thal, file = "/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3_Thal.Rdata")



