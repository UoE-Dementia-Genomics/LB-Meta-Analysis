#Setting up parallel processors
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
library(devtools)
library(qqman)
library(fastDummies)
library(CETYGO)
library(corrplot)

setwd("/lustre/projects/Research_Project-T112069/")
# Load in data
load("Meth/Methylation/QC/FinalData/betasPFC.Rdata")
load("Meth/Methylation/QC/Pheno.Rdata")
load("Meth/Methylation/QC/FinalData/betasCNG.Rdata")

betas <- cbind(betas.cng,betas.pfc)


#Remove phenotypic outliers (control sampels with evidence of brain atrophy)
rownames(pheno) <- pheno$Basename
pheno <- pheno[colnames(betas),]
neuropathOutliers <- c("C064","C077","PDC008","C016","C030","C032","C049","C072","PDC027","PDC030","PDC105")
pheno <- pheno[-which(pheno$ID %in% neuropathOutliers),]



betas <- betas[,rownames(pheno)]
pcAll <- prcomp(t(betas))
# Summarise the proportion of variance per principal component as a standardized percentage
PoV <- pcAll$sdev^2/sum(pcAll$sdev^2) * 100
cumpro <- cumsum(pcAll$sdev^2/sum(pcAll$sdev^2)) # summarise the cumulative proportion of variance explained 
# Create summary plots and export

# Summarise cumulative variance explained per PC

barplot(PoV[1:20],ylab = "Percentage of Variance",main = "PCA")





# Define subsetting variables 
testVars1 <- c("Sex","Age","Post.Mortem.Delay","Braak.LB.Numeric","Braak.Tangle","Brain_Region") # Biological
testVars2 <- c(c(paste("plateNumber",c(1:12),sep = "_"),"Brain.Bank_Imperial Brain Bank",            
                 "Brain.Bank_Newcastle Brain Tissue Resource", "Brain.Bank_Oxford Brain Bank",
                 "Brain.Bank_QSBB"  ))


pheno <- dummy_cols(pheno,select_columns = c("plateNumber","Brain.Bank"))# Make dummy_cols for technical effects

rownames(pheno) <- pheno$Basename

pheno$Sex <- as.numeric(as.factor(pheno$Gender)) # Make sex numeric
pheno[which(is.na(pheno$Post.Mortem.Delay)),"Post.Mortem.Delay"] <- mean(pheno$Post.Mortem.Delay,na.rm=TRUE) # mean impute PMI
pheno$Brain_Region <- as.numeric(as.factor(pheno$Brain_Region))

##### Run new celltype deconvolution #####
cellSet<-projectCellTypeWithError(betas, modelBrainCoef[["IDOL"]][[8]])
cellSet <- as.data.frame(cellSet)

pheno <- pheno[colnames(betas),]
pheno <- cbind(pheno,cellSet[rownames(pheno),])

# This chunk here makes the correlation matrix for plotting: First with Biological variables
princomps<- pcAll$x[pheno$Basename,c(1:20)] # Extract principal components into dataframe
testCor <- cor.mtest(data.frame(pheno[,testVars2],princomps),conf.level=0.95) # Create correlation matrix
testCorp <- testCor$p
testCorp <- testCorp[1:length(testVars2),(length(testVars2)+1):ncol(testCorp)]
corPCA <- cor(princomps,pheno[,testVars2]) # create matrix of pc's and test variables

pheno$PC4 <- pcAll$x[pheno$Basename,"PC4"]

# Outliers are evident on row C on plate 9 and on plate 7, these correspond to 
# known processing outliers and will be removed
ggplot(pheno, aes(x = as.factor(rowIndex), y = PC4))+
  geom_jitter()+geom_boxplot()+
  facet_grid(~as.factor(plateNumber))


# these correspond to the following numbers
Outliers <- rownames(pheno[which(pheno$PC4 > 15),])
# c("203991400162_R01C01",
#               "203991390061_R01C01",
#               "203991390103_R03C01",
#               "203991390103_R08C01",
#               "203991400112_R03C01",
#               "203973690002_R03C01",
#               "203973690024_R03C01",
#               "203977090074_R03C01",
#               "203977090067_R03C01",
#               "203991390089_R03C01",
#               "203991390053_R03C01",
#               "203991410028_R03C01",
#               "203991410111_R03C01",
#               "203991410053_R03C01")



pheno <- pheno[-which(rownames(pheno) %in% Outliers),]

pheno$Clinical.Cohort <- gsub('[[:digit:]]+', '', pheno$Cohort.Bin)# Subset clinical cohorts
pheno <- dummy_cols(pheno,select_columns = c("plateNumber","Brain.Bank"))# Make dummy_cols for technical effects
pheno <- dummy_cols(pheno,select_columns = c("Clinical.Cohort"))
pheno$Sex <- as.numeric(as.factor(pheno$Gender)) # Make sex numeric
pheno$Brain_Region <- as.numeric(as.factor(pheno$Brain_Region))
pheno[which(is.na(pheno$Post.Mortem.Delay)),"Post.Mortem.Delay"] <- mean(pheno$Post.Mortem.Delay,na.rm=TRUE) # mean impute PMI
# Subset cohort to remove samples marked as LB controls but with high braak NFT pathology
phenoP <- pheno[-which(pheno$Clinical.Cohort == "CON" & pheno$Braak.Tangle > 2),]

#remove redundant/outdated collumns
redundantCols <- c("DoubleN","NeuNP","Sox10P","error","nCGmissing","NeuN_neg","NeuN_pos","CellTypePredictionError","PC4")
phenoP <- phenoP[,-which(colnames(phenoP) %in% redundantCols)]
rownames(phenoP) <- phenoP$Basename

##### Run new celltype deconvolution #####
betasP <-  betas[,phenoP$Basename]

cellSet<-projectCellTypeWithError(betasP, modelBrainCoef[["IDOL"]][[8]])
phenoP <- cbind(phenoP,cellSet[rownames(phenoP),])

save(phenoP, file = "Meth/Methylation/QC/ewasPheno.Rdata")
