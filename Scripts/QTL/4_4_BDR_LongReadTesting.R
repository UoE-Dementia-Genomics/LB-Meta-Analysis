
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


# The below section takes the BDR data (available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197305)
# And extract the genotype for the lead SNP of interest (rs763443, chr4:90819961), as well as methylation in the SNCA region for testing
setwd("/lustre/projects/Research_Project-T112069/Meth/BDR")
bdrPhenoAll <- read.csv("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_pheno_nPath_PCs.csv",header = T,row.names = 1)

load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_combined_1221_QCd.rdat")
load("/lustre/projects/Research_Project-T112069/Meth/Methylation/QC/Pheno.Rdata")

# bim <- fread("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/bdrValidation/BDR_imputed_EUR_QCd.bim")
# plink --bfile BDR_imputed_EUR_QCd \
# --snp "4:90819961" \
# --recodeA \
# --out /lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/bdrValidation/pd_SNP_Extract

DLBleadSNP <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/bdrValidation/dlbLeadbdr.ped")
DLBleadSNP$geno <-paste0(DLBleadSNP$V7,DLBleadSNP$V8)

PDleadSNP <-  read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/bdrValidation/pd_SNP_Extract.raw",header = T)
DLBleadSNP$X4.90819961_T <- PDleadSNP$X4.90819961_T

rownames(DLBleadSNP) <- DLBleadSNP$V1
colnames(DLBleadSNP)[1] <- "DNA_ID"
bdrPhenoAll <- left_join(bdrPhenoAll, DLBleadSNP[,c("DNA_ID","geno","X4.90819961_T")])

rownames(bdrPhenoAll) <- bdrPhenoAll$X

bdrPhenoAll$cg01966878 <- betas["cg01966878",rownames(bdrPhenoAll)]
bdrPhenoAll$cg15133208 <- betas["cg15133208",rownames(bdrPhenoAll)]
bdrPhenoAll$cg20003494 <- betas["cg20003494",rownames(bdrPhenoAll)]
bdrPhenoAll$cg14346243 <- betas["cg14346243",rownames(bdrPhenoAll)]
save(bdrPhenoAll,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/bdrValidation/bdrValidation.Rdata")

### This section takes the BDR data processed above, merges it with Long read sequencing data from Leung 2024 (https://zenodo.org/records/10784918) for isoform testing 
load("C:/Users/jh1159/Dropbox/Josh Harvey PhD/Bioinformatics/MRC/bdr/bdrValidation.Rdata")
library(data.table)
longReadHits <- fread("C:/Users/jh1159/Dropbox/Josh Harvey PhD/Bioinformatics/MRC/bdr/ontBDR_collapsed_RulesFilter_result_classification.targetgenes_counts_filtered.txt")

library(edgeR)
library(dplyr)
longReadHits <- as.data.frame(longReadHits)
colSelect <- colnames(longReadHits)[grep("B1",colnames(longReadHits))] # select samples
longReadHits <- longReadHits[,c("associated_gene","structural_category","associated_transcript",colSelect)] # keep gene and trancript annotations / samples
longReadHits <- longReadHits[-which(longReadHits$associated_transcript == "novel"),] #only take known isoforms
longReadHits <- longReadHits[which(longReadHits$structural_category == "FSM"),] #only take full splice matches
remname <- names(which(colSums(longReadHits[,-c(1:3)]) == 0)) #subset out unwanted collumnes
longReadHits <- longReadHits[,-which(colnames(longReadHits) %in% remname)]
longReadHits <- longReadHits[-which(rowSums(longReadHits[,-c(1:3)]) == 0),]
rowData <- longReadHits[,c(1:3)] #Subset row information
colSelect <- colnames(longReadHits)[grep("B1",colnames(longReadHits))]

#Concatenate duplicate transcripts
concat <- longReadHits[,c("associated_transcript",colSelect)] %>%
  group_by(associated_transcript) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
concat <- as.data.frame(concat)
rownames(concat) <- concat$associated_transcript
concat <- concat[,-1]

#Format for DGE
dge <- DGEList(counts = concat)

# Compute log2CPM using edgeR
log2cpm_matrix <- cpm(dge, log = TRUE)

#Remove duplicates
rowDataSing <- rowData[-which(duplicated(rowData)),]
#Subset SNCA data
sncaMat <- log2cpm_matrix[rowDataSing[which(rowDataSing$associated_gene == "SNCA"),"associated_transcript"],]
#Format full BDR pheno file for merging
bdrPhenoAll$joinVar <- clean_vector <- gsub("[^a-zA-Z0-9]", "", bdrPhenoAll$BBNId)
#Subset collumns of interest
bdrSum <- bdrPhenoAll[which(bdrPhenoAll$BR == "Prefrontal"),c("joinVar","geno","X4.90819961_T","cg01966878","cg15133208","cg14346243","cg20003494","Gender","Age","NeuN","Double")]

#Format SNCA expression data for merging
colnames(sncaMat) <- gsub("B1.",replacement = "",colnames(sncaMat))
library(stringr)
colnames(sncaMat) <-str_sub(colnames(sncaMat),end = -3L)

# allign pheno and expression data
rownames(bdrSum) <- bdrSum$joinVar
bdrSum <- bdrSum[colnames(sncaMat),]

# Function for testing linear association of individual isoform expression and methylation or genotype 
run_lm <- function(expression_values,site1,site2) {
  model1 <- lm(expression_values ~ bdrSum[,site1]+as.factor(bdrSum$Gender)+bdrSum$Age+bdrSum$NeuN+bdrSum$Double)
  model2 <- lm(expression_values ~ bdrSum[,site2]+as.factor(bdrSum$Gender)+bdrSum$Age+bdrSum$NeuN+bdrSum$Double)
  
  # Extract p-values
  pval1 <- summary(model1)$coefficients[2, 4]  # p-value for geno
  pval2 <- summary(model2)$coefficients[2, 4]  # p-value for cg01966878
  
  return(c(pval1, pval2))
}

#format genotype for testing
bdrSum$geno_numeric <- ifelse(bdrSum$geno == "AA", 0,
                              ifelse(bdrSum$geno == "CA", 1, 
                                     ifelse(bdrSum$geno == "CC", 2, NA)))

# Run linear model test
lm_results <- t(apply(sncaMat, 1, run_lm,"geno_numeric","cg01966878"))

# This section tests the linear association
lm_results <- lm_results[order(lm_results[,1]),]
