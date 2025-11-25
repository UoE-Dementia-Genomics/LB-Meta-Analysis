# Browns Overlap analysis
setwd("/lustre/projects/Research_Project-T112069/")
.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(BiocManager)
library(EmpiricalBrownsMethod)
library(data.table)
library(BBmisc)

setwd("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/")

#Load in regions from the PD GWAS 2024 (Kim et al.)
genomicRegions <- read.csv("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/PD_Regions2024.csv")



load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
# Process meta summary statistics for DMR analysis
resMetaCC <- as.data.frame(allSumstat)
resMetaCC <- resMetaCC[order(resMetaCC$Fixed_P),]
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
row.names(epicMani) <- epicMani[,1]


##### Load packages and set working directory #####
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
library(GEOquery)
library(ggplot2)
library(dplyr)
library(data.table)
library(FlowSorted.DLPFC.450k)
library(CETYGO)

##### Load data: UKBBN #####
# Load in data
load("Meth/Methylation/QC/FinalData/betasPFC.Rdata")
load("Meth/Methylation/QC/ewasPheno.Rdata")
load("Meth/Methylation/QC/FinalData/betasCNG.Rdata")

# Merge data and index
betas <- cbind(betas.cng,betas.pfc)
betas <- betas[,rownames(phenoP)]
phenoP$count <- 1
betasP <-  betas[,phenoP$Basename]



##### Load data: NBB #####
# Load in data
NBB<-getGEO(filename='/lustre/projects/Research_Project-T112069/Meth/NBB/GSE203332_series_matrix.txt',
            GSEMatrix=TRUE,getGPL = FALSE)
NBBpheno <- pData(object = NBB)

# Process the pheno filename
NBBpheno <- NBBpheno[,c(49:65)]
colnames(NBBpheno) <-  gsub(pattern = ":ch1", replacement = "",colnames(NBBpheno)) # sort out the colnames
NBBpheno <- NBBpheno[-which(NBBpheno$neuropathological_diagnosis %in% c("AD","mixed_AD_LBD")),]
rownames(NBBpheno) <- paste(NBBpheno$`slide`,NBBpheno$`slide_position`,sep = "_")
NBBgeo <- NBBpheno
NBBpheno <- read.table("/lustre/projects/Research_Project-T112069/Meth/NBB/NBB_pheno.txt", header = T)
rownames(NBBpheno) <- NBBpheno$PlatePos_ID
NBBpheno <- cbind(NBBpheno, NBBgeo[rownames(NBBpheno),c("braak_nft_stage","cerad_stage")])
NBBpheno$braak_nft_stage <- as.numeric(NBBpheno$braak_nft_stage)
NBBpheno$PlateNum <- as.numeric(as.factor(NBBpheno$Plate))
NBBpheno$SexNum <- as.numeric(as.factor(NBBpheno$Sex))

#Filter beta matrix
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
filter <- read.table("/lustre/projects/Research_Project-T112069/Meth/reference/SNPProbes_McCartney.txt")
row.names(epicMani) <- epicMani[,1]
NBBbetas <- fread("/lustre/projects/Research_Project-T112069/Meth/NBB/GSE203332_Matrix_processed_dasen_betas.txt")
NBBbetas <- as.data.frame(NBBbetas)
rownames(NBBbetas) <- NBBbetas$ID_REF
NBBbetas$ID_REF <- NULL
NBBbetas <- NBBbetas[,NBBpheno$PlatePos_ID]
xyProbes <- rownames(epicMani[which(epicMani$CHR %in% c("X","Y")),])
NBBbetas <- as.data.frame(NBBbetas)
NBBbetas <- NBBbetas[-which(rownames(NBBbetas) %in% xyProbes),]

##### Load data: BDR #####
load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Pheno.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/BDR/BDR_LB_Betas.Rdata")



##### Process Methylation data #####
# Retain overlapping cpgs
overlap <- intersect(intersect(rownames(betasP), rownames(betasLB)),rownames(NBBbetas))


#z-score normalize betas
scaledUK <- t(scale(t(betasP)))
scaledBDR <- t(scale(t(betasLB)))
scaledNBB <- t(scale(t(NBBbetas)))

#merge
megaBetas <- cbind(scaledUK[overlap,],scaledBDR[overlap,],scaledNBB[overlap,])


##### PD Browns enrichment test
probeIDs <- c()
brownPDres <- genomicRegions
brownPDres[,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- NA

# Loop to run genomic enrichment on PD regions
for(x in 1:nrow(genomicRegions)){
  chr <- genomicRegions[x,"CHR"]
  start <- genomicRegions[x,"Locus.start"]
  end <- genomicRegions[x,"Locus.end"]
  print(c(x,chr,start,end))
  probeIDs <- rownames(epicMani[which(epicMani$CHR ==  chr & epicMani$MAPINFO >= start & epicMani$MAPINFO <= end),])
  
  nProbes <- length(probeIDs)
  
  if(nProbes == 0 ){
    brownPDres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(rep(NA,4),nProbes)
  }else{
  
  probeIDs <- intersect(probeIDs,rownames(resMetaCC))
  browns.Loc.i<-empiricalBrownsMethod(megaBetas[probeIDs,], resMetaCC[probeIDs,3], extra_info=T)
  
  brownPDres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(unlist(browns.Loc.i),nProbes)
  }
}


# AD genomic regions (Kunkle)
adSet <- matrix(c("chr1","207679307","207850539","CR1",
         "chr2","127882182","127894615","BIN1",
         "chr2","233976593","233981912","INPP5D",
         "chr6","32395036","32636434","HLA-DRB1",
         "chr6","40706366","41365821","TREM2",
         "chr6","47412916","47628558","CD2AP",
         "chr7","99932049","100190116","NYAP1",
         "chr7","143099107","143109208","EPHA1",
         "chr8","27195121","27238052","PTK2B",
         "chr8","27456253","27468503","CLU",
         "chr10","11703491","11723257","ECHDC3",
         "chr11","47372377","47466790","SPI1",
         "chr11","59856028","60097777","MS4A2",
         "chr11","85670385","85868640","PICALM",
         "chr11","121433926","121461593","SORL1",
         "chr14","53293307","53462216","FERMT2",
         "chr14","92926952","92957176","SLC24A4",
         "chr15","58873555","59120077","ADAM10",
         "chr16","19706199","19867021","IQCK",
         "chr16","79355847","79355847","WWOX",
         "chr17","61499732","61543566","ACE",
         "chr19","1050130","1075979","ABCA7",
         "chr20","54979828","55025377","CASS4",
         "chr21","28146668","28166355","ADAMTS1",
         "chr19","45116911","46318605","APOE"),ncol=4,byrow=T)

#Process data for browns
adSet <- as.data.frame(adSet)
colnames(adSet) <- c("CHR", "Locus.start", "Locus.end")
adSet$CHR <- gsub("chr","", adSet$CHR)
brownADres <- adSet
brownADres[,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- NA
brownADres$Locus.start <- as.numeric(brownADres$Locus.start)
brownADres$Locus.end <- as.numeric(brownADres$Locus.end)

#For loop to test enrichment
for(x in 1:nrow(brownADres)){
  chr <- brownADres[x,"CHR"]
  start <- as.numeric(brownADres[x,"Locus.start"])
  end <- as.numeric(brownADres[x,"Locus.end"])
  print(c(x,chr,start,end))
  probeIDs <- rownames(epicMani[which(epicMani$CHR ==  chr & epicMani$MAPINFO >= start & epicMani$MAPINFO <= end),])
  
  nProbes <- length(probeIDs)
  if(nProbes <= 1 ){
    brownADres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(rep(NA,4),nProbes)
  }else{
    
    probeIDs <- intersect(probeIDs,rownames(resMetaCC))
    browns.Loc.i<-empiricalBrownsMethod(megaBetas[probeIDs,], resMetaCC[probeIDs,3], extra_info=T)
    
    brownADres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(unlist(browns.Loc.i),nProbes)
  }
}

write.csv(brownADres, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/adResKunkle.csv")
write.csv(brownPDres, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/pdRes2024.csv")


# LBD genomic regions (Chia) 
dlbSet <- matrix(c("19","45324138","45429708","APOE",
                  "2","127891427","127892810","BIN1",
                  "1","155116151","156063880","GBA",
                  "4","939087","939087","TMEM175",
                  "4", "90743331","91005096","SNCA"),ncol=4,byrow=T)

dlbSet <- as.data.frame(dlbSet)
colnames(dlbSet) <- c("CHR", "Locus.start", "Locus.end")

brownDLBres <- dlbSet
brownDLBres[,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- NA
brownDLBres$Locus.start <- as.numeric(brownDLBres$Locus.start)
brownDLBres$Locus.end <- as.numeric(brownDLBres$Locus.end)

for(x in 1:nrow(brownDLBres)){
  chr <- brownDLBres[x,"CHR"]
  start <- as.numeric(brownDLBres[x,"Locus.start"])
  end <- as.numeric(brownDLBres[x,"Locus.end"])
  print(c(x,chr,start,end))
  probeIDs <- rownames(epicMani[which(epicMani$CHR ==  chr & epicMani$MAPINFO >= start & epicMani$MAPINFO <= end),])
  
  nProbes <- length(probeIDs)
  if(nProbes <= 1 ){
    brownDLBres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(rep(NA,4),nProbes)
  }else{
    
    probeIDs <- intersect(probeIDs,rownames(resMetaCC))
    browns.Loc.i<-empiricalBrownsMethod(megaBetas[probeIDs,], resMetaCC[probeIDs,3], extra_info=T)
    
    brownDLBres[x,c("P_test","P_Fisher", "Scale_Factor_C","DF","N")] <- c(unlist(browns.Loc.i),nProbes)
  }
}



write.csv(brownDLBres,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/dlbResChia.csv")

#The below code merges the results from all above analyses

colnames(brownADres)[4] <- "ID"
brownADres$Study <- "AD_Kunkle"
colnames(brownDLBres)[4] <- "ID"
brownPDres$Study <- "PD_Kim"
colnames(brownPDres)[2] <- "ID"
brownDLBres$Study <- "DLB_Chia"

testCols <- c("CHR" ,"Locus.start" ,"Locus.end" , "ID","P_test" ,     "Scale_Factor_C" ,"DF","N" , "Study")

#Merge and multiple testing correct all
allRes <- rbind(brownPDres[,testCols],brownADres[,testCols],brownDLBres[,testCols])
allRes <- allRes[order(allRes$P_test),]
allRes$P_adj <- p.adjust(allRes$P_test,method = "fdr")

# save results
write.csv(allRes, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/Browns/allRes.csv")
