
load("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/sncaAllQTL_fullcovar.Rdata")

.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(tidyverse)
library(magrittr)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(cowplot)
library(stringr)
setwd("/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/")
###### Read in summary statistics#]####
nalls <- fread("/lustre/projects/Research_Project-T112069/Genetics/sumstats/nallsEtAl2019_forPRSice2.txt")
chia <- fread("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/LBD_2021_formatted.txt") 
chia <- chia[which(chia$CHR == "4"),]
nalls <- nalls[which(nalls$CHR == "4",)]

###### Read in bim File and QTL data#####
load("/lustre/projects/Research_Project-T112069/Meth/QTL/Data/SNCA.matrixEQTL.chr4/SNCA.matrixEQTL.chr4.RData")
cisQTL <- as.data.frame(me$cis$eqtls)
# Format QTL data
cisQTL$CHR <- gsub(":","",str_sub(cisQTL$snp,end = 2L))
cisQTL$BP <-str_sub(cisQTL$snp,start = 3L,end = -3L)
cisQTL$BP <- as.numeric(cisQTL$BP)

# Format bim file for consistency
bim <- read.table("/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/filtered_chr_bp.bim")
bim$V2 <- paste0(bim$V1,":",bim$V4,"_",bim$V5)
write.table(bim,file = "/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/filtered_chr_bp.bim",quote = F,row.names = F,col.names = F)

## Define genomic locations for SNCA:
bim <- bim[which(bim$V1 == 4),]
bim <- bim[which((bim$V4 >= min(cisQTL$BP)) & (bim$V4 <= max(cisQTL$BP))),]
write.table(bim$V2, file ="/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/sncaTargets.txt",quote = F, row.names = F,col.names = F )

# Use plink to calculate minor allele frequencies for the 
# ## Plink
# module load PLINK/1.07-foss-2015b
# plink --bfile filtered_chr_bp --extract sncaTargets.txt --make-bed --out snca_extracted
# # Calculate allele frequencies (MAF)
# plink --bfile snca_extracted --freq --out snca_maf
# plink --bfile snca_extracted --r --matrix

#Subset the same region for the summary statistics
chia <- chia[which(chia$POS >= min(cisQTL$BP) & chia$POS <= max(cisQTL$BP)),]
nalls <- nalls[which(nalls$BP >= min(cisQTL$BP) & nalls$BP <= max(cisQTL$BP)),]

#The nalls allele frequency is in a differing format. Adjust to haromize
for(x in 1:nrow(nalls)){
  if(nalls[x,"freq"] <=0.5){
    nalls[x,"CHR_BP_A"] <- paste0(nalls[x,"CHR"],":",nalls[x,"BP"],"_",nalls[x,"A1"])
    nalls[x,"MAF"] <- nalls[x,"freq"]
  }else if(nalls[x,"freq"] >= 0.5){
    nalls[x,"CHR_BP_A"] <- paste0(nalls[x,"CHR"],":",nalls[x,"BP"],"_",nalls[x,"A2"])
    nalls[x,"MAF"] <- 1 - nalls[x,"freq"]
  }
}

#Final formatting for consistency (convert SE to Variance)
nallsFormatted <- nalls[,c("CHR_BP_A","b","se","p","MAF")]
nallsFormatted$se <- nallsFormatted$se^2
nallsFormatted$N <- nalls$N_cases + nalls$N_controls

chia$CHR_BP_A <- paste0(chia$CHR,":",chia$POS,"_",chia$A1)
chiaFormatted <- chia[,c("CHR_BP_A","BETA","SE","P","MAF","N")]
chiaFormatted$SE <- chiaFormatted$SE^2

colnames(nallsFormatted) <- c("SNP","BETA","VAR","P","MAF","N")
colnames(chiaFormatted) <- c("SNP","BETA","VAR","P","MAF","N")

# For cisQTLs, merge MAF and format
QTL_MAF <- read.table("/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/snca_maf.frq",header = T)
rownames(QTL_MAF) <- QTL_MAF$SNP
cisQTL$MAF <- QTL_MAF[cisQTL$snps,"MAF"]

cisQTL$se <- abs(cisQTL$beta)/abs(cisQTL$statistic)
cisQTL$VAR <- cisQTL$se^2
cisQTL <- cisQTL[-which(is.na(cisQTL$MAF)),]

###### Load and process the conditional results (Pihlstrom 2019)

#The data is annotated in the format:
# Unconditioned (raw/uncontrolled signal) = "Conditioned0"
# PD 5' Signal (rs356182 conditioned) = "Conditioned_1"
# PD 3' Signal (conditioned on rs763443) = "Conditioned_2"
# PD nominal signal (conditioned on rs356182 & rs763443) = "pd_Conditioned1_2"

pd_Conditioned0 <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/conditional/all_snca_0_cond_for_Josh.assoc.logistic",header = T)
pd_Conditioned1 <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/conditional/all_snca_1_cond_for_Josh.assoc.logistic",header = T)
pd_Conditioned2 <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/conditional/all_snca_2_cond_for_Josh.assoc.logistic",header = T)
pd_Conditioned1_2 <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/conditional/all_snca_1and2_cond_for_Josh.assoc.logistic",header = T)

pd_Conditioned0$SNP <- paste(pd_Conditioned0$SNP,pd_Conditioned0$A1,sep = "_")
pd_Conditioned1$SNP <- paste(pd_Conditioned1$SNP,pd_Conditioned1$A1,sep = "_")
pd_Conditioned2$SNP <- paste(pd_Conditioned2$SNP,pd_Conditioned2$A1,sep = "_")
pd_Conditioned1_2$SNP <- paste(pd_Conditioned1_2$SNP,pd_Conditioned1_2$A1,sep = "_")

pd_Conditioned0$BETA <- log(pd_Conditioned0$OR)
pd_Conditioned1$BETA <- log(pd_Conditioned1$OR)
pd_Conditioned2$BETA <- log(pd_Conditioned2$OR) 
pd_Conditioned1_2$BETA <- log(pd_Conditioned1_2$OR)

pd_Conditioned0$SE <- (log(pd_Conditioned0$U95) - log(pd_Conditioned0$L95)) / (2 * 1.96)
pd_Conditioned1$SE <- (log(pd_Conditioned1$U95) - log(pd_Conditioned1$L95)) / (2 * 1.96)
pd_Conditioned2$SE <- (log(pd_Conditioned2$U95) - log(pd_Conditioned2$L95)) / (2 * 1.96) 
pd_Conditioned1_2$SE <- (log(pd_Conditioned1_2$U95) - log(pd_Conditioned1_2$L95)) / (2 * 1.96)

pd_Conditioned0$VAR <- pd_Conditioned0$SE^2
pd_Conditioned1$VAR <- pd_Conditioned1$SE^2
pd_Conditioned2$VAR <- pd_Conditioned2$SE^2
pd_Conditioned1_2$VAR <- pd_Conditioned1_2$SE^2

pd_Conditioned0$SET <- "Unconditioned"
pd_Conditioned1$SET <- "Conditioned_1"
pd_Conditioned2$SET <- "Conditioned_2"
pd_Conditioned1_2$SET <- "Conditioned_1_2"

saveNames <- c("SNP","BETA","VAR","P","SET","NMISS")

allConditioned <- rbind(pd_Conditioned0[,saveNames],
                        pd_Conditioned1[,saveNames],
                        pd_Conditioned2[,saveNames],
                        pd_Conditioned1_2[,saveNames])

allConditioned <- allConditioned[-which(is.na(allConditioned$BETA)),]

# Add MAF
checkMAF <- read.table("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/coloc/conditional/all_snca_MAFs_for_Josh.frq",header = T)
checkMAF$snps <- paste(checkMAF$SNP,checkMAF$A1, sep = "_")
rownames(checkMAF) <- checkMAF$snps

allConditioned$MAF <- checkMAF[allConditioned$SNP,"MAF"]

#Save all outputs
save(chiaFormatted,nallsFormatted,cisQTL,allConditioned,file = "/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/inputsFormatted.Rdata")


############################################
# COLOC TESTING

load("/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/inputsFormatted.Rdata")
#Subset overlaps
nallsFormatted <- nallsFormatted[which(nallsFormatted$SNP %in% cisQTL$snps),]
chiaFormatted <- chiaFormatted[which(chiaFormatted$SNP %in% cisQTL$snps),]
allConditioned <- allConditioned[which(allConditioned$SNP %in% cisQTL$snps),]
#Annotate significant QTLs
sigCpGs <- unique(cisQTL[which(cisQTL$FDR < 0.05),"gene"])


library(coloc)
# Format summary statistics
chiaList <- list(beta = chiaFormatted$BETA,
                 varbeta = chiaFormatted$VAR,
                 type = "cc",
                 snp = chiaFormatted$SNP,
                 MAF = chiaFormatted$MAF,
                 N = chiaFormatted$N)

nallsList <- list(beta = nallsFormatted$BETA,
                  varbeta = nallsFormatted$VAR,
                  type = "cc",
                  snp = nallsFormatted$SNP,
                  MAF = nallsFormatted$MAF,
                  N = nallsFormatted$N)

allConditioned <- allConditioned[-which(is.na(allConditioned$MAF)),]

colocRes <- data.frame(nsnps=NA,PP.H0=NA,PP.H1=NA,PP.H2=NA,PP.H3=NA,PP.H4=NA,GWAS=NA,cpg = NA)

# subset SNCA region
cisQTL <- cisQTL[which(cisQTL$BP %in% c(90500000:91400000)),]

#The below for-loop tests colocalization
for(t in unique(cisQTL$gene)){
  print(t)
  qtlList <- list(beta = cisQTL[which(cisQTL$gene == t),"beta"],
                  varbeta = cisQTL[which(cisQTL$gene == t),"VAR"],
                  type = "quant",
                  snp = cisQTL[which(cisQTL$gene == t),"snps"],
                  MAF = cisQTL[which(cisQTL$gene == t),"MAF"],
                  N = 386)
  print("1")
  test.colocCHIA <- coloc.abf(qtlList,chiaList)
  print("2")
  test.colocNALLS <- coloc.abf(qtlList,nallsList)
  print("3")
  colocRes <- rbind(colocRes,c(test.colocCHIA$summary,"DLB (Chia)",t))
  colocRes <- rbind(colocRes,c(test.colocNALLS$summary,"PD (Nalls)",t))
  for(g in unique(allConditioned$SET)){
    print(paste(g))
    condList <- list(beta = allConditioned[which(allConditioned$SET == g),"BETA"],
                         varbeta = allConditioned[which(allConditioned$SET == g),"VAR"],
                         type = "cc",
                         snp = allConditioned[which(allConditioned$SET == g),"SNP"],
                         MAF = allConditioned[which(allConditioned$SET == g),"MAF"],
                         N = allConditioned[which(allConditioned$SET == g),"NMISS"])
                     test.colocCOND <- coloc.abf(qtlList,condList)
                     colocRes <- rbind(colocRes,c(test.colocCOND$summary,g,t))
  }
}


# The below code formats the coloc results to long for plotting
colocRes <- colocRes[-1,]
colocRes$nsnps <- as.numeric(colocRes$nsnps)
colocRes_long <- colocRes %>%
  pivot_longer(
    cols = starts_with("PP"),   
    names_to = "PP_type",     
    values_to = "PP_value"        
  )

# Format data
colocRes[,c("nsnps","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")] <- as.numeric(unlist(colocRes[,c("nsnps","PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")]))

write.csv(colocRes,file = "/lustre/projects/Research_Project-T112069/Meth/QTL/Coloc/finalRes_Conditional_Nalls_Chia.csv")



## The below code plots the main colocalisation findings

sigCpGs <- unique(colocRes[which(colocRes$PP.H4 >= 0.8),"cpg"])
cisQTL$plotAnnot <- NA
cisQTL[which(cisQTL$gene %in% sigCpGs),"plotAnnot"] <- cisQTL[which(cisQTL$gene %in% sigCpGs),"gene"] 
cisQTL <- cisQTL[rev(order(cisQTL$plotAnnot)),]

minBP <- min(cisQTL[-which(is.na(cisQTL$plotAnnot)),"BP"])
maxBP <- 90900000 # shorten region for plotting



library(stringr)
chiaFormatted$CHR <- gsub(":","",str_sub(chiaFormatted$SNP,end = 2L))
chiaFormatted$BP <-as.numeric(str_sub(chiaFormatted$SNP,start = 3L,end = -3L))

nallsFormatted$CHR <- gsub(":","",str_sub(nallsFormatted$SNP,end = 2L))
nallsFormatted$BP <-as.numeric(str_sub(nallsFormatted$SNP,start = 3L,end = -3L))

allConditioned$CHR <- gsub(":","",str_sub(allConditioned$SNP,end = 2L))
allConditioned$BP <-as.numeric(str_sub(allConditioned$SNP,start = 3L,end = -3L))

qtlPlot <- ggplot(cisQTL[-which(is.na(cisQTL$plotAnnot)),],aes(x = as.numeric(BP), y = -log10(pvalue),color = as.factor(plotAnnot)))+
  geom_point()+scale_color_brewer(palette = "Set2",na.value = "grey")+
  theme_bw()+labs(x = NULL, y = "-log10(P-mQTL)",color = "CpG Site")+
  theme(legend.position = c(0.8,0.7),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ylab("-log10(mQTL P-value)")+xlim(minBP,maxBP)+geom_vline(xintercept = 90763360,color = "darkred",linetype = "dashed")

chiaPlot <- ggplot(chiaFormatted[which(chiaFormatted$BP %in% c(minBP:maxBP)),],aes(x = as.numeric(BP), y = -log10(P)))+
  geom_point()+scale_color_brewer(palette = "Set2",na.value = "grey")+
  theme_bw()+labs(x = NULL, y = "-log10(P-meQTL)",color = "CpG Site")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ylab("-log10(LBD Status P-value)")+geom_vline(xintercept = 90763360,color = "darkred",linetype = "dashed")


nallsPlot <- ggplot(nallsFormatted[which(nallsFormatted$BP %in% c(minBP:maxBP)),],aes(x = as.numeric(BP), y = -log10(P)))+
  geom_point()+scale_color_brewer(palette = "Set2",na.value = "grey")+
  theme_bw()+
  theme(legend.position = c(0.8,0.6))+ylab("-log10(LBD Status P-value)")


Con1Plot <- ggplot(allConditioned[which((allConditioned$BP %in% c(minBP:maxBP)) & (allConditioned$SET == "Conditioned_1")),],aes(x = as.numeric(BP), y = -log10(P)))+
  geom_point()+scale_color_brewer(palette = "Set2",na.value = "grey")+
  theme_bw()+
  theme(legend.position = c(0.8,0.8),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ylab("-log10(PD Status P-value)")+geom_vline(xintercept = 90763360,color = "darkred",linetype = "dashed")

Con2Plot <- ggplot(allConditioned[which((allConditioned$BP %in% c(minBP:maxBP)) & (allConditioned$SET == "Conditioned_2")),],aes(x = as.numeric(BP), y = -log10(P)))+
  geom_point()+scale_color_brewer(palette = "Set2",na.value = "grey")+
  theme_bw()+ylab("-log10(PD Status P-value)")+xlab("Basepair Coordinate")+theme(plot.title = element_text(margin = margin(b = -20)))+geom_vline(xintercept = 90763360,color = "darkred",linetype = "dashed")

library(cowplot)

plot_grid(qtlPlot,chiaPlot,Con1Plot,Con2Plot,ncol = 1, align = "v",rel_heights = c(0.24,0.24,0.24,0.28))




plotQTL <- cisQTL[which(cisQTL$gene == "cg01966878"),]
rownames(plotQTL) <- plotQTL$snps

chiaPlot2 <- as.data.frame(chiaFormatted)
rownames(chiaPlot2) <- chiaPlot2$SNP
plotLabels <- intersect(rownames(plotQTL),rownames(chiaPlot2))

chiaFrame <- data.frame(qtlHit = -log10(plotQTL[plotLabels,"pvalue"]),chiaHit = -log10(chiaPlot2[plotLabels,"P"]))
chiaScatter <- ggplot(chiaFrame,aes(x = qtlHit, y = chiaHit))+
  geom_point(color = "#1B9E77")+theme_bw()+
  ylab("-log10(LBD Status P-value)")+
  theme(legend.position = c(0.8,0.8),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


condPlot1 <- as.data.frame(allConditioned[which((allConditioned$BP %in% c(minBP:maxBP)) & (allConditioned$SET == "Conditioned_1")),])
rownames(condPlot1) <- condPlot1$SNP
plotLabels <- intersect(rownames(plotQTL),rownames(condPlot1))
condFrame1 <- data.frame(qtlHit = -log10(plotQTL[plotLabels,"pvalue"]),condHit1 = -log10(condPlot1[plotLabels,"P"]))

condScatter1 <- ggplot(condFrame1,aes(x = qtlHit, y = condHit1))+
  geom_point(color = "#1B9E77")+theme_bw()+
  ylab("-log10(PD Status P-value)")+
  theme(legend.position = c(0.8,0.8),axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

condPlot2 <- as.data.frame(allConditioned[which((allConditioned$BP %in% c(minBP:maxBP)) & (allConditioned$SET == "Conditioned_2")),])
rownames(condPlot2) <- condPlot2$SNP
plotLabels <- intersect(rownames(plotQTL),rownames(condPlot2))
condFrame2 <- data.frame(qtlHit = -log10(plotQTL[plotLabels,"pvalue"]),condHit2 = -log10(condPlot2[plotLabels,"P"]))

condScatter2 <- ggplot(condFrame2,aes(x = qtlHit, y = condHit2))+
  geom_point(color = "#1B9E77")+theme_bw()+
  ylab("-log10(PD Status P-value)")+
  xlab("-log10(cg01966878 mQTL P-value)")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

placePlot <- ggplot()+theme_void()

plot_grid(qtlPlot,placePlot,chiaPlot,chiaScatter,Con1Plot,condScatter1,Con2Plot,condScatter2,
          ncol = 2,rel_heights = c(0.24,0.24,0.24,0.28),rel_widths = c(0.6,0.4), align = "v",axis = "lr")

save(cisQTL,chiaFormatted,allConditioned,condFrame1,condFrame2,chiaFrame, file = "/lustre/projects/Research_Project-T112069/Meth/QTL/ColocPlottingData.Rdata")
