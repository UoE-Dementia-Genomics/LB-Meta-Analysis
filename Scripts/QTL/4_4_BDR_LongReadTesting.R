
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

# Compute CPM using edgeR
cpm_matrix <- cpm(dge, log = FALSE)

# Apply inverse normal transformation
n <- ncol(cpm_matrix)
zvalues <- qnorm(ppoints(n))
z <- as.matrix(cpm_matrix)
for (i in 1:nrow(z)) z[i,] <- zvalues[order(order(z[i,]))]


#Remove duplicates
rowDataSing <- rowData[-which(duplicated(rowData)),]
#Subset SNCA data
sncaMat <- z[rowDataSing[which(rowDataSing$associated_gene == "SNCA"),"associated_transcript"],]
#Format full BDR pheno file for merging
cellDecon <- read.csv("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/bdr/bdr/BDR_fiveCellDecon.csv",header = T,row.names = 1)
bdrPhenoAll <- cbind(bdrPhenoAll,cellDecon[rownames(bdrPhenoAll),] )
# Add genotype PC's
genoPC <- read.table("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/bdr_pca.eigenvec")
bdrPhenoAll$gPC1 <- genoPC[match(bdrPhenoAll$DNA_ID,genoPC$V1),"V3"]
bdrPhenoAll$gPC2 <- genoPC[match(bdrPhenoAll$DNA_ID,genoPC$V1),"V4"]
bdrPhenoAll$gPC3 <- genoPC[match(bdrPhenoAll$DNA_ID,genoPC$V1),"V5"]


bdrSum <- bdrPhenoAll[which(bdrPhenoAll$BR == "Prefrontal"),c("joinVar","geno","cg01966878","Gender","Age","NeuNNeg_Sox10Neg_IRF8Pos",
                                                              "NeuNPos_SOX6Neg", "NeuNPos_SOX6Pos", "NeuNNeg_SOX10Pos" ,
                                                              "NeuNNeg_Sox10Neg_IRF8Neg","gPC1","gPC2","gPC3","Plate")]

#Format SNCA expression data for merging
colnames(sncaMat) <- gsub("B1.",replacement = "",colnames(sncaMat))
library(stringr)
colnames(sncaMat) <-str_sub(colnames(sncaMat),end = -3L)

# allign pheno and expression data
rownames(bdrSum) <- bdrSum$joinVar
bdrSum <- bdrSum[colnames(sncaMat),]

# Join RIN's
bdrRINS <- read.csv("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/bdr/bdr/BDR_RINS.csv",header = T)
bdrRINS$Name <- gsub("[^a-zA-Z0-9]", "", bdrRINS$BBN.ID)
rownames(bdrRINS) <- bdrRINS$Name

bdrSum$RIN <- bdrRINS[bdrSum$joinVar,"RINe"]

# Function for testing linear association of individual isoform expression and methylation or genotype 
run_lm <- function(expression_values) {
  model1 <- lm(expression_values ~ bdrSum$geno_numeric+as.factor(bdrSum$Gender)+as.factor(bdrSum$Plate)+bdrSum$Age+bdrSum$RIN+bdrSum$NeuNNeg_Sox10Neg_IRF8Pos+
                 bdrSum$NeuNPos_SOX6Neg+bdrSum$NeuNPos_SOX6Pos+bdrSum$NeuNNeg_Sox10Neg_IRF8Neg+bdrSum$gPC1+bdrSum$gPC2+bdrSum$gPC3)
  model2 <- lm(expression_values ~ bdrSum$cg01966878+as.factor(bdrSum$Gender)+as.factor(bdrSum$Plate)+bdrSum$Age+bdrSum$RIN+bdrSum$NeuNNeg_Sox10Neg_IRF8Pos+
                 bdrSum$NeuNPos_SOX6Neg+bdrSum$NeuNPos_SOX6Pos+bdrSum$NeuNNeg_Sox10Neg_IRF8Neg)
  
  # Extract p-values
  pval1 <- summary(model1)$coefficients[2, 4]  # p-value for geno
  pval2 <- summary(model2)$coefficients[2, 4]  # p-value for cg01966878
  
  return(c(pval1, pval2))
}

bdrSum$geno_numeric <- ifelse(bdrSum$geno == "AA", 0,
                              ifelse(bdrSum$geno == "CA", 1, 
                                     ifelse(bdrSum$geno == "CC", 2, NA)))

# Apply function to each gene (row-wise)
lm_results <- t(apply(sncaMat, 1, run_lm))
lm_results <- lm_results[order(lm_results[,1]),]

# Adjust P-values
lm_results_adj <- p.adjust(lm_results,method = "fdr")

# Order results by genotype associatin P-value
lm_results <- lm_results[order(lm_results[,1]),] 

# Extract expression data for plotting
bdrSum$ENST00000508895 <- sncaMat["ENST00000508895.5",]
bdrSum$ENST00000394991 <- sncaMat["ENST00000394991.8",]

# Association Plots
methENST00000508895 <- ggplot(bdrSum,aes(x = cg01966878, y = ENST00000508895))+geom_point(color = "#FC8D62")+
  geom_smooth(,method = "lm",color = "black")+
  labs(x = "cg01966878 Methylation (Beta Value)", y = "ENST00000508895 Expression (CPM)")+
  theme_bw()

genoENST00000508895 <- ggplot(bdrSum,aes(x = as.factor(geno), y = ENST00000508895))+geom_violin(fill = "#FC8D62",alpha = 0.6,color = NA)+
  geom_boxplot(width = 0.3,lwd = .8,outlier.shape = NA)+geom_jitter(width = 0.1)+
  labs(x = "rs7680557 genotype", y = "ENST00000508895 Expression (CPM)")+
  theme_bw()

methENST00000394991 <- ggplot(bdrSum,aes(x = cg01966878, y = ENST00000394991))+geom_point(color = "#E41A1C")+
  geom_smooth(,method = "lm",color = "black")+
  labs(x = "cg01966878 Methylation (Beta Value)", y = "ENST00000394991 Expression (CPM)")+
  theme_bw()

genoENST00000394991 <- ggplot(bdrSum,aes(x = as.factor(geno), y = ENST00000394991))+geom_violin(fill = "#E41A1C",alpha = 0.6,color = NA)+
  geom_boxplot(width = 0.3,lwd = .8,outlier.shape = NA)+geom_jitter(width = 0.1)+
  labs(x = "rs7680557 genotype", y = "ENST00000394991 Expression (CPM)")+
  theme_bw()


library(cowplot)
plot_grid(methENST00000394991,genoENST00000394991,methENST00000508895,genoENST00000508895)

# The below code plots the transcript structures
library(BiocManager)
library("ggtranscript")

library(biomaRt)
library(dplyr)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
snca_exons <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "exon_chrom_start",
    "exon_chrom_end"
  ),
  filters = "hgnc_symbol",
  values = "SNCA",
  mart = mart
)

# Rename columns to match ggtranscript input
snca_df <- snca_exons 
colnames(snca_df) <- c("transcript","start","end")
snca_df$type <- "exon"
snca_df$seqnames <- 4
snca_df$strand <- "-"
snca_df$transcript_name <- snca_df$transcript
# Compute introns
snca_df_full <- to_intron(snca_df)
snca_df_join <- full_join(snca_df_full,snca_df)

snca_exons <- snca_df_join %>% dplyr::filter(type == "exon")


snca_rescaled <- shorten_gaps(
  exons = snca_df, 
  introns = to_intron(snca_df, "transcript_name"), 
  group_var = "transcript_name"
)

calledTranscripts <- str_sub(rownames(lm_results),end = -3L)
table(unique(snca_rescaled$transcript_name) %in% calledTranscripts)
snca_rescaled <- snca_rescaled[which(snca_rescaled$transcript_name %in% calledTranscripts),]
transcriptOrder <- c("ENST00000394991","ENST00000508895")
transcriptOrder <- c(transcriptOrder,calledTranscripts[!calledTranscripts %in% transcriptOrder])

snca_rescaled$transcript_name <- factor(snca_rescaled$transcript_name,levels = rev(transcriptOrder))
snca_rescaled[which(snca_rescaled$transcript_name == "ENST00000394991"),"label"] <- "ENST00000394991"
snca_rescaled[which(snca_rescaled$transcript_name == "ENST00000508895"),"label"] <- "ENST00000508895"

# shorten_gaps() returns exons and introns all in one data.frame()
# let's split these for plotting 
snca_rescaled_exons <- snca_rescaled %>% dplyr::filter(type == "exon") 
snca_rescaled_introns <- snca_rescaled %>% dplyr::filter(type == "intron") 

snca_rescaled$label

transcriptStructure <- snca_rescaled_exons %>% 
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_name
  )) +
  geom_range(
    aes(fill = as.factor(label))
  ) +
  geom_intron(
    data = snca_rescaled_introns,
    aes(strand = strand), 
    arrow.min.intron.length = 300
  )+
  theme_cowplot()+ylab(NULL)+
  theme(legend.position = "none")+
  scale_fill_manual(values = c("#E41A1C","#FC8D62"),na.value = "grey80")

#Create a final grid plot of results
hitGrid <- plot_grid(methENST00000394991,genoENST00000394991,methENST00000508895,genoENST00000508895+ylim(0,5000),align = "v")
transcriptGrid <- plot_grid(transcriptStructure,hitGrid,rel_widths = c(0.43,0.6))

pdf("C:/Users/jh1159/Dropbox/Josh Harvey PhD/Bioinformatics/MRC/Coloc/transcriptPlots.pdf",height = 5.5,width = 11)
transcriptGrid
dev.off()

