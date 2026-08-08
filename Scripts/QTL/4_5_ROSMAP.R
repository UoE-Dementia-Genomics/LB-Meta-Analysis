################# ROSMAP analysis ####
# Scripts for analysis of the ROSMAP expression data vs. genotyping for DLB risk and DNA methylation

library(stringr)
library(cowplot)
library(ggplot2)


#Load and process Expression Data
rosmapCounts <- read.table("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/ROSMAP_RNAseq_FPKM_gene/ROSMAP_RNAseq_FPKM_gene.tsv",header = T)
rownames(rosmapCounts) <- rosmapCounts$gene_id
colnames(rosmapCounts) <-  gsub("X","",colnames(rosmapCounts))
colnames(rosmapCounts) <- str_sub(colnames(rosmapCounts),end = -3)
# 60230   640

# Variant data for rs7680557
rosmapVar <- read.table("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/dlb_SNP_Extract.raw",header = T)


# Load relevant metadata for ROSMAP phenotyping and cell deconvolution data
load("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/ROSMAP_metadata.Rdata")
load("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/ROSMAP_validation.Rdata")
rosmapCell <- read.csv("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/rosmap_CIBERSORT_results.txt", header = T, row.names = 1)


# Process ROSMAP ID's for merging
table(rosmapVar$FID == rosmapVar$IID)

for(i in 1:nrow(pheno)){
  print(i)
  if(any(grepl(pheno[i,"projid"],rosmapVar$FID))){
  pheno[i,"FID"] <- rosmapVar[grep(pheno[i,"projid"],rosmapVar$FID),"FID"]
  }
}

# Extract the rs7680557 variant and annotate to the pheno file
pheno$var <- rosmapVar[match(pheno$FID,rosmapVar$FID),"X4.90763360.C.A_A"]

# Load genotyping PC's and merge
pcs <- read.table("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/Data/LB/ROSMAPpcs.eigenvec",header = T)
rownames(pcs) <- pcs$FID
pheno$PC1 <- pcs[pheno$FID,"PC1"]
pheno$PC2 <- pcs[pheno$FID,"PC2"]
pheno$PC3 <- pcs[pheno$FID,"PC3"]

# Re-process ID's for integrating expression and sequencing data
colnames(rosmapCounts) <- paste("X",colnames(rosmapCounts),sep = "")
pheno$seqID <- paste("X",pheno$specimenID,sep = "")

# merge
pheno <- cbind(pheno,rosmapCell[match(rownames(rosmapCell),pheno$seqID),])

#Process rownames for indexing
rownames(pheno) <- pheno$individualID.x
rownames(ROSMAP_pheno) <- ROSMAP_pheno$individualID

#Mean impute missing PMIs
ROSMAP_pheno$pmi[which(is.na(ROSMAP_pheno$pmi))] <- mean(na.omit(ROSMAP_pheno$pmi))
lmCor <- lm(cg01966878~Gender + Age + pmi + batch + NeuNNeg_Sox10Neg_IRF8Pos +NeuNPos_SOX6Neg +NeuNPos_SOX6Pos + NeuNNeg_SOX10Pos +  NeuNNeg_Sox10Neg_IRF8Neg, data = ROSMAP_pheno)

#Residual correct DNA methylation data
ROSMAP_pheno$cg01966878_cor <- resid(lmCor) + lmCor$coefficients[1]

pheno$cg01966878 <- ROSMAP_pheno[rownames(pheno),"cg01966878"]
pheno$cg01966878_cor <- ROSMAP_pheno[rownames(pheno),"cg01966878_cor"]
pheno <- cbind(pheno,ROSMAP_pheno[rownames(pheno),c("NeuNNeg_Sox10Neg_IRF8Pos","NeuNPos_SOX6Neg","NeuNPos_SOX6Pos",
                                                    "NeuNNeg_SOX10Pos","NeuNNeg_Sox10Neg_IRF8Neg")])

pheno[which(pheno$age_death == "90+"),"age_death"] <- 90 # Deal with age ceiling
pheno$age_death <- as.numeric(pheno$age_death) #Force to numeric



# rownames(rosmapCounts) <- rosmapCounts$gene
# rosmapCounts <- rosmapCounts[,-1]

#Merge rawSNCA expression data
pheno$SNCA_TPM <- as.numeric(rosmapCounts[grep("ENSG00000145335",rownames(rosmapCounts)),pheno$seqID])
pheno$SNCA_AS1_TPM <- as.numeric(rosmapCounts[grep("ENSG00000247775",rownames(rosmapCounts)),pheno$seqID])


# Inverse normal transformation
n <- ncol(rosmapCounts)
zvalues <- qnorm(ppoints(n))
z <- as.matrix(rosmapCounts)
for (i in 1:nrow(z)) z[i,] <- zvalues[order(order(z[i,]))]

#Extract SNCA
SNCA_z <- as.numeric(z[grep("ENSG00000145335",rownames(z)),pheno$seqID])
SNCA_AS1_z <- as.numeric(z[grep("ENSG00000247775",rownames(z)),pheno$seqID])

#Joing INT normalised data
pheno$SNCA_AS1_z <- as.numeric(SNCA_AS1_z)
pheno$SNCA_z <- as.numeric(SNCA_z)

# Factorise variants for plotting
pheno$var_factor <- ifelse(pheno$var == 0,"CC",ifelse(pheno$var == 1, "CA","AA"))

#Plot associations
png("C:/Users/jh1159/OneDrive - University of Exeter/PostDoc/LB Meta Analysis Manuscript/Figures/supplementaryRosmapUpdated.png",width = 3000, height = 2100,res = 300)
plot_grid(
  ggplot(pheno[-which(is.na(pheno$var_factor)),], aes(x = as.factor(var_factor), y = SNCA_z))+
    geom_hline(yintercept =  0)+geom_violin(fill = "#E7298A",color= "#F0027F",alpha = 0.7)+
    geom_boxplot(width = 0.2,lwd = 1,outlier.shape = NA)+ geom_jitter(size = 0.2,width = 0.2)+theme_bw()+
    labs(y = "SNCA (FPKM + INT)", x = "rs7680557 Genotype")
  ,
  ggplot(pheno, aes(x = cg01966878_cor, y = SNCA_z))+
    geom_point(color = "#E7298A")+geom_smooth(method = "lm")+theme_bw()+
    labs(y = "SNCA (FPKM + INT)", x = "cg01966878 Methylation (Corrected beta)")
  ,
  ggplot(pheno[-which(is.na(pheno$var_factor)),], aes(x = as.factor(var_factor), y = SNCA_AS1_z))+
    geom_hline(yintercept = 0) + geom_violin(fill = "#E6AB02",color= "#E6AB02",alpha = 0.7)+
    geom_boxplot(width = 0.2,lwd = 1,outlier.shape = NA)+ geom_jitter(size = 0.2,width = 0.2)+theme_bw()+
    labs(y = "SNCA-AS1 (FPKM + INT)", x = "rs7680557 Genotype")
  ,
  ggplot(pheno, aes(x =cg01966878_cor, y = SNCA_AS1_z))+
    geom_point(color = "#E6AB02")+geom_smooth(method = "lm")+theme_bw()+
    labs(y = "SNCA-AS1 (FPKM + INT)", x = "cg01966878 Methylation (Corrected beta)")
,labels = c("A","B","C","D"))
dev.off()

pheno <- pheno[-which(is.na(pheno$var_factor) & is.na(pheno$cg01966878_cor)),]
write.csv(pheno[,c("var_factor","cg01966878_cor","SNCA_z","SNCA_AS1_z")],file = "SNCA_rosmap_plottingdata.csv") # save data

#Test linear associations
lm_SNCA_AS_var <- summary(lm(SNCA_AS1_z ~ var + age_death + msex + RIN + as.factor(libraryBatch) + pmi +PC1 + PC2 + PC3 + fpkm_astrocytes + fpkm_endothelial + fpkm_microglia + fpkm_neurons + fpkm_OPC + fpkm_oligodendrocytes, data = pheno))
lm_SNCA_AS_cpg <-summary(lm(SNCA_AS1_z ~ cg01966878_cor + age_death + msex + RIN + as.factor(libraryBatch) + pmi + fpkm_astrocytes + fpkm_endothelial + fpkm_microglia + fpkm_neurons + fpkm_OPC + fpkm_oligodendrocytes, data = pheno))
lm_SNCA_var <-summary(lm(SNCA_z ~ var + age_death + msex + RIN + as.factor(libraryBatch) + pmi +PC1 + PC2 + PC3 + fpkm_astrocytes + fpkm_endothelial + fpkm_microglia + fpkm_neurons + fpkm_OPC + fpkm_oligodendrocytes , data = pheno))
lm_SNCA_cpg <-summary(lm(SNCA_z ~ cg01966878_cor + age_death + msex + RIN + as.factor(libraryBatch) + pmi  + fpkm_astrocytes + fpkm_endothelial + fpkm_microglia + fpkm_neurons + fpkm_OPC + fpkm_oligodendrocytes , data = pheno))



# Store results
results <- as.data.frame(rbind(coefficients(lm_SNCA_AS_cpg)["cg01966878_cor",],
            coefficients(lm_SNCA_AS_var)["var",],
            coefficients(lm_SNCA_cpg)["cg01966878_cor",],
            coefficients(lm_SNCA_var)["var",]))

#Label comparisons
results$comp <- c("SNCA-AS1 vs cg01966878","SNCA-AS1 vs rs7680557",
                  "SNCAvs cg01966878","SNCA vs rs7680557")

# Multiple testing correction
results$FDR <- p.adjust(results$`Pr(>|t|)`)
