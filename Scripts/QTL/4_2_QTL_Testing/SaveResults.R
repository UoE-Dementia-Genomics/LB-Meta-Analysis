setwd("/lustre/projects/Research_Project-T112069/Meth/QTL/Data/")
.libPaths("/lustre/projects/Research_Project-T112069/packages")

OutPrefix="PFC"
FactCovar="Sex"
NumCovar=c("Age","PC1","PC2","PC3","NeuNNeg_Sox10Neg_IRF8Pos","NeuNPos_SOX6Neg","NeuNPos_SOX6Pos","NeuNNeg_Sox10Neg_IRF8Neg")
Dist=1e+6
cis.pval=1
trans.pval=1e-5
trans.cross.chr=T
# SavecsvCis=yes
# SavecsvTrans=yes
  
OutDir <- "/lustre/projects/Research_Project-T112069/Meth/QTL/Data/"
OutPrefix <- "PFC"

library(stringr)

setwd("/lustre/projects/Research_Project-T112069/Meth/QTL/Data")
OutPrefix <- "PFC"
cis_pval <- 1e5
trans_pval <- 1e5
save.csv.cis <- T
save.csv.trans <- T
dist_ <- 1e6

cis_file = paste0(OutPrefix,".eQTL.Cis.Pval.",cis_pval,".Dist.",dist_,".csv")
trans_file = paste0(OutPrefix,".eQTL.Trans",trans.cross.chr,".Pval.",trans_pval,".Dist.",dist_,".csv")
merge_file = paste0(OutPrefix,".eQTL",trans.cross.chr,".TransPvalue.",trans_pval,".CisPvalue.",cis_pval,".Dist.",dist_,".RData")

cat("\n")
print("Preparing merged QTL results...")
cat("\n")

if(!file.exists(merge_file)){
  for (i in 1:22) {
    file_ <- paste0(OutPrefix,".matrixEQTL.chr",i,"/",paste0(OutPrefix,".matrixEQTL.chr",i,".RData"))
    if(file.exists(file_)){
      print(paste("Reading QTL file Chr",i))
      load(file_)
      assign(x = paste0("eQTL.chr",i),value = me)
      remove(me)
    }
    
  }
  
  variables=ls()[str_detect(ls(),pattern ="eQTL.chr" )]
  mQTLs <- mget(variables)
  print("Writing merged rdat file...")
  save(mQTLs,file = merge_file)
}else{
  print("Loading merged rdat file...")
  load(merge_file)
}

mqtl.cis <- vector(mode = "list",length = length(mQTLs))
mqtl.trans <- vector(mode = "list",length = length(mQTLs))
for (i in 1:length(mQTLs)) {
  if(save.csv.cis)
    mqtl.cis[[i]] <- mQTLs[[i]]$cis$eqtls[mQTLs[[i]]$cis$eqtls$pvalue	< cis_pval,]
  if(save.csv.trans)
    mqtl.trans[[i]] <- mQTLs[[i]]$trans$eqtls[mQTLs[[i]]$trans$eqtls$pvalue	< trans_pval,]
}

if(save.csv.cis)
  mqtl.cis_out <- do.call(rbind.data.frame,mqtl.cis)
if(save.csv.trans)
  mqtl.trans_out <- do.call(rbind.data.frame,mqtl.trans)

if(save.csv.cis)
  if(ncol(mqtl.cis_out) > 0){
    names(mqtl.cis_out)[2] <- "gene"
    mqtl.cis_out <- mqtl.cis_out[order(mqtl.cis_out$pvalue , decreasing = F),]
    print("Writing Cis csv file...")
    write.csv(mqtl.cis_out,file = cis_file,row.names = F)
  }else{
    print("There is no Cis results!")
  }

if(save.csv.trans)
  if(ncol(mqtl.trans_out) > 0){
    names(mqtl.trans_out)[2] <- "gene"
    mqtl.trans_out <- mqtl.trans_out[order(mqtl.trans_out$pvalue , decreasing = F),]
    print("Writing Trans csv file...")
    write.csv(mqtl.trans_out,file =trans_file,row.names = F)
  }else{
    print("There is no Trans results!")
  }


setwd("/lustre/projects/Research_Project-T112069/Meth/QTL")
.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(data.table)

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
allSumstat <- as.data.frame(allSumstat)
allSumstat$FDR <- p.adjust(allSumstat$Fixed_P, method = "fdr")
fdrHits <- allSumstat[which(allSumstat$FDR < 0.05),]

cisQTLold <- read.csv("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/mQTL/PFC.eQTL.Cis.Pval.1e+05.Dist.1e+06.csv",header = T)
cisQTL <- read.csv("/lustre/projects/Research_Project-T112069/Meth/QTL/Data/PFC.eQTL.Cis.Pval.1e+05.Dist.1e+06.csv",header = T)
cisQTL <- cisQTL[which(cisQTL$gene %in% rownames(fdrHits)),]
cisQTL <- cisQTL[which(cisQTL$FDR < 0.05),]

cisQTLold <- cisQTLold[which(cisQTLold$gene %in% rownames(fdrHits)),]
cisQTLold <- cisQTLold[which(cisQTLold$FDR < 0.05),]
rownames(cisQTLold) <- paste(cisQTLold$snps,cisQTLold$gene, sep = "_")
rownames(cisQTL) <- paste(cisQTL$snps,cisQTL$gene, sep = "_")

plotNames <- intersect(rownames(cisQTLold),rownames(cisQTL))
plot(cisQTLold[plotNames,"beta"],cisQTL[plotNames,"beta"])

resMetaCC <- allSumstat
epicMani <- fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv", stringsAsFactors = F,skip = 6,data.table =F)
row.names(epicMani) <- epicMani[,1]
RESULTS <- cbind(resMetaCC,epicMani[rownames(resMetaCC),])

bim <- fread("/lustre/projects/Research_Project-T112069/Meth/QTL/Data/filtered_chr_bp.bim")
bim$SNP <- paste(bim$V2,bim$V5, sep = "_")
### 

bim <- bim[!duplicated(bim$SNP),]
bim <- as.data.frame(bim)
rownames(bim) <- bim$SNP
cisQTL$CHR <- bim[cisQTL$snps,"V1"]
cisQTL$BP <- bim[cisQTL$snps,"V4"]
cisQTL$A1 <- bim[cisQTL$snps,"V6"]
cisQTL$A2 <- bim[cisQTL$snps,"V5"]


qtlSummary <- data.frame(cpg = "NA",cpgCHR = "NA",cpgBP = "NA", snpID = "NA",snpPos = "NA",pvalue = "NA",FDR = "NA",effect = "NA",nQTL = "NA",minBP = "NA",maxBP = "NA")
for(i in unique(cisQTL$gene)){
  tempQTL <- cisQTL[which(cisQTL$gene == i),]
  nQTL <- nrow(tempQTL)
  minBP <- min(tempQTL$BP,na.rm = T)
  maxBP <- max(tempQTL$BP,na.rm = T)
  
  lead <- tempQTL[which(tempQTL$pvalue == min(tempQTL$pvalue)),]
  lead <- lead[1,]
  
  cpgPos <- epicMani[i,"MAPINFO"]
  cpgchr <- epicMani[i,"CHR"]
  
  qtlSummary <- rbind(qtlSummary,c(i,cpgchr,cpgPos,lead$snps,lead$BP,lead$pvalue,lead$FDR,lead$beta,nQTL,minBP,maxBP))
}
qtlSummary <- qtlSummary[-1,]
qtlSummary <- qtlSummary[order(as.numeric(qtlSummary$pvalue)),]

