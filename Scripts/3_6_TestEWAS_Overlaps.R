.libPaths("/lustre/projects/Research_Project-T112069/packages")
library(data.table)
library(stringr)
library(ggrepel)
library(dplyr)    
library(tidyverse)  
library(tidyr)
library(stringr)
library(plyr)

epicMani <-  fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv",skip = 7)
epicMani <- as.data.frame(epicMani)
rownames(epicMani) <- epicMani$IlmnID


load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
# Process meta summary statistics for DMR analysis
RESULTS <- as.data.frame(allSumstat)
colnames(RESULTS)[c(1,3)] <- c("Effect.Mixed_LB","P.Mixed_LB")

# Load in AD data (Smith, Pishva 2021)
adMeta <- fread("/lustre/projects/Research_Project-T112069/Reference/alldataresultswithCB.csv")
adMeta <- as.data.frame(adMeta)
rownames(adMeta) <- adMeta[,1]
adMeta <- adMeta[order(adMeta$meta.p),]

overlap <- intersect(rownames(adMeta),rownames(RESULTS))
RESULTS[overlap,"P.NFT"] <- adMeta[overlap,"meta.p"]
RESULTS[overlap,"Effect.NFT"] <- adMeta[overlap,"meta.es"]

# Load in HD regional specific (Wheildon 2025)
load("/lustre/projects/Research_Project-T112069/Meth/reference/EC_PAPER_model2_bacon_EPIC.RData")
load("/lustre/projects/Research_Project-T112069/Meth/reference/STR_PAPER_model2_bacon_EPIC.RData")
overlap <- intersect(rownames(modelEC),rownames(RESULTS))
RESULTS[overlap,"P.HD_EC"] <- modelEC[overlap,"pval.bc"]
RESULTS[overlap,"Effect.HD_EC"] <- modelEC[overlap,"es.b"]
overlap <- intersect(rownames(modelSTR),rownames(RESULTS))
RESULTS[overlap,"P.HD_STR"] <- modelSTR[overlap,"pval.bc"]
RESULTS[overlap,"Effect.HD_STR"] <- modelSTR[overlap,"es.b"]


# Generate all combinations of 2 elements
combinations <- t(combn(c("NFT", "Mixed_LB", "HD_EC", "HD_STR"), 2))

# Convert to a dataframe
combinations_df <- as.data.frame(combinations, stringsAsFactors = FALSE)

# Assign meaningful column names
colnames(combinations_df) <- c("Var1", "Var2")


# Create one dataframe for testing
overlapTests <- rbind(combinations_df,combinations_df,combinations_df,combinations_df)
overlapTests$Discovery <- c(rep("NFT",6),rep("Mixed_LB",6),rep("HD_EC",6),rep("HD_STR",6))
overlapTests$Threshold <- 1e-5
overlapTests[which(overlapTests$Discovery == "NFT"),"Threshold"] <- 1.238e-7
overlapTests[which(overlapTests$Discovery %in% c("HD_EC")),"Threshold"] <- 1e-4



for(x in 1:nrow(overlapTests)){
  print(x)
  discID <- overlapTests[x,"Discovery"]
  discLoci <- rownames(RESULTS[which(RESULTS[,paste("P",discID,sep = ".")] < overlapTests[x,"Threshold"]),])
  overlapTests[x,"nLoci"] <- length(discLoci)
  study1 <- paste("Effect",overlapTests[x,"Var1"],sep = ".")
  study2 <- paste("Effect",overlapTests[x,"Var2"],sep = ".")
  overlapTests[x,"overLoci"]<- length(which(sign(RESULTS[discLoci,study1]) == sign(RESULTS[discLoci,study2])))
  naTest <- length(which(is.na(sign(RESULTS[discLoci,study1]) == sign(RESULTS[discLoci,study2]))))
  overlapTests[x,"testedLoci"]<- overlapTests[x,"nLoci"] - naTest
  binomTest <- binom.test(overlapTests[x,"overLoci"],overlapTests[x,"testedLoci"],alternative = "greater")
  overlapTests[x,"Prob"] <- binomTest$estimate
  overlapTests[x,"P.val"] <- binomTest$p.value
  
  corTest <- cor.test(RESULTS[discLoci,study1],RESULTS[discLoci,study2],use="pairwise.complete.obs")
  overlapTests[x,"P.valCorr"] <- corTest$p.value
  overlapTests[x,"Cor"] <- corTest$estimate
}


saveInter <- overlapTests  
overlapTests <- overlapTests[which(overlapTests$Discover %in% c("NFT","Mixed_LB")),]

patterns <- c("NFT", "Mixed_LB","HD_EC","HD_STR")
replacements <- c("Braak-NFT", "Braak-LB","HD(Entorhinal Cortex)","HD(Striatum)")

# Replace the patterns in the entire dataframe
overlapTests <- as.data.frame(lapply(overlapTests, function(column) {
  if (is.character(column)) {
    for (i in seq_along(patterns)) {
      column <- gsub(patterns[i], replacements[i], column)
    }
  }
  column
}), stringsAsFactors = FALSE)

overlapTests$Var1 <- factor(overlapTests$Var1,levels = c("Braak-NFT", "Braak-LB","HD(Entorhinal Cortex)","HD(Striatum)"))
overlapTests$Var2 <- factor(overlapTests$Var2,levels = c("Braak-NFT", "Braak-LB","HD(Entorhinal Cortex)","HD(Striatum)"))
overlapTests$P.adjust <- p.adjust(overlapTests$P.val,method = "fdr")
overlapTests$Pcor.adjust <- p.adjust(overlapTests$P.valCorr,method = "fdr")

library(cowplot)
overlapHeatmaps <- ggplot(overlapTests, aes(x = Var1,y = Var2,fill = Prob, label = signif(Prob, 2)))+
  geom_tile(color = "black",size = .5)+
  geom_text(data = overlapTests[which(overlapTests$P.adjust < 0.05),],color = "black",fontface = "bold",size = 7)+
  geom_text(data = overlapTests[-which(overlapTests$P.adjust < 0.05),],color = "grey50",size = 5)+
  facet_wrap(~Discovery)+
  scale_fill_viridis_c(option = "A")+
  theme_classic()+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_color_manual(values = c("grey","black"))+
  theme(legend.position = "right",axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.background = element_rect(fill = NA), strip.text = element_text(face = "bold"))+
  labs(x = NULL, y = NULL, fill = "Overlap\nProbability")

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/OverlapHeatmapRerun.pdf",width = 8, height = 4.5)
overlapHeatmaps
dev.off()




write.csv(RESULTS[which(RESULTS$P.NFT < 1.238e-7),],file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/overlapADloci.csv")
write.csv(RESULTS[which(RESULTS$P.Mixed_LB < 1e-5),],file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/overlapLBloci.csv")



library(ggplot2)
discAD_testMIX <- ggplot(data = RESULTS[which(RESULTS$P.NFT < 1.238e-7),],aes(x = Effect.Mixed_LB,y = Effect.NFT))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  ylab("NFT Effect Estimate")+
  xlab("LB Full Cohort Estimate")+
  ylim(-0.005,0.005)+
  xlim(-0.0025,0.0025)+
  labs(color = "Nominal LB association (P < 1e-5)",
       shape = "Nominal LB association (P < 1e-5)")+
  theme(legend.position = "bottom")




# RESULTS <- RESULTS[order(RESULTS$P_Mixed_LB),]
discMIX_testAD <- ggplot(data = RESULTS[which(RESULTS$P.Mixed_LB < 1e-5),],aes(x = Effect.Mixed_LB,y = Effect.NFT))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  ylab("NFT Effect Estimate")+
  xlab("LB Full Cohort Estimate")+
  ylim(-0.002,0.002)+
  xlim(-0.007,0.007)

discMIX_testHD_EC <- ggplot(data = RESULTS[which(RESULTS$P.Mixed_LB < 1e-5),],aes(x = Effect.Mixed_LB,y = Effect.HD_EC))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  ylab("HD Effect Estimate")+
  xlab("LB Full Cohort Estimate")+
  ylim(-0.075,0.12)+
  xlim(-0.008,0.008)

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/OverlapScatters.pdf",width = 10, height = 2.7)
plot_grid(discMIX_testAD,discMIX_testHD_EC,discAD_testMIX,nrow = 1)
dev.off()

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/OverlapHeatmap_Scatters.pdf",width = 9, height = 7.5)
plot_grid(overlapHeatmaps,
          plot_grid(discMIX_testAD,discMIX_testHD_EC,discAD_testMIX,nrow = 1,labels = c("B)","C)","D)")),ncol = 1,rel_heights = c(0.65,0.35),labels = c("A)",""))
dev.off()


