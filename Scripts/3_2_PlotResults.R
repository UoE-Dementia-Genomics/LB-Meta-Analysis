setwd("/lustre/projects/Research_Project-T112069/")
# Script to plot qq plots, Manhattans and relevant comparisons of effects across analyses


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
library(sva)
library(reshape)
library(meta)
library(cowplot)
library(bacon)
library(data.table)
library(stringr)
library(ggrepel)



#########################
##### Data loading ######
#########################
# Load and merge cohort specific and meta analysis summary statistics
load("res_UKBBN_CrossCortex.Rdata")
load("res_UKBBN_cng.Rdata")
load("res_UKBBN_pfc.Rdata")
load("res_UKBBN_CrossCortex_Filtered.Rdata")
load("res_UKBBN_CrossCortex_Thal_SV3.Rdata")
res_UKBBN <- res_UKBBN[,-3]
colnames(res_UKBBN) <- c("Values","SE","T","P")
load("NBB_knownCovar_SV3.Rdata")
load("NBB_knownCovar_SV3_filtered.Rdata")
load("NBB_knownCovar_SV3_Thal.Rdata")
load("BDR_Res_SV1.Rdata")
load("BDR_Res_Thal_SV1.Rdata")


resUK <- res_UKBBN
resNBB <- res_NBB_SV3
resBDR <- BDR_res_LB_SV1

overlap <-Reduce(intersect, list(rownames(resUK),rownames(resNBB),rownames(resBDR)))
resUK <- resUK[overlap,]
resNBB <- resNBB[overlap,]
resBDR <- resBDR[overlap,]


resAll <- cbind(resUK, resNBB,resBDR)
resAll <- as.data.frame(resAll)

resUK_thal <- UKBBN_res_Thal_SV3
resUK_thal <- resUK_thal[overlap,-3]
resNBB_thal <- NBB_res_SV3_Thal[overlap,]
resBDR_thal <- BDR_res_LB_Thal_SV1[overlap,]

resUK_pure <- res_UKBBN_Filtered
resUK_pure <- resUK_pure[overlap,-3]
resNBB_pure <- NBB_res_SV3_filtered[overlap,]


resUK_pfc <- res_UKBBN_pfc[overlap,]
resUK_acc <- res_UKBBN_cng[overlap,]

colnames(resUK) <- c("Beta","SE","T","P")
colnames(resNBB) <- c("Beta","SE","T","P")
colnames(resBDR) <- c("Beta","SE","T","P")
colnames(resUK_pure) <- c("Beta","SE","T","P")
colnames(resNBB_pure) <- c("Beta","SE","T","P")
colnames(resUK_thal) <- c("Beta","SE","T","P")
colnames(resNBB_thal) <- c("Beta","SE","T","P")
colnames(resBDR_thal) <- c("Beta","SE","T","P")




load("MultiMeta_FullCohort.Rdata")
resMetaCC <- allSumstat[overlap,]
load("MultiMeta_PureLB.Rdata")
lowSumstat <- lowSumstat[overlap,]
load("MultiMeta_ThalControlled.Rdata")
thalSumstat <- thalSumstat[overlap,]
load("PFC_Meta_FullCohort.Rdata")
pfcSumstat <- pfcSumstat[overlap,]

#########################
##### QQ plots ##########
#########################
# Calculate Lambdas for each analysis
LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  return(lambda)
}

colnames(resMetaCC)[colnames(resMetaCC) == "Fixed_P"] <- "P"
colnames(lowSumstat)[colnames(lowSumstat) == "Fixed_P"] <- "P"
colnames(thalSumstat)[colnames(thalSumstat) == "Fixed_P"] <- "P"
colnames(pfcSumstat)[colnames(pfcSumstat) == "Fixed_P"] <- "P"


fullCohort <- list(Meta = as.data.frame(resMetaCC),
                   UK = as.data.frame(resUK), 
                   NBB = as.data.frame(resNBB), 
                   BDR = as.data.frame(resBDR))

pureCohort <- list(Meta = as.data.frame(lowSumstat),
                   UK = as.data.frame(resUK_pure), 
                   NBB = as.data.frame(resNBB_pure))

thalCohort <- list(Meta = as.data.frame(thalSumstat),
                   UK = as.data.frame(resUK_thal), 
                   NBB = as.data.frame(resNBB_thal),
                   BDR = as.data.frame(resBDR_thal))

pfcCohort <- list(Meta = as.data.frame(pfcSumstat),
                  UK_PFC = as.data.frame(resUK_pfc), 
                  UK_ACC = as.data.frame(resUK_acc))

process_pval_list <- function(datasets) {
  
  tSet.list <- lapply(names(datasets), function(nm) {
    
    df <- datasets[[nm]]
    
    # order by P (smallest → largest)
    df <- df[order(df$P), ]
    
    # observed p-values
    observed <- df$P
    
    # expected p-values
    n <- length(observed)
    expected <- seq(1/n, 1, by = 1/n)
    
    # dataset label
    dataset <- paste(
      nm,
      signif(LambdaInf(na.omit(observed)), 4)
    )
    
    # return ONLY required columns
    data.frame(
      dataset = dataset,
      observed.pvals = observed,
      expected.pvals = expected
    )
  })
  
  # bind all together
  tSet <- do.call(rbind, tSet.list)
  
  return(tSet)
}

tSet <- process_pval_list(fullCohort)
tSet_thal <- process_pval_list(thalCohort)
tSet_pure <- process_pval_list(pureCohort)
tSet_pfc <- process_pval_list(pfcCohort)

fullCohortPlot <- ggplot(tSet, aes(y = -log10(observed.pvals), x = -log10(expected.pvals), color = factor(dataset, levels = unique(tSet$dataset))))+
  geom_point()+
  geom_abline(slope = 1,linetype = "dashed")+
  scale_color_brewer(palette = "Dark2")+
  theme_cowplot()+
  scale_x_continuous(breaks = c(0:7),expand = expansion(mult = c(0, .1))) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2),expand = expansion(mult = c(0, .1))) +
  theme(legend.position = c(0.1, 0.8))+
  labs(color = "")+
  xlab("-log10(Expected p-values)")+
  ylab("-log10(Observed p-values)")

thalCohortPlot <- ggplot(tSet_thal, aes(y = -log10(observed.pvals), x = -log10(expected.pvals), color = factor(dataset, levels = unique(tSet_thal$dataset))))+
  geom_point()+
  geom_abline(slope = 1,linetype = "dashed")+
  scale_color_brewer(palette = "Dark2")+
  theme_cowplot()+
  scale_x_continuous(breaks = c(0:7),expand = expansion(mult = c(0, .1))) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2),expand = expansion(mult = c(0, .1))) +
  theme(legend.position = c(0.1, 0.8))+
  labs(color = "")+
  xlab("-log10(Expected p-values)")+
  ylab("-log10(Observed p-values)")

pureCohortPlot <- ggplot(tSet_pure, aes(y = -log10(observed.pvals), x = -log10(expected.pvals), color = factor(dataset, levels = unique(tSet_pure$dataset))))+
  geom_point()+
  geom_abline(slope = 1,linetype = "dashed")+
  scale_color_brewer(palette = "Dark2")+
  theme_cowplot()+
  scale_x_continuous(breaks = c(0:7),expand = expansion(mult = c(0, .1))) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2),expand = expansion(mult = c(0, .1))) +
  theme(legend.position = c(0.1, 0.8))+
  labs(color = "")+
  xlab("-log10(Expected p-values)")+
  ylab("-log10(Observed p-values)")

pfcCohortPlot <- ggplot(tSet_pfc, aes(y = -log10(observed.pvals), x = -log10(expected.pvals), color = factor(dataset, levels = unique(tSet_pfc$dataset))))+
  geom_point()+
  geom_abline(slope = 1,linetype = "dashed")+
  scale_color_brewer(palette = "Dark2")+
  theme_cowplot()+
  scale_x_continuous(breaks = c(0:7),expand = expansion(mult = c(0, .1))) +
  scale_y_continuous(breaks = seq(from = 0, to = 12, by = 2),expand = expansion(mult = c(0, .1))) +
  theme(legend.position = c(0.1, 0.8))+
  labs(color = "")+
  xlab("-log10(Expected p-values)")+
  ylab("-log10(Observed p-values)")

png("AllQQ.png",width = 1600, height = 1600,res = 200)
plot_grid(fullCohortPlot,pureCohortPlot,thalCohortPlot,pfcCohortPlot,ncol = 2,labels = c("A)","B)","C)","D)"))
dev.off()


################################
##### Manhattan plots ##########
################################
# Load manifest data and merge
epicMani <-  fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv",skip = 7)
epicMani <- as.data.frame(epicMani)
rownames(epicMani) <- epicMani$IlmnID
idx <- match(rownames(resMetaCC), rownames(epicMani))
RESULTS <- cbind(
  resMetaCC,
  epicMani[idx, , drop = FALSE]
)

sigRes <- RESULTS[which(RESULTS$Fixed_P <= 1e-5),]
sigRes <- cbind(sigRes,resAll[rownames(sigRes),])
sigRes <- sigRes[order(sigRes$Fixed_P),]

colnames(RESULTS)[which(colnames(RESULTS)== "Fixed_P")] <- "P"
colnames(RESULTS)[which(colnames(RESULTS)== "MAPINFO")] <- "BP"

testDF <- RESULTS[order(RESULTS$P),]
n <- length(which(p.adjust(testDF$P,method = "fdr") < 0.05))

# The below loop takes the top results table (testDF) and parses the UCSC refgene names to remove redundant gene annotations
# It does this for the number of FDR significant results
for(i in c(1:n)){
  cpgIndex <- rownames(testDF)[i]
  test1 <- as.character(testDF[cpgIndex,"UCSC_RefGene_Name"])
  
  if(str_detect(test1,";") == TRUE){
    strIndex <- as.data.frame(str_locate_all(test1,";"))$end
    length(strIndex)
    strIndex <- c(0,strIndex,str_length(test1))
    storage <- c()
    for(y in strIndex){
      if(which(strIndex == y) == 2){
        storage <- append(storage,str_sub(test1,
                                          start = strIndex[which(strIndex == y) -1],
                                          end = y - 1))
      }else if(y >0 & y != max(strIndex)){
        storage <-append(storage,str_sub(test1,
                                         start = strIndex[which(strIndex == y) -1]+ 1,
                                         end = y - 1))
      }else if(y == max(strIndex)){
        storage <-append(storage,str_sub(test1,
                                         start = strIndex[which(strIndex == y) -1]+ 1))
      } 
    }
    sumGene <- as.character(unique(storage))
    
    if(length(sumGene) == 1){
      testDF[cpgIndex,"UCSC1"] <- sumGene
    }else if(length(sumGene) == 2){
      testDF[cpgIndex,c("UCSC1","UCSC2")] <- sumGene[1:2]
    }else if(length(sumGene) > 2){
      testDF[cpgIndex,c("UCSC1","UCSC2","UCSC3")] <- sumGene[1:3]
    }
  } else if(str_detect(test1,";") == FALSE){
    testDF[cpgIndex,"UCSC1"] <- as.character(testDF[cpgIndex,"UCSC_RefGene_Name"])
  } 
}

#For samples with no UCSC annotation, use cpg identifier
testDF[which(testDF$UCSC1 == ""),"UCSC1"] <- rownames(testDF[which(testDF$UCSC1 == ""),])
testDF[which(is.na(testDF$UCSC1)),"UCSC1"] <- rownames(testDF[which(is.na(testDF$UCSC1)),])

# Reduce down to the relevant bits
RESULTS <- testDF[,c("CHR","BP","P","UCSC1")]
RESULTS$CHR <- as.numeric(RESULTS$CHR)
RESULTS$BP <- as.numeric(RESULTS$BP)
RESULTS$SNP <- RESULTS$UCSC1
RESULTS <- na.omit(RESULTS)

# Set FDR threshold for plotting
fdrThreshold <- max(RESULTS[which(p.adjust(RESULTS$P,method = "fdr") < 0.05),"P"])

#Transform results so they are plottable via ggplot, annotate genome wide + FDR significant loci
don <- RESULTS %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(RESULTS, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(-log10(P)>-log10(9e-8), "yes", "no")) %>%
  mutate( is_annotate=ifelse(-log10(P)>=-log10(fdrThreshold) , "yes", "no"))


# Transform
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Plot manhattan (this takes some tweaking to get the labels in the ideal location)
manhattan <- ggplot(don, aes(x=BPcum, y=-log10(P)))+
  geom_hline(yintercept = -log10(fdrThreshold),linetype = "dashed",color = "grey60")+
  geom_hline(yintercept = -log10(9e-8))+
  # Show all points
  geom_point( aes(color=as.factor(CHR)),  size=1.3) +
  scale_color_manual(values = rep(c("#A1D99B", "#006D2C"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ,expand = c(0.01,0)) +    # remove space between plot area and x axis
  # Custom the theme:
  theme_cowplot() +
  scale_y_continuous(expand = c(0, .1),limits = c(0,13) ) + 
  geom_text_repel(data=subset(don, is_highlight =="yes"), aes(label=SNP), size=3,
                  nudge_y       =  11 - -log10(subset(don, is_highlight =="yes")$P),
                  angle        = 90,
                  segment.size  = 0.2,
                  segment.color = "black",
                  direction     = "x",
                  hjust         = 0
  )+geom_text_repel(data=don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(1)),], aes(label=SNP), size=3,
                    nudge_y       =  9 - -log10(don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(1)),"P"]),
                    segment.linetype = "dotted",
                    nudge_x = 259903284 - don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(1)),"BPcum"],
                    angle        = 90,
                    segment.size  = 0.2,
                    segment.color = "grey60",
                    color = "grey40",
                    direction     = "x",
                    hjust         = 0,
                    vjust = 0,
                    segment.angle     = 120,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = TRUE,
                    segment.curvature = -1e-20
  )+geom_text_repel(data=don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(2:10)),], aes(label=SNP), size=3,
                    nudge_y       =  7.5 - -log10(don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(2:10)),"P"]),
                    segment.linetype = "dotted",
                    angle        = 90,
                    segment.size  = 0.2,
                    segment.color = "grey60",
                    color = "grey40",
                    direction     = "x",
                    hjust         = 0,
                    vjust = 0,
                    segment.angle     = 120,
                    segment.ncp       = 1,
                    segment.square    = TRUE,
                    segment.inflect   = TRUE,
                    segment.curvature = -1e-20
  )+
  geom_text_repel(data=don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(11)),], aes(label=SNP), size=3,
                  nudge_y       =  10 - -log10(don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(11)),"P"]),
                  nudge_x = 1816726044 - don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(11)),"BPcum"],
                  segment.linetype = "dotted",
                  angle        = 90,
                  segment.size  = 0.2,
                  segment.color = "grey60",
                  color = "grey40",
                  direction     = "x",
                  hjust         = 0,
                  vjust = 0,
                  segment.angle     = 120,
                  segment.ncp       = 1,
                  segment.square    = TRUE,
                  segment.inflect   = TRUE,
                  segment.curvature = -1e-20
  )+
  geom_text_repel(data=don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(12:22)),], aes(label=SNP), size=3,
                  nudge_y       =  7.5 - -log10(don[which(don$is_annotate == "yes" & don$is_highlight == "no"& don$CHR %in% c(12:22)),"P"]),
                  segment.linetype = "dotted",
                  angle        = 90,
                  segment.size  = 0.2,
                  segment.color = "grey60",
                  color = "grey40",
                  direction     = "x",
                  hjust         = 0,
                  vjust = 0,
                  segment.angle     = 120,
                  segment.ncp       = 1,
                  segment.square    = TRUE,
                  segment.inflect   = TRUE,
                  segment.curvature = -1e-20
  )+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    axis.text.x = element_text(color = "grey20", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain"),
    axis.text.y = element_text(color = "grey20", size = 8, hjust = .5, vjust = .5, face = "plain"),
    axis.title.x = element_text(color = "grey20", size = 10, face = "plain"),
    axis.title.y = element_text(color = "grey20", size = 10,angle = 90, face = "plain")
  )+xlab("Chromosome")


png("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Plots/Manhattan_Meta.png",width = 1700,height = 900,res = 200,type = "cairo-png")
manhattan
dev.off()


###################################
###### Plot effect size comp ######
###################################

load("MultiMeta_FullCohort.Rdata")
highNFTsub <- as.data.frame(allSumstat)
colnames(highNFTsub)[which(colnames(highNFTsub)%in% c("Fixed_Effect","Fixed_SE","Fixed_P"))] <- c("Effect","SE","P")
highNFTsub <- highNFTsub[which(p.adjust(highNFTsub$P,method = "fdr") < 0.05),]
load("MultiMeta_ThalControlled.Rdata")
thalSumstat <- as.data.frame(thalSumstat)


thalSig <- thalSumstat[which(thalSumstat$Fixed_P < 1e-5),]
idx <- match(rownames(thalSig), rownames(epicMani))
thalSig <- cbind(
  thalSig,
  epicMani[idx, , drop = FALSE]
)

# write.csv(thalSig, file = "Supplementary4_ThalMetaAnalysisResults.csv")

colnames(thalSumstat)[which(colnames(thalSumstat)%in% c("Fixed_Effect","Fixed_SE","Fixed_P"))] <- c("Effect_Sensitivity","SE_Sensitivity","P_Sensitivity")

load("MultiMeta_PureLB.Rdata")
lowSumstat <- as.data.frame(lowSumstat)

lowSig <- lowSumstat[which(lowSumstat$Fixed_P < 1e-5),]
idx <- match(rownames(lowSig), rownames(epicMani))
lowSig <- cbind(
  lowSig,
  epicMani[idx, , drop = FALSE]
)

# write.csv(lowSig, file = "Supplementary3_LowMetaAnalysisResults.csv")

colnames(lowSumstat)[which(colnames(lowSumstat)%in% c("Fixed_Effect","Fixed_SE","Fixed_P"))] <- c("Effect_Sensitivity","SE_Sensitivity","P_Sensitivity")
plotComp1 <- cbind(highNFTsub[,c("Effect","SE","P")],thalSumstat[rownames(highNFTsub),c("Effect_Sensitivity","SE_Sensitivity","P_Sensitivity")])
plotComp2 <- cbind(highNFTsub[,c("Effect","SE","P")],lowSumstat[rownames(highNFTsub),c("Effect_Sensitivity","SE_Sensitivity","P_Sensitivity")])
plotComp1$Analysis <- "Thal Controlled"
plotComp2$Analysis <- "Pure LB"
plotComp <- rbind(plotComp1,plotComp2)

write.csv(plotComp, file = "Supplementary5_OverlapSensitivity.csv")

effectComp <- ggplot(plotComp, aes(x = Effect, y = Effect_Sensitivity)) +
  # Error bars for Effect
  geom_errorbar(aes(ymin = Effect_Sensitivity - SE_Sensitivity, ymax = Effect_Sensitivity + SE_Sensitivity), width = 0, color = "gray50") +
  # Error bars for Effect_Low
  geom_errorbarh(aes(xmin = Effect - SE, xmax = Effect + SE), height = 0, color = "gray50") +
  # Points with color and size mappings
  geom_point(aes(fill = -log10(P_Sensitivity), size = -log10(P)), shape = 21, color = "black") +
  scale_fill_viridis_c() + # Color mapping for fill
  scale_size(range = c(2, 6)) + # Adjust point size range
  theme_bw() +
  labs(x = "Primary Meta Estimate (+/- SE)", y = "Sensitivity Meta Estimate (+/- SE)", fill = "Sensitivity \n P-value",size = "Primary \n P-value") +
  theme(legend.position = "left", 
        legend.box = "bottom",strip.text.x = element_text(size = 12))+
  ylim(-0.0081,0.0051)+
  xlim(-0.0081,0.0051)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  facet_grid(~Analysis)


load("PFC_Meta_FullCohort.Rdata")
colnames(pfcSumstat) <- paste("pfc_",colnames(pfcSumstat),sep = "")
pfcSummary <- cbind(pfcSumstat[rownames(highNFTsub),],highNFTsub)

colnames(res_UKBBN_pfc) <- paste("pfc_",colnames(res_UKBBN_pfc),sep = "")
colnames(res_UKBBN_cng) <- paste("acc_",colnames(res_UKBBN_cng),sep = "")
crossBrainSummary <- cbind(res_UKBBN_pfc[rownames(highNFTsub),],res_UKBBN_cng[rownames(highNFTsub),])

png("pfc_EffectSizeComp.png",width = 1200, height = 1800,res = 220)
plot_grid(
    ggplot(pfcSummary, aes(x = Effect, y = pfc_Fixed_Effect )) +
        # Error bars for Effect
        geom_errorbar(aes(ymin = pfc_Fixed_Effect  - pfc_Fixed_SE, ymax = pfc_Fixed_Effect + pfc_Fixed_SE), width = 0, color = "gray50") +
        # Error bars for Effect_Low
        geom_errorbarh(aes(xmin = Effect - SE, xmax = Effect + SE), height = 0, color = "gray50") +
        # Points with color and size mappings
        geom_point(aes(fill = -log10(pfc_Fixed_P), size = -log10(P)), shape = 21, color = "black") +
        scale_fill_viridis_c() + # Color mapping for fill
        scale_size(range = c(2, 6)) + # Adjust point size range
        theme_bw() +
        labs(x = "Primary Meta Estimate (+/- SE)", y = "PFC Specific Estimate (+/- SE)", fill = "PFC Specific \n P-value",size = "Primary \n P-value") +
        theme(legend.position = "right", 
              legend.box = "bottom",strip.text.x = element_text(size = 12))+
        ylim(-0.0082,0.0051)+
        xlim(-0.0082,0.0051)+
        geom_hline(yintercept = 0)+
        geom_vline(xintercept = 0)
        ,
    ggplot(crossBrainSummary, aes(x = pfc_Effect, y = acc_Effect)) +
      # Error bars for Effect
      geom_errorbar(aes(ymin = acc_Effect  - acc_SE, ymax = acc_Effect + acc_SE), width = 0, color = "gray50") +
      # Error bars for Effect_Low
      geom_errorbarh(aes(xmin = pfc_Effect - pfc_SE, xmax = pfc_Effect + pfc_SE), height = 0, color = "gray50") +
      # Points with color and size mappings
      geom_point(aes(fill = -log10(pfc_P), size = -log10(acc_P)), shape = 21, color = "black") +
      scale_fill_viridis_c() + # Color mapping for fill
      scale_size(range = c(2, 6)) + # Adjust point size range
      theme_bw() +
      labs(x = "UKBBN PFC Estimate (+/- SE)", y = "UKBBN ACC Estimate (+/- SE)", fill = "UKBBN PFC \n P-value",size = "UKBBN ACC \n P-value") +
      theme(legend.position = "right", 
            legend.box = "bottom",strip.text.x = element_text(size = 12))+
      ylim(-0.0082,0.0051)+
      xlim(-0.0082,0.0051)+
      geom_hline(yintercept = 0)+
      geom_vline(xintercept = 0),ncol = 1, labels = c("A)","B)"))
dev.off()

getwd()
write.csv(crossBrainSummary,file = "crossBrainEffectsSize.csv")
write.csv(pfcSummary,file = "pfcCompEffectsSize.csv")
###################################
###### Plot Forests ######
###################################
# Forest Plots
colnames(resAll) <- c("UKBBN_Effect", "UKBBN_SE"   , "UKBBN_T"   ,  "UKBBN_P"   ,  
                      "NBB_Effect" ,"NBB_SE","NBB_T"  ,  "NBB_P"  ,  
                      "BDR_Effect" ,"BDR_SE" ,  "BDR_T",    "BDR_P")


# Function to perform meta-analysis for each cohort (consistent with previous analysis)
meta_results <- data.frame()
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
resMetaCC <- allSumstat[overlap,]
resMetaCC <- as.data.frame(resMetaCC)
resMetaCC <- resMetaCC[order(resMetaCC$Fixed_P),]
#Subset to top 3 reaching stringent genome wide threshold
resAll <- resAll[rownames(resMetaCC[which(resMetaCC$Fixed_P < 9e-8),]),]

resAll$Variable <- rownames(resAll)

plot_data <- resAll %>%
  dplyr::select(Variable, UKBBN_Effect, UKBBN_SE, NBB_Effect, NBB_SE, BDR_Effect, BDR_SE) %>%
  pivot_longer(cols = -Variable, names_to = c("Cohort", ".value"), names_pattern = "(.*)_(Effect|SE)")


# Loop through each row (variable)
# Calculate relevant meta statistics for forest plots
for (i in 1:nrow(resAll)) {
  
  # Extract relevant values for the meta-analysis
  effects <- c(resAll$UKBBN_Effect[i], resAll$NBB_Effect[i], resAll$BDR_Effect[i])
  ses <- c(resAll$UKBBN_SE[i], resAll$NBB_SE[i], resAll$BDR_SE[i])
  labels <- c("UKBBN", "NBB", "BDR")
  
  # Run fixed-effect meta-analysis
  meta_res <- metagen(TE = effects, seTE = ses, studlab = labels)
  
  
  # Store results
  meta_results <- rbind(meta_results, data.frame(
    Variable = rownames(resAll)[i],
    FE_Effect = meta_res$TE.fixed,         # Fixed-effect estimate
    FE_SE = meta_res$seTE.fixed,           # Standard error
    FE_P = meta_res$pval.fixed,            # Fixed-effect P-value
    I2 = meta_res$I2,                      # Heterogeneity I²
    Q = meta_res$Q,                        # Cochran's Q
    p_Het = meta_res$pval.Q,               # Heterogeneity p-value
    UKBBN_Weight = meta_res$w.fixed[1] / sum(meta_res$w.fixed),   # Weight for UKBBN
    NBB_Weight = meta_res$w.fixed[2]/ sum(meta_res$w.fixed),     # Weight for NBB
    BDR_Weight = meta_res$w.fixed[3]/ sum(meta_res$w.fixed)
  )
  )
  plot_data[which(plot_data$Variable == rownames(resAll)[i] & plot_data$Cohort == "UKBBN"),"Weight"] <- meta_results$UKBBN_Weight[i]
  plot_data[which(plot_data$Variable == rownames(resAll)[i] & plot_data$Cohort == "NBB"),"Weight"] <- meta_results$NBB_Weight[i]
  plot_data[which(plot_data$Variable == rownames(resAll)[i] & plot_data$Cohort == "BDR"),"Weight"] <- meta_results$BDR_Weight[i]
  
}

#add collumns for heterogeneity statistics
plot_data[,c("I2","Q","p_Het")] <- NA

colnames(meta_results)[2:3] <- c("Effect","SE")
meta_results$Cohort <- "Fixed Effect"
meta_results$Weight <- NA
plot_data <- rbind(plot_data, meta_results[,colnames(plot_data)])
colnames(plot_data)[which(colnames(plot_data) == "I2")] <- "Isquare"

# Function to plot one variable at a time as a forest
plot_forest <- function(data, variable_name,minEf,maxEf) {
  # Filter data for the selected variable
  subplot_data <- data %>% filter(Variable == variable_name)
  
  # Extract heterogeneity stats for the Fixed Effect row
  I = round(subplot_data[which(subplot_data$Cohort == "Fixed Effect"),"Isquare"], 2)
  Q = round(subplot_data[which(subplot_data$Cohort == "Fixed Effect"),"Q"], 2)
  p_Het =  round(subplot_data[which(subplot_data$Cohort == "Fixed Effect"),"p_Het"], 2)
  
  
  heterogeneity_text <- paste("Heterogeneity: I² =", I,",p-value =",p_Het)
  
  # Create the forest plot
  ggplot(subplot_data, aes(x = Effect, y = factor(Cohort,levels = c("Fixed Effect","BDR","NBB","UKBBN")))) +
    geom_point(aes(size = Weight), shape = 22, fill = "grey", color = "grey", na.rm = TRUE) +
    geom_errorbarh(aes(xmin = Effect - 1.96*SE, xmax = Effect + 1.96*SE), height = 0, color = "black") +
    geom_point(data = subplot_data %>% filter(Cohort == "Fixed Effect"),
               aes(x = Effect, y = Cohort),
               shape = 18, size = 5, fill = "black") +
    geom_errorbarh(data = subplot_data %>% filter(Cohort == "Fixed Effect"),
                   aes(xmin = Effect - 1.96*SE, xmax = Effect + 1.96*SE, y = Cohort),
                   height = 0, color = "black") +
    labs(
      y= paste(variable_name),
      x = "Effect Estimate"
    ) +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dashed",color = "grey40")+
    theme(
      legend.position = "none",
      legend.box = "horizontal",
      plot.caption = element_text(hjust = 0),
      axis.ticks.x=element_blank(), 
      axis.ticks.y=element_blank() ,
      axis.line.y =element_blank()
    ) +
    annotate("text", x =-0.0059,
             y = "Fixed Effect", label = heterogeneity_text, hjust = 0, size = 3,vjust = 1.5)+
    xlim(minEf,maxEf)
}



plot1 <-   plot_forest(plot_data, "cg13847853",-0.006,0.008)
plot2 <-   plot_forest(plot_data, "cg08577553",-0.006,0.008)
plot2 <- plot2+ylab("cg08577553 \n (PTAFR)")
plot3 <-   plot_forest(plot_data, "cg27481153",-0.006,0.008)
plot3 <- plot3+ylab("cg27481153 \n (UBASH3B)")
forests <- plot_grid(plot1+xlab(NULL),plot3+xlab(NULL),plot2,ncol = 1,rel_heights = c(1,1,1.1),align = "v")
library(ggtranscript)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)

gtf_path <- file.path(tempdir(), "Homo_sapiens.GRCh37.87.chr.gtf.gz")

download.file(
  paste0(
    "https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/",
    "Homo_sapiens.GRCh37.87.chr.gtf.gz"
  ),
  destfile = gtf_path
)

gtf <- rtracklayer::import(gtf_path)

 
plotRegion <- function(chrnum,startBP,endBP,cpglocus){         
          region <- GRanges(
            seqnames = chrnum,
            ranges = IRanges(start = startBP-10000, end = endBP+ 10000),
            strand = "*"
          )
          #> 
          gtf <- gtf %>% dplyr::as_tibble()
          class(gtf)
          gtf_gr <- GRanges(
            seqnames = gtf$seqnames,
            ranges = IRanges(gtf$start, gtf$end),
            strand = gtf$strand
          )
          
          mcols(gtf_gr) <- gtf[, c(
            "type",
            "gene_name",
            "transcript_name",
            "transcript_biotype"
          )]
          
          seqlevelsStyle(gtf_gr) <- "UCSC"
          seqlevelsStyle(region) <- "UCSC"
          region_hits <- subsetByOverlaps(gtf_gr, region)
          
          region_hits_df <- as.data.frame(region_hits)
          region_hits_df <- region_hits_df %>%
            dplyr::filter(type %in% c(
              "gene",
              "transcript",
              "exon",
              "CDS",
              "UTR"
            ))
          
          region_hits_df <- region_hits_df %>%
            dplyr::select(
              seqnames,
              start,
              end,
              strand,
              type,
              gene_name,
              transcript_name,
              transcript_biotype
            )
          
          region_hits_df_exons <- region_hits_df %>% dplyr::filter(type == "exon")
          
          filtered_df <- region_hits_df_exons %>%
            group_by(transcript_name) %>%
            filter(
              max(end) >= startBP &  # not entirely upstream
                min(start) <= endBP    # not entirely downstream
            ) %>%
            ungroup()
          
          
          
          filtered_df %>%
            ggplot(aes(
              xstart = start,
              xend = end,
              y = transcript_name
            )) +
            geom_range(fill = "black"
            ) +
            geom_intron(
              data = to_intron(filtered_df, "transcript_name"),
              aes(strand = strand),arrow.min.intron.length = 200,
              arrow = grid::arrow(ends = "last", length = grid::unit(0.05, "inches"))
            )+theme_bw()+labs(y = NULL)+geom_vline(xintercept = cpglocus,color = "red")
}


miniman1 <- plotRegion("chr17",39696336-18000,39696336+18000,39696336)
miniman2 <- plotRegion("chr1",28520392-28000,28520392+28000,28520392)
miniman3 <- plotRegion("chr11",122587207-88000,122587207+88000,122587207)


regionPlots <- plot_grid(plot1,miniman1,plot3,miniman3,plot2,miniman2+xlab("Genomic Location"),rel_widths = c(0.4,0.6),align = "vh",rel_heights = c(1,1,1),ncol = 2, nrow = 3)

head(allSumstat)


pdf("Plot2_Forests_EffectComp.pdf",width = 10,height = 11.2)
plot_grid(regionPlots,effectComp,ncol = 1,rel_heights = c(0.6,0.4),labels = c("B)","C)"),label_size = 20)
dev.off()
