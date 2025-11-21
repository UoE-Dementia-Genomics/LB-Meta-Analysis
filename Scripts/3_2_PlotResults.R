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
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/UKBBN/res_UKBBN_CrossCortex.Rdata")
res_UKBBN <- res_UKBBN[,-3]
colnames(res_UKBBN) <- c("Values","SE","T","P")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/NBB/NBB_knownCovar_SV3.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/BDR/BDR_Res_SV1.Rdata")

resUK <- res_UKBBN
resNBB <- res_NBB_SV3
resBDR <- BDR_res_LB_SV1

overlap <-Reduce(intersect, list(rownames(resUK),rownames(resNBB),rownames(resBDR)))
resUK <- resUK[overlap,]
resNBB <- resNBB[overlap,]

resBDR <- resBDR[overlap,]

colnames(resUK) <- c("Beta_UK","SE_UK","T_UK","P_UK")
colnames(resNBB) <- c("Beta_NBB","SE_NBB","T_NBB","P_NBB")
colnames(resBDR) <- c("Beta_BDR","SE_BDR","T_BDR","P_BDR")

resAll <- cbind(resUK, resNBB,resBDR)
resAll <- as.data.frame(resAll)

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
resMetaCC <- allSumstat[overlap,]

#########################
##### QQ plots ##########
#########################
# Calculate Lambdas for each analysis
LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  return(lambda)
}

resUK <- as.data.frame(resUK)
resNBB <- as.data.frame(resNBB)
resBDR <- as.data.frame(resBDR)
infUK <- signif(LambdaInf(na.omit(resUK[order(resUK$P_UK),"P_UK"])),4)
infNBB <- signif(LambdaInf(resNBB[order(resNBB$P_NBB),"P_NBB"]),4)
infBDR <- signif(LambdaInf(resBDR[order(resBDR$P_BDR),"P_BDR"]),4)
infMeta <- signif(LambdaInf(resMetaCC[,3]),4)

# Make a mega data frame of all observed vs. expected p.values
# names are formatted in "Cohort" "Lambda Inflation"
tSet <- rbind(data.frame(dataset = paste("UKBBN",infUK),observed.pvals = resUK[order(resUK$P_UK),"P_UK"],expected.pvals = seq(from = 1/length(overlap),to = 1, by = 1/length(overlap))),
              data.frame(dataset = paste("NBB",infNBB),observed.pvals = resNBB[order(resNBB$P_NBB),"P_NBB"],expected.pvals = seq(from = 1/length(overlap),to = 1, by = 1/length(overlap))),
              data.frame(dataset = paste("BDR",infBDR),observed.pvals = resBDR[order(resBDR$P_BDR),"P_BDR"],expected.pvals = seq(from = 1/length(overlap),to = 1, by = 1/length(overlap))),
              data.frame(dataset = paste("Meta",infMeta),observed.pvals = resMetaCC[order(resMetaCC[,3]),3],expected.pvals = seq(from = 1/length(overlap),to = 1, by = 1/length(overlap))))

# Plot qq-plots
png("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Plots/quickQQ.png",width = 800,height = 800,res = 200,type = "cairo-png")
ggplot(tSet, aes(y = -log10(observed.pvals), x = -log10(expected.pvals), color = factor(dataset, levels = c(paste("Meta",infMeta),
                                                                                                            paste("UKBBN",infUK),
                                                                                                            paste("NBB",infNBB),
                                                                                                            paste("BDR",infBDR)))))+
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
dev.off()

################################
##### Manhattan plots ##########
################################
# Load manifest data and merge
epicMani <-  fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv",skip = 7)
epicMani <- as.data.frame(epicMani)
rownames(epicMani) <- epicMani$IlmnID
RESULTS <- cbind(resMetaCC,epicMani[rownames(resMetaCC),])
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

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")
highNFTsub <- as.data.frame(allSumstat)
colnames(highNFTsub)[which(colnames(highNFTsub)%in% c("Fixed_Effect","Fixed_SE","Fixed_P"))] <- c("Effect","SE","P")
highNFTsub <- highNFTsub[which(p.adjust(highNFTsub$P,method = "fdr") < 0.05),]
load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_ThalControlled.Rdata")
thalSumstat <- as.data.frame(thalSumstat)
colnames(thalSumstat)[which(colnames(thalSumstat)%in% c("Fixed_Effect","Fixed_SE","Fixed_P"))] <- c("Effect_Thal","SE_Thal","P_Thal")


plotComp <- cbind(highNFTsub[,c("Effect","SE","P")],thalSumstat[rownames(highNFTsub),c("Effect_Thal","SE_Thal","P_Thal")])

corPlot <- ggplot(plotComp, aes(x = Effect, y = Effect_Thal)) +
  # Error bars for Effect
  geom_errorbar(aes(ymin = Effect_Thal - SE_Thal, ymax = Effect_Thal + SE_Thal), width = 0, color = "gray50") +
  # Error bars for Effect_Low
  geom_errorbarh(aes(xmin = Effect - SE, xmax = Effect + SE), height = 0, color = "gray50") +
  # Points with color and size mappings
  geom_point(aes(fill = -log10(P_Thal), size = -log10(P)), shape = 21, color = "black") +
  scale_fill_viridis_c() + # Color mapping for fill
  scale_size(range = c(2, 6)) + # Adjust point size range
  theme_bw() +
  labs(x = "Primary Meta Estimate (+/- SE)", y = "Thal Controlled Meta Estimate (+/- SE)", fill = "-log10(Thal Controlled P-value)",size = "-log10(Primary P-value)") +
  theme(legend.position = "bottom", 
        legend.box = "vertical")+
  ylim(-0.0081,0.0051)+
  xlim(-0.0081,0.0051)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)


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
