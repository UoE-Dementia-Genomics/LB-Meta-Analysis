library(data.table)
library(stringr)
library(ggrepel)
library(dplyr)         
library(tidyverse)  
library(tidyr)
library(stringr)
library(plyr)
library(tibble)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)  
library(plyranges)
library(AnnotationHub)
library(reshape2)
library(ggplot2)
library(cowplot)
library(qqman)
library(parallel)
library(Haplin) 
library(gridExtra)
library(ggbio)  
library(EnsDb.Hsapiens.v75)
library(ggplotify)
library(ggrepel)

#load data
epicMani <-  fread("/lustre/projects/Research_Project-T112069/Meth/reference/MethylationEPIC_v-1-0_B4.csv",skip = 7)
epicMani <- as.data.frame(epicMani)
rownames(epicMani) <- epicMani$IlmnID

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")

# Process meta summary statistics for DMR analysis
resMetaCC <- as.data.frame(allSumstat)
resMetaCC <- resMetaCC[order(resMetaCC$Fixed_P),]

RESULTS <- cbind(resMetaCC,epicMani[rownames(resMetaCC),])

RESULTS <- cbind(resMetaCC,epicMani[rownames(resMetaCC),])
colnames(RESULTS)[which(colnames(RESULTS)== "Fixed_P")] <- "P"
colnames(RESULTS)[which(colnames(RESULTS)== "MAPINFO")] <- "BP"


# Format data for comb-p
dmr <- data.frame("chrom" = 	paste0("chr", RESULTS$CHR), 
                  "start" = 	as.numeric(RESULTS$BP),
                  "end" 	= 	as.numeric(RESULTS$BP),
                  "pvalue" = 	RESULTS$P)

dmr<-na.omit(dmr)
dmr<-dmr[order(dmr[,1],dmr[,2]),]
head(dmr)
write.table(dmr, file="/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/combp/highPath.txt",sep = "\t", row.names = FALSE,col.names = TRUE, quote = FALSE) # save to run combp
# The below line was run for combp analysis
# comb-p pipeline -c 4 --dist 500 --seed 1.0e-4 --anno hg19 -p  out highPath.txt

#######################################
# The below section was used for plotting up minimanhattan views of the DMR results and SNCA region




#SNCA PLOT
sub <- RESULTS[which(RESULTS$CHR == "4" & RESULTS$BP >=(90647041-10000) & RESULTS$BP <=(90760556+10000)),]
# sub <- sub[which(sub$CHR == "4"),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  # geom_rect(aes(xmin = 90647041,xmax = 90760556,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)

ensdb <- EnsDb.Hsapiens.v75
gr <- GRanges(seqnames = "4", IRanges((90647041),(90760556)), strand = "*")
a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "tx_id" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/SNCA_Plot.pdf",width = 6,height = 6)
plot_grid(a1+xlim((90647041-10000),(90760556+10000)), a2@ggplot+xlim((90647041-10000),(90760556+10000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.5,0.5))
dev.off()


#MYO16 PLOT
sub <- RESULTS
sub <- sub[which(sub$CHR == "13" & sub$BP > (109741030 - 100000)	& sub$BP < (109741117+100000)),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  geom_rect(aes(xmin = 109741030,xmax = 109741117,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)


ensdb <- EnsDb.Hsapiens.v75
gr <- GRanges(seqnames = "13", IRanges((109741030 - 100000),(109741117+100000)), strand = "*")

a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "tx_id" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/MYO16_Plot.pdf",width = 6,height = 6)
plot_grid(a1+xlim((109741030 - 100000),(109741117+100000)), a2@ggplot+xlim((109741030 - 100000),(109741117+100000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.5,0.5))
dev.off()


#NDRG4 PLOT
sub <- RESULTS
sub <- sub[which(sub$CHR == "16" & sub$BP > (58535416-10000)		& sub$BP < (58535556+10000)),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  geom_rect(aes(xmin = 58535416,xmax = 58535556,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)

ensdb <- EnsDb.Hsapiens.v75

gr <- GRanges(seqnames = "16", IRanges((58535416-10000),(58535556+10000)), strand = "*")

a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/NDRG4_Plot.pdf",width = 6,height = 6)
plot_grid(a1+xlim((58535416-10000),(58535556+10000)), a2@ggplot+xlim((58535416-10000),(58535556+10000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.4,0.7))
dev.off()


#LTF Plot
sub <- RESULTS
sub <- sub[which(sub$CHR == "3" & sub$BP > (46506206 - 20000)	& sub$BP < (46506519+ 10000)),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  geom_rect(aes(xmin = 46506206,xmax = 46506519,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)


ensdb <- EnsDb.Hsapiens.v75

gr <- GRanges(seqnames = "3", IRanges((46506206 - 20000),(46506519+ 10000)), strand = "*")

a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "tx_id" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)

pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/LTF_Plot.pdf",width = 6,height = 6)
plot_grid(a1+xlim((46506206 - 20000),(46506519+ 10000)), a2@ggplot+xlim((46506206 - 20000),(46506519+ 10000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.5,0.5))
dev.off()

#C3orf56 Plot
#chr3	126911830	126911953

sub <- RESULTS
sub <- sub[which(sub$CHR == "3" & sub$BP > (126911830 - 20000)	& sub$BP < (126911953+ 20000)),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  geom_rect(aes(xmin = 126911830,xmax = 126911953,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)

ensdb <- EnsDb.Hsapiens.v75

gr <- GRanges(seqnames = "3", IRanges((126911830 - 20000),(126911953+ 20000)), strand = "*")

a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "tx_id" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)


pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/C3orf56_Plot.pdf",width = 6,height = 4)
plot_grid(a1+xlim((126911830 - 20000),(126911953+ 20000)), a2@ggplot+xlim((126911830 - 20000),(126911953+ 20000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.7,0.3))
dev.off()

#CD244	
# chr1	160833741	160833863

sub <- RESULTS
sub <- sub[which(sub$CHR == "1" & sub$BP > (160833741 - 20000)	& sub$BP < (160833863+ 20000)),]
colnames(sub)[1] <- "Effect"

sub$Effect <- as.factor(sign(sub$Effect))
levels(sub$Effect) <- c("Hypo-","Hyper-")

a1 <- ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
  geom_rect(aes(xmin = 160833741,xmax = 160833863,ymin = 0, ymax = -log10(min(sub$P))),color = "lightblue",fill = "lightblue",alpha = 0.2)+
  geom_point()+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  theme_bw()+
  labs(color = "Direction")+
  scale_color_brewer(palette = "Set2",direction = -1)

ensdb <- EnsDb.Hsapiens.v75

gr <- GRanges(seqnames = "1", IRanges((160833741 - 20000),(160833863+ 20000)), strand = "*")

a2 <- autoplot(ensdb, GRangesFilter(gr),  names.expr = "tx_id" ,color = "black",fill ="black")+ 
  labs(x= "Genomic Location")+
  theme_bw()+
  scale_x_continuous(position="top")+
  theme(axis.title.x = element_text(size = 16),
        axis.line.y =element_blank())+xlab(NULL)


pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/CD244_Plot.pdf",width = 6,height = 6)
plot_grid(a1+xlim((160833741 - 20000),(160833863+ 20000)), a2@ggplot+xlim((160833741 - 20000),(160833863+ 20000)),
          ncol = 1, nrow = 2,
          rel_widths = c(1,1),
          align = "v",axis = "lr",rel_heights = c(0.5,0.5))
dev.off()


