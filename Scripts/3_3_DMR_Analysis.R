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
library(ggtranscript)

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






# Load DMR Data
dmrs <- read.table("out.regions-t.bed",header = T)
dmrs <- dmrs[order(dmrs$z_p),]
dmrs <- dmrs[which(dmrs$z_sidak_p < 0.05 & dmrs$n_probes > 2),]

for(t in 1:nrow(dmrs)){
  
  chr <- gsub("chr","",dmrs[t,"chrom"])
  minBP <- dmrs[t,"start"] - 1
  maxBP <- dmrs[t,"end"] + 1
  
  tempManiIlmnID <- unlist(c(epicMani[which(epicMani$CHR == chr &
                                     epicMani$MAPINFO > minBP &
                                     epicMani$MAPINFO < maxBP),"IlmnID"]))
  
  tempManiUCSC <- c(epicMani[which(epicMani$CHR == chr &
                               epicMani$MAPINFO > minBP &
                               epicMani$MAPINFO < maxBP),"UCSC_RefGene_Name"])
  
  dmrs[t,"minp"] <- min(allSumstat[match(tempManiIlmnID,rownames(allSumstat)),"Fixed_P"])
  
  dmrs[t,"Gene"] <- unique(unlist(strsplit(unlist(tempManiUCSC), ";")))
  
}

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
  scale_color_brewer(palette = "Set2",direction = -1)+theme(legend.position = "top")

write.csv(sub, file= "sncaPlotData.csv")

dir.create("~/R/win-library/4.4", recursive = TRUE)
.libPaths(c("~/R/win-library/4.4", .libPaths()))
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



plotRegion <- function(chrnum,startBP,endBP,dmrBP_min,dmrBP_max){         
  region <- GRanges(
    seqnames = chrnum,
    ranges = IRanges(start = startBP, end = endBP),
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
    "gene_name","transcript_id",
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
      transcript_name,transcript_id,
      transcript_biotype
    )
  
  region_hits_df_exons <- region_hits_df %>% dplyr::filter(type == "exon")
  
  filtered_df <- region_hits_df_exons %>%
    group_by(transcript_name) %>%
    filter(
      max(end) >= dmrBP_min &  # not entirely upstream
        min(start) <= dmrBP_max    # not entirely downstream
    ) %>%
    ungroup()
  
  
  
  filtered_df %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(fill = "black"
    ) +
    geom_intron(
      data = to_intron(filtered_df, "transcript_id"),
      aes(strand = strand),arrow.min.intron.length = 200,
      arrow = grid::arrow(ends = "last", length = grid::unit(0.05, "inches"))
    )+theme_bw()+labs(y = NULL)+scale_y_discrete(position = "right")
}


miniman1 <- plotRegion("chr4",90647041-10000,90760556+10000)
pdf("SNCA_Plot.pdf",height = 6, width = 6)
plot_grid(a1+xlim(90647041-10000,90760556+1000)+  theme(
  axis.text.x = element_blank(),   # Remove tick labels
  axis.title.x = element_blank()   # Remove axis title
        ),miniman1+xlim(90647041-10000,90760556+1000)+xlab("Basepair Coordinate"),ncol =1, align = "v")
dev.off()



      
plotMini <- function(chromNum, dmrStart,dmrEnd,buffer){      
      #SNCA PLOT
      sub <- RESULTS[which(RESULTS$CHR == gsub("chr","",chromNum) & RESULTS$BP >=(dmrStart-buffer) & RESULTS$BP <=(dmrEnd+buffer)),]
      # sub <- sub[which(sub$CHR == "4"),]
      colnames(sub)[1] <- "Effect"
      
      sub$effect_raw <- sub$Effect
      sub$Effect <- as.factor(sign(sub$Effect))
      levels(sub$Effect) <- c("Hypo-","Hyper-")
      
      
      
      ucsc <- sub[which(sub$BP >=(dmrStart) & sub$BP <=(dmrEnd)),"UCSC_RefGene_Name"] 
      
      ucsc<- unique(unlist(strsplit(ucsc, ";")))
      
      ymaxi <- (-log10(min(sub$P)))
      
      
      ggplot(sub,aes(x = BP, y = -log10(P), color = Effect))+
        geom_rect(aes(xmin = dmrStart,xmax = dmrEnd,ymin = 0, ymax = ymaxi+0.5),color = "lightblue",fill = "lightblue",alpha = 0.2)+
        geom_point()+
        geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
        theme_bw()+
        labs(color = "Direction")+
        scale_color_brewer(palette = "Set2",direction = -1)+theme(legend.position = "top")+
        ggtitle(paste(ucsc,"Region"))+scale_y_continuous(limits = c(0,ymaxi+0.5), expand = c(0,0))
      
}

miniman1 <- plotMini(dmrs[1,"chrom"],dmrs[1,"start"],dmrs[1,"end"],10000)+xlim(dmrs[1,"start"]-10000,dmrs[1,"end"]+10000)+  theme(
  axis.text.x = element_blank(),   # Remove tick labels
  axis.title.x = element_blank())
miniman1_data <- ggplot_build(miniman1)
miniman1_data <- cbind(miniman1_data$data[[2]],miniman1_data$data[[1]][,c("xmin","xmax")])
write.csv(miniman1_data,file = "miniman1_data.csv")

miniman2 <- plotMini(dmrs[2,"chrom"],dmrs[2,"start"],dmrs[2,"end"],10000)+xlim(dmrs[2,"start"]-7000,dmrs[2,"end"]+7000)+  theme(
  axis.text.x = element_blank(),   # Remove tick labels
  axis.title.x = element_blank())
miniman2_data <- ggplot_build(miniman2)
miniman2_data <- cbind(miniman2_data$data[[2]],miniman2_data$data[[1]][,c("xmin","xmax")])
write.csv(miniman2_data,file = "miniman2_data.csv")

miniman3 <- plotMini(dmrs[3,"chrom"],dmrs[3,"start"],dmrs[3,"end"],10000)+xlim(dmrs[3,"start"]-10000,dmrs[3,"end"]+10000)+  theme(
  axis.text.x = element_blank(),   # Remove tick labels
  axis.title.x = element_blank())
miniman3_data <- ggplot_build(miniman3)
miniman3_data <- cbind(miniman3_data$data[[2]],miniman3_data$data[[1]][,c("xmin","xmax")])
write.csv(miniman3_data,file = "miniman3_data.csv")

miniman4 <- plotMini(dmrs[4,"chrom"],dmrs[4,"start"],dmrs[4,"end"],10000)+xlim(dmrs[4,"start"]-10000,dmrs[4,"end"]+10000)+  theme(
  axis.text.x = element_blank(),   # Remove tick labels
  axis.title.x = element_blank())
miniman4_data <- ggplot_build(miniman4)
miniman4_data <- cbind(miniman4_data$data[[2]],miniman4_data$data[[1]][,c("xmin","xmax")])
write.csv(miniman4_data,file = "miniman4_data.csv")


region1 <- plotRegion(dmrs[1,"chrom"],dmrs[1,"start"]-50000,dmrs[1,"end"]+50000,dmrs[1,"start"]-10000,dmrs[1,"end"]+10000)
region2 <- plotRegion(dmrs[2,"chrom"],dmrs[2,"start"]-50000,dmrs[2,"end"]+50000,dmrs[2,"start"]-7000,dmrs[2,"end"]+7000)
region3<- plotRegion(dmrs[3,"chrom"],dmrs[3,"start"]-50000,dmrs[3,"end"]+50000,dmrs[3,"start"]-10000,dmrs[3,"end"]+10000)
region4 <- plotRegion(dmrs[4,"chrom"],dmrs[4,"start"]-50000,dmrs[4,"end"]+50000,dmrs[4,"start"]-10000,dmrs[4,"end"]+10000)

region1 <- region1+coord_cartesian(xlim = c(dmrs[1,"start"]-10000,dmrs[1,"end"]+10000))+xlab("Basepair Coordinate")

region2 <- region2+coord_cartesian(xlim = c(dmrs[2,"start"]-7000,dmrs[2,"end"]+7000))+xlab("Basepair Coordinate")

region3 <- region3+coord_cartesian(xlim = c(dmrs[3,"start"]-10000,dmrs[3,"end"]+10000))+xlab("Basepair Coordinate")

region4 <- region4+coord_cartesian(xlim = c(dmrs[4,"start"]-10000,dmrs[4,"end"]+10000))+xlab("Basepair Coordinate")

fillingPlot <- ggplot() + theme_void()
       miniman3,miniman4,region3,plot_grid(region4,fillingPlot,ncol = 1,align = "v",axis = "lr"),align = "v",axis = "lr",ncol = 2,nrow = 4)

legend<- get_legend(miniman4)
pdf("Plot3_DMRS.pdf",width = 10,height = 10)
plot_grid(plot_grid(miniman1+theme(legend.position = "none"),region1,
          miniman3+theme(legend.position = "none"),region3,align = "v",axis = "lr",ncol = 1,labels = c("C)","","E)","")),
          plot_grid(miniman2+theme(legend.position = "none"),region2,NULL,
          miniman4+theme(legend.position = "none"),region4,legend,labels = c("D)","","","F)","",""),align = "v",axis = "lr",ncol = 1,rel_heights = c(0.3,0.2,0.1)))
dev.off()

dmrPlots <- plot_grid(plot_grid(miniman1+theme(legend.position = "none"),
                                region1,
                                miniman3+theme(legend.position = "none"),
                                region3+ theme(
                                                axis.text.y = element_blank(),   # Remove tick labels
                                                axis.title.y = element_blank(),
                                                axis.ticks.y = element_blank()),
                                align = "v",axis = "lr",ncol = 1,labels = c("C)","","E)",""),rel_heights = c(0.2,0.2,0.2,0.3)),
                      plot_grid(miniman2+theme(legend.position = "none"),
                                region2,NULL,
                                miniman4+theme(legend.position = "none"),
                                region4,legend,
                                labels = c("D)","","","F)","",""),align = "v",axis = "lr",ncol = 1,rel_heights = c(0.2,0.1,0.1,0.2,0.1,0.2)))

pdf("Plot3_Merged.pdf",width = 10,height = 12)
plot_grid(goPlots,dmrPlots,ncol = 1,rel_heights = c(0.25,0.6))
dev.off()
