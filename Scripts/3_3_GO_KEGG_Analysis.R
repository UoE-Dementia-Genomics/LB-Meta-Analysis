######### Gene ontology (methylGSA) #####
########

.libPaths("/lustre/projects/Research_Project-T112069/packages")
options(rlib_downstream_check = FALSE)
library(BiocManager)
library(methylGSA)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ggplot2)
library(missMethyl)
library(stringr)
# Load in meta results:

load("/lustre/projects/Research_Project-T112069/Meth/EWAS/Meta/Meta/MultiMeta_FullCohort.Rdata")

resMetaCC <- as.data.frame(allSumstat)
resMetaCC <- resMetaCC[order(resMetaCC$Fixed_P),]


pvalMeta <- resMetaCC$Fixed_P
names(pvalMeta) <- rownames(resMetaCC)
sigCpg <- rownames(resMetaCC[which(resMetaCC$Fixed_P < 1e-5),])

# Run the ORA method in methylGSA


#### GO Analysis
# resORA_GO = methylRRA(cpg.pval = pvalMeta, method = "ORA",GS.type = "GO",array.type = "EPIC")
# save(resORA_GO,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_GO.Rdata")
# resORA_KEGG = methylRRA(cpg.pval = pvalMeta, method = "ORA",GS.type = "KEGG",array.type = "EPIC")
# save(resORA_KEGG,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_KEGG.Rdata")
# resORA_Re = methylRRA(cpg.pval = pvalMeta, method = "ORA",GS.type = "Reactome",array.type = "EPIC")
# save(resORA_Re,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_Re.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_GO.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_KEGG.Rdata")
load("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_Re.Rdata")

resORA_GO <- resORA_GO[which(resORA_GO$padj < 0.05 & resORA_GO$Count >= 3),]
head(resORA_GO)
resORA_GO$Description <- factor(resORA_GO$Description, levels = resORA_GO$Description) 
oraGOplot <- ggplot(resORA_GO, aes(x = as.factor(Description), y = -log10(padj)))+
  geom_bar(stat = "identity",fill = "lightgreen",color = "black")+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed")+
  coord_flip()+
  ylab("-log10(FDR p-value)")+
  xlab("Ontology Term")+
  theme_bw()

png("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/goPlot.png",res = 200, height  = 1400, width = 1400)
oraGOplot
dev.off()
write.csv(resORA_GO,file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/goORAres.csv")


#### Reactome Analysis
resORA_Re <- resORA_Re[which(resORA_Re$pvalue < 0.05 & resORA_Re$Count >= 3),]
resORA_Re$Description <- factor(resORA_Re$Description, levels = resORA_Re$Description) 

oraREplot <- ggplot(resORA_Re, aes(x = as.factor(Description), y = -log10(pvalue)))+
  geom_bar(stat = "identity")+
  coord_flip()+
  ylab("-log10(p-value)")+
  xlab("Reactome Term")

#### KEGG Analysis
resORA_KEGG <- resORA_KEGG[which(resORA_KEGG$pvalue < 0.05 & resORA_KEGG$Count >= 3),]
resORA_KEGG$Description <- factor(resORA_KEGG$Description, levels = resORA_KEGG$Description) 

oraKEGGplot <- ggplot(resORA_KEGG, aes(x = factor(Description,levels = rev(resORA_KEGG$Description)), y = -log10(pvalue),fill = padj < 0.05,color = padj < 0.05))+
  geom_bar(stat = "identity")+
  coord_flip()+
  ylab("-log10(p-value)")+
  xlab("KEGG Term")+
  scale_fill_manual(values = c("grey","#1B9E77"))+
  scale_color_manual(values = c("grey","black"))+
  theme_bw()+labs(fill = "FDR Significant",color = "FDR Significant")

#### For the GO analysis, use revigo to reduce down redundant GO terms
library(rrvgo)

scores <- setNames(-log10(resORA_GO$padj), resORA_GO$ID)
simMatrixBP <- calculateSimMatrix(rownames(resORA_GO),
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")


reducedTermsBP <- reduceSimMatrix(simMatrixBP,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")


simMatrix <- calculateSimMatrix(rownames(resORA_GO),
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")


reducedTermsCC <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")

simMatrix <- calculateSimMatrix(rownames(resORA_GO),
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")


reducedTermsMF <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
# Note, no significant MF terms so excluded at this point


# Merge reduced terms
reducedTerms <- rbind(reducedTermsCC,reducedTermsBP)

# Annotate reduced terms
resORA_GO$parentTerm <- resORA_GO$Description
resORA_GO[rownames(reducedTerms),"parentTerm"] <- reducedTerms$parentTerm
resORA_GO$Sort <- as.numeric(!resORA_GO$Description == resORA_GO$parentTerm)
resORA_GO <- resORA_GO[order(resORA_GO$parentTerm,resORA_GO$Sort),]
resORA_GO$plotIndex <- 1:nrow(resORA_GO)

#Only retain relevant terms
resORA_GO <- resORA_GO[which(resORA_GO$Description %in% unique(resORA_GO$parentTerm)),]


# Save GO results
write.csv(resORA_GO, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_GOrrvgoRerun.csv")
write.csv(resORA_KEGG, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_KEGGRerun.csv")
write.csv(resORA_Re, file = "/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/resORA_ReRerun.csv")

resORA_GO <- resORA_GO[order(resORA_GO$pvalue),]

#Plot the relevant terms
library(cowplot)
pdf("/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/goRRVGO_Plot_Rerun.pdf",width = 12,height = 4.2)
plot_grid(
  ggplot(resORA_GO, aes(y = factor(Description, levels = rev(resORA_GO$Description)),x = -log10(pvalue)))+
    geom_point(shape = 21,fill = "#1B9E77",size = 6)+
    theme_bw()+
    scale_x_continuous(limits = c(3,5), expand = c(0,0))+
    scale_size_continuous(limits = c(0,0.06))+
    labs(x = "-log10(P-value)", y = "Ontology Term"),
  oraKEGGplot
)
dev.off()


#This section of code is adapted from the methylGSA package to extract the RRA gene ID's used for
# testing

FullAnnot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]

## get the first gene in each USCS_RefGene_Name
temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
              '[', 1, FUN.VALUE=character(1))
FullAnnot$UCSC_RefGene_Name = temp
## get the first gene group in each UCSC_RefGene_Group
temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Group,split=";"),
              '[', 1, FUN.VALUE=character(1))
FullAnnot$UCSC_RefGene_Group = temp


cpg.pval <- pvalMeta
cpg.intersect = intersect(names(cpg.pval), rownames(FullAnnot))
cpg.pval = cpg.pval[cpg.intersect]
## match user input CpG to our FullAnnot database
FullAnnot.sub = FullAnnot[names(cpg.pval), ]
## change CpG ids to gene symbols
names(cpg.pval) = FullAnnot.sub$UCSC_RefGene_Name
## convert cpg.pval to a list, each element of this
# list is a gene and it's corresponding cpg pvalue
geneID.list = split(cpg.pval, names(cpg.pval))
library(RobustRankAggreg)
## for each gene, compute its rho score
rho = vapply(geneID.list, rhoScores, FUN.VALUE = 1)

GS.list = getGS(names(geneID.list), GS.type = "GO")
GS.list = lapply(GS.list, na.omit)

## calculate gene set size
GS.sizes = vapply(GS.list, length, FUN.VALUE = 0)
## filter gene sets by their sizes
GS.list.sub = GS.list[GS.sizes>=100 & GS.sizes<=500]
message(length(GS.list.sub), " gene sets are being tested...")
size = vapply(GS.list.sub, length, FUN.VALUE = 0)
ID = names(GS.list.sub)
gs.pval = rep(NA, length(GS.list.sub))
Count = rep(NA, length(GS.list.sub))

rhoadj = p.adjust(rho, method = "BH")
DEgenes = names(rhoadj)[rhoadj<0.05]
# countTable <- table(GS.list.sub)

deSummary <- data.frame(DEgenes = DEgenes,Count = NA, rhoP = NA, adjP = NA, minP = NA)
rownames(deSummary) <- deSummary$DEgenes
for(x in rownames(deSummary)){
  countSet <- unlist(geneID.list[x])
  
  deSummary[x,"Count"] <- length(countSet)
  deSummary[x,"rhoP"] <- rho[x]
  deSummary[x,"adjP"] <- rhoadj[x]
  deSummary[x,"minP"] <- min(countSet)
  
}

write.csv(deSummary,"/lustre/projects/Research_Project-T112069/Meth/Meta/FiveCell/methylGSA/DEgenesRerun.csv")
