
args<-commandArgs(trailingOnly = TRUE)

OutPrefix <- args[1]
chr <- args[2]
Dist <- as.numeric(trimws(args[3]))
cis.pval <- as.numeric(trimws(args[4]))
trans.pval <- as.numeric(trimws(args[5]))
trans.cross.chr <- ifelse(tolower(args[6])=="yes",T,F)

OutDir <- dirname(OutPrefix)
OutPrefix <- basename(OutPrefix)
covariates.file <- paste0(OutDir,"/",OutPrefix,".covariates.txt")
SNP.file <- paste0(OutDir,"/",OutPrefix,".snps.chr",chr,".txt")
expression.file <- paste0(OutDir,"/",OutPrefix,".exp.txt")
SNP.location.file <- paste0(OutDir,"/",OutPrefix,".snps.loc.chr",chr,".txt")
gene.location.file <- paste0(OutDir,"/",OutPrefix,".gene.loc.txt")

cat("############################### Main.R script ###########################################\n")
cat("Arguments:\n")
cat("  Output directory: ",OutDir,"\n")
cat("  Output files prefix: ",OutPrefix,"\n")
cat("  Covariates file: ",covariates.file,"\n")
cat("  SNP data file: ",SNP.file,"\n")
cat("  SNP location file: ",SNP.location.file,"\n")
cat("  Expression data file: ",expression.file,"\n")
cat("  Gene location file: ",gene.location.file,"\n")
cat("  cis distance: ",Dist,"\n")
cat("  cis Pvalue: ",cis.pval,"\n")
cat("  trans Pvalue: ",trans.pval,"\n")
cat("  Get the trans results in cross chromosom mode? ",ifelse(trans.cross.chr,"Yes","No"),"\n")
cat("##########################################################################################\n")

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.
.libPaths("/lustre/projects/Research_Project-T112069/packages")
# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Covariates file name
# Set to character() for no covariates

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# Only associations significant at this level will be saved
pvOutputThreshold_cis = cis.pval;
pvOutputThreshold_tra = trans.pval;

# Distance for local gene-SNP pairs
cisDist = Dist;

cat("Reading gene/probe expression/methylation data...\n")
tryCatch(expr = {
  expression = SlicedData$new();
  expression$fileDelimiter = "\t";      # the TAB character
  expression$fileOmitCharacters = "NA"; # denote missing values;
  expression$fileSkipRows = 1;          # one row of column labels
  expression$fileSkipColumns = 1;       # one column of row labels
  expression$LoadFile(expression.file);
},
error=function(e){
  stop("The following error occured during reading the expression file:\n",conditionMessage(e))
})
cat("Reading gene/probe location data...\n")
genepos = read.table(gene.location.file, header = TRUE, stringsAsFactors = FALSE);
if(!trans.cross.chr){
  index <- which(genepos$chr==chr)
  genepos <- genepos[index,]
  expression$RowReorder(index)
}
if(nrow(expression) == 0){
  stop("There is no expression data for chr ",chr)
}

cat("Reading covariate data...\n")
tryCatch(expr = {
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  cvrt$LoadFile(covariates.file);
},
error=function(e){
  stop("The following error occured during reading the covariates file:\n",conditionMessage(e))
})

cat("Reading genotype data...\n")
tryCatch(expr = {
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$LoadFile(SNP.file);
},
error=function(e){
  stop("The following error occured during reading the genotype file related to chr ",chr,":\n",conditionMessage(e))
})
cat("Reading SNP location data...\n")
snpspos = read.table(SNP.location.file, header = TRUE, stringsAsFactors = FALSE);

dir.create(paste0(OutDir,"/",OutPrefix,".matrixEQTL.chr",chr))

output_file_name_cis =paste0(OutDir,"/",OutPrefix,".matrixEQTL.chr",chr,paste0("/",OutPrefix,".matrixEQTL.cis.chr",chr,".out"))
output_file_name_tra = paste0(OutDir,"/",OutPrefix,".matrixEQTL.chr",chr,paste0("/",OutPrefix,".matrixEQTL.trans.chr",chr,".out"))


snpspos$pos <- sub("_.*", "", snpspos$pos)
snpspos$pos <- as.numeric(snpspos$pos)
cat("Running QTL analysis...\n")
me = Matrix_eQTL_main(
  snps = snps,
  gene = expression,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE );

unlink(output_file_name_tra);
unlink(output_file_name_cis);

cat("Saving the results on: ",OutDir,"/",OutPrefix,".matrixEQTL.chr",chr,"/",OutPrefix,".matrixEQTL.chr",chr,".RData\n")
save(me,file=paste0(OutDir,"/",OutPrefix,".matrixEQTL.chr",chr,paste0("/",OutPrefix,".matrixEQTL.chr",chr,".RData")))

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
