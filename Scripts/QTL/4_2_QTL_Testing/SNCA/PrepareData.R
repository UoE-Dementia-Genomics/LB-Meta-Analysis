args<-commandArgs(trailingOnly = TRUE)

suppressMessages(library(stringr))
suppressMessages(library(data.table))

fam.file <- args[1]
eigenvec.file <- args[2]
exp.rds.file <- args[3]
exp.pheno.file <- args[4]
geneLocation.csv.file <- args[5]
covar.fact <- args[6] 
covar.num <- c("Age","NeuNNeg_Sox10Neg_IRF8Pos","NeuNPos_SOX6Neg","NeuNPos_SOX6Pos","NeuNNeg_Sox10Neg_IRF8Neg")
OutPrefix <- args[8] 
chr <- args[9]

cat("############################### PrepareData.R script ###########################################\n")
cat("Input arguments:\n")
cat("     fam file: ",fam.file,"\n")
cat("     eigenvec file: ",eigenvec.file,"\n")
cat("     Methylation/Expression matrix file: ",exp.rds.file,"\n")
cat("     Phenotype file: ",exp.pheno.file,"\n")
cat("     CpG/Gene location file file: ",geneLocation.csv.file,"\n")
cat("     Factor covariates: ",covar.fact,"\n")
cat("     Numeric covariates: ",covar.num,"\n")
cat("     Chromosome: ",chr,"\n")
cat("     Output files prefix: ",OutPrefix,"\n")
cat("#################################################################################################\n")

if(chr=="all"){
  chr <- seq(1,22,1)
}else{
  chr <- as.numeric(str_split(trimws(chr), pattern = ",", simplify = T)[1,])
}

exp.txt.file = paste0(OutPrefix,".exp.txt")
geneLocation.txt.file = paste0(OutPrefix,".gene.loc.txt")
covariat.file = paste0(OutPrefix,".covariates.txt")

if(file.exists(fam.file)){
  cat("Reading fam file...\n")
  samples <- read.table(file=fam.file,sep="\t",header=F,stringsAsFactors=F)
  rownames(samples) <- samples$V2
}else{
  stop("Unable to find ",fam.file)
}

if(file.exists(eigenvec.file)){
  cat("Reading genotype eigenvectors...\n")
  eigenvec <- read.table(file=eigenvec.file,header=F,row.names = 2)[,c(-1)]
}else{
  stop("Unable to find ",eigenvec.file)
}

if(file.exists(exp.rds.file)){
  cat("Reading expression data...\n")
  exp_all <- as.data.frame(readRDS(exp.rds.file))
}else{
  stop("Unable to find ",exp.rds.file)
}

if(!all(sapply(list(colnames(exp_all),rownames(eigenvec)), FUN = identical, rownames(samples)))){
  warning("Sample names in expression matrix and genotype data are not exactly matched! The intersection will be consider.")
  shared_names <- Reduce(intersect, list(colnames(exp_all),rownames(samples),rownames(eigenvec)))
  exp_all <- exp_all[,colnames(exp_all) %in% shared_names]
  samples <- samples[rownames(samples) %in% shared_names,]
  
  exp_all <- exp_all[shared_names]
  samples <- samples[match(shared_names,rownames(samples)),]
  eigenvec <- eigenvec[match(shared_names,rownames(eigenvec)),]
}

if(file.exists(geneLocation.csv.file)){
  cat("Reading Gene location data...\n")
  gene_loc_all <- read.csv(file=geneLocation.csv.file,stringsAsFactors = F,row.names = 1)
  gene_loc_all <- cbind.data.frame(rownames(gene_loc_all),gene_loc_all$chr,gene_loc_all$start,gene_loc_all$end)
  names(gene_loc_all) <- c("geneid","chr","left","right")
}else{
  stop("Unable to find ",geneLocation.csv.file)
}
if(!identical(row.names(exp_all) , gene_loc_all$geneid)){
  warning("Gene/CpG names in expression matrix and genotype data are not exactly matched! Trying to match them...")
  index <- match(row.names(exp_all) , gene_loc_all$geneid)
  gene_loc_all <- gene_loc_all[index,]
}

exp_all <- cbind.data.frame(gene_loc_all$geneid,exp_all)
names(exp_all)[1] <- "geneid"

if(!file.exists(exp.txt.file)){
  cat("Are the gene ids matched in gene expression and location data? ")
  cat(ifelse(identical(exp_all$geneid,gene_loc_all$geneid),"Yes","NO"),"\n")
  cat("Saving Expression data...\n")
  write.table(exp_all,file = exp.txt.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(exp.txt.file," exists\n")
}
if(!file.exists(geneLocation.txt.file)){
  cat("Saving gene location data...\n")
  write.table(gene_loc_all,file = geneLocation.txt.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(geneLocation.txt.file," exists\n")
}

if(file.exists(exp.pheno.file)){
  exp_sample <- read.csv(exp.pheno.file,row.names=1,stringsAsFactors = F) 
  
  if(!identical(rownames(exp_sample) , colnames(exp_all))){
    warning("Sample names in expression matrix and phenotype data are not exactly matched! The intersection will be consider.")
    shared_names <- intersect(colnames(exp_all),rownames(exp_sample))
    exp_sample <- exp_sample[rownames(exp_sample) %in% shared_names,]
    exp_sample <- exp_sample[match(shared_names,rownames(exp_sample)),]
  }
}
for (c in covar.fact) {
  exp_sample[,c] <- as.numeric(as.factor(exp_sample[,c]))
}

for (c in covar.num) {
  exp_sample[,c] <- as.numeric(exp_sample[,c])
}

for (i in chr) {
  cat("************************\n")
  cat("Working on chr ",i,":\n")

  genotype.txt.file = paste0(OutPrefix,".chr",i,".genotype.txt")
  genotype.raw.file = paste0(OutPrefix,".chr",i,".raw")
  snp.file = paste0(OutPrefix,".snps",".chr",i,".txt")
  snp.loc.file=paste0(OutPrefix,".snps.loc",".chr",i,".txt")

  if((!file.exists(snp.file)) | (!file.exists(snp.loc.file))){
    if(file.exists(genotype.txt.file)){
      cat(paste("Reading Genotype text data chr",i,"...\n"))
      snps<-fread(genotype.txt.file, data.table = FALSE,nThread=16)
    } else
      if(file.exists(genotype.raw.file)){

        cat(paste("Reading Genotype raw data chr",i,"...\n"))
        snps<-fread(genotype.raw.file, data.table = FALSE,nThread=16)

        snps <- t(snps)
        snps <- cbind.data.frame(rownames(snps),snps)
        names(snps) <- c("snpid",snps[2,c(-1)])
        snps <- snps[c(-1:-6),]
        write.table(snps,file=genotype.txt.file,sep="\t",row.names = F,col.names = T,quote = F)
      }else{
        stop("Unable to find genotype files for chr ",i,":\n",genotype.txt.file,"\n",genotype.raw.file)
      }

    if(!identical(colnames(exp_all),colnames(snps)[-1])){
      index <- match(colnames(exp_all)[-1],colnames(snps)[-1])
      snps <- snps[,c(1,index+1)]
    }
    snps_loc <- str_split(snps$snpid,pattern=':',simplify=T)[,c(1,2)]
    snps_loc <- cbind.data.frame(snps$snpid,snps_loc)
    names(snps_loc) <- c("snpid","chr","pos")

    cat(paste("Are the SNP ids matched in in chr",i," data? "))
    cat(ifelse(identical(snps$snpid,snps_loc$snpid),"Yes","NO"),"\n")

    cat(paste("Are the sample ids matched in expression data and chr",i," SNP data? "))
    cat(ifelse(identical(colnames(exp_all)[-1], colnames(snps)[-1]),"Yes","NO"),"\n")

    cat(paste("Saving Chr",i,"...\n"))
    write.table(snps,file = snp.file,quote = F,col.names = T,row.names = F,sep = '\t')
    write.table(snps_loc,file = snp.loc.file,quote = F,col.names = T,row.names = F,sep = '\t')


  }else
    cat("Both files already exist:\n",snp.file,"\n",snp.loc.file,"\n")

}
cat("******************************************\n")
if(!file.exists(covariat.file)){
  cat("Generating covariates file...\n")
  
  covariates <- cbind.data.frame(exp_sample[,covar.fact],exp_sample[,covar.num],eigenvec$V3,eigenvec$V4,eigenvec$V5)
  names(covariates) <- c(covar.fact,covar.num,"PC1","PC2","PC3")
  
  rownames(covariates) <- colnames(exp_all)[-1]
  covariates <- t(covariates)
  covariates <- cbind.data.frame(rownames(covariates),covariates)
  names(covariates)[1] <- "id"
  
  cat("Sample ids matched in covariate file with the others? ")
  if(file.exists(exp.pheno.file)){
    cat(ifelse(all(sapply(list(colnames(exp_all)[-1], rownames(exp_sample),rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"),"\n")
  } else{
    cat(ifelse(all(sapply(list(colnames(exp_all)[-1],rownames(eigenvec),colnames(covariates)[-1]), FUN = identical, rownames(samples))),"Yes","NO"),"\n")
  }
  
  cat("Saving covaraites file...\n")
  write.table(covariates,file = covariat.file,quote = F,col.names = T,row.names = F,sep = '\t')
}else{
  cat(covariat.file," exist\n")
}

cat("Preparing data is done!\n")

