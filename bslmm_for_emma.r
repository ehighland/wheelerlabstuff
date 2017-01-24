####by Heather E. Wheeler 20150108####
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c('22','YRI') #uncomment for testing in RStudio
"%&%" = function(a,b) paste(a,b,sep="")
library(dplyr)

###############################################
### Directories & Variables
exp.dir <- "/home/aly/Matrix_eQTL/Expression/"
snp.dir <- "/home/aly/Matrix_eQTL/SNPs/"
snp.annot.dir <- "/home/aly/Matrix_eQTL/SNPs_Location/"
#out.dir <- "/home/wheelerlab1/bslmm_scripts/testoutput/" #edit to your own output directory
out.dir <- "/home/emma/gemmaResults/"
wk.dir <- getwd() %&% "/" #working directory to retrieve direct bslmm output, bslmm automatically puts output in ./output

tis <- "HM3_LCLs"  
chromosome <- as.numeric(args[1]) 
chrname <- "chr" %&% chromosome
pop <- args[2]

getquant <- function(x) quantile(x,c(0.5,0.025,0.975)) ##pulls the median and the 95% credible sets

################################################
exp.file <- exp.dir %&% pop %&% "_Expression.txt"
exp.annot.file <- "/home/wheelerlab1/elastic_net/GRCh37_hg19_ILMN_Human-6_v2_gene_annotation_for_elasticNet.txt"
snp.file <- snp.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.txt.gz"
snp.annot.file <- snp.annot.dir %&% pop %&% "_" %&% chromosome %&% ".SNP.Location.txt.gz"

##get gene pos info
gencode <- read.table(exp.annot.file,header=TRUE)
##get snp pos info
snpcode <- read.table(snp.annot.file,header=TRUE)
##get snp allele info (needed for weight output)
allelecode <- read.table("/home/wheelerlab1/elastic_net/allele_annot_files/chr" %&% chromosome %&% "_" %&% pop %&% ".txt.gz")
colnames(allelecode) <- c("CHR","POS","SNP","refAllele","effectAllele")
rownames(allelecode) <- allelecode$POS #name by position b/c we found duplicate rsids

##read exp and chr gt dosages
exp <- read.table(exp.file, header=TRUE)
gt <- read.table(snp.file,header=TRUE)

##order by id to ensure exp id's match gt id's
exp <- exp[,order(colnames(exp))]
gt <- gt[,order(colnames(gt))]
stopifnot(colnames(exp)==colnames(gt))

##join pos info 
popgt <- left_join(snpcode,gt,by=c("snp"="id"))
popgt <- popgt[duplicated(popgt$snp)==FALSE,] #remove duplicated rsids with incorrect pos
popgt <- popgt[duplicated(popgt$pos)==FALSE,] #remove duplicated pos 
rownames(popgt) <- popgt[,3] #name by position b/c we found duplicate rsids
popgt <- popgt[popgt[,3] %in% allelecode$POS,] #only keep SNPs in allelecode file (removes esv SNPs)

##join gene info
popexp <- left_join(gencode,exp,by=c("geneid"="id"))

popsamplelist <- colnames(exp)[-1]

#pull gene info & expression from pop of interest
popexp <- dplyr::filter(popexp,chrom==chromosome)
explist <- as.character(popexp$geneid)

resultsarray <- array(0,c(length(explist),19))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975")
dimnames(resultsarray)[[2]] <- resultscol

working100K <- out.dir %&% "working_" %&% pop %&% "_" %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chromosome %&% "_" %&% date %&% ".txt"
write(resultscol,file=working100K,ncolumns=19,sep="\t")

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  start <- popexp$s1[i] - 1e6 ### 1Mb gene lower bound for cis-eQTLS
  end <- popexp$s2[i] + 1e6 ###  1Mb gene upper bound for cis-eQTLs
  cisgenos <- subset(popgt,popgt[,3]>=start & popgt[,3]<=end) ### pull cis-SNP genotypes
  rownames(cisgenos) <- cisgenos$pos #carry positions along
  cismat <- as.matrix(cisgenos[,4:dim(cisgenos)[2]]) #get dosages only in matrix format for glmnet
  expmat <- as.matrix(popexp[,9:dim(popexp)[2]]) #make exp only in matrix format for glmnet
  expmat <- t(expmat) #transpose to match previous code
  colnames(expmat) <- popexp$geneid #carry gene IDs along
  if(dim(cisgenos)[1] > 0){
    cisgenos <- mutate(data.frame(cisgenos),snp=rownames(cisgenos))
    cisbim <- inner_join(allelecode,cisgenos,by=c('POS'='pos')) %>% mutate(SNP=as.character(SNP))
    annotfile <- cbind(cisbim[,3],cisbim[,2],cisbim[,1]) #annot file is rs, pos, chr#
    genofile <- cbind(cisbim[,3:5],cisbim[,8:dim(cisbim)[2]])
    phenofile <- data.frame(expmat[,gene])
    stopifnot(rownames(phenofile)==colnames(genofile)[4:dim(genofile)[2]]) #check for correct sample id order

    write.table(annotfile, file=out.dir %&% "tmp2.annot." %&% chromosome %&% ".s." %&% pop, quote=F, row.names=F, col.names=F, sep=",")
    write.table(genofile, file=out.dir %&% "tmp2.geno." %&% chromosome %&% ".s." %&% pop, quote=F, row.names=F, col.names=F, sep=",")
    write.table(phenofile, file=out.dir %&% "tmp2.pheno." %&% chromosome %&% ".s." %&% pop, quote=F, row.names=F, col.names=F, sep=",")

    runBSLMM <- "gemma -g " %&% out.dir %&% "tmp2.geno." %&% chromosome %&% ".s." %&% pop %&% " -p " %&% out.dir %&% 
      "tmp2.pheno." %&% chromosome %&% ".s." %&% pop %&% " -a " %&% out.dir %&% "tmp2.annot." %&% chromosome %&% ".s." %&% 
      pop %&% " -bslmm 1 -seed 12345 -s 100000 -o tmp2." %&% chromosome %&% ".s." %&% pop
    system(runBSLMM)

    hyp <- read.table(wk.dir %&% "output/tmp2." %&% chromosome %&% ".s." %&% pop %&% ".hyp.txt",header=T)
    hyp50 <- hyp[(dim(hyp)[1]/2+1):dim(hyp)[1],] #take second half of sampling iterations
    quantres <- apply(hyp50,2,getquant)
    res <- c(gene,quantres[1,],quantres[2,],quantres[3,])

  }else{
    res <- c(gene,rep(NA,18))
  }
  names(res) <- c("gene","h50","pve50","rho50","pge50","pi50","n_gamma50","h025","pve025","rho025","pge025","pi025","n_gamma025","h975","pve975","rho975","pge975","pi975","n_gamma975") 
  resultsarray[gene,] <- res
  write(res,file=working100K,ncolumns=19,append=T,sep="\t")
}

write.table(resultsarray,file=out.dir %&% pop %&% "_" %&% tis %&% "_exp_BSLMM-s100K_iterations_chr" %&% chromosome %&% "_" %&% date %&% ".txt",quote=F,row.names=F,sep="\t")
