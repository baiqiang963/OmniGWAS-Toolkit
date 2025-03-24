library(MungeSumstats)
args=commandArgs(T)

GWAS_input=args[1]
out_dir=args[2]
setwd(out_dir)
dat <- readRDS(GWAS_input)
colnames(dat)[which(colnames(dat)=="Pvalue")]="P"
dat1<-dat[,c('SNP','CHR','BP')]
write.table(dat1,paste0(basename(GWAS_input),'.snp_loc'),row.names=F,col.names=F,quote=F,sep='\t')
PVAL_FILE<-dat[,c('SNP','P')]
write.table(PVAL_FILE,paste0(basename(GWAS_input),'.PVAL_FILE'),row.names=F,col.names=F,quote=F,sep='\t')

