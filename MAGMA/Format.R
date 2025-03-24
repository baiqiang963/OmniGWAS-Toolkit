library(MungeSumstats)
args=commandArgs(T)

GWAS_input=args[1]
out_dir=args[2]
setwd(out_dir)
dat <- MungeSumstats::read_sumstats(path = GWAS_input, 
                                    standardise_headers = TRUE)
MungeSumstats::write_sumstats(sumstats_dt = dat,
                                                save_path = basename(GWAS_input),
                                                tabix_index = TRUE,
                                                write_vcf = FALSE,
                                                return_path = TRUE)
dat1<-dat[,c('SNP','CHR','BP')]
write.table(dat1,paste0(basename(GWAS_input),'.snp_loc'),row.names=F,col.names=F,quote=F,sep='\t')
PVAL_FILE<-dat[,c('SNP','P')]
write.table(PVAL_FILE,paste0(basename(GWAS_input),'.PVAL_FILE'),row.names=F,col.names=F,quote=F,sep='\t')

