require(data.table)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
trait1=args[1]
trait2=args[2]
ASSET_casecontrol_table_path=args[3]
save_folder=args[4]
if (!file.exists(file.path(save_folder,basename(trait1)))&!file.exists(file.path(save_folder,basename(trait2)))){
cc=fread(ASSET_casecontrol_table_path)
cc$file=basename(cc$file)
CASE1=cc$case[cc$file==basename(trait1)]
CNTL1=cc$control[cc$file==basename(trait1)]
CASE2=cc$case[cc$file==basename(trait2)]
CNTL2=cc$control[cc$file==basename(trait2)]

case_trait1=CASE1
control_trait1=CNTL1
case_trait2=CASE2
control_trait2=CNTL2


setwd(save_folder)

dat1=fread(trait1)
dat1$N = as.double(dat1$N)
dat1$N=4/(1/case_trait1 + 1/control_trait1)
dat1$Z=dat1$BETA/dat1$SE
dat1 <- dat1[!(dat1$CHR == 6 & dat1$BP >= 26000000 & dat1$BP <= 34000000), ]
dat1=dat1[,c('SNP', 'CHR', 'BP', 'A1', 'A2', 'N', 'Z')]
write.table(dat1,paste0(file_path_sans_ext(basename(trait1)),'.gz'),sep = "\t",quote = FALSE,row.names = FALSE)

dat2=fread(trait2)
dat2$N = as.double(dat2$N)
dat2$N=4/(1/case_trait2 + 1/control_trait2)
dat2$Z=dat2$BETA/dat2$SE
dat2 <- dat2[!(dat2$CHR == 6 & dat2$BP >= 26000000 & dat2$BP <= 34000000), ]
dat2=dat2[,c('SNP', 'CHR', 'BP', 'A1', 'A2', 'N', 'Z')]
write.table(dat2,paste0(file_path_sans_ext(basename(trait2)),'.gz'),sep = "\t",quote = FALSE,row.names = FALSE)}else{
  print("已存在针对MixeR格式化完成的GWAS")
}

