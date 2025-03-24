#===================================================
#要求两个性状的GWAS数据拥有'SNP','CHR','BP','BETA','SE'列
#ASSET包引用：Bhattacharjee S, Qi G, Chatterjee N, Wheeler W (2023). ASSET: An R package for subset-based association analysis of heterogeneous traits and subtypes. doi:10.18129/B9.bioc.ASSET
#qqman包引用：Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software, 3(25), 731, https://doi.org/10.21105/joss.00731.
#Usage:module load R-4.2.0
#Usage:Rscript ASSET_easy.R GWAS1的绝对路径 GWAS1的CASE数量 GWAS1的control数量 GWAS2的绝对路径 GWAS2的CASE数量 GWAS2的control数量 输出路径
#===================================================
print('=============================')
print('这是一个将ASSET_meta分析与qqman曼哈顿图结合的ASSET分析脚本,配合sh脚本使用')
print('Usage:Rscript ASSET_easy.R GWAS1的绝对路径 GWAS2的绝对路径 ASSET_casecontrol_table_path 输出路径')
print('=============================')
#清空当前环境
rm(list=ls())
#加载R包
suppressPackageStartupMessages({
  library(ASSET)
  library(data.table)
  library(tools)
  library(cli)
})

#命令行传递
args <- commandArgs(trailingOnly = TRUE)
GWAS1=args[1]
GWAS2=args[2]
ASSET_casecontrol_table_path=args[3]
OutDir=args[4]

cc=fread(ASSET_casecontrol_table_path)
cc$file=basename(cc$file)
CASE1=cc$case[cc$file==basename(GWAS1)]
CNTL1=cc$control[cc$file==basename(GWAS1)]
CASE2=cc$case[cc$file==basename(GWAS2)]
CNTL2=cc$control[cc$file==basename(GWAS2)]
  

EX=GWAS1
OT=GWAS2
#在OutDir目录下创建'GWAS1_GWAS2'文件夹，用以保存结果
setwd(OutDir)
outname=paste0(file_path_sans_ext(basename(EX)),'_',file_path_sans_ext(basename(OT)),"_ASSETres")
if(! dir.exists(outname)){
  dir.create(outname)
}
setwd(outname)
#输入信息输出
cli::cli_h1("ASSET分析")
cli::cli_alert_info('GWAS1:{args[1]}\nGWAS1_CASE:{CASE1}\nGWAS1_CNTL:{CNTL1}\nGWAS2:{args[4]}\nGWAS2_CASE:{CASE2}\nGWAS2_CNTL:{CNTL2}\n输出地址:{OutDir}')


#读取GWAS1
data1=fread(EX)
colnames(data1)=paste('GWAS1',colnames(data1),sep="_")
colnames(data1)[colnames(data1) == "GWAS1_SNP"] <- "SNP"
#读取GWAS2
data2=fread(OT)
colnames(data2)=paste('GWAS2',colnames(data2),sep="_")
colnames(data2)[colnames(data2) == "GWAS2_SNP"] <- "SNP"
#匹配GWAS1和GWAS2的SNP
data=merge(data1,data2,by='SNP')
snps <- as.vector(data[,data$SNP])
nSNP <- length(snps)
print(nSNP)
#匹配GWAS1和GWAS2的性状名，在这里统一为GWAS1和GWAS2
traits.lab <- c('GWAS1','GWAS2')
nTraits<-length(traits.lab)
print(nTraits)
#匹配GWAS1和GWAS2的BETA
beta.hat <- cbind(data$GWAS1_BETA,data$GWAS2_BETA)
colnames(beta.hat)<-c('GWAS1','GWAS2')
#匹配GWAS1和GWAS2的SE
sigma.hat <- cbind(data$GWAS1_SE,data$GWAS2_SE)
colnames(sigma.hat)<-c('GWAS1','GWAS2')
#匹配GWAS1和GWAS2的CASE数量
data$GWAS1_CASE=CASE1
data$GWAS2_CASE=CASE2
ncase<- cbind(data$GWAS1_CASE,data$GWAS2_CASE)
colnames(ncase)<-c('GWAS1','GWAS2')
#部分GWAS内各个SNP的CASE数量不一致，因此默认选择第一行的CASE数量进行匹配
ncase<-ncase[1,1:2]
#匹配GWAS1和GWAS2的Control数量
data$GWAS1_CNTL=CNTL1
data$GWAS2_CNTL=CNTL2
ncntl<- cbind(data$GWAS1_CNTL,data$GWAS2_CNTL)
colnames(ncntl)<-c('GWAS1','GWAS2')
#部分GWAS内各个SNP的Control数量不一致，因此默认选择第一行的CASE数量进行匹配
ncntl=ncntl[1,1:2]
#得到同时包含Concordant和Discordant的结果
res<-h.traits(snp.vars=snps, traits.lab, beta.hat, sigma.hat, ncase, ncntl,
              cor=NULL, cor.numr=FALSE, search=NULL, side=2, meta=FALSE,
              zmax.args=NULL, meth.pval="DLM")
summary=h.summary(res)
Subset.2sided=summary$Subset.2sided
CHR_BP_table=fread(EX)
CHR_BP_table=CHR_BP_table[,c('SNP','CHR','BP')]
Subset.2sided=merge(Subset.2sided, CHR_BP_table, by = "SNP", all.x = TRUE)
saveRDS(Subset.2sided, paste(outname,"_all.rds"))

name33=paste(outname,"_discordant_all.rds",sep="")
name333=paste(outname,"_discordant_all.txt",sep="")
name44=paste(outname,"_concordant_all.rds",sep="")
name444=paste(outname,"_concordant_all.txt",sep="")
#得到包含Discordant的结果
discordant<-Subset.2sided[complete.cases(Subset.2sided$OR.1, Subset.2sided$OR.2), ]
discordant<- discordant[, c("CHR","BP","SNP", "Pvalue",'OR.1','OR.2')]
saveRDS(discordant,name33)
write.table(discordant,name333,quote=F,row.names = F)
#得到包含Concordant的结果
concordant<- Subset.2sided[grepl("GWAS1,GWAS2",Subset.2sided$Pheno.1) | grepl("GWAS1,GWAS2",Subset.2sided$Pheno.2), ]
concordant<- concordant[,c("CHR","BP","SNP", "Pvalue","OR.1","OR.2")]
concordant$OR_concordant <- ifelse(is.na(concordant$OR.1), concordant$OR.2, concordant$OR.1)
concordant$OR=concordant$OR_concordant
write.table(concordant,name444,quote=F,row.names = F)

#从此开始，将利用得到的Con和Disc结果，画出两张曼哈顿图，仅画p>1e-100的位点
require(qqman)

concordant_order<-concordant[order(concordant$CHR,concordant$BP),]
discordant_order<-discordant[order(discordant$CHR,discordant$BP),]
concordant_order=concordant_order[concordant_order$P>1e-100,]
discordant_order=discordant_order[discordant_order$P>1e-100,]
colnames(concordant_order)=c('CHR','BP','SNP','P')
colnames(discordant_order)=c('CHR','BP','SNP','P')
p1name='Concrodant'
p2name='Discrodant'
tiffname=paste0(outname,'_Miami.tiff')
tiff(filename =tiffname)
layout(matrix(c(2, 1), nrow = 2, byrow = TRUE))
p1=manhattan(concordant_order,chr = "CHR",
             bp = "BP",
             p = "P",
             snp = "SNP",
             col = c("gray10", "gray60"),
             chrlabs = NULL,
             suggestiveline =  -log10(1e-5),
             genomewideline = -log10(5e-08),
             highlight = NULL,
             logp = TRUE,
             annotatePval = NULL,
             annotateTop = F)
title(main = p1name)
p2=manhattan(discordant_order,chr = "CHR",
             bp = "BP",
             p = "P",
             snp = "SNP",
             col = c("gray10", "gray60"),
             chrlabs = NULL,
             suggestiveline =  -log10(1e-5),
             genomewideline = -log10(5e-08),
             highlight = NULL,
             logp = TRUE,
             annotatePval = NULL,
             annotateTop = F)
title(main = p2name)
dev.off()

