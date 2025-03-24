#===============================================
#Automatically add proxy SNPs, automate MR-PRESSO, external command line written by Baiqiang Xue.
#usage: module load R-4.2.0
#usage: Rscript TwoSampleMR_easy.R 暴露数据绝对路径 结局数据绝对路径 人群:EUR/EAS/... 工具变量阈值,例如5e-8/1e-7/1e-6/1e-5 是否关闭去除连锁不平衡,t/f(关/开) 去除连锁不平衡的r2,例如0.01 F统计值公式选择,1/2 是否查找并补充代理SNP,t/f 保存路径
#F统计值公式1:F=R^2(N−2)/(1−R^2),参考文献：Shim H, et al. A multivariate genome-wide association analysis of 10 LDL subfractions, and their response to statin treatment, in 1868 caucasians. PLoS ONE. 2015;10:e0120758. doi: 10.1371/journal.pone.0120758
#F统计值公式2:F=(beta/se)^2,参考文献：Hemani G., Zheng J., Elsworth B., Wade K.H., Haberland V., Baird D., Laurin C., Burgess S., Bowden J., Langdon R., et al. The MR-Base platform supports systematic causal inference across the human phenome. eLife.2018;7:e34408. doi: 10.7554/eLife.34408
#代理SNP基于hg19，如需更改，请在代码中更改。
#
#如果需要将此脚本上传fat节点运行，请保证proxy="f"
#===============================================
print('===============================================')
print('欢迎使用TwoSampleMR_easy')
print('usage: module load R-4.2.0')
print('usage: Rscript TwoSampleMR_easy.R 暴露数据绝对路径 结局数据绝对路径 LD参考人群:EUR/EAS/... LDr2例如0.01 F统计值公式选择,1/2 保存路径')
print('请保证暴露和结局的GWAS中snp_col ,beta_col,se_col,effect_allele_col,other_allele_col,pval_col,samplesize_col的列名必须是SNP，BETA，SE，A1，A2，P，N')
print('eaf是非必要的。如果存在eaf列，请保证它的列名为EAF。')
print('如果运行中发生错误，极大概率是工具变量不足，此时可以考虑将工具变量阈值放宽')
print('===============================================')
#加载R包
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(MRPRESSO)
  library(dplyr)
  library(LDlinkR)
  library(cli)
})

#命令行传递
args <- commandArgs(trailingOnly = TRUE)
EXP=args[1]
OT=args[2]
GWAS_population=args[3]
pval=as.numeric(args[4])
clump_r2=as.numeric(args[5])
F_statistic=as.numeric(args[6])
out_dir=args[7]

proxy="f"
print(paste('暴露数据绝对路径',EXP,sep=":"))
print(paste('结局数据绝对路径',OT,sep=":"))
print(paste('人群参考面板',GWAS_population,sep=":"))
print(paste('工具变量阈值',pval,":"))
print(paste('去除连锁不平衡r2<',clump_r2,sep=""))
if (F_statistic==1){
  print(paste('F统计值公式','F=R^2(N−2)/(1−R^2)',sep=":"))}else{
  print(paste('F统计值公式','F=(beta/se)^2',sep=":"))
}
print(paste('是否查找并补充代理SNP',proxy,sep=":"))
print(paste('保存路径',out_dir,":"))

setwd(out_dir) # 结果保存路径
EXP_NAME=tools::file_path_sans_ext(basename(EXP))
OT_NAME=tools::file_path_sans_ext(basename(OT))
if (!dir.exists(paste0(EXP_NAME,'_',OT_NAME))) {
  pval1=as.character(pval)
  pval2=gsub("-","",pval1)
  pval2=gsub("0","",pval2)
  XXXX=paste0(EXP_NAME,'_',OT_NAME,"_",pval2)
  dir.create(XXXX)
}
setwd(XXXX)
tryCatch({
#输出参数文件
sink('args.txt')
print('暴露数据绝对路径')
print(args[1])
print('结局数据绝对路径')
print(args[2])
print('GWAS研究/LD面板选择的人群')
print(args[3])
print('工具变量阈值')
print(args[4])
print('去除连锁不平衡r2<')
print(args[5])
print('F统计值公式')
print(args[6])
cli_alert_info("输出路径:{XXXX}")
sink()

#读取暴露文件绝对路径，并格式化
exp=fread(EXP)
if ("OR" %in% colnames(exp)){
exp$BETA=log(as.numeric(exp$OR))}
#SNP-性状关联分析，找到性状的差异性SNP。5e-8值最好，1e-5宽泛   5e-8做不出来，再改成1e-5。
exp=exp[which(exp$P<pval),] 
#格式化暴露数据
exp_dat <- format_data(
  exp,
  type='exposure',
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="A1",
  other_allele_col = "A2",
  pval_col = "P",
  samplesize_col = "N")
#自定义暴露性状的名字，默认去除后缀的文件名
exp_dat$exposure<-EXP_NAME
#筛选SNP列为rs开头的SNP或行
index <- grepl("^rs", exp_dat$SNP)
exp_dat<-exp_dat[index,]
#去除连锁不平衡
  exposure_all_dat <- exp_dat %>%
    rename(rsid = SNP,
           pval = pval.exposure) %>%
    ieugwasr::ld_clump(dat = .,
                       clump_kb = 10000,
                       clump_r2 = clump_r2,
                       clump_p = 1,
                       plink_bin = '/home/tengli/R/x86_64-pc-linux-gnu-library/4.2/plinkbinr/bin/plink_Linux', 
                       bfile = GWAS_population)%>%
    rename(SNP = rsid,
           pval.exposure = pval)
#格式化结局数据
outcome_all_dat<-read_outcome_data(snps=exposure_all_dat$SNP,
                                   filename=OT,
                                   snp_col ="SNP",  
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col ="A1",
                                   other_allele_col = "A2",
                                   pval_col = "P",
                                   chr_col = "CHR",
                                   pos_col ="BP")
#自定义结局性状的名字
outcome_all_dat$outcome<-OT_NAME
#寻找存在于暴露但结局丢失的工具变量
missing_snp<-anti_join(exposure_all_dat, outcome_all_dat, by = "SNP")
#启动代理查找时，如果发现结局丢失的工具变量>0,则自动联网补充代理工具变量/SNP。
if (proxy=='t') {
  if (nrow(missing_snp)>0) {
    outcome_dat=read_outcome_data(
      filename=OT,
      snp_col ="SNP",  
      beta_col = "BETA",
      se_col = "SE",
      effect_allele_col ="A1",
      other_allele_col = "A2",
      pval_col = "P",	
      samplesize_col = "N")
    outcome_dat$outcome<-OT_NAME
    outcome_dat$id.outcome<-outcome_all_dat$id.outcome[1]
    missing_snp_number=nrow(missing_snp)
    print(paste('Number of IVs missing from outcome GWAS:',missing_snp_number))
    LDproxy_batch(missing_snp$SNP, 
                  pop = GWAS_population, 
                  r2d = "r2", 
                  token = 'a7bf3ccf48cf', 
                  append = TRUE,
                  genome_build = "grch37")
    proxies_list <- fread('combined_query_snp_list_grch37.txt',header=F)
    proxies_list<- proxies_list[,c(2:3,9)]
    colnames(proxies_list)<-c('query_snp','SNP','R2')
    proxies_list<-proxies_list[proxies_list$R2>0.8,]
    proxies_list=merge(proxies_list,outcome_dat,by = c("SNP"),all = FALSE)
    grouped_proxies_list <- proxies_list %>%
      group_by(query_snp)%>%
      slice_max(order_by = R2, n = 1)%>%
      slice_head(n=1)%>%
      ungroup()
    proxies_list=grouped_proxies_list[,-1]
    colnames(proxies_list[1])='SNP'
    write.table(proxies_list,'proxies_list.txt',quote=F,row.names=F)
    proxies_list=proxies_list[,-2]
    outcome_all_dat<-rbind(outcome_all_dat,proxies_list)
    #打印成功添加代理SNP的数量
    print(paste0('Proxy SNP added successfully:',nrow(proxies_list)))
    #打印结局仍然丢失SNP的数量
    print(paste0('Still missing:',missing_snp_number-nrow(proxies_list)))}
}
#暴露与结局数据协调
dat <- harmonise_data(exposure_all_dat, outcome_all_dat)
#计算F统计值
if (F_statistic==1){
  dat=cbind(dat,Rsq=dat$beta.exposure^2/(dat$beta.exposure^2+dat$samplesize.exposure*dat$se.exposure^2),  F=(dat$beta.exposure^2/(dat$beta.exposure^2+dat$samplesize.exposure*dat$se.exposure^2))*(dat$samplesize.exposure-2)/(1-(dat$beta.exposure^2/(dat$beta.exposure^2+dat$samplesize.exposure*dat$se.exposure^2))))}else{
    dat=cbind(dat,F=(dat$beta.exposure/dat$se.exposure)^2)
  }
#保留F>10的SNP
dat=dat[dat$F>10,]
#提取去除回文序列SNP后的数据
dat=subset(dat,dat$mr_keep!=FALSE)
#打印MRpresso之前的工具变量数量
print(paste0('number of IVs without MR-presso:',length(dat$SNP)))
#如果工具变量为0，则跳过当前迭代，进入下一次迭代
if (nrow(dat)==0){
  sink('Error.txt')
  print(paste0('Number of SNPs below the p-threshold:',nrow(exp_dat)))
  print(paste0('Number of SNPs after clumping:',nrow(exposure_all_dat)))
  if (proxy=='t'){
    print(paste('Number of IVs missing from outcome GWAS:',missing_snp_number))
    print(paste0('Proxy SNP added successfully:',nrow(proxies_list)))
    print(paste0('Still missing:',missing_snp_number-nrow(proxies_list)))
  }
  cli_alert_danger('Entry into MR analysis without instrumental variable')
  print(paste0('Automatically proceed to the MR analysis of the next exposure with ',OT_NAME))
  sink()
  next
}
#输出工具变量列表
write.table(dat,"harmonise.txt",quote=FALSE,sep="\t",row.names=FALSE)
#进行MR-PRESSO分析
Presso=mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure",OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 3000,SignifThreshold = 0.05) 
sink("MRpresso.txt")
print(Presso)
sink()
#判断MR-PRESSO是否发现水平多效性/异质性工具变量,如果未发现,则直接进行MR分析并输出结果；如果发现，则剔除异常工具变量后，进行MR分析并输出结果。
if(any(is.na(Presso$`Main MR results`[2, 3:6]))){
  print("MR-Presso NO outlier snp")
  
  mr_data_input <- mr_input(
    bx = dat$beta.exposure, 
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome, 
    snps = dat$SNP, 
    effect_allele = dat$effect_allele.exposure, 
    other_allele = dat$other_allele.exposure,
    eaf = dat$eaf.exposure
  )
  
  
  mr_feivw_mendrand <- MendelianRandomization::mr_ivw(mr_data_input, model = "fixed")
  mr_reivw_mendrand <- MendelianRandomization::mr_ivw(mr_data_input, model = "random")
  mr_egger_mendrand <- MendelianRandomization::mr_egger(mr_data_input)
  mr_mode_mendrand <- MendelianRandomization::mr_mbe(mr_data_input)
  mr_median_mendrand <- MendelianRandomization::mr_median(mr_data_input)
  res <- mr(dat)
  Presso_result=Presso$`Main MR results`[1, c(3:4,6)]
  colnames(Presso_result)=c('b','se','pval')
  Presso_result$method='MRPRESSO'
  Presso_result=cbind(res[1,c(1:4,6)],Presso_result)
  res=rbind(Presso_result,res)
  res=res[c('id.exposure','id.outcome','exposure','outcome','method','nsnp','b','se','pval')]
  heterogeneity<-mr_heterogeneity(dat)
  pleiotropy<-mr_pleiotropy_test(dat)
  
  sink("result.txt")
  print(mr_feivw_mendrand)
  print(mr_reivw_mendrand)
  print(mr_median_mendrand)
  print(mr_egger_mendrand)
  print(mr_mode_mendrand)
  print(heterogeneity)
  print(pleiotropy)
  print(res)
  sink()
  
  pdf(file="MRplot.pdf")
  print(mr_scatter_plot(res, dat))
  print(mr_forest_plot(mr_singlesnp(dat)))
  print(mr_funnel_plot(mr_singlesnp(dat)))
  print(mr_leaveoneout_plot(mr_leaveoneout(dat)))
  print(mr_plot(mr_data_input, line = "ivw", interactive = FALSE, labels = FALSE,orientate = TRUE))
  dev.off()
}else{
  print("MR-Presso found outlier snp and removing")
  OutliersIndices <- Presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  dat=dat[-OutliersIndices,]
  write.table(dat,"harmonise_RMoutliers.txt",quote=F,sep="\t",row.names=F)
  print(paste0('MR-Presso find outlier snps and remove them!!!:',OutliersIndices))
  mr_data_input <- mr_input(
    bx = dat$beta.exposure, 
    bxse = dat$se.exposure,
    by = dat$beta.outcome,
    byse = dat$se.outcome, 
    snps = dat$SNP, 
    effect_allele = dat$effect_allele.exposure, 
    other_allele = dat$other_allele.exposure,
    eaf = dat$eaf.exposure
  )
  
  
  mr_feivw_mendrand <- MendelianRandomization::mr_ivw(mr_data_input, model = "fixed")
  mr_reivw_mendrand <- MendelianRandomization::mr_ivw(mr_data_input, model = "random")
  mr_egger_mendrand <- MendelianRandomization::mr_egger(mr_data_input)
  mr_mode_mendrand <- MendelianRandomization::mr_mbe(mr_data_input)
  mr_median_mendrand <- MendelianRandomization::mr_median(mr_data_input)
  res <- mr(dat)
  Presso_result=Presso$`Main MR results`[2, c(3:4,6)]
  colnames(Presso_result)=c('b','se','pval')
  Presso_result$method='MRPRESSO'
  Presso_result=cbind(res[1,c(1:4,6)],Presso_result)
  res=rbind(Presso_result,res)
  res=res[c('id.exposure','id.outcome','exposure','outcome','method','nsnp','b','se','pval')]
  heterogeneity<-mr_heterogeneity(dat)
  pleiotropy<-mr_pleiotropy_test(dat)
  
  sink("result.txt")
  print(mr_feivw_mendrand)
  print(mr_reivw_mendrand)
  print(mr_median_mendrand)
  print(mr_egger_mendrand)
  print(mr_mode_mendrand)
  print(heterogeneity)
  print(pleiotropy)
  print(res)
  sink()
  
  pdf(file="MRplot.pdf")
  print(mr_scatter_plot(res, dat))
  print(mr_forest_plot(mr_singlesnp(dat)))
  print(mr_funnel_plot(mr_singlesnp(dat)))
  print(mr_leaveoneout_plot(mr_leaveoneout(dat)))
  print(mr_plot(mr_data_input, line = "ivw", interactive = FALSE, labels = FALSE,orientate = TRUE))
  dev.off()
}
#以excel形式输出结果，便于整理
res$OR=exp(res$b)
res$'OR_lower_95%CI'=exp(res$b-1.96*res$se)
res$'OR_upper_95%CI'=exp(res$b+1.96*res$se)
res$egger_intercept_p=pleiotropy$pval
writexl::write_xlsx(res,'result.xlsx')}, error = function(e) {
  message("An error occurred: ", e$message)
  cli_alert_danger("以工具变量阈值{pval}运行{EXP_NAME}与{OT_NAME}出现错误")
  #清空当前环境
  rm(list=ls())
  })
