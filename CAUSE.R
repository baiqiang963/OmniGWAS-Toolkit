args<-commandArgs(T)
EX=args[1]
OT=args[2]
r2_thresh=as.numeric(args[3])
pval=as.numeric(args[4])
bfile=args[5]
out_dir=args[6]

suppressPackageStartupMessages({
  library(data.table)
  library(cause)
  library(dplyr)
  library(readr)
  library(genetics.binaRies)
  library(plinkbinr)
  library(cli)
})

setwd(out_dir)
EXP_NAME=tools::file_path_sans_ext(basename(EX))
OT_NAME=tools::file_path_sans_ext(basename(OT))
if (!dir.exists(paste0(EXP_NAME,'_',OT_NAME))) {
  pval1=as.character(pval)
  pval2=gsub("-","",pval1)
  pval2=gsub("0","",pval2)
  folder_path=paste0(EXP_NAME,'_',OT_NAME,"_",pval2)
  dir.create(folder_path)
}

cause_script <- function(folder_path,EX,OT) {
X1<-fread(EX)
X2<-fread(OT)
X <- gwas_merge(X1, X2, snp_name_cols = c("SNP", "SNP"), 
                       beta_hat_cols = c("BETA", "BETA"), 
                       se_cols = c("SE", "SE"), 
                       A1_cols = c("A1", "A1"), 
                       A2_cols = c("A2", "A2"), 
                       pval_cols = c("P", "P"))
set.seed(100)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)
X_clump <- X %>%
           rename(rsid = snp,
                  pval = p1) %>%
           ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval,
                     plink_bin = '/home/tengli/R/x86_64-pc-linux-gnu-library/4.2/plinkbinr/bin/plink_Linux', 
                     bfile = bfile)
top_vars <- X_clump$rsid
res <- cause(X=X,variants= top_vars,param_ests= params)
setwd(folder_path)
sink('MRcause.txt')
remaining_snps <- length(top_vars)
cat("cause分析的SNP数量：", remaining_snps, "\n")
cat('sharing model:\n')
print(res$loos[[2]])
print(loo::pareto_k_table(res$loos[[2]]))
print(res$loos[[3]])
print(loo::pareto_k_table(res$loos[[3]]))
cat('causal model:\n')
print(res$elpd)
print(summary(res, ci_size=0.95))
sink()
pdf(file='CausePlot.pdf')
plot(res)
plot(res, type="data")
dev.off()
}
tryCatch({
  cause_script(folder_path,EX,OT)
},error = function(e) {
  message("An error occurred: ", e$message)
  cli_alert_danger("执行{EXP_NAME}与{OT_NAME}的CAUSE分析发生错误")
}
)


