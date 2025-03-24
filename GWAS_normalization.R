#require environment
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
#BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
#BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")
suppressPackageStartupMessages(
  #因为很多时候Allele1—— the Effect Allele；Allele2——the other Allele，这与此包相反。
#最终增加一步A1和A2调换过程，使A1重现变回fect Allele
    {library(MungeSumstats)
   library(cli)
    library(optparse)
    library(data.table)
    library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
    }
)
option_list <- list(
  make_option(c("-i", "--GWAS_input"), type = "character",
              help = "GWAS_input_path"),
  make_option(c("-c", "--casecontrol_table"), type = "character",
              help = "casecontrol_table_path"),
  make_option(c("-o", "--formatted_GWAS_output"), type = "character",
              help = "formatted_GWAS_output_dir")
)
parser <- OptionParser(usage = " -i GWAS_input_path1,GWAS_input_path2,... -o formatted_GWAS_output_dir", 
                       option_list = option_list)
args <- parse_args(parser)

#args$GWAS_input=""
#args$formatted_GWAS_output="/home/tengli/BQ_GWAS/hg19/test"

# 获取所有 GWAS 文件路径
gwas_files <- unlist(strsplit(args$GWAS_input, ","))
cli_alert_info("发现{length(gwas_files)}个GWAS输入文件")
# 检查是否提供了文件路径
if (length(gwas_files) == 0) {
  stop("未提供 GWAS 文件路径。")
}

for (gwas_file in gwas_files) {
  cli_alert_info("正在处理文件: {gwas_file}")
#检查GWAS汇总数据版本
cli_alert_info("识别GWAS汇总数据版本。。。")
dat <- MungeSumstats::read_sumstats(path = gwas_file, 
                                    standardise_headers = TRUE)
#找到内容中有以"rs"开头的列
rs_cols <- sapply(dat,function(col) {
  any(grepl("^rs", col))
})
rs_col_names <- names(rs_cols[rs_cols])
names(dat)[names(dat) %in% rs_col_names] <- "SNP"

rs_snps <- grep("^rs", dat$SNP, value = TRUE)
if (length(rs_snps) == 0) {
  cli_alert_warning("文件 {gwas_file} 中不存在以 'rs' 开头的 SNP，跳过处理。")
  next
}

# 读取 dpSNP155 GRCh37 参考数据
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37
max_attempts <- 100  # 防止无限循环
attempt <- 1
success <- FALSE

while (!success && attempt <= max_attempts) {
  # 随机挑选一个 SNP
  random_snp <- sample(rs_snps, 1)
  
  # 尝试从参考数据库查询
  snp_dpSNP155 <- tryCatch(
    {
      snpsById(snps, random_snp)  # 主查询
    },
    error = function(e) {
      return(NULL)  # 如果报错返回NULL
    }
  )
  
  # 判断是否查询成功（根据返回类型调整）
  if (!is.null(snp_dpSNP155) && length(snp_dpSNP155) > 0) {
    success <- TRUE
  } else {
    cli_alert_warning("SNP {random_snp} 在dpSNP155中未找到，尝试重新抽取 (已尝试 {attempt} 次)")
    attempt <- attempt + 1
  }
}

if (!success) {
  cli_alert_danger("无法在 {max_attempts} 次尝试内找到有效SNP，请检查数据一致性")
  next
}
snp_row <- dat[dat$SNP == random_snp, ]
# 模糊查找 CHR 所在的列
chr_col <- grep("chr|chrom", colnames(dat), ignore.case = TRUE, value = TRUE)
if (length(chr_col) == 0) {
  cli_alert_warning("文件 {gwas_file} 中未找到包含 'chr' 的列，跳过处理。")
  next
}
CHR_random_snp <- snp_row[, ..chr_col][[1]]
# 模糊查找 BP 所在的列
bp_col <- grep("bp|pos", colnames(dat), ignore.case = TRUE, value = TRUE)
if (length(bp_col) == 0) {
  cli_alert_warning("文件 {gwas_file} 中未找到包含 'bp' 的列，跳过处理。")
  next
}
BP_random_snp <- snp_row[, ..bp_col][[1]]

# 对比 GWAS 与 dpSNP155 GRCh37 参考数据中的 random_snp 的 CHR 和 BP 是否一致
if (CHR_random_snp == as.character(seqnames(snp_dpSNP155)) & BP_random_snp == pos(snp_dpSNP155)) {
  ref_genome <- "GRCh37"
} else {
  ref_genome <- "GRCh38"
}
cli_alert_success("识别 GWAS 汇总数据版本: {ref_genome}")
#添加N列
cc=fread(args$casecontrol_table)
dat$N=cc$total[cc$file==args$GWAS_input] 
cli_alert_info("添加N列:{dat$N[1]}")
cli_alert_info("预处理前：\n{head(dat)}")
#如果数据自带A1,A2。那么将A1定为alt，A2定为ref 放错机制
if ("A1" %in% colnames(dat)){
  colnames(dat)[which(names(dat)=="A1")]="alt"
}
if ("A2" %in% colnames(dat)){
  colnames(dat)[which(names(dat)=="A2")]="ref"
}
#检测可使用的最大处理器线程
max_threads <- parallel::detectCores()
# 自动进行版本转化和预处理
if (ref_genome == "GRCh38") {
  cli_alert_info("正在将 GRCh38 转换为 GRCh37")
  dat <- format_sumstats(dat,
                         impute_beta = TRUE,
                         ref_genome = ref_genome,
                         convert_ref_genome = "GRCh37",
                         compute_z = TRUE,
                         nThread = max_threads,
                         return_data = TRUE)
} else {
  dat <- format_sumstats(dat,
                         impute_beta = TRUE,
                         ref_genome = ref_genome,
                         compute_z = TRUE,
                         nThread = max_threads,
                         return_data = TRUE)
}
if (!dir.exists(args$formatted_GWAS_output)) {
  dir.create(args$formatted_GWAS_output, recursive = TRUE)
}
#将A1，A2转换回来。A1为effect，A2为reference
dat <- dat[, c("A2", "A1") := .(A1, A2)]
cli_alert_info("A1，A2转换后：A1为effect allele，A2为other allele")
print(head(dat))
# 输出格式化后的 GWAS 数据
if (!"BETA" %in% colnames(dat)){
  dat$BETA=log(dat$OR)
}
head(dat)
dat=dat[,c("SNP","CHR","BP","A1","A2","BETA","SE","P","N","Z")]
fwrite(dat,file.path(args$formatted_GWAS_output,basename(gwas_file)),row.names = F,quote = F,sep=" ")
cli_alert_success("GWAS 汇总数据输出成功: {file.path(args$formatted_GWAS_output,basename(gwas_file))}")}

