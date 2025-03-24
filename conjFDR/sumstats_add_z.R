require(data.table)
args <- commandArgs(trailingOnly = TRUE)
GWAS_path=args[1]
out_folder=args[2]
setwd(out_folder)

dat=fread(GWAS_path)
if (! "Z" %in% colnames(dat) {
# Checking if the 'BETA' and 'SE' columns exist in the data.
if (!all(c("BETA", "SE") %in% colnames(dat))) {
  stop("Error: The data frame is missing either the 'BETA' column or the 'SE' column.\n")
}
dat$Z=as.numeric(data$BETA)/as.numeric(data$SE)
write.table(dat,basename(GWAS_path),quote=F,row.names=F)}

