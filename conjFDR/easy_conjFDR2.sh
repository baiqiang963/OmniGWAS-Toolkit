#========================================================================
#Usage: sh easy_conj2FDR   (FUMA progress)
#========================================================================



cd $save_folder
unzip 'FUMA_*.zip' -d .
if [ "$2" == "conjfdr" ]; then
module load R-4.2.0
Rscript /home/tengli/BQ_GWAS/pleiofdr-master/fuma/conj_fuma_combined_snps.R \
result.clump.lead.csv \
snps.txt \
$trait1_file \
$trait2_file \
$(basename "$trait1_file" | sed 's/\..*$//') \
$(basename "$trait2_file" | sed 's/\..*$//')
fi
