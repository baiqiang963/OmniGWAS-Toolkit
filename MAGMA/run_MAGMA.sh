#============================================
#usage: sh /home/baiqiang/MAGMA/run_MAGMA.sh GWAS outdir GWAS_build population N txt
#$1 GWAS
#$2 outdir
#$3 GWAS_build：hg19/hg38
#$4 population: EUR/EAS
#$5 GWAS N
#$6 txt/rds
#This script is designed for cell type analysis in FUMA.
#============================================
module load R-4.2.0
if [ "$6" == "txt" ]; then
        Rscript $(dirname "$0")/Format.R $1 $2
elif [ "$6" == "rds" ]; then
        Rscript $(dirname "$0")/Format_rds.R $1 $2
else
        echo "请选择输入文件的类型txt/rds"
    fi

cd $2
#module unload R-4.2.0

module load magma_v1.10
if [ "$3" == "hg19" ]; then
  /home/pub/softwares/magma_v1.10/magma \
    --annotate \
    --snp-loc $(basename "$1").snp_loc \
    --gene-loc /home/baiqiang/MAGMA/NCBI37.3.gene.loc \
    --out $2/$(basename "$1")

elif [ "$3" == "hg38" ]; then
  /home/pub/softwares/magma_v1.10/magma \
    --annotate \
    --snp-loc $(basename "$1").snp_loc \
    --gene-loc /home/baiqiang/MAGMA/NCBI38.gene.loc \
    --out $2/$(basename "$1")
else
  echo "Error: Unsupported genome version '$3'. Supported versions are hg19 and hg38."
  exit 1
fi

if [ "$4" == "EUR" ]; then
  /home/pub/softwares/magma_v1.10/magma \
    --bfile $(dirname "$(dirname "$0")")/ref/EUR \
    --pval $(basename "$1").PVAL_FILE N=$5 \
    --gene-annot $2/$(basename "$1").genes.annot \
    --out $2/$(basename "$1")

elif [ "$4" == "EAS" ]; then
  /home/pub/softwares/magma_v1.10/magma \
    --bfile $(dirname "$(dirname "$0")")/ref/EAS \
    --pval $(basename "$1").PVAL_FILE N=$5 \
    --gene-annot $2/$(basename "$1").genes.annot \
    --out $2/$(basename "$1")
else
  echo "Error: Unsupported genome version $4. Supported versions are EUR and EAS."
  exit 1
fi
cd