#!/bin/bash
# 显示帮助信息
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  --GWAS1 PATH                GWAS1文件路径（多个用逗号分隔）"
    echo "  --GWAS2 PATH                GWAS2文件路径（多个用逗号分隔）"
    echo "  --formattedGWAS_path PATH   预处理数据存储路径（必需）"
    echo "  --output_path PATH         输出路径（必需）"
    echo "  --casecontrol_table_path PATH 样本量表路径"
    echo "  --Ref [EUR|EAS|AFR|AMR|SAS] 参考人群（默认：EUR）"
    echo "  --IV_P VALUE               IV P值阈值（默认：5e-8）"
    echo "  --LD_R2 VALUE              LD R2阈值（默认：0.001）"
    echo "  --analysis NAME1 NAME2 ... 指定要执行的分析（可多选）"
    echo "                              可选分析：${!ANALYSIS[*]}"
    echo "  --help                     显示帮助信息"
    exit 0
}
##-----------------------------------------------------------------------------------------------------------------------
#加载模块
module load R-4.2.0
module load anaconda/Anaconda2
#LDSC软件所在路径
s_LDSC=/home/pub/softwares/ldsc
#pleiofdr-master所在路径
s_pleiofdr=/home/tengli/BQ_GWAS/pleiofdr-master
#python_convert-master所在路径
s_python_convert=/home/tengli/BQ_GWAS/python_convert-master
#R脚本所在的目录相同
SCRIPT_DIR="$(dirname "$0")"
##各个分析对应的脚本
Formatted_RSCRIPT=$SCRIPT_DIR/GWAS_normalization.R
HDL_RSCRIPT=$SCRIPT_DIR/HDL-master/HDL.run.R
LDSC_PYSCRIPT1=$s_LDSC/munge_sumstats.py
LDSC_PYSCRIPT2=$s_LDSC/ldsc.py
ASSET_RSCRIPT=$SCRIPT_DIR/ASSET_easy.R
conjFDR_SCRIPT=$SCRIPT_DIR/conjFDR/easy_conjFDR.sh
GSMR_RSCRIPT=$SCRIPT_DIR/GSMR.R
TwoSampleMR_RSCRIPT=$SCRIPT_DIR/TwoSampleMR_easy_10000.R
CAUSE_RSCRIPT=$SCRIPT_DIR/CAUSE.R
MiXeR_SCRIPT=$SCRIPT_DIR/MixeR/mixeR
MiXeR_PAR_SCRIPT=$SCRIPT_DIR/MixeR/mixeR_GWAS_convert.R
MAGMA_SCRIPT=$SCRIPT_DIR/MAGMA/run_MAGMA.sh

##参数
#GWAS1和GWAS2均允许输入单个或多个以","间隔的GWAS文件的路径。在因果推断中，GWAS1代表暴露，GWAS2代表结局。
#Example:GWAS1=./Genetic/Example/Epilepsy.all_2022ILAE.gz,./Genetic/Example/XBIO.gz
#Example:GWAS2=./Genetic/Example/Schizophrenia_PGC3.ALL.gz,./Genetic/Example/BBA.gz
#formattedGWAS_path代表流程预处理的GWAS保存结果路径
#Example:formattedGWAS_path=./Genetic/Example/formattedGWAS/
#output_path分析结果输出路径
#casecontrol_table_path:ASSET和MixeR分析的GWAS样本量分组文件路径(自建四列的csv或txt文件，file列填写来自GWAS1和GWAS2输入路径，对应的case列和control列以及total列分别填写它们的案例数和对照数以及样本总数)
#Example:./Genetic/Example/ASSET_CC.csv
#Ref:LDSC，TwoSampleMR，GSMR分析的参考面板人群(EUR/EAS/AFR/AMR/SAS)[欧洲人/东亚人/..........]
#TwoSampleMR,GSMR,CAUSE工具变量参数:IV_P,LD_R2

# 参数默认值
SCRIPT_DIR="$(dirname "$0")"
GWAS1=""
GWAS2=""
formattedGWAS_path=""
output_path=""
casecontrol_table_path=""
Ref="EUR"
IV_P="5e-8"
LD_R2="0.001"

##-----------------------------------------------------------------------------------------------------------------------
##选择需要执行的分析（"Analysis"=T，执行对应分析）
#HDL(使用的UKB的参考面板，因此仅支持欧洲人)
#LDSC
#ASSET
#TwoSampleMR
#GSMR(Generalised Summary-data-based Mendelian Randomisation v2)
#CAUSE
#MiXeR（分析时，在$MiXeR_PAR_SCRIPT脚本中自动移除chr6：26000000-34000000的MHC区域）
#MAGMA
# 初始化所有分析开关为 F（默认不执行）
declare -A ANALYSIS=(
    [HDL]=F
    [LDSC]=F
    [ASSET]=F
    [TwoSampleMR]=F
    [GSMR]=F
    [CAUSE]=F
    [MiXeR]=F
    [MAGMA]=F
	[conjFDR]=F
)

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        --GWAS1) GWAS1="$2"; shift 2 ;;
        --GWAS2) GWAS2="$2"; shift 2 ;;
        --formattedGWAS_path) formattedGWAS_path="$2"; shift 2 ;;
        --output_path) output_path="$2"; shift 2 ;;
        --casecontrol_table_path) casecontrol_table_path="$2"; shift 2 ;;
        --Ref) Ref="$2"; shift 2 ;;
        --IV_P) IV_P="$2"; shift 2 ;;
        --LD_R2) LD_R2="$2"; shift 2 ;;
        --analysis)
            shift
            while [[ $# -gt 0 && ! "$1" == --* ]]; do
                if [[ -n "${ANALYSIS[$1]:-}" ]]; then
                    ANALYSIS["$1"]=T
                else
                    echo "错误：未知分析名称 '$1'，有效选项：${!ANALYSIS[*]}" >&2
                    exit 1
                fi
                shift
            done
            ;;
        --help) show_help ;;
        *)
            echo "错误：未知参数 $1" >&2
            exit 1
            ;;
    esac
done

# 检查必需参数
[[ -z "$formattedGWAS_path" ]] && { echo "错误：必须指定 --formattedGWAS_path" >&2; exit 1; }
[[ -z "$output_path" ]] && { echo "错误：必须指定 --output_path" >&2; exit 1; }
##----------------------------------------------------------------------------------------------------
HDL=${ANALYSIS[HDL]}
LDSC=${ANALYSIS[LDSC]}
ASSET=${ANALYSIS[ASSET]}
TwoSampleMR=${ANALYSIS[TwoSampleMR]}
GSMR=${ANALYSIS[GSMR]}
CAUSE=${ANALYSIS[CAUSE]}
MVMR=${ANALYSIS[MVMR]}
MiXeR=${ANALYSIS[MiXeR]}
MAGMA=${ANALYSIS[MAGMA]}
conjFDR=${ANALYSIS[conjFDR]}

#输出路径系列
mkdir $formattedGWAS_path
mkdir $output_path
#将原始输入路径转化为格式化后的路径
# 将 GWAS1 拆分为数组
gwas1_array=()
IFS=',' read -r -a gwas1_array <<< "$GWAS1"

# 初始化 GWAS1F
GWAS1F=""
# 循环遍历 GWAS1 数组
for gwas1 in "${gwas1_array[@]}"; do
  # 提取文件名
  filename=$(basename "$gwas1")
  
  # 组合新路径
  new_path="$formattedGWAS_path/$filename"
  
  # 将新路径添加到 GWAS1F
  if [[ -z "$GWAS1F" ]]; then
    GWAS1F="$new_path"
  else
    GWAS1F="$GWAS1F,$new_path"
  fi
done

# 将 GWAS2 拆分为数组
gwas2_array=()
IFS=',' read -r -a gwas2_array <<< "$GWAS2"

# 初始化 GWAS2F
GWAS2F=""
# 循环遍历 GWAS1 数组
for gwas2 in "${gwas2_array[@]}"; do
  # 提取文件名
  filename=$(basename "$gwas2")
  
  # 组合新路径
  new_path="$formattedGWAS_path/$filename"
  
  # 将新路径添加到 GWAS2F
  if [[ -z "$GWAS2F" ]]; then
    GWAS2F="$new_path"
  else
    GWAS2F="$GWAS2F,$new_path"
  fi
done

# 将 GWAS1F 和 GWAS2F 拆分为数组
IFS=',' read -r -a gwas1_array <<< "$GWAS1F"
IFS=',' read -r -a gwas2_array <<< "$GWAS2F"
##预处理
#检查预处理结果文件夹，判断GWAS数据是否存在预处理过的版本
#针对没有发现预处理过的版本，执行原始GWAS数据的预处理
if [[ -z "$GWAS1" && -z "$GWAS2" ]]; then
  echo "错误：未发现GWAS1和GWAS2中文件！"
    exit 1  # 非零退出码表示错误
elif [[ -z "$GWAS1" ]]; then
  GWASs="$GWAS2"
elif [[ -z "$GWAS2" ]]; then
  GWASs="$GWAS1"
else
  GWASs="$GWAS1,$GWAS2"
fi
gwas_array=()
IFS=',' read -r -a gwas_array <<< "$GWASs"
declare -A gwas_files
for gwas in "${gwas_array[@]}"; do
  filename=$(basename "$gwas")
  gwas_files["$filename"]="$gwas"
done

# 检查 formattedGWAS_path 中是否存在这些文件
missing_files=()  # 用于存储不存在的文件
missing_paths=()  # 用于存储不存在的文件的原路径

for filename in "${!gwas_files[@]}"; do
  if [[ ! -f "$formattedGWAS_path/$filename" ]]; then
    missing_files+=("$filename")
    missing_paths+=("${gwas_files["$filename"]}")
  fi
done

# 根据检查结果执行操作
if [[ ${#missing_files[@]} -eq 0 ]]; then
  echo "所有文件已存在，跳过操作。"
else
  echo "以下文件不存在，执行预处理操作:"
  for ((i=0; i<${#missing_files[@]}; i++)); do
    filename="${missing_files[$i]}"
    original_path="${missing_paths[$i]}"
	echo "- 文件名: $filename, 原路径: $original_path"
    Rscript $Formatted_RSCRIPT \
	-i $original_path \
	-c $casecontrol_table_path \
	-o $formattedGWAS_path
  done
fi

##分析
Ref_path=$SCRIPT_DIR/ref/$Ref
#嵌套循环遍历数组：外层GWAS1数组，内层GWAS2数组，适用于循环进行MVMR分析之外的分析
for gwas1 in "${gwas1_array[@]}"; do
  filename1=$(basename "$gwas1" | sed 's/\.[^.]*$//')
  for gwas2 in "${gwas2_array[@]}"; do
    filename2=$(basename "$gwas2" | sed 's/\.[^.]*$//')
	echo "正在分析: trait1=$filename1, trait2=$filename2"
	#HDL分析	
    if [[ "$HDL" == "T" && -f "$HDL_RSCRIPT" ]]; then
      echo "条件满足，执行HDL：$gwas1 -> $gwas2"
	  module load R-4.2.0
      Rscript "$HDL_RSCRIPT" \
      gwas1.df=$gwas1 \
      gwas2.df=$gwas2 \
	  LD.path=$SCRIPT_DIR/HDL-master/Reference_panels/HapMap3 \
      output.file="$output_path/${filename1}_${filename2}_HDL.txt"
    else
      echo "$HDL_RSCRIPT条件不满足，跳过执行。"
	fi
	#LDSC分析
	if [[ "$LDSC" == "T" && -f "$LDSC_PYSCRIPT1" && -f "$LDSC_PYSCRIPT2" ]]; then
      echo "条件满足，正在执行 $LDSC：$gwas1 -> $gwas2"
	  if [[ "$Ref" == "EUR" ]]; then
        ld_chr=$SCRIPT_DIR/ref/LDSC/eur_w_ld_chr/
      elif [[ "$Ref" == "EAS" ]]; then
        ld_chr=$SCRIPT_DIR/ref/LDSC/eas_ldscores/
	  fi
	  python $LDSC_PYSCRIPT1 \
	  --sumstats $gwas1 \
	  --signed-sumstats BETA,0 \
	  --N-col N \
	  --a1 A1 \
      --a2 A2 \
	  --p P \
	  --out $output_path/$filename1.ldsc \
	  --merge-alleles $SCRIPT_DIR/ref/LDSC/w_hm3.snplist 
	  python $LDSC_PYSCRIPT1 \
	  --sumstats $gwas2 \
	  --signed-sumstats BETA,0 \
	  --N-col N \
	  --a1 A1 \
      --a2 A2 \
	  --p P \
	  --out $output_path/$filename2.ldsc \
	  --merge-alleles $SCRIPT_DIR/ref/LDSC/w_hm3.snplist
      python $LDSC_PYSCRIPT2 \
      --ref-ld-chr $ld_chr \
	  --w-ld-chr $ld_chr \
      --rg $output_path/$filename1.ldsc.sumstats.gz,$output_path/$filename2.ldsc.sumstats.gz \
      --out ${output_path}/${filename1}_${filename2}_LDSC
    else
      echo "$LDSC_PYSCRIPT条件不满足，跳过执行。"
    fi
	#TwoSampleMR
	if [[ "$TwoSampleMR" == "T" && -f "$TwoSampleMR_RSCRIPT" ]]; then
	  echo "条件满足，正在执行TwoSampleMR：$gwas1 -> $gwas2"
	  module load R-4.2.0
	  Rscript "$TwoSampleMR_RSCRIPT" \
	  $gwas1 \
	  $gwas2 \
	  $SCRIPT_DIR/ref/$Ref \
	  $IV_P \
	  $LD_R2 \
	  2 \
	  $output_path
	else
      echo "$TwoSampleMR_RSCRIPT条件不满足，跳过执行。"
	fi
	#GSMR
	IV_PP=$(echo "$IV_P" | tr -d '-' | tr '[:upper:]' '[:lower:]')
	HM_RMoutliers="$output_path/${filename1}_${filename2}_${IV_PP}/harmonise_RMoutliers.txt"
    HM_harmonise="$output_path/${filename1}_${filename2}_${IV_PP}/harmonise.txt"
	if [[ "$GSMR" == "T" && ( -f "$HM_RMoutliers" || -f "$HM_harmonise" ) ]]; then
	  echo "条件满足，正在执行GSMR：$gwas1 -> $gwas2"
	  if [[ -f "$HM_RMoutliers" ]]; then
        HM="$HM_RMoutliers"
      elif [[ -f "$HM_harmonise" ]]; then
        HM="$HM_harmonise"
	  fi
	  Rscript $SCRIPT_DIR/get_eaf_from_1000G_for_GSMR.R \
	  $SCRIPT_DIR \
	  $HM
	  Rscript "$GSMR_RSCRIPT" \
	  $HM \
	  $IV_P \
	  $LD_R2
	  module load gcta_1.93
	  gcta64 \
	  --bfile $Ref_path \
	  --extract $output_path/${filename1}_${filename2}_${IV_PP}/gsmr_example_snps.allele \
	  --update-ref-allele $output_path/${filename1}_${filename2}_${IV_PP}/gsmr_example_snps.allele \
	  --recode \
	  --out $output_path/${filename1}_${filename2}_${IV_PP}/gsmr_example
	  module unload gcta_1.93
	  Rscript "$GSMR_RSCRIPT" \
	  $HM \
	  $IV_P \
	  $LD_R2
	else
      echo "GSMR分析条件不满足，跳过执行。请检查是否已执行TwoSampleMR"
	fi
	#CAUSE
	if [[ "$CAUSE" == "T" && -f "$CAUSE_RSCRIPT" ]]; then
	  echo "条件满足，正在执行CAUSE：$gwas1 -> $gwas2"
	  module load R-4.2.0
	  Rscript "$CAUSE_RSCRIPT" \
	  $gwas1 \
	  $gwas2 \
	  $LD_R2 \
	  $IV_P \
	  $SCRIPT_DIR/ref/$Ref \
	  $output_path
	else
      echo "CAUSE分析条件不满足，跳过执行。"
	fi
	#ASSET分析
	if [[ "$ASSET" == "T" && -f "$ASSET_RSCRIPT" ]]; then
      echo "条件满足，正在执行ASSET_RSCRIPT：$gwas1 -> $gwas2"
	  module load R-4.2.0
      Rscript "$ASSET_RSCRIPT" \
      $gwas1 \
	  $gwas2 \
	  $casecontrol_table_path \
      $output_path
    else
      echo "$ASSET_RSCRIPT条件不满足，跳过执行。"
    fi
	#MiXeR分析
	if [[ "$MiXeR" == "T" && -f "$MiXeR_SCRIPT" && -f "$MiXeR_PAR_SCRIPT" ]]; then
	  echo "条件满足，正在执行MiXeR_SCRIPT：$gwas1 -> $gwas2"
	  mixerout=$output_path/${filename1}_${filename2}_MiXeRres
	  mkdir $mixerout
	  module load R-4.2.0
	  Rscript $MiXeR_PAR_SCRIPT \
	  $gwas1 \
	  $gwas2 \
	  $casecontrol_table_path \
	  $mixerout
	  cp $MiXeR_SCRIPT $mixerout
	  cd $mixerout
      sed -i "4a\Genetic=$(dirname "$0")" mixeR
	  sed -i "12a\cd $mixerout" mixeR
      sed -i "13a\export trait1=$(basename "$gwas1" | sed 's/\..*$//').gz" mixeR
      sed -i "14a\export trait2=$(basename "$gwas2" | sed 's/\..*$//').gz" mixeR
      sed -i "15a\export trait1_basename_noextension=$(basename "$gwas1" | sed 's/\..*$//')" mixeR
      sed -i "16a\export trait2_basename_noextension=$(basename "$gwas2" | sed 's/\..*$//')" mixeR
      sed -i "17a\export MPLCONFIGDIR=$mixerout" mixeR
	  # 添加 PBS 路径
      export PATH="/opt/tsce4/torque6/bin:$PATH"
      qsub \
	  -q fat \
	  -N mixeR_job \
      -t 1-20 \
      $mixerout/mixeR
	else
      echo "$MiXeR_SCRIPT条件不满足，跳过执行。"
    fi
	#conjFDR分析
	if [[ "$conjFDR" == "T" && -f "$conjFDR_SCRIPT" ]]; then
	  echo "条件满足，正在执行conjFDR：$gwas1 -> $gwas2"
	  module load R-4.2.0
	  N1=$(Rscript -e "library(data.table); data <- fread('$gwas1'); cat(data[['N']][1])")
	  N2=$(Rscript -e "library(data.table); data <- fread('$gwas2'); cat(data[['N']][1])")
	  $conjFDR_SCRIPT \
	  $gwas1 \
	  $N1 \
	  $gwas2 \
	  $N2 \
	  $output_path \
	  conjfdr \
	  $s_pleiofdr \
	  $s_python_convert
	else
      echo "$conjFDR_SCRIPT条件不满足，跳过执行。"
	fi
  done
done

#MAGMA分析    
if [[ "$MAGMA" == "T" && -f "$MAGMA_SCRIPT" ]]; then
  for gwas in "${gwas1_array[@]}" "${gwas2_array[@]}"; do
    echo "执行 MAGMA 分析：$gwas"
	module load R-4.2.0
	N3=$(Rscript -e "library(data.table); data <- fread('$gwas'); cat(data[['N']][1])")
    $MAGMA_SCRIPT \
	$gwas \
	$output_path/${filename1}_MAGMA \
	hg19 \
	$Ref \
	$N3 \
	txt
  done
else
  echo "$MAGMA_SCRIPT条件不满足，跳过执行。"
fi
