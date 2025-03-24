# 介绍 #

# 环境配置 #
第一次运行前所需
## 模块配置 ##

## 原始分析软件配置 ##
### (1)conjFDR ###
github:https://github.com/precimed/pleiofdr
```
cd yourpath
wget https://github.com/precimed/pleiofdr/archive/refs/heads/master.zip yourpath
unzip pleiofdr-master.zip 
cd yourpath/pleiofdr-master
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/CTG_COG_2018.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/SSGAC_EDU_2016.mat
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_bfile.tar.gz
wget https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/9545380.ref
```
以上步骤操作完成后，	
打开`config_default.txt`将`reffile=`改为`reffile=/yourpath/pleiofdr-master/ref9545380_1kgPhase3eur_LDr2p1.mat`并保存	
打开`config_default.txt`将`traitfolder=`改为`traitfolder=/yourpath/pleiofdr-master/traitfolder`并保存	
### (2)MiXeR ###
github:https://github.com/comorment/mixer	
### (3)ldsc ###
github:https://github.com/bulik/ldsc
```
cd yourpath
git clone https://github.com/bulik/ldsc.git
```
### (3)python_convert ###
github:https://github.com/precimed/python_convert
```
cd yourpath
wget https://github.com/precimed/python_convert/archive/refs/heads/master.zip
unzip python_convert-master.zip
```
### 修改Genetic.sh中的脚本路径指向原始分析软件路径 ###
```
cd yourpath/OmniGWAS-Toolkit
```
打开`Genetic.sh`将`s_pleiofdr=/home/tengli/BQ_GWAS/pleiofdr-master`改为`s_pleiofdr=/yourpath/pleiofdr-master`	
打开`Genetic.sh`将`s_python_convert=/home/tengli/BQ_GWAS/python_convert-master`改为`s_python_convert=/yourpath/python_convert-master`	
打开`Genetic.sh`将`s_LDSC=/home/pub/softwares/ldsc`改为`s_LDSC=/yourpath/ldsc`	
保存`Genetic.sh`	
## 人群参考面板配置 ##
(1)Reference for causal analysis
```
cd OmniGWAS-Toolkit
wget http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz ./ref
cd OmniGWAS-Toolkit/ref
tar -zxvf 1kg.v3.tgz
```
(2)Reference for HDL (just EUR)
```
cd OmniGWAS-Toolkit/HDL-master/Reference_panels/HapMap3
wget -c -t 1 \
https://www.dropbox.com/s/6js1dzy4tkc3gac/UKB_imputed_SVD_eigen99_extraction.tar.gz?dl=0 \
--no-check-certificate -O UKB_imputed_SVD_eigen99_extraction.tar.gz
tar -xzvf UKB_imputed_SVD_eigen99_extraction.tar.gz
```
(3)Reference for MiXeR
通过网盘分享的文件：1000G_EUR_Phase3_plink.rar	
链接: https://pan.baidu.com/s/1HXVgjjajzRuoWmBndQKSyw?pwd=uxek 提取码: uxek	
```
cp 1000G_EUR_Phase3_plink.rar OmniGWAS-Toolkit/MixeR/reference/ldsc
cd OmniGWAS-Toolkit/MixeR/reference/ldsc
mkdir -p 1000G_EUR_Phase3_plink
unrar e 1000G_EUR_Phase3_plink.rar OmniGWAS-Toolkit/MixeR/reference/ldsc/1000G_EUR_Phase3_plink/
```
(4)info reference for GSMR
cite from:https://github.com/HaobinZhou/Get_MR/

通过网盘分享的文件：fileFrequency.frq
链接: https://pan.baidu.com/s/1eZYexlzdUPdSEbiY5_yoCg?pwd=jrr5 提取码: jrr5
```
cp fileFrequency.frq OmniGWAS-Toolkit
```
#仅对gwas1.txt和gwas2.txt执行数据预处理(自动识别hg38或hg19坐标轴并填补rs号，自动将hg38数据转为hg19，自动标准化表头，添加N列，自动计算并添加Z分数列)
## 方式1：
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt \
    --GWAS2 gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results
```
## 方式2：
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt，gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results
```
#执行基本遗传分析（启用 HDL 和 LDSC）​
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt \
    --GWAS2 gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis HDL LDSC
```	
#执行基本因果推断(因->果：gwas1a->gwas2a)（启用TwoSampleMR和GSMR以及CAUSE）​
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt \
    --GWAS2 gwas2a.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis TwoSampleMR GSMR CAUSE
```	
#多GWAS文件输入执行因果推断(因->果：gwas1a->gwas2a,gwas1a->gwas2b,gwas1b->gwas2a,gwas1b->gwas2b)
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis TwoSampleMR GSMR CAUSE
```	
#注意事项：执行GSMR分析的时候，必须连同TwoSampleMR一起执行
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis TwoSampleMR GSMR
```
#指定参考人群（非指定情况下，默认为EUR）和 IV 参数（非指定情况下，默认为5E-8）以及LD R2（默认0.001）
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --Ref EAS \
    --IV_P 1e-6 \
	--LD_R2 0.1 \
    --analysis TwoSampleMR GSMR CAUSE
```


#对4个GWAS汇总数据挨个进行MAGMA分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis MAGMA
```	
#对单个GWAS汇总数据进行MAGMA分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis MAGMA
```
#对多个东亚人群GWAS汇总数据执行所有分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1aEAS.txt,gwas1bEAS.txt \
    --GWAS2 gwas2aEAS.txt,gwas2bEAS.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --Ref EAS \
    --IV_P 1e-6 \
	--LD_R2 0.1 \
	--ASSET_casecontrol_table_path ./results/cc.csv
    --analysis HDL LDSC ASSET TwoSampleMR GSMR CAUSE MiXeR MAGMA
```
# 注意事项 #
conjFDR HDL MiXeR仅能执行欧洲人群的post-GWAS分析