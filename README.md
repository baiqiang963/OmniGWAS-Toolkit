# OmniGWAS工具包 (OGT)

## 概述  
一套全面的模块化流程，用于全基因组关联研究（GWAS）下游整合分析，旨在简化基因组学研究的后续探索。OGT整合了多种稳健框架：  
- ​**遗传相关性与遗传力**：HDL（LD分数回归）和LDSC，用于性状相关性分析和遗传力估计。  
- ​**精细定位与跨性状分析**：ASSET和conjFDR，用于共定位和多效性检测。  
- ​**因果推断**：GSMR、TwoSampleMR和CAUSE，支持孟德尔随机化及多变量因果建模。  
- ​**多基因结构解析**：MiXeR量化多基因重叠，MAGMA实现基因集/通路富集分析。  
- ​**标准化与格式转换**：自动化预处理（`GWAS_normalization.R`）和格式统一（`MixeR_GWAS_convert.R`）。  

OGT以互操作性为核心设计，支持跨方法验证（如CAUSE与TwoSampleMR对比）和批量化处理，通过用户友好的封装脚本实现复杂性状遗传架构的可扩展分析。其模块化结构允许灵活定制，适用于从功能注释（MAGMA）到多祖先荟萃分析等多样化研究目标。  

## 持续进化与协作  
- ​**持续扩展**：OGT正在积极开发中，计划整合新兴方法（如非线性孟德尔随机化工具DeepNull、跨祖先多基因评分工具PRS-CSx）并响应社区工具需求。  
- ​**用户优先优化**：定期更新将优化工作流（如并行化MAGMA运行、容器化依赖环境），以减少运行时间和手动操作。  
- ​**文档与教程**：持续完善的指南和单命令演示数据集将降低新用户的学习门槛。  

## 选择OGT的理由  
- ​**节省时间**：消除50%以上的方法整合冗余代码。  
- ​**一致性保障**：通过版本控制的流程确保跨项目可重复性。  
- ​**前沿技术触达**：无缝对接最新算法，保持研究领先性。  
---
# OmniGWAS-Toolkit (OGT)
## Introduction
A comprehensive, modular pipeline for post-GWAS integrative analysis, designed to streamline downstream genomic investigations.  
The OGT integrates multiple robust frameworks:  
​**Genetic Correlation & Heritability:** HDL (LD score regression) and LDSC for trait correlation and heritability estimation.  
​​**Fine-Mapping & Cross-Trait Analysis:** ASSET and conjFDR for colocalization and pleiotropy detection.	  
​**Causal Inference:​** GSMR, TwoSampleMR, and CAUSE for Mendelian randomization and multivariable causal modeling.	  
​​**Polygenic Architecture:** MiXeR for polygenic overlap quantification and MAGMA for gene-set/pathway enrichment.	  
​​**Standardization & Conversion:​** Automated preprocessing (GWAS_normalization.R) and format harmonization (MixeR_GWAS_convert.R).	  
Built with interoperability in mind, OGT supports cross-method validation (e.g., CAUSE vs. TwoSampleMR) and batch processing via user-friendly wrappers, enabling scalable exploration of complex trait architectures. Its modular design allows flexible customization for diverse study goals, from functional annotation (MAGMA) to multi-ancestry meta-analysis.  

## Future-Proof & Collaborative
​​​**Continuous Expansion​​​**: The OGT is under active development, with planned integration of emerging methods (e.g., DeepNull for nonlinear MR, PRS-CSx for cross-ancestry polygenic scoring) and community-driven tool requests.  
​​​**User-Centric Optimization​​​**: Regular updates will refine workflows (e.g., parallelized MAGMA runs, containerized environments for dependency control) to reduce runtime and manual intervention.  
​​​**Documentation & Tutorials​​​**: Evolving guides and one-command demo datasets will lower the learning curve for new users.  

## Why choose OGT?​​

​​​** Time-Saving **:​​ Eliminate 50%+ of boilerplate code for method integration.  
​​** ​Consistency **: Ensure reproducibility across projects with version-controlled pipelines.  
​​​** Cutting-Edge Access **: Stay ahead with seamless adoption of the latest algorithms.  
---
# 环境配置 #
install
## 模块配置 ##
![](https://github.com/baiqiang963/OmniGWAS-Toolkit/blob/main/pictures/R-4.2.0.jpg)
R version > 4.2.0
```
vim ~/my_modules/R-4.2.0   # 创建模块文件
```
```
#在~/my_modules/R-4.2.0中设置环境变量
module-whatis "R-4.2.0"
prepend-path    PATH           /yourpath/software/R-4.2.0/bin
prepend-path    LD_LIBRARY_PATH /yourpath/software/R-4.2.0/lib
...
```
```
#检查
module avai
```
![](https://github.com/baiqiang963/OmniGWAS-Toolkit/blob/main/pictures/Anaconda2.jpg)
```
mkdir -p ~/my_modules/anaconda  # 创建模块目录
vim ~/my_modules/anaconda/Anaconda2   # 创建模块文件
```
```
#在~/my_modules/anaconda/Anaconda2中设置环境变量
prepend-path    PATH           /yourpath/anaconda2//bin
prepend-path    LD_LIBRARY_PATH /yourpath/anaconda2//lib
...
```
```
#检查
module avai
```
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
#运行示例
##仅对gwas1.txt和gwas2.txt执行数据预处理  
**自动识别hg38或hg19坐标轴并填补rs号，自动将hg38数据转为hg19，自动标准化表头，添加N列，自动计算并添加Z分数列**
### 方式1：
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt \
    --GWAS2 gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results
```
### 方式2：
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt，gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results
```
---
**以下分析过程中，会自动预处理GWAS数据后进行**
## 执行基本遗传分析  
启用 HDL 和 LDSC​
```
sh ./Genetic.sh \
    --GWAS1 gwas1.txt \
    --GWAS2 gwas2.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis HDL LDSC
```	
## 执行基本因果推断  
因->果：gwas1a->gwas2a（启用TwoSampleMR和GSMR以及CAUSE）​
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt \
    --GWAS2 gwas2a.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis TwoSampleMR
```	
### 多GWAS文件输入执行GSMR因果推断  
因->果：gwas1a->gwas2a,gwas1a->gwas2b,gwas1b->gwas2a,gwas1b->gwas2b
```
sh ./Genetic.sh \
    --GWAS1 gwas1EAS.txt,gwas1EAS.txt \
    --GWAS2 gwas2EAS.txt,gwas2EAS.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis TwoSampleMR GSMR
```	
**注意事项:在集成脚本中，GSMR与CAUSE分析需要TwoSampleMR中的部分结果作为前置文件,因此需要连同TwoSampleMR仪器分析**
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --Ref EAS \
    --IV_P 1e-6 \
	--LD_R2 0.1 \
    --analysis TwoSampleMR GSMR
```
**指定参考人群（非指定情况下，默认为EUR）和 IV 参数（非指定情况下，默认为5E-8）以及LD R2（默认0.001）**
## MAGMA分析
### 对多个GWAS汇总数据挨个进行MAGMA分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis MAGMA
```	
### 对单个GWAS汇总数据进行MAGMA分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --analysis MAGMA
```
## 完整执行
### 对东亚人群的多个GWAS汇总数据执行所有分析
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
    --analysis LDSC ASSET TwoSampleMR GSMR CAUSE MAGMA
```
### 对欧洲人群的多个GWAS汇总数据执行所有分析
```
sh ./Genetic.sh \
    --GWAS1 gwas1a.txt,gwas1b.txt \
    --GWAS2 gwas2a.txt,gwas2b.txt \
    --formattedGWAS_path ./formatted \
    --output_path ./results \
    --Ref EAS \
    --IV_P 1e-6 \
	--LD_R2 0.1 \
	--ASSET_casecontrol_table_path ./results/cc.csv
    --analysis HDL LDSC ASSET TwoSampleMR GSMR CAUSE MAGMA conjFDR MiXeR
```
**注意事项：conjFDR HDL MiXeR仅能执行欧洲人群的post-GWAS分析**