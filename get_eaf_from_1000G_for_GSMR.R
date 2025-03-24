library(data.table)
args=commandArgs(T)
Genetic_path=args[1]
dat_path=args[2]
get_eaf_from_1000G<-function(dat,path,type="exposure"){
  corrected_eaf_expo<-function(data_MAF){
    effect=data_MAF$effect_allele.exposure
    other=data_MAF$other_allele.exposure
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    
    cor1<-which(data_MAF$effect_allele.exposure !=data_MAF$A1)
    
    
    data_MAF$eaf.exposure=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.exposure[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.exposure ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.exposure[error]=NA
    data_MAF$type[error]="error"
    
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    
    return(data_MAF)
    
  }
  
  corrected_eaf_out<-function(data_MAF){
    effect=data_MAF$effect_allele.outcome
    other=data_MAF$other_allele.outcome
    A1=data_MAF$A1
    A2=data_MAF$A2
    MAF_num=data_MAF$MAF
    EAF_num=1-MAF_num
    
    harna<-is.na(data_MAF$A1)
    harna<-data_MAF$SNP[which(harna==T)]
    
    cor1<-which(data_MAF$effect_allele.outcome !=data_MAF$A1)
    
    
    data_MAF$eaf.outcome=data_MAF$MAF
    data_MAF$type="raw"
    data_MAF$eaf.outcome[cor1]=EAF_num[cor1]
    data_MAF$type[cor1]="corrected"
    cor2<-which(data_MAF$other_allele.outcome ==data_MAF$A1)
    cor21<-setdiff(cor2,cor1)
    cor12<-setdiff(cor1,cor2)
    error<-c(cor12,cor21)
    data_MAF$eaf.outcome[error]=NA
    data_MAF$type[error]="error"
    
    data_MAF<-list(data_MAF=data_MAF,cor1=cor1,harna=harna,error=error)
    
    return(data_MAF)
    
  }
  
  if(type=="exposure" && (is.na(dat$eaf.exposure[1])==T || is.null(dat$eaf.exposure)==T)){
    r<-nrow(dat)
    
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    
    dat<-corrected_eaf_expo(dat)
    
    cor1<-dat$cor1
    
    harna<-dat$harna
    
    error<-dat$error
    
    dat<-dat$data_MAF
    
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    
    print(paste0("一共有",length(harna),"个SNP在1000G中找不到，占比",length(harna)/r*100,"%"))
    
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    
    return(dat)
  }
  
  if(type=="outcome" && (is.na(dat$eaf.outcome[1])==T || is.null(dat$eaf.outcome)==T)){
    r<-nrow(dat)
    
    setwd(path)
    MAF<-fread("fileFrequency.frq",header = T)
    
    dat<-merge(dat,MAF,by.x = "SNP",by.y = "SNP",all.x = T)
    
    dat<-corrected_eaf_out(dat)
    
    cor1<-dat$cor1
    
    harna<-dat$harna
    
    error<-dat$error
    
    dat<-dat$data_MAF
    
    print(paste0("一共有",(r-length(harna)-length(error)),"个SNP成功匹配EAF,占比",(r-length(harna)-length(error))/r*100,"%"))
    
    print(paste0("一共有",length(cor1),"个SNP是major allele，EAF被计算为1-MAF,在成功匹配数目中占比",length(cor1)/(r-length(harna)-length(error))*100,"%"))
    
    print(paste0("一共有",length(harna),"个SNP在1000G找不到，占比",length(harna)/r*100,"%"))
    
    print(paste0("一共有",length(error),"个SNP在输入数据与1000G中效应列与参照列，将剔除eaf，占比",length(error)/r*100,"%"))
    
    print("输出数据中的type列说明：")
    print("raw：EAF直接等于1000G里的MAF数值，因为效应列是minor allele")
    print('corrected：EAF等于1000G中1-MAF，因为效应列是major allele')
    print("error：输入数据与1000G里面提供的数据完全不一致，比如这个SNP输入的效应列是C，参照列是G，但是1000G提供的是A-T，这种情况下，EAF会被清空（NA），当成匹配失败")
    
    return(dat)
  }
  else{return(dat)}
}
dat=fread(dat_path)
dat=get_eaf_from_1000G(dat,Genetic_path,type="exposure")
fwrite(dat,dat_path,row.names = F,quote = F)