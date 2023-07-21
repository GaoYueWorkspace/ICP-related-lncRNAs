# 寻找模式多变的lnc-icp
path=getwd()

COAD=read.table("./COAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
ESCA=read.table("./ESCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
HNSC=read.table("./HNSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
KIRC=read.table("./KIRC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LIHC=read.table("./LIHC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LUAD=read.table("./LUAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LUSC=read.table("./LUSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
PRAD=read.table("./PRAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
READ=read.table("./READ/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
STAD=read.table("./STAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
THCA=read.table("./THCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
UCEC=read.table("./UCEC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
BRCA=read.table("./BRCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
BLCA=read.table("./BLCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
SKCM=read.table("./SKCM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
CHOL=read.table("./CHOL/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
GBM=read.table("./GBM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
OV=read.table("./OV/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)

library(tidyverse)
############################# 函数编译(genepair) ##############################
count_model=function(cancer){
  cancer_num=paste(cancer$lnc,cancer$icp,sep="-")
  cancer_num=as.data.frame(cancer_num)
  cancer_num[,"model"]=cancer$model
  colnames(cancer_num)[1]="ID"
  cancer_num=cancer_num %>% count(ID, model) %>% group_by(ID) %>%
    summarise(n = n_distinct(model)) 
  colnames(cancer_num)=c("lnc-icp","model_num")
  return(cancer_num)
}
################################################################################

############################ 统计各个癌症(genepair) ############################
COAD=count_model(COAD)
BLCA=count_model(BLCA)
BRCA=count_model(BRCA)
CHOL=count_model(CHOL)
ESCA=count_model(ESCA)
GBM=count_model(GBM)
HNSC=count_model(HNSC)
KIRC=count_model(KIRC)
LIHC=count_model(LIHC)
LUAD=count_model(LUAD)
LUSC=count_model(LUSC)
OV=count_model(OV)
PRAD=count_model(PRAD)
READ=count_model(READ)
SKCM=count_model(SKCM)
STAD=count_model(STAD)
THCA=count_model(THCA)
UCEC=count_model(UCEC)

lnc_icp_pair=list("COAD"=COAD,"BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"ESCA"=ESCA,
                  "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
                  "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                  "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)
################################################################################

library(dplyr)
library(purrr)
######  genepair ########
lst <- map(lnc_icp_pair, ~ .x[["lnc-icp"]])
lst <- reduce(lst, c)
df <- as.data.frame(table(lst))
colnames(df)=c("lnc-icp","count")
remove(lst)
cl_lncicp=df
remove(df)

############## 函数编译(过滤初步的模式多变的lnc-icp对) #########################
cancer_filter=function(cancer){
  cancer_filter=paste(cancer$lnc,cancer$icp,sep="-")
  cancer_filter=as.data.frame(cancer_filter)
  cancer_filter[,"imm"]=cancer$imm
  cancer_filter[,"model"]=cancer$model
  colnames(cancer_filter)[1]="ID"
  cancer_filter=cancer_filter[cl_lncicp$`lnc-icp`,]
  cancer_filter=cancer_filter%>%
    group_by(ID) %>%
    summarise(n_model = n_distinct(model), model_list = paste(sort(unique(model)), collapse = ",")) %>%
    ungroup()
  cancer_filter=as.data.frame(cancer_filter)
  return(cancer_filter)
}
################################################################################

COAD=read.table("./COAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
ESCA=read.table("./ESCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
HNSC=read.table("./HNSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
KIRC=read.table("./KIRC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LIHC=read.table("./LIHC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LUAD=read.table("./LUAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
LUSC=read.table("./LUSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
PRAD=read.table("./PRAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
READ=read.table("./READ/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
STAD=read.table("./STAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
THCA=read.table("./THCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
UCEC=read.table("./UCEC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
BRCA=read.table("./BRCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
BLCA=read.table("./BLCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
SKCM=read.table("./SKCM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
CHOL=read.table("./CHOL/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
GBM=read.table("./GBM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)
OV=read.table("./OV/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)


COAD=cancer_filter(COAD)
BLCA=cancer_filter(BLCA)
BRCA=cancer_filter(BRCA)
CHOL=cancer_filter(CHOL)
ESCA=cancer_filter(ESCA)
GBM=cancer_filter(GBM)
HNSC=cancer_filter(HNSC)
KIRC=cancer_filter(KIRC)
LIHC=cancer_filter(LIHC)
LUAD=cancer_filter(LUAD)
LUSC=cancer_filter(LUSC)
OV=cancer_filter(OV)
PRAD=cancer_filter(PRAD)
READ=cancer_filter(READ)
SKCM=cancer_filter(SKCM)
STAD=cancer_filter(STAD)
THCA=cancer_filter(THCA)
UCEC=cancer_filter(UCEC)

COAD=na.omit(COAD)
BLCA=na.omit(BLCA)
BRCA=na.omit(BRCA)
CHOL=na.omit(CHOL)
ESCA=na.omit(ESCA)
GBM=na.omit(GBM)
HNSC=na.omit(HNSC)
KIRC=na.omit(KIRC)
LIHC=na.omit(LIHC)
LUAD=na.omit(LUAD)
LUSC=na.omit(LUSC)
OV=na.omit(OV)
PRAD=na.omit(PRAD)
READ=na.omit(READ)
SKCM=na.omit(SKCM)
STAD=na.omit(STAD)
THCA=na.omit(THCA)
UCEC=na.omit(UCEC)

cl_lncicp_list=list("COAD"=COAD,"BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"ESCA"=ESCA,
                  "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
                  "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                  "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)


result <- bind_rows(cl_lncicp_list)
cl_lncicp2=as.data.frame(table(result$ID))
cl_lncicp2=cl_lncicp2[which(cl_lncicp2$Freq>=5),] # 选择在大于等于5个癌症中出现的lnc-icp对
result=result[which(result$ID%in%cl_lncicp2$Var1),]

result <- result %>%
  group_by(ID) %>%
  summarize(count = n_distinct(model_list),model_list = paste(sort(unique(model_list)), collapse = "\t、\t"))

result2=as.data.frame(t(as.data.frame(strsplit(result$ID,split = "-"))))
colnames(result2)=c("lnc_ENSEMBL","icp_ENSEMBL")
rownames(result2)=1:nrow(result2)

lnc=read.table("lncRNA_ENSG.txt",sep = "\t")
icp=read.table("ICP.txt",sep = "\t",header = T)

for(i in 1:nrow(result2)){
  result2[i,"lnc_SYMBOL"]=lnc[which(lnc$V1%in%result2[i,1]),2]
  result2[i,"icp_SYMBOL"]=icp[which(icp$ENSG%in%result2[i,2]),2]
}

result2[,"ID"]=result$ID

setwd(path)
write.table(result2,"cl_lncicp.txt",sep="\t",quote = F,row.names = F)

