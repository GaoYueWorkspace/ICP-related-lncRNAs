path=getwd()
result_lnc=read.table("result_lncRNA.txt",header = T)

########################## 获取表达谱 ##########################################
library(data.table)
library(readr)
# cancer=function(exp){
#   exp=as.data.frame(exp)
#   exp[,1]=substr(exp[,1],start = 1,stop = 15)
#   rownames(exp)=exp[,1]
#   exp=exp[,-1]
#   index=which(substr(colnames(exp),start = 14,stop = 15)>10)
#   exp=exp[,-index]
#   return(exp)
# }
# 
# BLCA=read_tsv("BLCA/TCGA-BLCA.htseq_fpkm.tsv")
# BLCA=cancer(BLCA)
# BRCA=read_tsv("BRCA/TCGA-BRCA.htseq_fpkm.tsv")
# BRCA=cancer(BRCA)
# CHOL=read_tsv("CHOL/TCGA-CHOL.htseq_fpkm.tsv")
# CHOL=cancer(CHOL)
# COAD=read_tsv("COAD/TCGA-COAD.htseq_fpkm.tsv")
# COAD=cancer(COAD)
# ESCA=read_tsv("ESCA/TCGA-ESCA.htseq_fpkm.tsv")
# ESCA=cancer(ESCA)
# GBM=read_tsv("GBM/TCGA-GBM.htseq_fpkm.tsv.gz")
# GBM=cancer(GBM)
# HNSC=read_tsv("HNSC/TCGA-HNSC.htseq_fpkm.tsv")
# HNSC=cancer(HNSC)
# KIRC=read_tsv("KIRC/TCGA-KIRC.htseq_fpkm.tsv")
# KIRC=cancer(KIRC)
# LIHC=read_tsv("LIHC/TCGA-LIHC.htseq_fpkm.tsv")
# LIHC=cancer(LIHC)
# LUAD=read_tsv("LUAD/TCGA-LUAD.htseq_fpkm.tsv")
# LUAD=cancer(LUAD)
# LUSC=read_tsv("LUSC/TCGA-LUSC.htseq_fpkm.tsv")
# LUSC=cancer(LUSC)
# OV=read_tsv("OV/TCGA-OV.htseq_fpkm.tsv.gz")
# OV=as.data.frame(OV)
# OV[,1]=substr(OV[,1],start = 1,stop = 15)
# rownames(OV)=OV[,1]
# OV=OV[,-1]
# PRAD=read_tsv("PRAD/TCGA-PRAD.htseq_fpkm.tsv")
# PRAD=cancer(PRAD)
# READ=read_tsv("READ/TCGA-READ.htseq_fpkm.tsv")
# READ=cancer(READ)
# SKCM=read_tsv("SKCM/TCGA-SKCM.htseq_fpkm.tsv")
# SKCM=cancer(SKCM)
# STAD=read_tsv("STAD/TCGA-STAD.htseq_fpkm.tsv")
# STAD=cancer(STAD)
# THCA=read_tsv("THCA/TCGA-THCA.htseq_fpkm.tsv")
# THCA=cancer(THCA)
# UCEC=read_tsv("UCEC/TCGA-UCEC.htseq_fpkm.tsv")
# UCEC=cancer(UCEC)
# 
# exp_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
#              "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
#              "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
#              "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
# remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
#        LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)
# remove(cancer)
# setwd(path)
# 
# save(exp_all,file = "exp_all.RData")
load("exp_all.RData")

## 筛选表达谱，只含有结果lnc
filter_cancer=function(df){
  return(df[which(rownames(df)%in%result_lnc$ENSG),])
}

exp_all=lapply(exp_all,filter_cancer)
remove(filter_cancer)

# 创建一个函数，用于去掉列名的最后一个字符
remove_last_char <- function(x) {
  # 使用substr函数截取列名的第一个字符到倒数第二个字符
  substr(x, 1, nchar(x) - 1)
}

# 使用lapply函数对列表中的每个数据框应用该函数
exp_all <- lapply(exp_all, function(df) {
  # 使用colnames函数获取和修改列名
  colnames(df) <- remove_last_char(colnames(df))
  # 返回修改后的数据框
  return(df)
})
remove(remove_last_char)

############################# 读入临床数据 #####################################
COAD=read.table("临床数据/COAD.txt",sep="\t",header = T)
ESCA=read.table("临床数据/ESCA.txt",sep="\t",header = T)
HNSC=read.table("临床数据/HNSC.txt",sep="\t",header = T)
KIRC=read.table("临床数据/KIRC.txt",sep="\t",header = T)
LIHC=read.table("临床数据/LIHC.txt",sep="\t",header = T)
LUAD=read.table("临床数据/LUAD.txt",sep="\t",header = T)
LUSC=read.table("临床数据/LUSC.txt",sep="\t",header = T)
PRAD=read.table("临床数据/PRAD.txt",sep="\t",header = T)
READ=read.table("临床数据/READ.txt",sep="\t",header = T)
STAD=read.table("临床数据/STAD.txt",sep="\t",header = T)
THCA=read.table("临床数据/THCA.txt",sep="\t",header = T)
UCEC=read.table("临床数据/UCEC.txt",sep="\t",header = T)
BRCA=read.table("临床数据/BRCA.txt",sep="\t",header = T)
BLCA=read.table("临床数据/BLCA.txt",sep="\t",header = T)
SKCM=read.table("临床数据/SKCM.txt",sep="\t",header = T)
CHOL=read.table("临床数据/CHOL.txt",sep="\t",header = T)
GBM=read.table("临床数据/GBM.txt",sep="\t",header = T)
OV=read.table("临床数据/OV.txt",sep="\t",header = T)
clinical_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
             "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
             "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
             "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)

## 过滤样本
# BLCA
tmp=intersect(colnames(exp_all$BLCA),clinical_all$BLCA$sample)
clinical_all$BLCA=clinical_all$BLCA[which(clinical_all$BLCA$sample%in%tmp),]
exp_all$BLCA=exp_all$BLCA[,clinical_all$BLCA$sample]

# BRCA
tmp=intersect(colnames(exp_all$BRCA),clinical_all$BRCA$sample)
clinical_all$BRCA=clinical_all$BRCA[which(clinical_all$BRCA$sample%in%tmp),]
exp_all$BRCA=exp_all$BRCA[,clinical_all$BRCA$sample]

# CHOL
tmp=intersect(colnames(exp_all$CHOL),clinical_all$CHOL$sample)
clinical_all$CHOL=clinical_all$CHOL[which(clinical_all$CHOL$sample%in%tmp),]
exp_all$CHOL=exp_all$CHOL[,clinical_all$CHOL$sample]

# COAD
tmp=intersect(colnames(exp_all$COAD),clinical_all$COAD$sample)
clinical_all$COAD=clinical_all$COAD[which(clinical_all$COAD$sample%in%tmp),]
exp_all$COAD=exp_all$COAD[,clinical_all$COAD$sample]

# ESCA
tmp=intersect(colnames(exp_all$ESCA),clinical_all$ESCA$sample)
clinical_all$ESCA=clinical_all$ESCA[which(clinical_all$ESCA$sample%in%tmp),]
exp_all$ESCA=exp_all$ESCA[,clinical_all$ESCA$sample]

# GBM
tmp=intersect(colnames(exp_all$GBM),clinical_all$GBM$sample)
clinical_all$GBM=clinical_all$GBM[which(clinical_all$GBM$sample%in%tmp),]
exp_all$GBM=exp_all$GBM[,clinical_all$GBM$sample]

# HNSC
tmp=intersect(colnames(exp_all$HNSC),clinical_all$HNSC$sample)
clinical_all$HNSC=clinical_all$HNSC[which(clinical_all$HNSC$sample%in%tmp),]
exp_all$HNSC=exp_all$HNSC[,clinical_all$HNSC$sample]

# KIRC
tmp=intersect(colnames(exp_all$KIRC),clinical_all$KIRC$sample)
clinical_all$KIRC=clinical_all$KIRC[which(clinical_all$KIRC$sample%in%tmp),]
exp_all$KIRC=exp_all$KIRC[,clinical_all$KIRC$sample]

# LIHC
tmp=intersect(colnames(exp_all$LIHC),clinical_all$LIHC$sample)
clinical_all$LIHC=clinical_all$LIHC[which(clinical_all$LIHC$sample%in%tmp),]
exp_all$LIHC=exp_all$LIHC[,clinical_all$LIHC$sample]

# LUAD
tmp=intersect(colnames(exp_all$LUAD),clinical_all$LUAD$sample)
clinical_all$LUAD=clinical_all$LUAD[which(clinical_all$LUAD$sample%in%tmp),]
exp_all$LUAD=exp_all$LUAD[,clinical_all$LUAD$sample]

# LUSC
tmp=intersect(colnames(exp_all$LUSC),clinical_all$LUSC$xena_sample)
clinical_all$LUSC=clinical_all$LUSC[which(clinical_all$LUSC$xena_sample%in%tmp),]
exp_all$LUSC=exp_all$LUSC[,clinical_all$LUSC$xena_sample]

# OV
tmp=intersect(colnames(exp_all$OV),clinical_all$OV$sample)
clinical_all$OV=clinical_all$OV[which(clinical_all$OV$sample%in%tmp),]
exp_all$OV=exp_all$OV[,clinical_all$OV$sample]

# PRAD
tmp=intersect(colnames(exp_all$PRAD),clinical_all$PRAD$sample)
clinical_all$PRAD=clinical_all$PRAD[which(clinical_all$PRAD$sample%in%tmp),]
exp_all$PRAD=exp_all$PRAD[,clinical_all$PRAD$sample]

# READ
tmp=intersect(colnames(exp_all$READ),clinical_all$READ$sample)
clinical_all$READ=clinical_all$READ[which(clinical_all$READ$sample%in%tmp),]
exp_all$READ=exp_all$READ[,clinical_all$READ$sample]

# SKCM
tmp=intersect(colnames(exp_all$SKCM),clinical_all$SKCM$sample)
clinical_all$SKCM=clinical_all$SKCM[which(clinical_all$SKCM$sample%in%tmp),]
exp_all$SKCM=exp_all$SKCM[,clinical_all$SKCM$sample]

# STAD
tmp=intersect(colnames(exp_all$STAD),clinical_all$STAD$sample)
clinical_all$STAD=clinical_all$STAD[which(clinical_all$STAD$sample%in%tmp),]
exp_all$STAD=exp_all$STAD[,clinical_all$STAD$sample]

# THCA
tmp=intersect(colnames(exp_all$THCA),clinical_all$THCA$sample)
clinical_all$THCA=clinical_all$THCA[which(clinical_all$THCA$sample%in%tmp),]
exp_all$THCA=exp_all$THCA[,clinical_all$THCA$sample]

# UCEC
tmp=intersect(colnames(exp_all$UCEC),clinical_all$UCEC$sample)
clinical_all$UCEC=clinical_all$UCEC[which(clinical_all$UCEC$sample%in%tmp),]
exp_all$UCEC=exp_all$UCEC[,clinical_all$UCEC$sample]


remove(tmp)
########################### 单cox分析 ##########################################
library(survival)
### BLCA
# 定义一个函数，用于对一个基因进行cox回归分析，返回回归系数、风险比、p值等信息
cox_gene <- function(x) {
  # 将基因表达水平和临床数据合并为一个数据框
  data <- cbind(clinical_all$BLCA, x)
  # 对合并后的数据进行cox回归分析，以生存时间和生存状态为因变量，以基因表达水平为自变量
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  # 提取cox回归分析的结果，包括回归系数、风险比、p值等
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  # 将结果组合为一个向量，返回该向量
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}

# 使用apply函数，对基因表达谱矩阵的每一行（每一个基因）应用上述函数，得到一个包含所有基因的cox分析结果的矩阵
cox_result <- apply(exp_all$BLCA, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
BLCA1=cox_result


logrank_gene=function(x){
  train_time=clinical_all$BLCA$OS.time
  train_status=clinical_all$BLCA$OS
  train_group=matrix(0,2,ncol(exp_all$BLCA))
  train_group[1,]=as.numeric(as.character(exp_all$BLCA[which(rownames(exp_all$BLCA)%in%x),]))
  colnames(train_group)=colnames(exp_all$BLCA)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
BLCA=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### BRCA
cox_gene <- function(x) {
  data <- cbind(clinical_all$BRCA, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$BRCA, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
BRCA1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$BRCA$OS.time
  train_status=clinical_all$BRCA$OS
  train_group=matrix(0,2,ncol(exp_all$BRCA))
  train_group[1,]=as.numeric(as.character(exp_all$BRCA[which(rownames(exp_all$BRCA)%in%x),]))
  colnames(train_group)=colnames(exp_all$BRCA)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
BRCA=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### CHOL
cox_gene <- function(x) {
  data <- cbind(clinical_all$CHOL, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$CHOL, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
CHOL1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$CHOL$OS.time
  train_status=clinical_all$CHOL$OS
  train_group=matrix(0,2,ncol(exp_all$CHOL))
  train_group[1,]=as.numeric(as.character(exp_all$CHOL[which(rownames(exp_all$CHOL)%in%x),]))
  colnames(train_group)=colnames(exp_all$CHOL)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
CHOL=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### COAD
cox_gene <- function(x) {
  data <- cbind(clinical_all$COAD, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$COAD, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
COAD1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$COAD$OS.time
  train_status=clinical_all$COAD$OS
  train_group=matrix(0,2,ncol(exp_all$COAD))
  train_group[1,]=as.numeric(as.character(exp_all$COAD[which(rownames(exp_all$COAD)%in%x),]))
  colnames(train_group)=colnames(exp_all$COAD)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
COAD=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### ESCA
cox_gene <- function(x) {
  data <- cbind(clinical_all$ESCA, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$ESCA, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
ESCA1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$ESCA$OS.time
  train_status=clinical_all$ESCA$OS
  train_group=matrix(0,2,ncol(exp_all$ESCA))
  train_group[1,]=as.numeric(as.character(exp_all$ESCA[which(rownames(exp_all$ESCA)%in%x),]))
  colnames(train_group)=colnames(exp_all$ESCA)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
ESCA=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### GBM
cox_gene <- function(x) {
  data <- cbind(clinical_all$GBM, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$GBM, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
GBM1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$GBM$OS.time
  train_status=clinical_all$GBM$OS
  train_group=matrix(0,2,ncol(exp_all$GBM))
  train_group[1,]=as.numeric(as.character(exp_all$GBM[which(rownames(exp_all$GBM)%in%x),]))
  colnames(train_group)=colnames(exp_all$GBM)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
GBM=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### HNSC
cox_gene <- function(x) {
  data <- cbind(clinical_all$HNSC, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$HNSC, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
HNSC1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$HNSC$OS.time
  train_status=clinical_all$HNSC$OS
  train_group=matrix(0,2,ncol(exp_all$HNSC))
  train_group[1,]=as.numeric(as.character(exp_all$HNSC[which(rownames(exp_all$HNSC)%in%x),]))
  colnames(train_group)=colnames(exp_all$HNSC)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
HNSC=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### KIRC
cox_gene <- function(x) {
  data <- cbind(clinical_all$KIRC, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$KIRC, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
KIRC1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$KIRC$OS.time
  train_status=clinical_all$KIRC$OS
  train_group=matrix(0,2,ncol(exp_all$KIRC))
  train_group[1,]=as.numeric(as.character(exp_all$KIRC[which(rownames(exp_all$KIRC)%in%x),]))
  colnames(train_group)=colnames(exp_all$KIRC)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
KIRC=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### LIHC
cox_gene <- function(x) {
  data <- cbind(clinical_all$LIHC, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$LIHC, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
LIHC1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$LIHC$OS.time
  train_status=clinical_all$LIHC$OS
  train_group=matrix(0,2,ncol(exp_all$LIHC))
  train_group[1,]=as.numeric(as.character(exp_all$LIHC[which(rownames(exp_all$LIHC)%in%x),]))
  colnames(train_group)=colnames(exp_all$LIHC)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
LIHC=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### LUAD
cox_gene <- function(x) {
  data <- cbind(clinical_all$LUAD, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$LUAD, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
LUAD1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$LUAD$OS.time
  train_status=clinical_all$LUAD$OS
  train_group=matrix(0,2,ncol(exp_all$LUAD))
  train_group[1,]=as.numeric(as.character(exp_all$LUAD[which(rownames(exp_all$LUAD)%in%x),]))
  colnames(train_group)=colnames(exp_all$LUAD)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
LUAD=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### LUSC
cox_gene <- function(x) {
  data <- cbind(clinical_all$LUSC, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$LUSC, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
LUSC1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$LUSC$OS.time
  train_status=clinical_all$LUSC$OS
  train_group=matrix(0,2,ncol(exp_all$LUSC))
  train_group[1,]=as.numeric(as.character(exp_all$LUSC[which(rownames(exp_all$LUSC)%in%x),]))
  colnames(train_group)=colnames(exp_all$LUSC)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
LUSC=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### OV
cox_gene <- function(x) {
  data <- cbind(clinical_all$OV, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$OV, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
OV1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$OV$OS.time
  train_status=clinical_all$OV$OS
  train_group=matrix(0,2,ncol(exp_all$OV))
  train_group[1,]=as.numeric(as.character(exp_all$OV[which(rownames(exp_all$OV)%in%x),]))
  colnames(train_group)=colnames(exp_all$OV)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
OV=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

#### PRAD
cox_gene <- function(x) {
  data <- cbind(clinical_all$PRAD, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$PRAD, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
PRAD1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$PRAD$OS.time
  train_status=clinical_all$PRAD$OS
  train_group=matrix(0,2,ncol(exp_all$PRAD))
  train_group[1,]=as.numeric(as.character(exp_all$PRAD[which(rownames(exp_all$PRAD)%in%x),]))
  colnames(train_group)=colnames(exp_all$PRAD)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
PRAD=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### READ
cox_gene <- function(x) {
  data <- cbind(clinical_all$READ, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$READ, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
READ1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$READ$OS.time
  train_status=clinical_all$READ$OS
  train_group=matrix(0,2,ncol(exp_all$READ))
  train_group[1,]=as.numeric(as.character(exp_all$READ[which(rownames(exp_all$READ)%in%x),]))
  colnames(train_group)=colnames(exp_all$READ)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
READ=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### SKCM
cox_gene <- function(x) {
  data <- cbind(clinical_all$SKCM, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$SKCM, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
SKCM1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$SKCM$OS.time
  train_status=clinical_all$SKCM$OS
  train_group=matrix(0,2,ncol(exp_all$SKCM))
  train_group[1,]=as.numeric(as.character(exp_all$SKCM[which(rownames(exp_all$SKCM)%in%x),]))
  colnames(train_group)=colnames(exp_all$SKCM)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
SKCM=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### STAD
cox_gene <- function(x) {
  data <- cbind(clinical_all$STAD, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$STAD, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
STAD1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$STAD$OS.time
  train_status=clinical_all$STAD$OS
  train_group=matrix(0,2,ncol(exp_all$STAD))
  train_group[1,]=as.numeric(as.character(exp_all$STAD[which(rownames(exp_all$STAD)%in%x),]))
  colnames(train_group)=colnames(exp_all$STAD)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
STAD=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### THCA
cox_gene <- function(x) {
  data <- cbind(clinical_all$THCA, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$THCA, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
THCA1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$THCA$OS.time
  train_status=clinical_all$THCA$OS
  train_group=matrix(0,2,ncol(exp_all$THCA))
  train_group[1,]=as.numeric(as.character(exp_all$THCA[which(rownames(exp_all$THCA)%in%x),]))
  colnames(train_group)=colnames(exp_all$THCA)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
THCA=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

### UCEC
cox_gene <- function(x) {
  data <- cbind(clinical_all$UCEC, x)
  fit <- coxph(Surv(OS.time, OS) ~ x, data = data)
  coef <- fit$coef # 回归系数
  hr <- exp(coef) # 风险比
  ci <- exp(confint(fit)) # 风险比的置信区间
  pval <- summary(fit)$coef[5] # p值
  result <- c(coef, hr, ci[1], ci[2], pval)
  return(result)
}
cox_result <- apply(exp_all$UCEC, 1, cox_gene)
cox_result=as.data.frame(t(cox_result))
colnames(cox_result)=c("回归系数","风险比","置信区间1","置信区间2","p值")
cox_result=cox_result[which(cox_result$p值<0.05),]
UCEC1=cox_result
logrank_gene=function(x){
  train_time=clinical_all$UCEC$OS.time
  train_status=clinical_all$UCEC$OS
  train_group=matrix(0,2,ncol(exp_all$UCEC))
  train_group[1,]=as.numeric(as.character(exp_all$UCEC[which(rownames(exp_all$UCEC)%in%x),]))
  colnames(train_group)=colnames(exp_all$UCEC)
  cutoff=median(train_group[1,])
  index_low=train_group[1,]<=cutoff
  index_high=!(index_low)
  train_group[2,index_low]=1
  train_group[2,index_high]=2
  train=list(train_time,train_status,train_group[2,])
  y=Surv(train_time,train_status)
  dif=survfit(y~train_group[2,])
  r=survdiff(y~train_group[2,],train)
  p=1-pchisq(r$chisq,length(r$n)-1)
  return(p)
}
surv_gene=as.data.frame(rownames(cox_result))
logrank_result=apply(surv_gene,1,logrank_gene)
logrank_result=data.frame("gene"=rownames(cox_result),"p"=logrank_result)
logrank_result=logrank_result[which(logrank_result$p<0.05),]
UCEC=logrank_result$gene
remove(cox_gene,cox_result,logrank_gene,logrank_result,surv_gene)

surv_lnc_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
             "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
             "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
             "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)

surv_lnc_cox_all=list("BLCA"=BLCA1,"BRCA"=BRCA1,"CHOL"=CHOL1,"COAD"=COAD1,"ESCA"=ESCA1,
                  "GBM"=GBM1,"HNSC"=HNSC1,"KIRC"=KIRC1,"LIHC"=LIHC1,"LUAD"=LUAD1,
                  "LUSC"=LUSC1,"OV"=OV1,"PRAD"=PRAD1,"READ"=READ1,"SKCM"=SKCM1,
                  "STAD"=STAD1,"THCA"=THCA1,"UCEC"=UCEC1)
remove(COAD1,BLCA1,BRCA1,CHOL1,ESCA1,GBM1,HNSC1,KIRC1,LIHC1,
       LUAD1,LUSC1,OV1,PRAD1,READ1,SKCM1,STAD1,THCA1,UCEC1)


gene=surv_lnc_all$BLCA
tmp=surv_lnc_all$BRCA
gene=c(gene,tmp)
tmp=surv_lnc_all$CHOL
gene=c(gene,tmp)
tmp=surv_lnc_all$COAD
gene=c(gene,tmp)
tmp=surv_lnc_all$ESCA
gene=c(gene,tmp)
tmp=surv_lnc_all$GBM
gene=c(gene,tmp)
tmp=surv_lnc_all$HNSC
gene=c(gene,tmp)
tmp=surv_lnc_all$KIRC
gene=c(gene,tmp)
tmp=surv_lnc_all$LIHC
gene=c(gene,tmp)
tmp=surv_lnc_all$LUAD
gene=c(gene,tmp)
tmp=surv_lnc_all$LUSC
gene=c(gene,tmp)
tmp=surv_lnc_all$OV
gene=c(gene,tmp)
tmp=surv_lnc_all$PRAD
gene=c(gene,tmp)
tmp=surv_lnc_all$READ
gene=c(gene,tmp)
tmp=surv_lnc_all$SKCM
gene=c(gene,tmp)
tmp=surv_lnc_all$STAD
gene=c(gene,tmp)
tmp=surv_lnc_all$THCA
gene=c(gene,tmp)
tmp=surv_lnc_all$UCEC
gene=c(gene,tmp)


gene=as.data.frame(table(gene))
# surv_lnc=as.character(gene$gene)
surv_lnc=as.character(gene[which(gene$Freq>3),1])  #在至少三个癌症中都出现
anno=read.table("注释文件/SYMBOL_ENSEMBL.txt",header = T)
surv_lnc_cox=do.call(rbind,surv_lnc_cox_all)
library(stringr)
library(tidyverse)
name=rownames(surv_lnc_cox)
name=as.data.frame(strsplit(name,"\\."))
name=as.data.frame(t(name))
surv_lnc_cox[,"cancer"]=name$V1
surv_lnc_cox[,"gene"]=name$V2
rownames(surv_lnc_cox)=1:nrow(surv_lnc_cox)
surv_lnc_cox=surv_lnc_cox[which(surv_lnc_cox$gene%in%surv_lnc),]
surv_lnc_cox$风险比=log(surv_lnc_cox$风险比)
surv_lnc_cox=surv_lnc_cox[,c(2,6,7)]
for(i in 1:nrow(surv_lnc_cox)){
  surv_lnc_cox[i,"symbol"]=anno[which(anno$ENSG%in%surv_lnc_cox[i,3]),2]
}
library(ggplot2)
library(reshape2)
surv_lnc_cox=surv_lnc_cox[which(abs(surv_lnc_cox$风险比)<2),]
ggplot(data = surv_lnc_cox, aes(cancer,symbol ,fill = 风险比))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),limits=c(-2,2)) +  
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()

surv_lnc=surv_lnc_cox[,3:4]
surv_lnc=unique(surv_lnc)
write.table(surv_lnc,"surv_lnc.txt",quote = F,sep = "\t",row.names = F)
