path=getwd()
load("exp_all.RData")
library(data.table)
library(readr)
library(survival)
############################# 获取所有的结果三元组 #############################

COAD=read.table("./COAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
ESCA=read.table("./ESCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
HNSC=read.table("./HNSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
KIRC=read.table("./KIRC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LIHC=read.table("./LIHC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LUAD=read.table("./LUAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LUSC=read.table("./LUSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
PRAD=read.table("./PRAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
READ=read.table("./READ/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
STAD=read.table("./STAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
THCA=read.table("./THCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
UCEC=read.table("./UCEC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
BRCA=read.table("./BRCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
BLCA=read.table("./BLCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
SKCM=read.table("./SKCM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
CHOL=read.table("./CHOL/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
GBM=read.table("./GBM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
OV=read.table("./OV/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]

triplet_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
                 "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
                 "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                 "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)

############################# 读入临床数据 #####################################
setwd(path)
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

########################### 过滤样本 ##########################################
## 过滤样本

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

################################ 转置矩阵 #####################################
exp_all=lapply(exp_all,function(x){
  x=as.data.frame(t(x))
})

############################### 多cox分析 #####################################

### BLCA
exp=exp_all$BLCA
triplet=triplet_all$BLCA
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$BLCA$OS.time,clinical_all$BLCA$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$BLCA$OS.time
  train_status=clinical_all$BLCA$OS
  train_group=matrix(0,2,nrow(exp_all$BLCA))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$BLCA)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
BLCA=result
remove(exp,multi_cox,result,triplet)

### BRCA
exp=exp_all$BRCA
triplet=triplet_all$BRCA
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$BRCA$OS.time,clinical_all$BRCA$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$BRCA$OS.time
  train_status=clinical_all$BRCA$OS
  train_group=matrix(0,2,nrow(exp_all$BRCA))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$BRCA)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
BRCA=result
remove(exp,multi_cox,result,triplet)

### CHOL
gc()
exp=exp_all$CHOL
triplet=triplet_all$CHOL
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$CHOL$OS.time,clinical_all$CHOL$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$CHOL$OS.time
  train_status=clinical_all$CHOL$OS
  train_group=matrix(0,2,nrow(exp_all$CHOL))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$CHOL)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
CHOL=result
remove(exp,multi_cox,result,triplet)

### COAD
gc()
exp=exp_all$COAD
triplet=triplet_all$COAD
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$COAD$OS.time,clinical_all$COAD$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$COAD$OS.time
  train_status=clinical_all$COAD$OS
  train_group=matrix(0,2,nrow(exp_all$COAD))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$COAD)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
COAD=result
remove(exp,multi_cox,result,triplet)

### ESCA
gc()
exp=exp_all$ESCA
triplet=triplet_all$ESCA
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$ESCA$OS.time,clinical_all$ESCA$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$ESCA$OS.time
  train_status=clinical_all$ESCA$OS
  train_group=matrix(0,2,nrow(exp_all$ESCA))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$ESCA)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
ESCA=result
remove(exp,multi_cox,result,triplet)

### GBM
gc()
exp=exp_all$GBM
triplet=triplet_all$GBM
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$GBM$OS.time,clinical_all$GBM$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$GBM$OS.time
  train_status=clinical_all$GBM$OS
  train_group=matrix(0,2,nrow(exp_all$GBM))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$GBM)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
GBM=result
remove(exp,multi_cox,result,triplet)

### HNSC
gc()
exp=exp_all$HNSC
triplet=triplet_all$HNSC
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$HNSC$OS.time,clinical_all$HNSC$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$HNSC$OS.time
  train_status=clinical_all$HNSC$OS
  train_group=matrix(0,2,nrow(exp_all$HNSC))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$HNSC)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
HNSC=result
remove(exp,multi_cox,result,triplet)

### KIRC
gc()
exp=exp_all$KIRC
triplet=triplet_all$KIRC
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$KIRC$OS.time,clinical_all$KIRC$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$KIRC$OS.time
  train_status=clinical_all$KIRC$OS
  train_group=matrix(0,2,nrow(exp_all$KIRC))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$KIRC)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
KIRC=result
remove(exp,multi_cox,result,triplet)

### LIHC
gc()
exp=exp_all$LIHC
triplet=triplet_all$LIHC
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$LIHC$OS.time,clinical_all$LIHC$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$LIHC$OS.time
  train_status=clinical_all$LIHC$OS
  train_group=matrix(0,2,nrow(exp_all$LIHC))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$LIHC)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
LIHC=result
remove(exp,multi_cox,result,triplet)

### LUAD
gc()
exp=exp_all$LUAD
triplet=triplet_all$LUAD
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$LUAD$OS.time,clinical_all$LUAD$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$LUAD$OS.time
  train_status=clinical_all$LUAD$OS
  train_group=matrix(0,2,nrow(exp_all$LUAD))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$LUAD)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
LUAD=result
remove(exp,multi_cox,result,triplet)

### LUSC
gc()
exp=exp_all$LUSC
triplet=triplet_all$LUSC
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$LUSC$OS.time,clinical_all$LUSC$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$LUSC$OS.time
  train_status=clinical_all$LUSC$OS
  train_group=matrix(0,2,nrow(exp_all$LUSC))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$LUSC)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
LUSC=result
remove(exp,multi_cox,result,triplet)

### OV
gc()
exp=exp_all$OV
triplet=triplet_all$OV
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$OV$OS.time,clinical_all$OV$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$OV$OS.time
  train_status=clinical_all$OV$OS
  train_group=matrix(0,2,nrow(exp_all$OV))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$OV)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
OV=result
remove(exp,multi_cox,result,triplet)

### PRAD
gc()
exp=exp_all$PRAD
triplet=triplet_all$PRAD
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$PRAD$OS.time,clinical_all$PRAD$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$PRAD$OS.time
  train_status=clinical_all$PRAD$OS
  train_group=matrix(0,2,nrow(exp_all$PRAD))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$PRAD)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
PRAD=result
remove(exp,multi_cox,result,triplet)

### READ
gc()
exp=exp_all$READ
triplet=triplet_all$READ
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$READ$OS.time,clinical_all$READ$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$READ$OS.time
  train_status=clinical_all$READ$OS
  train_group=matrix(0,2,nrow(exp_all$READ))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$READ)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
READ=result
remove(exp,multi_cox,result,triplet)

### SKCM
gc()
exp=exp_all$SKCM
triplet=triplet_all$SKCM
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$SKCM$OS.time,clinical_all$SKCM$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$SKCM$OS.time
  train_status=clinical_all$SKCM$OS
  train_group=matrix(0,2,nrow(exp_all$SKCM))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$SKCM)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
SKCM=result
remove(exp,multi_cox,result,triplet)

### STAD
gc()
exp=exp_all$STAD
triplet=triplet_all$STAD
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$STAD$OS.time,clinical_all$STAD$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$STAD$OS.time
  train_status=clinical_all$STAD$OS
  train_group=matrix(0,2,nrow(exp_all$STAD))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$STAD)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
STAD=result
remove(exp,multi_cox,result,triplet)

### THCA
gc()
exp=exp_all$THCA
triplet=triplet_all$THCA
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$THCA$OS.time,clinical_all$THCA$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$THCA$OS.time
  train_status=clinical_all$THCA$OS
  train_group=matrix(0,2,nrow(exp_all$THCA))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$THCA)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
THCA=result
remove(exp,multi_cox,result,triplet)

### UCEC
gc()
exp=exp_all$UCEC
triplet=triplet_all$UCEC
tmp=c(unique(triplet$icp),unique(triplet$lnc),unique(triplet$imm))
tmp=unique(tmp)
exp=exp[,tmp]
exp=cbind(clinical_all$UCEC$OS.time,clinical_all$UCEC$OS,exp)
colnames(exp)[1:2]=c("OS.time","OS")

multi_cox=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  train_time=clinical_all$UCEC$OS.time
  train_status=clinical_all$UCEC$OS
  train_group=matrix(0,2,nrow(exp_all$UCEC))
  train_group[1,]=score
  colnames(train_group)=rownames(exp_all$UCEC)
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

result=apply(triplet,1,multi_cox)
result=cbind(triplet,result)
result=result[which(result$result<0.05),]
UCEC=result
remove(exp,multi_cox,result,triplet)

surv_triplet_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
                  "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
                  "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                  "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)

surv_triplet=do.call(rbind,surv_triplet_all)
name=rownames(surv_triplet)
name=as.data.frame(strsplit(name,"\\."))
name=as.data.frame(t(name))
surv_triplet[,"cancer"]=name$V1
remove(name)
##### 统计在每种癌症中的预后相关三元组的数量
tmp=as.data.frame(table(surv_triplet$cancer))
static_in_cancer=as.numeric(tmp$Freq)
static_in_cancer=as.data.frame(static_in_cancer)
static_in_cancer=as.data.frame(t(static_in_cancer))
colnames(static_in_cancer)=tmp$Var1
rownames(static_in_cancer)="Freq"
write.table(static_in_cancer,"雷达图数据.txt",quote = F,row.names = F,sep="\t")

# static_in_cancer=read.table("雷达图数据.txt",header = T)
library(ggplot2)
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)
ggradar(static_in_cancer,
        grid.max = max(static_in_cancer[1,]),
        grid.min = 0,
        grid.mid = max(static_in_cancer[1,])/2,
        axis.label.size = 6,
        group.point.size = 3,
        group.colours = "#5b9bd5")
# static_in_cancer=as.data.frame(t(static_in_cancer))
# static_in_cancer[,"cancer"]=rownames(static_in_cancer)
# ggplot(static_in_cancer, aes(x = cancer, y = Freq, fill = cancer)) +
#   geom_col(width = 1) +
#   coord_polar() +
#   scale_fill_manual(values = c("#e3a646","#cc762b","#d68158","#b3524e","#d0655f",
#                                "#b36383","#d378a0","#edbecc","#f5e6f0","#e5c4dc",
#                                "#d1bdd7","#c89ec7","#7393c4","#8cabd5","#a3b9ce",
#                                "#cfe1f2","#c9bb93","#e9d9b0")) +
#   theme_void()

surv_triplet[,"id"]=paste(surv_triplet$lnc,surv_triplet$icp,surv_triplet$imm,sep = "+")
static=as.data.frame(table(surv_triplet$id))
table(static$Freq)
barplot(c(77322,5911,683,67,9)/83992)
static=static[which(static$Freq>4),]
static[,1]=as.character(static[,1])
surv_triplet_plot=surv_triplet[which(surv_triplet$id%in%static$Var1),]
anno=read.table("注释文件/SYMBOL_ENSEMBL.txt",header = T)
for(i in 1:nrow(surv_triplet_plot)){
  surv_triplet_plot[i,"icp_symbol"]=anno[which(anno$ENSG%in%surv_triplet_plot[i,1]),2]
}
for(i in 1:nrow(surv_triplet_plot)){
  surv_triplet_plot[i,"lnc_symbol"]=anno[which(anno$ENSG%in%surv_triplet_plot[i,2]),2]
}
for(i in 1:nrow(surv_triplet_plot)){
  surv_triplet_plot[i,"imm_symbol"]=anno[which(anno$ENSG%in%surv_triplet_plot[i,3]),2]
}
colnames(surv_triplet_plot)[1:3]=c("icp_ensembl","lnc_ensembl","imm_ensembl")
surv_triplet_plot[,"id"]=paste(surv_triplet_plot$lnc_symbol,
                               surv_triplet_plot$icp_symbol,
                               surv_triplet_plot$imm_symbol,
                               sep = ";")
plot.data=surv_triplet_plot[,c(6,5,4)]
rownames(plot.data)=1:nrow(plot.data)
library(ggplot2)
library(reshape2)
ggplot(data = plot.data, aes(cancer,id ,fill = result))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("white","#d94c4c"),limits=c(0,0.05)) +  
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()


remove(i,anno,static,static_in_cancer,tmp)
############################ 绘制森林图 #######################################
gc()
exp=exp_all$KIRC
clinical_KIRC=read.csv("临床数据/clinical_use.tsv",header = T,sep = "\t")
clinical_KIRC=clinical_KIRC[,c(2,4,12)]
clinical_KIRC=unique(clinical_KIRC)
rownames(clinical_KIRC)=clinical_KIRC$case_submitter_id
name=substr(rownames(exp),start = 1,stop = 12)
clinical_KIRC=clinical_KIRC[name,]

exp=cbind(clinical_all$KIRC$OS.time,clinical_all$KIRC$OS,clinical_KIRC$age_at_index,clinical_KIRC$gender,exp)
colnames(exp)[1:4]=c("OS.time","OS","age","gender")
exp_pre=exp
x=c("ENSG00000234883","ENSG00000169245","ENSG00000163600")
formula <- as.formula(paste("Surv(OS.time, OS) ~", 
                            paste(x, collapse = "+")))
fit <- coxph(formula, data = exp)

coef <- as.numeric(fit$coef) # 回归系数
score1=coef[1]*(as.numeric(exp[,x[1]]))
score2=coef[2]*(as.numeric(exp[,x[2]]))
score3=coef[3]*(as.numeric(exp[,x[3]]))
score=score1+score2+score3
train_time=clinical_all$KIRC$OS.time
train_status=clinical_all$KIRC$OS
train_group=matrix(0,2,nrow(exp_all$KIRC))
train_group[1,]=score
colnames(train_group)=rownames(exp_all$KIRC)
cutoff=median(train_group[1,])
index_low=train_group[1,]<=cutoff
index_high=!(index_low)
train_group[2,index_low]=1
train_group[2,index_high]=2
train=list(train_time,train_status,train_group[2,])
y=Surv(train_time,train_status)
dif=survfit(y~train_group[2,])
r=survdiff(y~train_group[2,],train)
1-pchisq(r$chisq,length(r$n)-1)
remove(fit,formula,index_high,index_low,r,score1,score2,score3,train,train_time
       ,train_status,x,y,cutoff,coef,dif)

risk=exp[,1:4]
train_group=as.data.frame(t(train_group))
risk=cbind(risk,train_group)
colnames(risk)[5:6]=c("risk","group")
risk$group=factor(risk$group,levels = c(1,2),labels = c("low","high"))

fit <- coxph(Surv(OS.time, OS) ~age+gender+risk, data = risk)
ggforest(fit, data = risk)

train_group=cbind(clinical_all$KIRC$OS.time,clinical_all$KIRC$OS,train_group)
colnames(train_group)=c("OS.time","OS","Risk","group")
fit <- survfit(Surv(OS.time, OS)~group,data = train_group)
## 绘制logrank的KM曲线
library(survminer)
library(survival)
ggsurvplot(fit, data = train_group,
           pval = T,conf.int = T,
           surv.median.line = "hv",
           palette = c("#E7B800","#2E9FDF"))




################################## 绘制校准曲线 ################################
library(rms)
############## 三元组 ############
# 打包数据
dd<-datadist(risk)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~risk,data = risk, x=T, y=T, surv=T)
# 3.生成1/3/5年的死亡概率
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)

# 4.绘制nomogram
Cox_nomo1<-cph(Surv(OS.time, OS) ~age+group,data = risk, x=T, y=T, surv=T)
nomo_2a<-nomogram(Cox_nomo1, fun=list(surv1,surv3,surv5), lp=F,
                  funlabel =c("1-year Death", "3-year Death", "5-year Death"),
                  fun.at =c(0.05, seq(0.1,0.9, by=0.1), 0.95)
)
plot(nomo_2a,
     #col.grid=c("pink","cyan"),
     xfrac =0.3, #设置变量名与线段的横向占比
     
     cex.var =1, # 加粗变量字体
     
     cex.axis =1, #设置数据轴字体的代销
     
     lmgp =0.3) # 设置文字与数据轴刻度的距离
## 校准曲线
# 一年校准曲线
coxm_triplet_1 <- cph(Surv(OS.time, OS) ~risk,data=risk,
              surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_triplet_1<-calibrate(coxm_triplet_1,u=365,cmethod='KM',m=mmm,B=1000)

# 三年校准曲线
coxm_triplet_2 <- cph(Surv(OS.time, OS) ~risk,data=risk,
              surv=T,x=T,y=T,time.inc = 3*365)
cal_triplet_2<-calibrate(coxm_triplet_2,u=3*365,cmethod='KM',m=mmm,B=1000)

# 五年校准曲线
coxm_triplet_3 <- cph(Surv(OS.time, OS) ~risk+group,data=risk,
              surv=T,x=T,y=T,time.inc = 5*365)
cal_triplet_3<-calibrate(coxm_triplet_3,u=5*365,cmethod='KM',m=mmm,B=1000)

############# lnc #############
exp=exp_pre[,c("OS.time", "OS","ENSG00000234883")]
dd<-datadist(exp)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp, x=T, y=T, surv=T)
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)
coxm_lnc_1 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
              surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_lnc_1<-calibrate(coxm_lnc_1,u=365,cmethod='KM',m=mmm,B=1000)
coxm_lnc_2 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
              surv=T,x=T,y=T,time.inc = 3*365)
cal_lnc_2<-calibrate(coxm_lnc_2,u=3*365,cmethod='KM',m=mmm,B=1000)
coxm_lnc_3 <- cph(Surv(OS.time, OS) ~ENSG00000234883,data = exp,
              surv=T,x=T,y=T,time.inc = 5*365)
cal_lnc_3<-calibrate(coxm_lnc_3,u=5*365,cmethod='KM',m=mmm,B=1000)

############## icp #############
exp=exp_pre[,c("OS.time", "OS","ENSG00000169245")]
dd<-datadist(exp)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~ENSG00000169245,data = exp, x=T, y=T, surv=T)
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)
coxm_icp_1 <- cph(Surv(OS.time, OS) ~ENSG00000169245,data = exp,
              surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_icp_1<-calibrate(coxm_icp_1,u=365,cmethod='KM',m=mmm,B=1000)
coxm_icp_2 <- cph(Surv(OS.time, OS) ~ENSG00000169245,data = exp,
              surv=T,x=T,y=T,time.inc = 3*365)
cal_icp_2<-calibrate(coxm_icp_2,u=3*365,cmethod='KM',m=mmm,B=1000)
coxm_icp_3 <- cph(Surv(OS.time, OS) ~ENSG00000169245,data = exp,
              surv=T,x=T,y=T,time.inc = 5*365)
cal_icp_3<-calibrate(coxm_icp_3,u=5*365,cmethod='KM',m=mmm,B=1000)

################# imm #################
exp=exp_pre[,c("OS.time", "OS","ENSG00000163600")]
dd<-datadist(exp)
options(datadist='dd')
Cox_nomo1<-cph(Surv(OS.time, OS) ~ENSG00000163600,data = exp, x=T, y=T, surv=T)
surv<-Survival(Cox_nomo1)
surv1 <-function(x)surv(1*365,lp=1-x) #“lp=1-x”计算死亡率，“lp=x”计算生存率
surv3 <-function(x)surv(1*1095,lp=1-x)
surv5 <-function(x)surv(1*1825,lp=1-x)
coxm_imm_1 <- cph(Surv(OS.time, OS) ~ENSG00000163600,data = exp,
              surv=T,x=T,y=T,time.inc = 365)
mmm<-floor(nrow(risk)/3)
cal_imm_1<-calibrate(coxm_imm_1,u=365,cmethod='KM',m=mmm,B=1000)
coxm_imm_2 <- cph(Surv(OS.time, OS) ~ENSG00000163600,data = exp,
              surv=T,x=T,y=T,time.inc = 3*365)
cal_imm_2<-calibrate(coxm_imm_2,u=3*365,cmethod='KM',m=mmm,B=1000)
coxm_imm_3 <- cph(Surv(OS.time, OS) ~ENSG00000163600,data = exp,
              surv=T,x=T,y=T,time.inc = 5*365)
cal_imm_3<-calibrate(coxm_imm_3,u=5*365,cmethod='KM',m=mmm,B=1000)

################ 绘图 ###############
# 1年
plot(cal_triplet_1,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 1-year OS',#便签
     ylab='Actual 1-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.85,0.95),ylim = c(0.85,0.95),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)

plot(cal_lnc_1,lwd=2,lty=0,
     errbar.col=c("#5ebd9a"),
     col=c("#5ebd9a"),#设置一个颜色
     xlim = c(0.85,0.95),ylim = c(0.85,0.95),
     riskdist = F,subtitles = F,
     add = T
     )
lines(cal_lnc_1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5ebd9a"), pch = 16)

plot(cal_icp_1,lwd=2,lty=0,
     errbar.col=c("#6e9bc5"),
     col=c("#6e9bc5"),#设置一个颜色
     xlim = c(0.85,0.95),ylim = c(0.85,0.95),
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_icp_1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#6e9bc5"), pch = 16)

plot(cal_imm_1,lwd=2,lty=0,
     errbar.col=c("#f8ce6a"),
     col=c("#f8ce6a"),#设置一个颜色
     xlim = c(0.85,0.95),ylim = c(0.85,0.95),
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_imm_1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#f8ce6a"), pch = 16)

legend("bottomright", #图例的位置
       legend = c("MIR155HG;CXCL10;ICOS","MIR155HG","CXCL10","ICOS"), #图例文字
       col =c("#eb4d4b","#5ebd9a","#6e9bc5","#f8ce6a"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框


# 3年
plot(cal_triplet_2,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 3-year OS',#便签
     ylab='Actual 3-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.65,0.9),ylim = c(0.65,0.9),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)

plot(cal_lnc_2,lwd=2,lty=0,
     errbar.col=c("#5ebd9a"),
     col=c("#5ebd9a"),#设置一个颜色
     xlim = c(0.65,0.9),ylim = c(0.65,0.9),
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_lnc_2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5ebd9a"), pch = 16)

plot(cal_icp_2,lwd=2,lty=0,
     errbar.col=c("#6e9bc5"),
     col=c("#6e9bc5"),#设置一个颜色
     xlim = c(0.65,0.9),ylim = c(0.65,0.9),
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_icp_2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#6e9bc5"), pch = 16)

plot(cal_imm_2,lwd=2,lty=0,
     errbar.col=c("#f8ce6a"),
     col=c("#f8ce6a"),#设置一个颜色
     xlim = c(0.65,0.9),ylim = c(0.65,0.9),
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_imm_2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#f8ce6a"), pch = 16)

legend("bottomright", #图例的位置
       legend = c("MIR155HG;CXCL10;ICOS","MIR155HG","CXCL10","ICOS"), #图例文字
       col =c("#eb4d4b","#5ebd9a","#6e9bc5","#f8ce6a"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框


# 5年
plot(cal_triplet_3,lwd=2,lty=0, ##设置线条形状和尺寸
     errbar.col=c("#eb4d4b"), ##设置一个颜色
     xlab='Predicted Probability of 5-year OS',#便签
     ylab='Actual 5-year OS(proportion)',#标签
     col=c("#eb4d4b"),#设置一个颜色
     xlim = c(0.5,0.8),ylim = c(0.5,0.8),
     riskdist = F,subtitles = F) ##x轴和y轴范围

lines(cal_triplet_3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#eb4d4b"), pch = 16)

plot(cal_lnc_3,lwd=2,lty=0,
     errbar.col=c("#5ebd9a"),
     col=c("#5ebd9a"),#设置一个颜色
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_lnc_3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#5ebd9a"), pch = 16)

plot(cal_icp_3,lwd=2,lty=0,
     errbar.col=c("#6e9bc5"),
     col=c("#6e9bc5"),#设置一个颜色
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_icp_3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#6e9bc5"), pch = 16)

plot(cal_imm_3,lwd=2,lty=0,
     errbar.col=c("#f8ce6a"),
     col=c("#f8ce6a"),#设置一个颜色
     riskdist = F,subtitles = F,
     add = T
)
lines(cal_imm_3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#f8ce6a"), pch = 16)

legend("bottomright", #图例的位置
       legend = c("MIR155HG;CXCL10;ICOS","MIR155HG","CXCL10","ICOS"), #图例文字
       col =c("#eb4d4b","#5ebd9a","#6e9bc5","#f8ce6a"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 0.7,#图例字体大小
       bty = "n")#不显示图例边框


