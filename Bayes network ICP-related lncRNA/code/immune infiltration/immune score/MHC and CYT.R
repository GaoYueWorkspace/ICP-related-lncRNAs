library(readxl)
library(readr)
load("exp_all.RData")
load("sample.RData")
anno=read_csv("33个癌症的MHC和CYT.csv")
anno=as.data.frame(anno)
rownames(anno)=substr(anno$Sample,1,15)
anno$Sample=substr(anno$Sample,1,15)
hub_lnc=read.table("hub_lnc.txt",header = T)
hub_lnc=hub_lnc[which(hub_lnc$pubmed>18),]
hub_lnc=hub_lnc[which(hub_lnc$pubmed!=21),]

other_lnc=read.table("other_lnc.txt",header = T)
intersect(hub_lnc$lnc,other_lnc$lnc)
# hub_lnc=other_lnc
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

##################################### MHC ####################################
## BLCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"BLCA"),1])
dat1=exp_all$BLCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "spearman")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}

result=as.data.frame(t(apply(dat1,1,s_cor)))
dat_cor=result$R.rho
dat_cor=as.data.frame(dat_cor)
dat_cor=as.data.frame(t(dat_cor))
rownames(dat_cor)="BLCA"
colnames(dat_cor)=rownames(result)

dat_p=result$p
dat_p=as.data.frame(dat_p)
dat_p=as.data.frame(t(dat_p))
rownames(dat_p)="BLCA"
colnames(dat_p)=rownames(result)

## BRCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"BRCA"),1])
dat1=exp_all$BRCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="BRCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="BRCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## CHOL
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"CHOL"),1])
dat1=exp_all$CHOL[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="CHOL"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="CHOL"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## COAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"COAD"),1])
dat1=exp_all$COAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="COAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="COAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## ESCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"ESCA"),1])
dat1=exp_all$ESCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="ESCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="ESCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## GBM
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"GBM"),1])
dat1=exp_all$GBM[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="GBM"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="GBM"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## HNSC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"HNSC"),1])
dat1=exp_all$HNSC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="HNSC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="HNSC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## KIRC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"KIRC"),1])
dat1=exp_all$KIRC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="KIRC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="KIRC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LIHC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LIHC"),1])
dat1=exp_all$LIHC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LIHC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LIHC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LUAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LUAD"),1])
dat1=exp_all$LUAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LUSC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LUSC"),1])
dat1=exp_all$LUSC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUSC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUSC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## OV
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"OV"),1])
dat1=exp_all$OV[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="OV"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="OV"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## PRAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"PRAD"),1])
dat1=exp_all$PRAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="PRAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="PRAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## READ
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"READ"),1])
dat1=exp_all$READ[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="READ"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="READ"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## SKCM
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"SKCM"),1])
dat1=exp_all$SKCM[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="SKCM"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="SKCM"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## STAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"STAD"),1])
dat1=exp_all$STAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="STAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="STAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## THCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"THCA"),1])
dat1=exp_all$THCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="THCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="THCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## UCEC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"UCEC"),1])
dat1=exp_all$UCEC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="UCEC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="UCEC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

plot.data.cor=list("MHC"=dat_cor)
plot.data.p=list("MHC"=dat_p)

############################### CYT #########################################
## BLCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"BLCA"),1])
dat1=exp_all$BLCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "spearman")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}

result=as.data.frame(t(apply(dat1,1,s_cor)))
dat_cor=result$R.rho
dat_cor=as.data.frame(dat_cor)
dat_cor=as.data.frame(t(dat_cor))
rownames(dat_cor)="BLCA"
colnames(dat_cor)=rownames(result)

dat_p=result$p
dat_p=as.data.frame(dat_p)
dat_p=as.data.frame(t(dat_p))
rownames(dat_p)="BLCA"
colnames(dat_p)=rownames(result)

## BRCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"BRCA"),1])
dat1=exp_all$BRCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="BRCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="BRCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## CHOL
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"CHOL"),1])
dat1=exp_all$CHOL[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="CHOL"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="CHOL"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## COAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"COAD"),1])
dat1=exp_all$COAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="COAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="COAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## ESCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"ESCA"),1])
dat1=exp_all$ESCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="ESCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="ESCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## GBM
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"GBM"),1])
dat1=exp_all$GBM[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="GBM"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="GBM"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## HNSC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"HNSC"),1])
dat1=exp_all$HNSC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="HNSC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="HNSC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## KIRC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"KIRC"),1])
dat1=exp_all$KIRC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="KIRC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="KIRC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LIHC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LIHC"),1])
dat1=exp_all$LIHC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LIHC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LIHC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LUAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LUAD"),1])
dat1=exp_all$LUAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## LUSC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"LUSC"),1])
dat1=exp_all$LUSC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUSC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="LUSC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## OV
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"OV"),1])
dat1=exp_all$OV[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="OV"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="OV"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## PRAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"PRAD"),1])
dat1=exp_all$PRAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="PRAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="PRAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## READ
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"READ"),1])
dat1=exp_all$READ[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="READ"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="READ"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## SKCM
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"SKCM"),1])
dat1=exp_all$SKCM[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="SKCM"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="SKCM"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## STAD
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"STAD"),1])
dat1=exp_all$STAD[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="STAD"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="STAD"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## THCA
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"THCA"),1])
dat1=exp_all$THCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="THCA"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="THCA"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)

## UCEC
coo_sample=intersect(anno$Sample,sample[which(sample$label%in%"UCEC"),1])
dat1=exp_all$UCEC[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]

result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp=result$R.rho
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="UCEC"
colnames(tmp)=rownames(result)
dat_cor=rbind(dat_cor,tmp)

tmp=result$p
tmp=as.data.frame(tmp)
tmp=as.data.frame(t(tmp))
rownames(tmp)="UCEC"
colnames(tmp)=rownames(result)
dat_p=rbind(dat_p,tmp)


tmp.cor=list("CYT"=dat_cor)
tmp.p=list("CYT"=dat_p)
plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)


heatdata=as.matrix(plot.data.cor$MHC)
anno_melt=as.matrix(plot.data.p$MHC)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

heatdata=as.matrix(plot.data.cor$CYT)
anno_melt=as.matrix(plot.data.p$CYT)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)
