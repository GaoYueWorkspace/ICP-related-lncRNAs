library(data.table)
library(readr)
anno=fread("infiltration_estimation_for_tcga.csv.gz")
anno=as.data.frame(anno)
rownames(anno)=anno$cell_type
path=getwd()

# 构建样本及分类标签数据
cancer=function(exp){
  exp=as.data.frame(exp)
  exp[,1]=substr(exp[,1],start = 1,stop = 15)
  rownames(exp)=exp[,1]
  exp=exp[,-1]
  index=which(substr(colnames(exp),start = 14,stop = 15)>10)
  exp=exp[,-index]
  return(exp)
}

BLCA=read_tsv("BLCA/TCGA-BLCA.htseq_fpkm.tsv")
BLCA=cancer(BLCA)
sample=colnames(BLCA)
sample=as.data.frame(sample)
sample[,"label"]="BLCA"

BRCA=read_tsv("BRCA/TCGA-BRCA.htseq_fpkm.tsv")
BRCA=cancer(BRCA)
tmp=colnames(BRCA)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="BRCA"
sample=rbind(sample,tmp)

CHOL=read_tsv("CHOL/TCGA-CHOL.htseq_fpkm.tsv")
CHOL=cancer(CHOL)
tmp=colnames(CHOL)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="CHOL"
sample=rbind(sample,tmp)

COAD=read_tsv("COAD/TCGA-COAD.htseq_fpkm.tsv")
COAD=cancer(COAD)
tmp=colnames(COAD)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="COAD"
sample=rbind(sample,tmp)

ESCA=read_tsv("ESCA/TCGA-ESCA.htseq_fpkm.tsv")
ESCA=cancer(ESCA)
tmp=colnames(ESCA)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="ESCA"
sample=rbind(sample,tmp)

GBM=read_tsv("GBM/TCGA-GBM.htseq_fpkm.tsv.gz")
GBM=cancer(GBM)
tmp=colnames(GBM)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="GBM"
sample=rbind(sample,tmp)

HNSC=read_tsv("HNSC/TCGA-HNSC.htseq_fpkm.tsv")
HNSC=cancer(HNSC)
tmp=colnames(HNSC)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="HNSC"
sample=rbind(sample,tmp)

KIRC=read_tsv("KIRC/TCGA-KIRC.htseq_fpkm.tsv")
KIRC=cancer(KIRC)
tmp=colnames(KIRC)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="KIRC"
sample=rbind(sample,tmp)

LIHC=read_tsv("LIHC/TCGA-LIHC.htseq_fpkm.tsv")
LIHC=cancer(LIHC)
tmp=colnames(LIHC)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="LIHC"
sample=rbind(sample,tmp)

LUAD=read_tsv("LUAD/TCGA-LUAD.htseq_fpkm.tsv")
LUAD=cancer(LUAD)
tmp=colnames(LUAD)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="LUAD"
sample=rbind(sample,tmp)

LUSC=read_tsv("LUSC/TCGA-LUSC.htseq_fpkm.tsv")
LUSC=cancer(LUSC)
tmp=colnames(LUSC)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="LUSC"
sample=rbind(sample,tmp)

OV=read_tsv("OV/TCGA-OV.htseq_fpkm.tsv.gz")
OV=as.data.frame(OV)
OV[,1]=substr(OV[,1],start = 1,stop = 15)
rownames(OV)=OV[,1]
OV=OV[,-1]
tmp=colnames(OV)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="OV"
sample=rbind(sample,tmp)

PRAD=read_tsv("PRAD/TCGA-PRAD.htseq_fpkm.tsv")
PRAD=cancer(PRAD)
tmp=colnames(PRAD)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="PRAD"
sample=rbind(sample,tmp)

READ=read_tsv("READ/TCGA-READ.htseq_fpkm.tsv")
READ=cancer(READ)
tmp=colnames(READ)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="READ"
sample=rbind(sample,tmp)

SKCM=read_tsv("SKCM/TCGA-SKCM.htseq_fpkm.tsv")
SKCM=cancer(SKCM)
tmp=colnames(SKCM)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="SKCM"
sample=rbind(sample,tmp)

STAD=read_tsv("STAD/TCGA-STAD.htseq_fpkm.tsv")
STAD=cancer(STAD)
tmp=colnames(STAD)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="STAD"
sample=rbind(sample,tmp)

THCA=read_tsv("THCA/TCGA-THCA.htseq_fpkm.tsv")
THCA=cancer(THCA)
tmp=colnames(THCA)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="THCA"
sample=rbind(sample,tmp)

UCEC=read_tsv("UCEC/TCGA-UCEC.htseq_fpkm.tsv")
UCEC=cancer(UCEC)
tmp=colnames(UCEC)
tmp=as.data.frame(tmp)
colnames(tmp)="sample"
tmp[,"label"]="UCEC"
sample=rbind(sample,tmp)

exp_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
             "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
             "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
             "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC,tmp)
remove(cancer)
setwd(path)
sample$sample=substr(sample$sample,1,15)
save(sample,file="sample.RData")

###############################################################################

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
save(exp_all,file="exp_all.RData")
remove(remove_last_char)
###############################################################################
cell=colnames(anno)
cell=as.data.frame(cell)


hub_lnc=read.table("hub_lnc.txt",header = T)
hub_lnc=hub_lnc[which(hub_lnc$pubmed>18),]
hub_lnc=hub_lnc[which(hub_lnc$pubmed!=21),]

anno=anno[,c(1,2,3,4,5,6)]

other_lnc=read.table("other_lnc.txt",header = T)
intersect(hub_lnc$lnc,other_lnc$lnc)
# hub_lnc=read.table("other_lnc.txt",header = T)
# gene=rownames(exp_all$BLCA)
# hub_lnc=setdiff(gene,hub_lnc$lnc)
# hub_lnc=hub_lnc[sample(1:60438,100,replace=F)]
# hub_lnc=as.data.frame(hub_lnc)
# colnames(hub_lnc)="lnc"

### B Cell

## BLCA
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BLCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BRCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"CHOL"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"COAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"ESCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"GBM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"HNSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"KIRC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LIHC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"OV"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"PRAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"READ"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"SKCM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"STAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"THCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"UCEC"),1])
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

plot.data.cor=list("B Cell"=dat_cor)
plot.data.p=list("B Cell"=dat_p)

##############################################################################
### CD4+ Cell

## BLCA
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BLCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BRCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"CHOL"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"COAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"ESCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"GBM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"HNSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"KIRC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LIHC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"OV"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"PRAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"READ"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"SKCM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"STAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"THCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"UCEC"),1])
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


tmp.cor=list("CD4+ Cell"=dat_cor)
tmp.p=list("CD4+ Cell"=dat_p)

plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)

##############################################################################

### CD8+ Cell

## BLCA
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BLCA"),1])
dat1=exp_all$BLCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[3,]), method = "spearman")
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BRCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"CHOL"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"COAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"ESCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"GBM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"HNSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"KIRC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LIHC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"OV"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"PRAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"READ"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"SKCM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"STAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"THCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"UCEC"),1])
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


tmp.cor=list("CD8+ Cell"=dat_cor)
tmp.p=list("CD8+ Cell"=dat_p)

plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)

##############################################################################
### Neutrophil Cell

## BLCA
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BLCA"),1])
dat1=exp_all$BLCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[4,]), method = "spearman")
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BRCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"CHOL"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"COAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"ESCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"GBM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"HNSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"KIRC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LIHC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"OV"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"PRAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"READ"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"SKCM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"STAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"THCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"UCEC"),1])
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


tmp.cor=list("Neutrophil"=dat_cor)
tmp.p=list("Neutrophil"=dat_p)

plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)

###############################################################################

### Macrophage

## BLCA
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BLCA"),1])
dat1=exp_all$BLCA[,coo_sample]
dat1=dat1[hub_lnc$lnc,]
dat2=anno[coo_sample,]
dat2=as.data.frame(t(dat2))
dat2=dat2[-1,]
s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[5,]), method = "spearman")
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"BRCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"CHOL"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"COAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"ESCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"GBM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"HNSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"KIRC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LIHC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"LUSC"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"OV"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"PRAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"READ"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"SKCM"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"STAD"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"THCA"),1])
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
coo_sample=intersect(anno$cell_type,sample[which(sample$label%in%"UCEC"),1])
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


tmp.cor=list("Macrophage"=dat_cor)
tmp.p=list("Macrophage"=dat_p)

plot.data.cor=append(plot.data.cor,tmp.cor)
plot.data.p=append(plot.data.p,tmp.p)


heatdata=as.matrix(plot.data.cor$`B Cell`)
anno_melt=as.matrix(plot.data.p$`B Cell`)

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


heatdata=as.matrix(plot.data.cor$`CD4+ Cell`)
anno_melt=as.matrix(plot.data.p$`CD4+ Cell`)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

heatdata=as.matrix(plot.data.cor$`CD8+ Cell`)
anno_melt=as.matrix(plot.data.p$`CD8+ Cell`)
library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

heatdata=as.matrix(plot.data.cor$Neutrophil)
anno_melt=as.matrix(plot.data.p$Neutrophil)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

heatdata=as.matrix(plot.data.cor$Macrophage)
anno_melt=as.matrix(plot.data.p$Macrophage)

library(ggplot2)
library(reshape2)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")

ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)

