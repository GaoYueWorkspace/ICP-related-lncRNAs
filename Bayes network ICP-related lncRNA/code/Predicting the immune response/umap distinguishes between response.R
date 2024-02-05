exp=read.table("Gide_PD1.txt",header = T)
exp_resp=read.table("Gide_PD1_clinical.txt",header = T)

path=getwd()

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
setwd(path)
triplet=do.call(rbind,triplet_all)
remove(triplet_all)
triplet=unique(triplet)

name=rownames(exp)
result_lnc=read.table("result_lncRNA.txt",header = T)
result_imm=read.table("result_imm.txt",header = T)
result_icp=read.table("result_icp.txt",header = T)
result_lnc=result_lnc[which(result_lnc$SYMBOL%in%name),]
result_icp=result_icp[which(result_icp$SYMBOL%in%name),]
result_imm=result_imm[which(result_imm$SYMBOL%in%name),]

triplet=triplet[which(triplet$lnc%in%result_lnc$ENSG),]
triplet=triplet[which(triplet$icp%in%result_icp$ENSG),]
triplet=triplet[which(triplet$imm%in%result_imm$ENSG),]

anno=read.table("SYMBOL_ENSEMBL.txt",header = T)
for(i in 1:nrow(triplet)){
  triplet[i,"lnc_anno"]=anno[which(anno$ENSG%in%triplet[i,2]),2]
  triplet[i,"icp_anno"]=anno[which(anno$ENSG%in%triplet[i,1]),2]
  triplet[i,"imm_anno"]=anno[which(anno$ENSG%in%triplet[i,3]),2]
  
}
triplet=triplet[,4:6]
rownames(triplet)=1:nrow(triplet)
triplet[,"id"]=paste(triplet$lnc_anno,triplet$icp_anno,triplet$imm_anno,sep = ";")
remove(result_icp,result_imm,result_lnc,anno,i,name)

library(survival)
library(umap)
library(ggplot2)
library(ggforce)
library(tidyverse)



gene=union(union(triplet$lnc_anno,triplet$icp_anno),triplet$imm_anno)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)


exp=as.data.frame(t(exp))


os=rbind(exp_resp[,2:3])

exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

triplet2=apply(triplet,1,function(x){gsub("-","\\.",x)})
triplet2=as.data.frame(t(triplet2))
triplet2=triplet2[,-4]

risk=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  return(score)
}

result=apply(triplet2,1,risk)
colnames(result)=triplet$id
result=as.data.frame(result)

exp_resp=c(exp_resp$Response)

set.seed(110)
umap=umap(result,method = 'naive', n_neighbors = 5)
# 提取umap值作图用
plot.data<- data.frame(umap$layout)
plot.data$label <- as.factor(exp_resp) # 加入label列

colnames(plot.data)=c("umap1","umap2","label")
head(plot.data)
plot.labels=plot.data$label

plot.data %>%
  ggplot(aes(x = umap1, 
             y = umap2, 
             color = plot.labels))+
  geom_point(size=2,pch=16)+  # 16、19、20、21
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")+
  stat_ellipse(level = 0.9)


