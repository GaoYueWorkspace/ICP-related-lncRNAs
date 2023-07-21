path=getwd()
setwd(path)
library(umap)
library(Rtsne)
library(readr)
library(tidyverse)

##############################
# 所有的结果lnc作为输入特征
# 所有的hublnc作为输入特征
# 所有的specific作为输入特征
#############################

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

# # 1、所有的lnc作为输入特征，不scale
# lnc_all=read.table("lncRNA.txt",header = T)
# result_df <- as.data.frame(lnc_all[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_all$ENSG,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))


# 2、使用speci_lnc,只在一个癌症中出现，不scale
# lnc_specific=read.table("specific_lnc1.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 3、使用speci_lnc,在≤2个癌症中出现，不scale
lnc_specific=read.table("specific_lnc2.txt",header = T)
result_df <- as.data.frame(lnc_specific[,1])
for (df in exp_all) {
  sub_df <- df[lnc_specific$x,]
  result_df <- cbind(result_df, sub_df)
}
result_df=result_df[,-1]
result_df=as.data.frame(t(result_df))
result_umap=umap(result_df)

plot.data = as.data.frame(result_umap$layout)
colnames(plot.data)=c("umap1","umap2")
head(plot.data)
plot.labels=as.factor(sample$label)
plot.data=cbind(plot.data,plot.labels)

color=read.table("颜色组.txt",comment.char = "")

plot.data %>%
  ggplot(aes(x = umap1, 
             y = umap2, 
             color = plot.labels))+
  geom_point(size=0.7,pch=16)+  # 16、19、20、21
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")+
  scale_color_manual(values = color$V2)

# 4、使用speci_lnc,在≤3个癌症中出现，不scale
# lnc_specific=read.table("specific_lnc3.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 5、使用speci_lnc,在≤4个癌症中出现，不scale
# lnc_specific=read.table("specific_lnc4.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 6、所有的lnc作为输入特征,scale
# lnc_all=read.table("lncRNA.txt",header = T)
# result_df <- as.data.frame(lnc_all[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_all$ENSG,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=scale(result_df)
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 7、使用speci_lnc,只在一个癌症中出现，scale
# lnc_specific=read.table("specific_lnc1.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=scale(result_df)
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 8、使用speci_lnc,在≤2个癌症中出现，scale
# lnc_specific=read.table("specific_lnc2.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=scale(result_df)
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 9、使用speci_lnc,在≤3个癌症中出现，scale
# lnc_specific=read.table("specific_lnc3.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=scale(result_df)
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

# 10、使用speci_lnc,在≤4个癌症中出现，scale
# lnc_specific=read.table("specific_lnc4.txt",header = T)
# result_df <- as.data.frame(lnc_specific[,1])
# for (df in exp_all) {
#   sub_df <- df[lnc_specific$x,]
#   result_df <- cbind(result_df, sub_df)
# }
# result_df=result_df[,-1]
# result_df=scale(result_df)
# result_df=as.data.frame(t(result_df))
# result_umap=umap(result_df)
# 
# plot.data = as.data.frame(result_umap$layout)
# colnames(plot.data)=c("umap1","umap2")
# head(plot.data)
# plot.labels=as.factor(sample$label)
# plot.data=cbind(plot.data,plot.labels)
# plot.data %>%
#   ggplot(aes(x = umap1, 
#              y = umap2, 
#              color = plot.labels))+
#   geom_point(size=0.7)+
#   labs(x = "UMAP1",
#        y = "UMAP2",
#        subtitle = "UMAP plot")+
#   scale_color_manual(values = rainbow(18))

############################### TSNE ###########################################
lnc_specific=read.table("specific_lnc2.txt",header = T)
result_df <- as.data.frame(lnc_specific[,1])
for (df in exp_all) {
  sub_df <- df[lnc_specific$x,]
  result_df <- cbind(result_df, sub_df)
}
result_df=result_df[,-1]
result_df=as.data.frame(t(result_df))
tsne_test = Rtsne(
  result_df,
  dims = 2,
  theta = 0,
  pca=T,  ## 调用pca结果，计算快。
  pca_center = TRUE,
  pca_scale = F,
  perplexity = 30, ##  30
  verbose = F) # 进行t-SNE降维分析
plot.data = as.data.frame(tsne_test$Y)
colnames(plot.data) = c("tSNE1","tSNE2")

plot.labels=as.factor(sample$label)
plot.data=cbind(plot.data,plot.labels)

color=read.table("颜色组.txt",comment.char = "")

plot.data %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2, 
             color = plot.labels))+
  geom_point(size=0.7,pch=16)+  # 16、19、20、21
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "t-SNE plot")+
  scale_color_manual(values = color$V2)
