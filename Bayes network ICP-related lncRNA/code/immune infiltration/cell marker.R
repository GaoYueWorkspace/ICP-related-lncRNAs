path=getwd()
anno=read.table("SYMBOL_ENSEMBL.txt",sep = "\t",header = T)
lung_TCell=read.table("Lung_TCell.txt",header = T,sep="\t")
lung_TCell=lung_TCell[which(lung_TCell$Cell.marker%in%anno$SYMBOL),]
for(i in 1:nrow(lung_TCell)){
  lung_TCell[i,"Ensembl"]=anno[which(anno$SYMBOL%in%lung_TCell[i,4]),1]
}
write.table(lung_TCell,"Lung_TCell.txt",row.names = F,quote = F,sep = "\t")
###############################################################################
lung_CD4=read.table("Lung_CD4.txt",header = T,sep="\t")
lung_CD4=lung_CD4[which(lung_CD4$Cell.marker%in%anno$SYMBOL),]
for(i in 1:nrow(lung_CD4)){
  lung_CD4[i,"Ensembl"]=anno[which(anno$SYMBOL%in%lung_CD4[i,4]),1]
}
write.table(lung_CD4,"Lung_CD4.txt",row.names = F,quote = F,sep = "\t")
###############################################################################
lung_CD8=read.table("Lung_CD8.txt",header = T,sep="\t")
lung_CD8=lung_CD8[which(lung_CD8$Cell.marker%in%anno$SYMBOL),]
for(i in 1:nrow(lung_CD8)){
  lung_CD8[i,"Ensembl"]=anno[which(anno$SYMBOL%in%lung_CD8[i,4]),1]
}
write.table(lung_CD8,"Lung_CD8.txt",row.names = F,quote = F,sep = "\t")
################################################################################
lung_BCell=read.table("Lung_BCell.txt",header = T,sep="\t")
lung_BCell=lung_BCell[which(lung_BCell$Cell.marker%in%anno$SYMBOL),]
for(i in 1:nrow(lung_BCell)){
  lung_BCell[i,"Ensembl"]=anno[which(anno$SYMBOL%in%lung_BCell[i,4]),1]
}
write.table(lung_BCell,"Lung_BCell.txt",row.names = F,quote = F,sep = "\t")
################################################################################
load("exp_all.RData")
hub_lnc=read.table("hub_lnc.txt",header = T)
hub_lnc=hub_lnc[which(hub_lnc$pubmed>10),]
lnc_anno=read.table("lncRNA_ENSG.txt")
lnc_anno=lnc_anno[which(lnc_anno$V1%in%rownames(exp_all$BLCA)),]

a=sample(setdiff(lnc_anno$V1,hub_lnc$lnc),4)
other_lnc=data.frame("lnc"=a)
for(i in 1:4){
  other_lnc[i,"SYMBOL"]=lnc_anno[which(lnc_anno$V1%in%other_lnc[i,1]),2]
}
hub_lnc=hub_lnc[,c(1,3)]
hub_lnc=rbind(hub_lnc,other_lnc)
  
lung_TCell=lung_TCell[which(lung_TCell$Cell.marker%in%c("CD28","CTLA4","STAT1")),]
lung_BCell=lung_BCell[which(lung_BCell$Cell.marker%in%c("CD19","TCL1A","BLK")),]
lung_CD8=lung_CD8[-2,]

### Tcell
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_TCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD28"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD28"
colnames(tmp_p)=hub_lnc$SYMBOL
plot.data.cor=tmp_cor
plot.data.p=tmp_p

###########################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_TCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CTLA4"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CTLA4"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)
#####################################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_TCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[3,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="STAT1"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="STAT1"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)

### BCell
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_BCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="BLK"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="BLK"
colnames(tmp_p)=hub_lnc$SYMBOL
plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)

###########################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_BCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD19"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD19"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)
#####################################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_BCell$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[3,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="TCL1A"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="TCL1A"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)
####################################

### CD4
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_CD4$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD3D"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD3D"
colnames(tmp_p)=hub_lnc$SYMBOL
plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)

###########################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_CD4$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD4"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD4"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)
###########################################

### CD8
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_CD8$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[1,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CCL5"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CCL5"
colnames(tmp_p)=hub_lnc$SYMBOL
plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)

###########################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_CD8$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[2,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD3E"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD3E"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)
###########################
dat1=exp_all$LUAD[hub_lnc$lnc,]
dat2=exp_all$LUAD[lung_CD8$Ensembl,]

s_cor=function(x){
  test=cor.test(as.numeric(x), as.numeric(dat2[3,]), method = "pearson")
  estimate=test$estimate
  pvalue=test$p.value
  r<-c(R=estimate,p=pvalue)
  return(r)
}
result=as.data.frame(t(apply(dat1,1,s_cor)))
tmp_cor=result$R.cor
tmp_cor=as.data.frame(tmp_cor)
tmp_cor=as.data.frame(t(tmp_cor))
rownames(tmp_cor)="CD8A"
colnames(tmp_cor)=hub_lnc$SYMBOL

tmp_p=result$p
tmp_p=as.data.frame(tmp_p)
tmp_p=as.data.frame(t(tmp_p))
rownames(tmp_p)="CD8A"
colnames(tmp_p)=hub_lnc$SYMBOL

plot.data.cor=rbind(plot.data.cor,tmp_cor)
plot.data.p=rbind(plot.data.p,tmp_p)

heatdata=as.matrix(plot.data.cor)
anno_melt=as.matrix(plot.data.p)
heatdata_melt=melt(heatdata)
anno_melt=melt(anno_melt)
heatdata_melt$star <- ifelse(anno_melt$value< 0.05, "*", "")
library(ggplot2)
library(reshape2)

ggplot(data = heatdata_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colors = c("#322b80","white","#d94c4c"),
                       limits=c(-1,1)) +  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1))+
  coord_fixed()+
  geom_text(aes(label = star), size = 4)



ggplot()+
  geom_point(data = heatdata_melt,aes(Var2, Var1, size = value,fill=value),
             shape=21,color="#c4bcba")+
  scale_fill_gradient2(low="#4178a7",mid="#ffffff",high = "#e2201c",limits=c(-1,1))+
  scale_x_discrete(expand = c(0.025,0.025))+
  scale_y_discrete(expand = c(0.025,0.025))+
  geom_vline(aes(xintercept=seq(0.5,20.5,1)),color="#bbbbbb")+
  geom_hline(aes(yintercept=seq(0.5,20.5,1)),color="#bbbbbb")+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 6, hjust = 1))
  