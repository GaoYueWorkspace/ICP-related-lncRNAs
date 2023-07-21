result_lnc=read.table("result_lncRNA.txt",header = T)
result_imm=read.table("result_imm.txt",header = T)
result_icp=read.table("result_icp.txt",header = T)

### GSE127471
GSE127471=read.table("GSE127471.txt",header = T)
gene=intersect(result_lnc$SYMBOL,rownames(GSE127471))
GSE127471=GSE127471[gene,]
GSE127471_anno=read.table("GSE127471_anno.txt",header = T,sep = "\t")
group=as.factor(GSE127471_anno$cell)
GSE127471=as.data.frame(t(GSE127471))
exp=cbind(GSE127471,group)
test=function(x){
  a=with(GSE127471,kruskal.test(x , g=group, data = GSE127471,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE127471,2,test)
result=as.data.frame(result)
write.table(result,"GSE127471_lnc.txt",quote = F,sep = "\t")


exp$group=as.character(exp$group)
gene=aggregate(.~group,data = exp,mean)

max_cell=function(x){
  max_c=max(x)
  index=which(x%in%max_c)
  return(gene$group[index])
}
cell_high=apply(gene,2,max_cell)
cell_high=cbind(colnames(gene),cell_high)
cell_high=cell_high[-1,]
cell_high=as.data.frame(cell_high)
cell_high["MIR155HG",]

###############################################################################
GSE72056=read.table("GSE72056.txt",header = T)
gene=intersect(result_lnc$SYMBOL,rownames(GSE72056))
GSE72056=GSE72056[gene,]
GSE72056_anno=read.table("GSE72056_anno.txt",header = T,sep = "\t")
group=as.factor(GSE72056_anno$cell)
GSE72056=as.data.frame(t(GSE72056))
exp=cbind(GSE72056,group)
test=function(x){
  a=with(GSE72056,kruskal.test(x , g=group, data = GSE72056,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE72056,2,test)
result=as.data.frame(result)
write.table(result,"GSE72056_lnc.txt",quote = F,sep = "\t")

exp$group=as.character(exp$group)
gene=aggregate(.~group,data = exp,mean)

max_cell=function(x){
  max_c=max(x)
  index=which(x%in%max_c)
  return(gene$group[index])
}
cell_high=apply(gene,2,max_cell)
cell_high=cbind(colnames(gene),cell_high)
cell_high=cell_high[-1,]
cell_high=as.data.frame(cell_high)
cell_high["MIR155HG",]

plot.data=exp[,"MIR155HG"]
plot.data=cbind(plot.data,exp$group)
colnames(plot.data)=c("value","group")
plot.data=as.data.frame(plot.data)
plot.data$value=as.numeric(plot.data$value)
tmp=aggregate(value~group,plot.data,mean)
tmp=tmp[order(tmp$value),]
plot.data$group=factor(plot.data$group,levels = tmp$group)
library(vioplot)
plot.data=plot.data[which(plot.data$value!=0),]
vioplot(value~group, data = plot.data,
        col=c("#d67d79","#9fb9cd","#c5d832","#ba86b4","#d57732","#e4b3b3"),las=2)
result["MIR155HG",]

################################ GSE117570 #####################################
library(tidyverse)
GSE117570=read.table("GSE117570.txt",header = T)
anno=read.table("单细胞数据/SYMBOL_ENSEMBL.txt",header = T)
load("genelength.RData")
gene=gene[which(gene$gene_id%in%anno$ENSG),]
for(i in 1:nrow(gene)){
  gene[i,"gene_name"]=anno[which(anno$ENSG%in%gene[i,1]),2]
}

GSE117570[,"gene_name"]=rownames(GSE117570)
tmp <- inner_join(GSE117570,gene,by = 'gene_name') # 匹配gene_count与gene_length
temple=tmp[,1:1931]
temple%>% mutate(across(where(is.character), as.numeric))  -> temple
tmp=cbind(tmp[,1932],temple,tmp[,1934])
colnames(tmp)[1]="gene_name"
colnames(tmp)[1933]="length"
tmp=tmp[!duplicated(tmp$gene_name),]
fpkm <- data.frame(row.names = tmp$gene_name)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}
GSE117570=log2(1+fpkm)

gene=intersect(result_lnc$SYMBOL,rownames(GSE117570))
GSE117570=GSE117570[gene,]
GSE117570_anno=read.table("GSE117570_anno.txt",header = T,sep = "\t")
group=as.factor(GSE117570_anno$cell)
GSE117570=as.data.frame(t(GSE117570))
exp=cbind(GSE117570,group)
test=function(x){
  a=with(GSE117570,kruskal.test(x , g=group, data = GSE117570,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE117570,2,test)
result=as.data.frame(result)
write.table(result,"GSE117570_lnc.txt",quote = F,sep = "\t")

exp$group=as.character(exp$group)
gene=aggregate(.~group,data = exp,mean)

max_cell=function(x){
  max_c=max(x)
  index=which(x%in%max_c)
  return(gene$group[index])
}
cell_high=apply(gene,2,max_cell)
cell_high=cbind(colnames(gene),cell_high)
cell_high=cell_high[-1,]
cell_high=as.data.frame(cell_high)
cell_high["MIR155HG",]
plot.data=exp[,"MIR155HG"]
plot.data=cbind(plot.data,exp$group)
colnames(plot.data)=c("value","group")
plot.data=as.data.frame(plot.data)
plot.data$value=as.numeric(plot.data$value)
tmp=aggregate(value~group,plot.data,mean)
tmp=tmp[order(tmp$value),]
plot.data$group=factor(plot.data$group,levels = tmp$group)
library(vioplot)
plot.data=plot.data[which(plot.data$value!=0),]
vioplot(value~group, data = plot.data,
        col=c("#d67d79","#9fb9cd","#c5d832","#ba86b4","#d57732","#e4b3b3"),las=2)
result["MIR155HG",]



#####################################################################################

### GSE127471
GSE127471=read.table("GSE127471.txt",header = T)
gene=intersect(result_imm$SYMBOL,rownames(GSE127471))
GSE127471=GSE127471[gene,]
GSE127471_anno=read.table("GSE127471_anno.txt",header = T,sep = "\t")
group=as.factor(GSE127471_anno$cell)
GSE127471=as.data.frame(t(GSE127471))
exp=cbind(GSE127471,group)
test=function(x){
  a=with(GSE127471,kruskal.test(x , g=group, data = GSE127471,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE127471,2,test)
result=as.data.frame(result)
write.table(result,"GSE127471_imm.txt",quote = F,sep = "\t")

###############################################################################
GSE72056=read.table("GSE72056.txt",header = T)
gene=intersect(result_imm$SYMBOL,rownames(GSE72056))
GSE72056=GSE72056[gene,]
GSE72056_anno=read.table("GSE72056_anno.txt",header = T,sep = "\t")
group=as.factor(GSE72056_anno$cell)
GSE72056=as.data.frame(t(GSE72056))
exp=cbind(GSE72056,group)
test=function(x){
  a=with(GSE72056,kruskal.test(x , g=group, data = GSE72056,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE72056,2,test)
result=as.data.frame(result)
write.table(result,"GSE72056_imm.txt",quote = F,sep = "\t")

################################ GSE117570 #####################################
library(tidyverse)
GSE117570=read.table("GSE117570.txt",header = T)
anno=read.table("单细胞数据/SYMBOL_ENSEMBL.txt",header = T)
load("genelength.RData")
gene=gene[which(gene$gene_id%in%anno$ENSG),]
for(i in 1:nrow(gene)){
  gene[i,"gene_name"]=anno[which(anno$ENSG%in%gene[i,1]),2]
}

GSE117570[,"gene_name"]=rownames(GSE117570)
tmp <- inner_join(GSE117570,gene,by = 'gene_name') # 匹配gene_count与gene_length
temple=tmp[,1:1931]
temple%>% mutate(across(where(is.character), as.numeric))  -> temple
tmp=cbind(tmp[,1932],temple,tmp[,1934])
colnames(tmp)[1]="gene_name"
colnames(tmp)[1933]="length"
tmp=tmp[!duplicated(tmp$gene_name),]
fpkm <- data.frame(row.names = tmp$gene_name)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}
GSE117570=log2(1+fpkm)

gene=intersect(result_imm$SYMBOL,rownames(GSE117570))
GSE117570=GSE117570[gene,]
GSE117570_anno=read.table("GSE117570_anno.txt",header = T,sep = "\t")
group=as.factor(GSE117570_anno$cell)
GSE117570=as.data.frame(t(GSE117570))
exp=cbind(GSE117570,group)
test=function(x){
  a=with(GSE117570,kruskal.test(x , g=group, data = GSE117570,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE117570,2,test)
result=as.data.frame(result)
write.table(result,"GSE117570_imm.txt",quote = F,sep = "\t")


###################################################################################

### GSE127471
GSE127471=read.table("GSE127471.txt",header = T)
gene=intersect(result_icp$SYMBOL,rownames(GSE127471))
GSE127471=GSE127471[gene,]
GSE127471_anno=read.table("GSE127471_anno.txt",header = T,sep = "\t")
group=as.factor(GSE127471_anno$cell)
GSE127471=as.data.frame(t(GSE127471))
exp=cbind(GSE127471,group)
test=function(x){
  a=with(GSE127471,kruskal.test(x , g=group, data = GSE127471,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE127471,2,test)
result=as.data.frame(result)
write.table(result,"GSE127471_icp.txt",quote = F,sep = "\t")

###############################################################################
GSE72056=read.table("GSE72056.txt",header = T)
gene=intersect(result_icp$SYMBOL,rownames(GSE72056))
GSE72056=GSE72056[gene,]
GSE72056_anno=read.table("GSE72056_anno.txt",header = T,sep = "\t")
group=as.factor(GSE72056_anno$cell)
GSE72056=as.data.frame(t(GSE72056))
exp=cbind(GSE72056,group)
test=function(x){
  a=with(GSE72056,kruskal.test(x , g=group, data = GSE72056,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE72056,2,test)
result=as.data.frame(result)
write.table(result,"GSE72056_icp.txt",quote = F,sep = "\t")

################################ GSE117570 #####################################
library(tidyverse)
GSE117570=read.table("GSE117570.txt",header = T)
anno=read.table("单细胞数据/SYMBOL_ENSEMBL.txt",header = T)
load("genelength.RData")
gene=gene[which(gene$gene_id%in%anno$ENSG),]
for(i in 1:nrow(gene)){
  gene[i,"gene_name"]=anno[which(anno$ENSG%in%gene[i,1]),2]
}

GSE117570[,"gene_name"]=rownames(GSE117570)
tmp <- inner_join(GSE117570,gene,by = 'gene_name') # 匹配gene_count与gene_length
temple=tmp[,1:1931]
temple%>% mutate(across(where(is.character), as.numeric))  -> temple
tmp=cbind(tmp[,1932],temple,tmp[,1934])
colnames(tmp)[1]="gene_name"
colnames(tmp)[1933]="length"
tmp=tmp[!duplicated(tmp$gene_name),]
fpkm <- data.frame(row.names = tmp$gene_name)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}
GSE117570=log2(1+fpkm)

gene=intersect(result_icp$SYMBOL,rownames(GSE117570))
GSE117570=GSE117570[gene,]
GSE117570_anno=read.table("GSE117570_anno.txt",header = T,sep = "\t")
group=as.factor(GSE117570_anno$cell)
GSE117570=as.data.frame(t(GSE117570))
exp=cbind(GSE117570,group)
test=function(x){
  a=with(GSE117570,kruskal.test(x , g=group, data = GSE117570,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(GSE117570,2,test)
result=as.data.frame(result)
write.table(result,"GSE117570_icp.txt",quote = F,sep = "\t")
