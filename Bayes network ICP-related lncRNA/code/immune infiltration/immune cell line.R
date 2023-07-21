# 加载包
library(dplyr)
library(tidyverse)
# 读取基因列表和8个细胞系表达谱数据，创建数据框
gene_list = read.table("result_lncRNA.txt",header = T)[,2]
cell_lines=list.files("./免疫细胞系数据")
cell_lines=sub(".txt", "", cell_lines)
exp_data <- lapply(cell_lines, function(x) read.table(paste0("./免疫细胞系数据/", x,".txt"), header = TRUE)) # 读取每个细胞系的表达谱数据

group=lapply(exp_data,ncol)
group=as.data.frame(group)
group=group-1
group=as.data.frame(t(group))
rownames(group)=1:19
group=rep(cell_lines,group$V1)
group=as.factor(group)
# group <- as.factor(group)
# write.table(levels(group),"anno.txt",row.names = F,quote = F)

exp_data <- lapply(exp_data, function(x) {
  rownames(x)=x$ID_REF
  x=x[,-1]
})

exp_data <- Reduce(cbind, exp_data) # 合并数据
exp_data <- subset(exp_data, rownames(exp_data) %in% gene_list) # 筛选出基因列表中的基因


library(ggplot2) 
library(tidyr)
exp_data=as.data.frame(t(exp_data))
exp_data=cbind(exp_data,group)

# a=exp_data[,c(1,232)]
# ggplot(exp_data,aes(x=group,y=exp_data$ARIH2OS))+geom_boxplot(outlier.shape = NA)

test=function(x){
  a=with(exp_data,kruskal.test(x , g=group, data = exp_data,p.adj="bon",group=T))
  return(a$p.value)
}
result=apply(exp_data,2,test)
result=as.data.frame(result)

boxplot(exp_data[,"PSMB8-AS1"]~group,data = exp_data)
boxplot(exp_data[,"FGD5-AS1"]~group,data = exp_data)
boxplot(exp_data[,"LINC01614"]~group,data = exp_data)
boxplot(exp_data[,"MAGI2-AS3"]~group,data = exp_data)
boxplot(exp_data[,"MIR155HG"]~group,data = exp_data)
boxplot(exp_data[,"SOCAR"]~group,data = exp_data)
boxplot(exp_data[,"ADAMTS9-AS2"]~group,data = exp_data)
boxplot(exp_data[,"LINC01094"]~group,data = exp_data)
boxplot(exp_data[,"MIR100HG"]~group,data = exp_data)
boxplot(exp_data[,"USP30-AS1"]~group,data = exp_data)



exp_data$group=as.character(exp_data$group)
gene=aggregate(.~group,data = exp_data,mean)

max_cell=function(x){
  max_c=max(x)
  index=which(x%in%max_c)
  return(gene$group[index])
}
cell_high=apply(gene,2,max_cell)
cell_high=cbind(colnames(gene),cell_high)
cell_high=cell_high[-1,]
select_cell=c("CD4 T cell activated","CD4 T cell resting","CD8 T cell activated",
              "CD8 T cell resting","NKT activated","T gamma delta","T helper 17")
cell_high=as.data.frame(cell_high)
T_cell_lnc=cell_high[which(cell_high$cell_high%in%select_cell),1]
lnc_anno=read.table("lncRNA_ENSG.txt")
T_cell_lnc=as.data.frame(T_cell_lnc)
for (i in 1:69) {
  T_cell_lnc[i,"ENSEMBL"]=lnc_anno[which(lnc_anno$V2%in%T_cell_lnc[i,1]),1]
}
T_cell_lnc=T_cell_lnc[,c(2,1)]
colnames(T_cell_lnc)=c("ENSEMBL","SYMBOL")
write.table(T_cell_lnc,"T_cell_lnc.txt",quote = F,sep = "\t",row.names = F)


select_gene=c("PSMB8-AS1","FGD5-AS1","MAGI2-AS3","MIR155HG","ADAMTS9-AS2",
              "LINC01094","MIR100HG","USP30-AS1")
cell_high[which(cell_high$V1%in%select_gene),]

plot.data=exp_data[,"MIR155HG"]
plot.data=cbind(plot.data,exp_data$group)
colnames(plot.data)=c("value","group")
plot.data=as.data.frame(plot.data)
plot.data$value=as.numeric(plot.data$value)
tmp=aggregate(value~group,plot.data,mean)
tmp=tmp[order(tmp$value),]
plot.data$group=factor(plot.data$group,levels = tmp$group)
library(vioplot)
vioplot(value~group, data = plot.data,
        col=c("#d67d79","#9fb9cd","#c5d832","#98afd9","#4d79b6","#ba86b4",
              "#c2dabd","#d57732","#b4c950","#97c7cb","#8cb0b0","#e2b719",
              "#dc736b","#7a9ed5","#c5c37b","#94a8c5","#db7aa1","#dfcc7d",
              "#e4b3b3"),las=2)

plot.data=exp_data[,"ADAMTS9-AS2"]
plot.data=cbind(plot.data,exp_data$group)
colnames(plot.data)=c("value","group")
plot.data=as.data.frame(plot.data)
plot.data$value=as.numeric(plot.data$value)
tmp=aggregate(value~group,plot.data,mean)
tmp=tmp[order(tmp$value),]
plot.data$group=factor(plot.data$group,levels = tmp$group)
library(vioplot)
vioplot(value~group, data = plot.data,
        col=c("#d67d79","#9fb9cd","#c5d832","#98afd9","#4d79b6","#ba86b4",
              "#c2dabd","#d57732","#b4c950","#97c7cb","#8cb0b0","#e2b719",
              "#dc736b","#7a9ed5","#c5c37b","#94a8c5","#db7aa1","#dfcc7d",
              "#e4b3b3"),las=2)

data=read.table("雷达图数据.txt",sep = "\t",header = T)
library(ggplot2)
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)
ggradar(data,
        grid.max = max(data[1,]),
        grid.min = 0,
        grid.mid = max(data[1,])/2,
        axis.label.size = 6,
        group.point.size = 3,
        group.colours = "#5b9bd5")
data[1,]/231

