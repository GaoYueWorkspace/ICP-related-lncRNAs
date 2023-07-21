path=getwd()
library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(edgeR)# 差异分析时使用
library(limma)# 可做差异分析
library(purrr)
library(ppcor)
library(progress)
library(data.table)
library(tidyr)


#### 一、读入表达谱，并规范化数据 ####

### 1、读入表达谱
exp=read_tsv("TCGA-BLCA.htseq_fpkm.tsv")
exp=as.data.frame(exp)

### 2、去除ENSG的版本号，保留15位
exp[,1]=substr(exp[,1],start = 1,stop = 15)
rownames(exp)=exp[,1]
exp=exp[,-1]

### 3、去除超过70% 0值的基因
## 计算每行0值的个数
fc_zero<-function(f){
  zero_num<-sum(f==0)
  return(zero_num)
}
zero<-as.data.frame(apply(exp,1,fc_zero))
colnames(zero)="num"
zero[,"ratio"]<-zero$num/length(colnames(exp))
exp=exp[which(zero$ratio<0.7),]


#### 二、区分正常样本，并完成排序，正常在前 ####
index=which(substr(colnames(exp),start = 14,stop = 15)>10)
normal_exp=exp[,index]
exp=exp[,-index]
# 排序样本
nt=data.frame(normal_exp,exp)
exp=nt
remove(nt,fc_zero,index,normal_exp,zero)
dir.create("./表达谱")
setwd("./表达谱")
save(exp,file="after_group_exp.RData")
setwd(path)


#### 三、差异基因筛选 ####
load(".\\表达谱\\after_group_exp.RData")
index=which(substr(colnames(exp),start = 14,stop = 15)>10)
group<-c(rep("normal",length(index)),rep("disease",length(colnames(exp))-length(index)))
targets<-cbind(colnames(exp),group)#第一列GSM号，第二列为样本分组
targets<-as.data.frame(targets)
colnames(targets)=c("FileName","Target")
lev<-unique(targets$Target)
design <- model.matrix(~0+factor(targets$Target, levels=lev)) #样本矩阵
colnames(design) <- lev #更改列名为levels名
# 两两之间的比较，求差异基因使用topTable函数
cont.wt <- makeContrasts(disease-normal,
                         levels=design) 
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2,adjust="fdr",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
tT=tT[which(tT$FDR<0.01),] # 14960个差异表达的基因
exp=exp[rownames(tT),]
exp=exp[,-c(1:length(index))]
setwd("./表达谱")
save(exp,file="disease_exp.RData") # 保存仅差异表达基因的表达谱
setwd(path)
remove(cont.wt,design,fit,fit2,group,lev,targets,tT,index)


#### 四、拆分三种基因的表达谱 ####

### 1、取出表达谱中的基因
load("./表达谱/disease_exp.RData")
name=as.data.frame(rownames(exp))
colnames(name)="ENSG"

### 2、读入注释文件
SYMBOL_ENSEMBL=read.table("./注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("./注释文件/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("./注释文件/ICP.txt",sep="\t",header=T)
IMM=read.table("./注释文件/imm_gene.txt",sep="\t",header=T)

### 3、划分为三种基因
lncRNA_list=name[which(name$ENSG %in% lncRNA$ENSG),1]
ICP_list=name[which(name$ENSG %in% ICP$ENSG),1]
IMM_list=name[which(name$ENSG %in% IMM$ENSG),1]

### 4、划分出三种表达谱
name[which(name$ENSG %in% lncRNA_list),"class"]="lnc"
name[which(name$ENSG %in% ICP_list),"class"]="ICP"
name[which(name$ENSG %in% IMM_list),"class"]="IMM"
table(name$class)
lnc_exp=exp[lncRNA_list,]
ICP_exp=exp[ICP_list,]
IMM_exp=exp[IMM_list,]
# 表达谱转置一下，用于后续的map2，变成行是样本，列是基因
ICP_exp=as.data.frame(t(ICP_exp))
lnc_exp=as.data.frame(t(lnc_exp))
IMM_exp=as.data.frame(t(IMM_exp))
# 保存三种表达谱
setwd("./表达谱")
save(lnc_exp,file="lnc_exp.RData")
save(ICP_exp,file="ICP_exp.RData")
save(IMM_exp,file="IMM_exp.RData")
setwd(path)


#### 五、计算皮尔森偏相关 #### 
dir.create("./皮尔森相关结果")
setwd("./皮尔森相关结果")
### 1、计算lnc2icp
k=1
lnc2icp=matrix(ncol = 2,nrow = length(ICP_list)*length(lncRNA_list))
for(i in 1:length(lncRNA_list)){
  for(j in 1:length(ICP_list)){
    lnc2icp[k,]=c(lncRNA_list[i],ICP_list[j])
    k=k+1
  }
}
lnc2icp=as.data.frame(lnc2icp)
lnc2icp_data1=lnc_exp[,lnc2icp[,1]]
lnc2icp_data2=ICP_exp[,lnc2icp[,2]]

# 计算皮尔逊偏相关系数
PC_lnc2icp=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
} 
system.time({lnc2icp_PCS=map2_df(lnc2icp_data1,lnc2icp_data2,PC_lnc2icp)})
lnc2icp=cbind(lnc2icp,lnc2icp_PCS)
colnames(lnc2icp)=c("lnc","icp","R","p")
lnc2icp2=lnc2icp
lnc2icp2=lnc2icp2[which(lnc2icp2$p<0.05),]
lnc2icp2=lnc2icp2[which(abs(lnc2icp2$R)>0.5),]
write.table(lnc2icp2,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_lnc2icp,lnc2icp_PCS,lnc2icp_data1,lnc2icp_data2,lnc2icp,k)


### 2、计算icp2imm
k=1
icp2imm=matrix(ncol = 2,nrow = length(ICP_list)*length(IMM_list))
for(i in 1:length(ICP_list)){
  for(j in 1:length(IMM_list)){
    icp2imm[k,]=c(ICP_list[i],IMM_list[j])
    k=k+1
  }
}
icp2imm=as.data.frame(icp2imm)
icp2imm_data1=ICP_exp[,icp2imm[,1]]
icp2imm_data2=IMM_exp[,icp2imm[,2]]

# 计算皮尔逊相关系数
PC_icp2imm=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
} 

system.time({
  icp2imm_PCS=map2_df(icp2imm_data1,icp2imm_data2,PC_icp2imm)
})
icp2imm=cbind(icp2imm,icp2imm_PCS)
colnames(icp2imm)=c("icp","imm","R","p")
icp2imm2=icp2imm
icp2imm2=icp2imm2[which(icp2imm2$p<0.05),]
icp2imm2=icp2imm2[which(abs(icp2imm2$R)>0.5),]
write.table(icp2imm2,"icp2imm_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_icp2imm,icp2imm_PCS,icp2imm_data1,icp2imm_data2,icp2imm,k)


### 3、计算lnc2imm
k=1
lnc2imm=matrix(ncol = 2,nrow = length(lncRNA_list)*length(IMM_list))
for(i in 1:length(lncRNA_list)){
  for(j in 1:length(IMM_list)){
    lnc2imm[k,]=c(lncRNA_list[i],IMM_list[j])
    k=k+1
  }
}
lnc2imm=as.data.frame(lnc2imm)
lnc2imm_data1=lnc_exp[,lnc2imm[,1]]
lnc2imm_data2=IMM_exp[,lnc2imm[,2]]

# 计算皮尔逊相关
PC_lnc2imm=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
} 
system.time({
  lnc2imm_PCS=map2_df(lnc2imm_data1,lnc2imm_data2,PC_lnc2imm)
})
lnc2imm=cbind(lnc2imm,lnc2imm_PCS)
colnames(lnc2imm)=c("lnc","imm","R","p")
lnc2imm2=lnc2imm
lnc2imm2=lnc2imm2[which(lnc2imm2$p<0.05),]
lnc2imm2=lnc2imm2[which(abs(lnc2imm2$R)>0.5),]
write.table(lnc2imm2,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_lnc2imm,lnc2imm_PCS,lnc2imm_data1,lnc2imm_data2,lnc2imm,k,i,j)

# lnc2icp2=read.table("./皮尔森相关结果/lnc2icp_result.txt",header=T)
# lnc2imm2=read.table("./皮尔森相关结果/lnc2imm_result.txt",header=T)
# icp2imm2=read.table("./皮尔森相关结果/icp2imm_result.txt",header=T)

setwd(path)


#### 六、以皮尔森筛选后的基因互作构建共表达网络
lnc2icp=lnc2icp2[,-c(4)]
icp2imm=icp2imm2[,-c(4)]
lnc2imm=lnc2imm2[,-c(4)]

net=data.frame()
net=rbind(net,lnc2icp)
colnames(net)=colnames(icp2imm)
net=rbind(net,icp2imm)
colnames(net)=colnames(lnc2imm)
net=rbind(net,lnc2imm) # net为共表达网络
colnames(net)=c("gene1","gene2","cor")
dir.create("./共表达网络")
setwd("./共表达网络")
write.table(net,"net.txt",row.names = F,quote = F,sep = "\t")
setwd(path)
remove(lnc2icp,icp2imm,lnc2imm,net)

#### 七、python进行page rank算法，生成result文件 ####

