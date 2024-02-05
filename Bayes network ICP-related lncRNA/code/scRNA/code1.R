library(Seurat)
library(dplyr)
library(edgeR)
library(progress)
library(purrr)

path=getwd()
pbmc=readRDS("pbmc.rds")
pbmc.markers=read.table("markers.txt",sep = "\t")
# 提取counts矩阵
count <- GetAssayData(pbmc, slot = "counts")
# 转换为数据框
count <- as.data.frame(count)
count=count[unique(pbmc.markers$gene),]
name=as.data.frame(rownames(count))
colnames(name)="SYMBOL"

# 读入注释文件
SYMBOL_ENSEMBL=read.table("./注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("./注释文件/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("./注释文件/ICP.txt",sep="\t",header=T)
IMM=read.table("./注释文件/imm_gene.txt",sep="\t",header=T)

# 划分为三种基因
lncRNA_list=name[which(name$SYMBOL %in% lncRNA$SYMBOL),1]
ICP_list=name[which(name$SYMBOL %in% ICP$SYMBOL),1]
IMM_list=name[which(name$SYMBOL %in% IMM$SYMBOL),1]

SYMBOL_ENSEMBL=SYMBOL_ENSEMBL[which(SYMBOL_ENSEMBL$SYMBOL%in%rownames(count)),]
count=count[SYMBOL_ENSEMBL$SYMBOL,]
SYMBOL_ENSEMBL=SYMBOL_ENSEMBL[!duplicated(SYMBOL_ENSEMBL$SYMBOL),]
count=count[SYMBOL_ENSEMBL$SYMBOL,]
max(apply(count,1,max))
min(apply(count,1,min))
exp=count
dir.create("./表达谱")
save(exp,file="./表达谱/exp.RData") 

max(apply(exp,1,max))
min(apply(exp,1,min))

rownames(exp)=SYMBOL_ENSEMBL$ENSG

# 表达谱按细胞类型拆分
name=as.data.frame(rownames(exp))
colnames(name)="ENSG"
lncRNA_list=name[which(name$ENSG %in% lncRNA$ENSG),1]
ICP_list=name[which(name$ENSG %in% ICP$ENSG),1]
IMM_list=name[which(name$ENSG %in% IMM$ENSG),1]
cell_list=readRDS("./细胞列表/cell_list.rds")
cell_name=names(cell_list)
for(i in 1:length(cell_name)){
  dir.create(paste("./",cell_name[i],sep = ""))
  setwd(paste("./",cell_name[i],sep = ""))
  cell_exp=exp[,cell_list[[i]]]
  # 拆分成三种表达谱
  lnc_exp=cell_exp[lncRNA_list,]
  ICP_exp=cell_exp[ICP_list,]
  IMM_exp=cell_exp[IMM_list,]
  # 表达谱转置一下，用于后续的map2，变成行是样本，列是基因
  ICP_exp=as.data.frame(t(ICP_exp))
  lnc_exp=as.data.frame(t(lnc_exp))
  IMM_exp=as.data.frame(t(IMM_exp))
  # 保存三种表达谱
  dir.create("./表达谱")
  setwd("./表达谱")
  save(lnc_exp,file="lnc_exp.RData")
  save(ICP_exp,file="ICP_exp.RData")
  save(IMM_exp,file="IMM_exp.RData")
  setwd(path)
}


for(i in 1:length(cell_name)){
  #### 计算皮尔森偏相关 #### 
  setwd(paste("./",cell_name[i],sep = ""))
  load("表达谱/ICP_exp.RData")
  load("表达谱/IMM_exp.RData")
  load("表达谱/lnc_exp.RData")
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
  # lnc2icp2=lnc2icp2[which(abs(lnc2icp2$R)>0.3),]
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
  # icp2imm2=icp2imm2[which(abs(icp2imm2$R)>0.3),]
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
  # lnc2imm2=lnc2imm2[which(abs(lnc2imm2$R)>0.3),]
  write.table(lnc2imm2,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
  remove(PC_lnc2imm,lnc2imm_PCS,lnc2imm_data1,lnc2imm_data2,lnc2imm,k,i,j)
  setwd("../")
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
}



