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
dir.create("./互信息结果")
setwd("./互信息结果")

library(bnlearn)
library(purrr)
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

# 计算互信息系数
MI_lnc2icp=function(x,y){
  MIC=ci.test(as.numeric(x),as.numeric(y),test = "mi-g")
  MI_estimate=MIC$statistic
  MI_pvalue=MIC$p.value
  r<-c(MI=MI_estimate,p=MI_pvalue)
  return(r)
} 
system.time({lnc2icp_MIS=map2_df(lnc2icp_data1,lnc2icp_data2,MI_lnc2icp)})
lnc2icp=cbind(lnc2icp,lnc2icp_MIS)
colnames(lnc2icp)=c("lnc","icp","MI","p")
lnc2icp2=lnc2icp
lnc2icp2=lnc2icp2[which(lnc2icp2$p<0.05),]
# lnc2icp2=lnc2icp2[which(abs(lnc2icp2$R)>0.45),]
write.table(lnc2icp2,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")
remove(MI_lnc2icp,lnc2icp_MIS,lnc2icp_data1,lnc2icp_data2,lnc2icp,k)

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

# 计算互信息系数
MI_icp2imm=function(x,y){
  MIC=ci.test(as.numeric(x),as.numeric(y),test = "mi-g")
  MI_estimate=MIC$statistic
  MI_pvalue=MIC$p.value
  r<-c(MI=MI_estimate,p=MI_pvalue)
  return(r)
} 
system.time({icp2imm_MIS=map2_df(icp2imm_data1,icp2imm_data2,MI_icp2imm)})
icp2imm=cbind(icp2imm,icp2imm_PCS)
colnames(icp2imm)=c("icp","imm","MI","p")
icp2imm2=icp2imm
icp2imm2=icp2imm2[which(icp2imm2$p<0.05),]
# icp2imm2=icp2imm2[which(abs(icp2imm2$R)>0.45),]
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

# 计算互信息系数
MI_lnc2imm=function(x,y){
  MIC=ci.test(as.numeric(x),as.numeric(y),test = "mi-g")
  MI_estimate=MIC$statistic
  MI_pvalue=MIC$p.value
  r<-c(MI=MI_estimate,p=MI_pvalue)
  return(r)
} 
system.time({lnc2imm_MIS=map2_df(lnc2imm_data1,lnc2imm_data2,MI_lnc2imm)})
lnc2imm=cbind(lnc2imm,lnc2imm_PCS)
colnames(lnc2imm)=c("lnc","imm","R","p")
lnc2imm2=lnc2imm
lnc2imm2=lnc2imm2[which(lnc2imm2$p<0.05),]
# lnc2imm2=lnc2imm2[which(abs(lnc2imm2$R)>0.45),]
write.table(lnc2imm2,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_lnc2imm,lnc2imm_PCS,lnc2imm_data1,lnc2imm_data2,lnc2imm,k,i,j)

