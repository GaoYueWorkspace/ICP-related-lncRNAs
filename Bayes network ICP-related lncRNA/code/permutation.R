path=getwd()
library(readr)
#library(clusterProfiler)
#library(org.Hs.eg.db)
library(stringr)
library(edgeR)# 差异分析时使用
library(limma)# 可做差异分析
library(purrr)
library(ppcor)
library(progress)
library(data.table)
library(tidyr)
library(dplyr)
set.seed(12345)

#### 一、读入biaodapu，并规范化数据 ####

### 1、读入biaodapu
exp=read_tsv("E:\\曲课题\\数据\\data\\TCGA\\TCGA-GBM.htseq_fpkm.tsv\\TCGA-GBM.htseq_fpkm.tsv")
exp=as.data.frame(exp)
 
### 2、去除ENSG的版本号，保留15位
exp[,1]=substr(exp[,1],start = 1,stop = 15)
rownames(exp)=exp[,1]
exp=exp[,-1]
#exp[,"name"]=substr(rownames(exp),start = 1,stop = 15)
#rownames(exp)=exp[,"name"]
#exp=exp[,-1285]
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
remove(fc_zero,index,normal_exp,zero)
dir.create("./biaodapu")
setwd("./biaodapu")
save(exp,file="after_group_exp.RData")
setwd(path)


# N 30 ,C 81    #100    #25
#### 三、差异基因筛选 ####
load("./biaodapu/after_group_exp.RData")
index=which(substr(colnames(exp),start = 14,stop = 15)>10)
group<-c(rep("normal",length(index)),rep("disease",length(colnames(exp))-length(index)))
#group<-c(rep("normal",100),rep("disease",25))
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
tT = subset(tT, select="P.Value")
pb <- progress_bar$new(total=100)
system.time({
  # 创建一个空的结果数据框
  result_df <- data.frame(matrix(ncol = 0, nrow = length(rownames(exp))))
  for(i in 1:100){
    pb$tick()
  
    # 1、随机扰动样本标签
    all_sample <- colnames(exp)
    N_sample <- sample(all_sample, 100, replace = FALSE)
    T_sample  <- setdiff(all_sample, N_sample)
    final_sample_order <- c(N_sample, T_sample)
    exp1 <- exp[, final_sample_order]
    
    # 2、计算p值
    group1 <- c(rep("normal", length(N_sample)), rep("disease", length(T_sample)))
    targets1 <- data.frame(FileName = c(N_sample, T_sample), Target = group1)
    lev1 <- unique(targets1$Target)
    design1 <- model.matrix(~0 + factor(targets1$Target, levels = lev1))
    colnames(design1) <- lev1
    
    cont.wt1 <- makeContrasts(disease-normal, levels = design1)
    fit1 <- lmFit(exp1, design1)
    fit21 <- contrasts.fit(fit1, cont.wt1)
    fit21 <- eBayes(fit21)
    rtT <- topTable(fit21, adjust = "fdr", n = Inf)
    rtT <- subset(rtT, select = "P.Value")
    
    # 调整rtT的列名为"P.Value"后面的值随着i的变化而变化
    colname <- paste0("P.Value", i)
    colnames(rtT) <- colname
    
    # 将rtT按照tT的行名重新排序，并保持列名不变
    rtT_order <- rtT[rownames(tT), , drop = FALSE]
    
    # 合并两个数据框的列
    result_df <- cbind(result_df, rtT_order)
  }
})
result1 <- data.frame(matrix(ncol = 100, nrow = nrow(result_df)))
rownames(result1) <- rownames(result_df)
colnames(result1) <- colnames(result_df)
# 循环遍历第二个数据框的每一列
for (i in 1:ncol(result_df)) {
  # 将第一个数据框的第一列与第二个数据框的当前列进行比较
  # 并将结果存储在结果数据框相应的列中
  result1[, i] <- ifelse(tT[, 1] > result_df[, i], 1, 0)
}
count<- as.data.frame(rowSums(result1 == 1))
count100 = count %>% mutate(p = count$`rowSums(result1 == 1)`/100 ) 
count100[,"FDR"]=p.adjust(count100$p,method = "BH")
count_filter= count100 %>%filter(count100$FDR<0.01)#25990
expT=exp[rownames(count_filter),]
expT=expT[,-c(1:100)]
setwd("./biaodapu")
save(expT,file="disease_expT.RData") # 保存仅差异表达基因的biaodapu
setwd(path)
remove(cont.wt,design,fit,fit2,group,lev,targets,tT,index)
remove(cont.wt1,design1,fit1,fit21,group1,lev1,targets1)
#去除正常样本后的差异基因原先expT,后来expTT

#### 四、拆分三种基因的biaodapu ####

### 1、取出biaodapu中的基因
load("./biaodapu/disease_expT.RData")
name1=as.data.frame(rownames(expT))
colnames(name1)="ENSG"

### 2、读入注释文件
SYMBOL_ENSEMBL=read.table("E:\\曲课题\\数据\\注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("E:\\曲课题\\数据\\注释文件/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("E:\\曲课题\\数据\\注释文件/ICP.txt",sep="\t",header=T)
IMM=read.table("E:\\曲课题\\数据\\注释文件/imm_gene.txt",sep="\t",header=T)

### 3、划分为三种基因
lncRNA_list1=name1[which(name1$ENSG %in% lncRNA$ENSG),1]
ICP_list1=name1[which(name1$ENSG %in% ICP$ENSG),1]
IMM_list1=name1[which(name1$ENSG %in% IMM$ENSG),1]

### 4、划分出三种biaodapu
name1[which(name1$ENSG %in% lncRNA_list1),"class"]="lnc"
name1[which(name1$ENSG %in% ICP_list1),"class"]="ICP"
name1[which(name1$ENSG %in% IMM_list1),"class"]="IMM"
table(name1$class)#icp67,imm1503,lnc1260
lnc_exp1=expT[lncRNA_list1,]
ICP_exp1=expT[ICP_list1,]
IMM_exp1=expT[IMM_list1,]
# biaodapu转置一下，用于后续的map2，变成行是样本，列是基因
ICP_exp1=as.data.frame(t(ICP_exp1))
lnc_exp1=as.data.frame(t(lnc_exp1))
IMM_exp1=as.data.frame(t(IMM_exp1))
# 保存三种biaodapu
setwd("./biaodapu")
save(lnc_exp1,file="lnc_exp1.RData")
save(ICP_exp1,file="ICP_exp1.RData")
save(IMM_exp1,file="IMM_exp1.RData")
setwd(path)


#### 五、计算皮尔森偏相关 #### 
dir.create("./piersheng")
setwd("./piersheng")
### 1、计算lnc2icp
k=1
lnc2icp=matrix(ncol = 2,nrow = length(ICP_list1)*length(lncRNA_list1))
for(i in 1:length(lncRNA_list1)){
  for(j in 1:length(ICP_list1)){
    lnc2icp[k,]=c(lncRNA_list1[i],ICP_list1[j])
    k=k+1
  }
}
lnc2icp=as.data.frame(lnc2icp)
lnc2icp_data1=lnc_exp1[,lnc2icp[,1]]
lnc2icp_data2=ICP_exp1[,lnc2icp[,2]]
PC_lnc2icp=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
} 
system.time({lnc2icp_PCS=map2_df(lnc2icp_data1,lnc2icp_data2,PC_lnc2icp)})
#随机扰动
PC_lnc2icp1=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_pvalue=PC$p.value
  r<-c(p=PC_pvalue)
  return(r)
} 
pb <- progress_bar$new(total=100)
system.time({
  # 创建一个空的结果数据框
  lnc2icp_PCS100=data.frame(matrix(ncol =0,nrow = length(ICP_list1)*length(lncRNA_list1)))
  for(i in 1:100){
    pb$tick()
    
    # 1、随机扰动样本标签
    lnc2icp_data11 <- lnc2icp_data1[sample(nrow(lnc2icp_data1)), ]
    lnc2icp_data22 <- lnc2icp_data2[sample(nrow(lnc2icp_data2)), ]
    # 2、计算p值
    
    system.time({lnc2icp_PCS1=map2_df(lnc2icp_data11,lnc2icp_data22,PC_lnc2icp1)})
    
    # 合并两个数据框的列
    lnc2icp_PCS100=cbind(lnc2icp_PCS100,lnc2icp_PCS1)
  }
})
result_lnc2icp <- data.frame(matrix(ncol  = 100, nrow= nrow(lnc2icp_PCS100)))

# 循环遍历第二个数据框的每一列
for (i in 1:ncol(lnc2icp_PCS100)) {
  # 将第一个数据框的第一列与第二个数据框的当前列进行比较
  # 并将结果存储在结果数据框相应的列中
  result_lnc2icp[, i] <- ifelse(lnc2icp_PCS[, 2] > lnc2icp_PCS100[, i], 1, 0)
}
count_lnc2icp <- as.data.frame(rowSums(result_lnc2icp == 1))
count100_lnc2icp = count_lnc2icp %>% mutate(p = count_lnc2icp$`rowSums(result_lnc2icp == 1)`/100 ) 
pvalue=as.data.frame(count100_lnc2icp[,-1])
cor=as.data.frame(lnc2icp_PCS$R.cor)
lnc2icp_new=cbind(lnc2icp,cor,pvalue)
colnames(lnc2icp_new)=c("lnc","icp","R","p")
lnc2icp_filter=lnc2icp_new %>% filter(lnc2icp_new$p<0.05) #59267  #45542
lnc2icp_filter=lnc2icp_filter[which(abs(lnc2icp_filter$R)>0.5),]
write.table(lnc2icp_filter,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_lnc2icp,lnc2icp_PCS,lnc2icp_data1,lnc2icp_data2,lnc2icp,k,count_lnc2icp)
remove(PC_lnc2icp1,lnc2icp_PCS1,lnc2icp_data11,lnc2icp_data22,result_lnc2icp,lnc2icp_PCS100,pvalue,lnc2icp_new,count100_lnc2icp)

##
#lnc2icp1=cbind(lnc2icp,lnc2icp_PCS$R.cor)
#colnames(lnc2icp1)=c("lnc","icp","R")
#lnc2icp_filter <- merge(lnc2icp1, lnc2icp2, by = c("icp", "lnc"))
#lnc2icp_filter=lnc2icp_filter[which(abs(lnc2icp_filter$R)>0.5),]
#write.table(lnc2icp_filter,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")




### 2、计算icp2imm
k=1
icp2imm=matrix(ncol = 2,nrow = length(ICP_list1)*length(IMM_list1))
for(i in 1:length(ICP_list1)){
  for(j in 1:length(IMM_list1)){
    icp2imm[k,]=c(ICP_list1[i],IMM_list1[j])
    k=k+1
  }
}
icp2imm=as.data.frame(icp2imm)
icp2imm_data1=ICP_exp1[,icp2imm[,1]]
icp2imm_data2=IMM_exp1[,icp2imm[,2]]
PC_icp2imm=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
} 
PC_icp2imm1=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_pvalue=PC$p.value
  r<-c(p=PC_pvalue)
  return(r)
} 
system.time({icp2imm_PCS=map2_df(icp2imm_data1,icp2imm_data2,PC_icp2imm)})
pb <- progress_bar$new(total=100)
system.time({
  # 创建一个空的结果数据框
  icp2imm_PCS100=data.frame(matrix(ncol =0,nrow = length(ICP_list1)*length(IMM_list1)))
  for(i in 1:100){
    pb$tick()
    
    # 1、随机扰动样本标签
    icp2imm_data11 <- icp2imm_data1[sample(nrow(icp2imm_data1)), ]
    icp2imm_data22 <- icp2imm_data2[sample(nrow(icp2imm_data2)), ]
    # 2、计算p值
    
    system.time({icp2imm_PCS1=map2_df(icp2imm_data11,icp2imm_data22,PC_icp2imm1)})
    
    # 合并两个数据框的列
    icp2imm_PCS100=cbind(icp2imm_PCS100,icp2imm_PCS1)
  }
})
result_icp2imm <- data.frame(matrix(ncol  = 100, nrow= nrow(icp2imm_PCS100)))

# 循环遍历第二个数据框的每一列
for (i in 1:ncol(icp2imm_PCS100)) {
  # 将第一个数据框的第一列与第二个数据框的当前列进行比较
  # 并将结果存储在结果数据框相应的列中
  result_icp2imm[, i] <- ifelse(icp2imm_PCS[, 2] > icp2imm_PCS100[, i], 1, 0)
}
count_icp2imm <- as.data.frame(rowSums(result_icp2imm == 1))
count100_icp2imm = count_icp2imm %>% mutate(p = count_icp2imm$`rowSums(result_icp2imm == 1)`/100 ) 
pvalue=as.data.frame(count100_icp2imm[,-1])
cor=as.data.frame(icp2imm_PCS$R.cor)
icp2imm_new=cbind(icp2imm,cor,pvalue)
colnames(icp2imm_new)=c("icp","imm","R","p")
icp2imm_filter=icp2imm_new %>% filter(icp2imm_new$p<0.05) #59267  #45542
icp2imm_filter=icp2imm_filter[which(abs(icp2imm_filter$R)>0.5),]
write.table(icp2imm_filter,"icp2imm_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_icp2imm,icp2imm_PCS,icp2imm_data1,icp2imm_data2,icp2imm,k,count_icp2imm)
remove(PC_icp2imm1,icp2imm_PCS1,icp2imm_data11,icp2imm_data22,result_icp2imm,icp2imm_PCS100,pvalue,icp2imm_new,count100_icp2imm)
###
#icp2imm1=cbind(icp2imm,icp2imm_PCS$R.cor)
#colnames(icp2imm1)=c("icp","imm","R")
#icp2imm_filter <- merge(icp2imm1, icp2imm2, by = c("icp", "imm"))
#icp2imm_filter=icp2imm_filter[which(abs(icp2imm_filter$R)>0.5),]
#write.table(icp2imm_filter,"icp2imm_result.txt",quote = F,row.names = F,sep="\t")




### 3、计算lnc2imm
k=1
lnc2imm=matrix(ncol = 2,nrow = length(lncRNA_list1)*length(IMM_list1))
for(i in 1:length(lncRNA_list1)){
  for(j in 1:length(IMM_list1)){
    lnc2imm[k,]=c(lncRNA_list1[i],IMM_list1[j])
    k=k+1
  }
}
lnc2imm=as.data.frame(lnc2imm)
lnc2imm_data1=lnc_exp1[,lnc2imm[,1]]
lnc2imm_data2=IMM_exp1[,lnc2imm[,2]]

# 计算皮尔逊相关
PC_lnc2imm=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_estimate=PC$estimate
  PC_pvalue=PC$p.value
  r<-c(R=PC_estimate,p=PC_pvalue)
  return(r)
}
PC_lnc2imm1=function(x,y){
  PC=cor.test(as.numeric(x),as.numeric(y),method = "pearson")
  PC_pvalue=PC$p.value
  r<-c(p=PC_pvalue)
  return(r)
} 
system.time({
  lnc2imm_PCS=map2_df(lnc2imm_data1,lnc2imm_data2,PC_lnc2imm)
})
pb <- progress_bar$new(total=100)
system.time({
  # 创建一个空的结果数据框
  lnc2imm_PCS100=data.frame(matrix(ncol =0,nrow = length(lncRNA_list1)*length(IMM_list1)))
  for(i in 1:100){
    pb$tick()
    
    # 1、随机扰动样本标签
    lnc2imm_data11 <- lnc2imm_data1[sample(nrow(lnc2imm_data1)), ]
    lnc2imm_data22 <- lnc2imm_data2[sample(nrow(lnc2imm_data2)), ]
    # 2、计算p值
    
    system.time({lnc2imm_PCS1=map2_df(lnc2imm_data11,lnc2imm_data22,PC_lnc2imm1)})
    
    # 合并两个数据框的列
    lnc2imm_PCS100=cbind(lnc2imm_PCS100,lnc2imm_PCS1)
  }
})
 result_lnc2imm <- data.frame(matrix(ncol  = 100, nrow= nrow(lnc2imm_PCS100)))

# 循环遍历第二个数据框的每一列
for (i in 1:ncol(lnc2imm_PCS100)) {
  # 将第一个数据框的第一列与第二个数据框的当前列进行比较
  # 并将结果存储在结果数据框相应的列中
  result_lnc2imm[, i] <- ifelse(lnc2imm_PCS[, 2] > lnc2imm_PCS100[, i], 1, 0)
}
count_lnc2imm <- as.data.frame(rowSums(result_lnc2imm == 1))
count100_lnc2imm = count_lnc2imm %>% mutate(p = count_lnc2imm$`rowSums(result_lnc2imm == 1)`/100 ) 
pvalue=as.data.frame(count100_lnc2imm[,-1])
cor=as.data.frame(lnc2imm_PCS$R.cor)
lnc2imm_new=cbind(lnc2imm,cor,pvalue)
colnames(lnc2imm_new)=c("lnc","imm","R","p")
lnc2imm_filter=lnc2imm_new %>% filter(lnc2imm_new$p<0.05) #59267  #45542
lnc2imm_filter=lnc2imm_filter[which(abs(lnc2imm_filter$R)>0.5),]
write.table(lnc2imm_filter,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
remove(PC_lnc2imm,lnc2imm_PCS,lnc2imm_data1,lnc2imm_data2,lnc2imm,k,i,j)
remove(PC_lnc2imm1,lnc2imm_PCS1,lnc2imm_data11,lnc2imm_data22,result_lnc2imm,lnc2imm_PCS100,pvalue,lnc2imm_new,count100_lnc2imm)
# lnc2icp2=read.table("./piersheng/lnc2icp_result.txt",header=T)
# lnc2imm2=read.table("./piersheng/lnc2imm_result.txt",header=T)
# icp2imm2=read.table("./piersheng/icp2imm_result.txt",header=T)

setwd(path)


#### 六、以皮尔森筛选后的基因互作构建gongbiaoda_net
lnc2icp2=read.table("./piersheng/lnc2icp_result.txt",header=T)
lnc2imm2=read.table("./piersheng/lnc2imm_result.txt",header=T)
icp2imm2=read.table("./piersheng/icp2imm_result.txt",header=T)

lnc2icp=lnc2icp2[,-c(4)]
icp2imm=icp2imm2[,-c(4)]
lnc2imm=lnc2imm2[,-c(4)]

net=data.frame()
net=rbind(net,lnc2icp)
colnames(net)=colnames(icp2imm)
net=rbind(net,icp2imm)
colnames(net)=colnames(lnc2imm)
net=rbind(net,lnc2imm) # net为gongbiaoda_net
colnames(net)=c("gene1","gene2","cor")
dir.create("./gongbiaoda_net")
setwd("./gongbiaoda_net")
write.table(net,"net.txt",row.names = F,quote = F,sep = "\t")
setwd(path)
remove(lnc2icp,icp2imm,lnc2imm,net)

#### 七、python进行page rank算法，生成result文件 ####

