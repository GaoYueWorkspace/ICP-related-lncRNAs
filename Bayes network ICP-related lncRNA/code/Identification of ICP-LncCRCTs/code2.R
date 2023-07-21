path=getwd()
#library(readr)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(stringr)
library(edgeR)# 差异分析时使用
library(limma)# 可做差异分析
#library(purrr)
library(ppcor)
library(progress)
library(data.table)
library(tidyr)

#### 文件准备 ###
SYMBOL_ENSEMBL=read.table("./注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("./注释文件/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("./注释文件/ICP.txt",sep="\t",header=T)
IMM=read.table("./注释文件/imm_gene.txt",sep="\t",header=T)
lnc2icp2=read.table("./皮尔森相关结果/lnc2icp_result.txt",header=T)
lnc2imm2=read.table("./皮尔森相关结果/lnc2imm_result.txt",header=T)
icp2imm2=read.table("./皮尔森相关结果/icp2imm_result.txt",header=T)
load("./表达谱/ICP_exp.RData")
load("./表达谱/lnc_exp.RData")
load("./表达谱/IMM_exp.RData")


#### 八、选择影响前30%的基因作为与ICP基因联系密切的基因 ####
result=read.table("./page-rank结果/result.txt",header = T)
selection=result[which(!(result$ENSG %in% ICP$ENSG)),]
selection=selection[order(-selection$pagerank),]
gene_300=selection[1:200,]

gene_300[which(gene_300$ENSG %in% lncRNA$ENSG),"class"]="lnc"
gene_300[which(gene_300$ENSG %in% IMM$ENSG),"class"]="IMM"

table(gene_300$class)

imm_gene=gene_300[which(gene_300$class=="IMM"),]
lnc_gene=gene_300[which(gene_300$class=="lnc"),]
ICP_gene=result[which(result$ENSG %in% ICP$ENSG),]
remove(result,selection)

# 更新互作对
lnc2icp2=lnc2icp2[which(lnc2icp2$lnc %in% lnc_gene$ENSG),]
lnc2icp2=lnc2icp2[which(lnc2icp2$icp %in% ICP_gene$ENSG),]

lnc2imm2=lnc2imm2[which(lnc2imm2$lnc %in% lnc_gene$ENSG),]
lnc2imm2=lnc2imm2[which(lnc2imm2$imm %in% imm_gene$ENSG),]

icp2imm2=icp2imm2[which(icp2imm2$icp %in% ICP_gene$ENSG),]
icp2imm2=icp2imm2[which(icp2imm2$imm %in% imm_gene$ENSG),]

#### 九、将肿瘤纯度得分作为共变量，计算lnc与ICP的偏相关系数 ####

gene_pair=lnc2icp2[,c(1,2)]


### 1、读入肿瘤纯度得分文件
purity_score=read.table("./计算偏相关系数/inte_BLCA.txt",header = T,sep = "\t")
purity_score=purity_score[,c(1,3)]

### 2、对肿瘤纯度得分文件进行预处理

# 过滤掉正常样本
purity_score=purity_score[-which(substr(purity_score$V1,start = 14,stop = 15)>10),]

# 样本名与表达谱样本名格式保持一致
barcodes=purity_score$V1
barcodes=as.data.frame(strsplit(barcodes,"-"))
barcodes=apply(barcodes,2,function(x) return(paste(x[1],x[2],x[3],x[4],sep = ".")))
purity_score[,1]=barcodes
remove(barcodes)

### 3、使样本顺序保持一致
lnc_exp=as.data.frame(t(lnc_exp))
lnc_exp=lnc_exp[,purity_score$V1]

ICP_exp=as.data.frame(t(ICP_exp))
ICP_exp=ICP_exp[,purity_score$V1]

IMM_exp=as.data.frame(t(IMM_exp))
IMM_exp=IMM_exp[,purity_score$V1]


### 4、计算PCC系数及P值

PCOR_LI=function(x){
  PCC=pcor.test(as.numeric(lnc_exp[x[1],]),as.numeric(ICP_exp[x[2],]),as.numeric(purity_score$TumorPurity),method = "pearson")
  r=c(estimate=PCC$estimate,p=PCC$p.value)
  return(r)
}
PCC_LI=apply(gene_pair,1,PCOR_LI)
PCC_LI=as.data.frame(t(PCC_LI))
PCC_LI=cbind(gene_pair,PCC_LI)
PCC_LI[,"RS"]=(-log10(PCC_LI$p))*(sign(PCC_LI$estimate))
remove(PCOR_LI)
index=which(PCC_LI$estimate>0.25 & PCC_LI$p<0.05)
PCC_LI=PCC_LI[index,]
lnc2icp2=lnc2icp2[index,] # 更新lnc2icp关系对
remove(index,gene_pair)

#### 十、进行GSEA富集分析，对lnc-icp对进行打分 ####

### 1、读入免疫相关基因的相关通路文件
gmt<-read.table("./注释文件/imm_genelist.txt",sep = "\t",header = T)
colnames(gmt)<-c("ont","gene")# ont是生物学通路
length(table(gmt$ont))

### 2、读入IMM基因的表达谱文件
load("./表达谱/after_group_exp.RData")
gmtlast=data.frame()

### 3、将表达谱中基因的Ensembl ID转换为Symbol ID
name=as.data.frame(rownames(exp))
colnames(name)="name"
SYMBOL_ENSEMBL=read.table("./注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
exp=exp[which(name$name %in% SYMBOL_ENSEMBL$ENSG),] # 筛选出有注释的基因
name=as.data.frame(rownames(exp))
colnames(name)="ensg"

for (i in 1:nrow(name)) {
  name[i,"symbol"]=SYMBOL_ENSEMBL[(which(SYMBOL_ENSEMBL$ENSG %in% name[i,1])),2]
}

freq=as.data.frame(table(name$symbol))
freq=freq[which(freq$Freq==1),]
name=name[which(name$symbol %in% freq$Var1),]
exp=exp[name$ensg,]
rownames(exp)=name$symbol
remove(freq)

### 4、筛选表达谱中含有的IMM基因通路
for(i in 1:nrow(gmt)){
  if(gmt$gene[i] %in% name$symbol){
    gmtlast<-rbind(gmtlast,gmt[i,])
  } 
}
GO_term=gmtlast[which(gmtlast$ont=="GO"),]
GO_term=GO_term[which(GO_term$gene %in% gene_300$SYMBOL),]
gmtlast=gmtlast[which(!(gmtlast$ont=="GO")),]
gmtlast=rbind(gmtlast,GO_term)
gmt=gmtlast # 表达谱中存在于GeneList的基因
remove(gmtlast,GO_term)


### 5、差异基因筛选
index=which(substr(colnames(exp),start = 14,stop = 15)>10)
group<-c(rep("normal",length(index)),rep("disease",length(colnames(exp))-length(index)))
targets<-cbind(colnames(exp),group)#第一列GSM号，第二列为样本分组
targets<-as.data.frame(targets)
colnames(targets)=c("FileName","Target")
lev<-unique(targets$Target)
design <- model.matrix(~0+factor(targets$Target, levels=lev)) #样本矩阵
colnames(design) <- lev #更改列名为levels名
###两两之间的比较，求差异基因使用topTable函数
cont.wt <- makeContrasts(disease-normal,
                         levels=design) 
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2,adjust="fdr",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")


### 6、制作基因集
GSEAdata<-data.frame(gene=rownames(tT),logFC=tT$logFC)
geneList<-GSEAdata$logFC#第二列可以是folodchange，也可以是logFC
names(geneList)=GSEAdata$gene 
geneList=sort(geneList,decreasing = T) #log2FC从高到低排序

### 7、GSEA富集
gene_enrichment_w=data.frame()

for(n in 1:length(table(gmt$ont))){
  gmt1<-gmt[which(gmt$ont==unique(gmt$ont)[n]),]# 取出富集到第n个通路上的基因集
  gmts<-gmt1
  for(m in 1:nrow(gmt1)){
    gmt2<-data.frame(ont=m,gene=gmt1[-m,2])
    gmts<-rbind(gmts,gmt2)
  }
  set.seed(1234)
  GSEAresult2<-GSEA(geneList,TERM2GENE = gmts,pvalueCutoff=1)
  GSEAresult2<-as.data.frame(GSEAresult2)
  if(nrow(GSEAresult2)==0){
    next
  }
  ES<-GSEAresult2$enrichmentScore[which(GSEAresult2$ID==unique(gmt$ont)[n])]
  if(length(ES)==0){
    next
  }
  GSEAresult2<-GSEAresult2[-which(GSEAresult2$ID==unique(gmt$ont)[n]),]
  GSEAresult2$ID<-as.numeric(GSEAresult2$ID)
  GSEAresult2<-GSEAresult2[order(GSEAresult2$ID),]
  allxishu<-ES/GSEAresult2$enrichmentScore
  w=data.frame("gene"=gmt1$gene,"w"=allxishu)
  gene_enrichment_w=rbind(gene_enrichment_w,w)
}

colnames(gene_enrichment_w)=c("SYMBOL","w")
for(i in 1:length(imm_gene$ENSG)){
  imm_gene[i,"SYMBOL"]=SYMBOL_ENSEMBL[which(SYMBOL_ENSEMBL$ENSG %in% imm_gene[i,1]),2]
}
imm_p_w=merge(imm_gene,gene_enrichment_w,by="SYMBOL")# 其中包括了一个基因可能对应多个通路的情况
keys = colnames(imm_p_w)[!grepl('w',colnames(imm_p_w))]
imm_p_w=as.data.table(imm_p_w)
imm_p_w=imm_p_w[,list(w= mean(w)),keys]# 有多个的取均值
imm_p_w=as.data.frame(imm_p_w)
colnames(imm_p_w)[c(1,2,5)]=c("SYMBOL","ENSEMBL","w")
dir.create("./富集分析")
setwd("./富集分析")
write.table(imm_p_w,"imm_pagerank_and_enrichment.txt",quote = F,sep = "\t",row.names = F)
setwd(path)
remove(allxishu,cont.wt,design,ES,fit,fit2,gene_enrichment_w,geneList,gmt1,gmt2,
       gmts,group,GSEAdata,GSEAresult2,keys,lev,m,n,i,tT,w,targets,index,gmt)


## 8、过滤imm基因
lnc2imm2=lnc2imm2[which(lnc2imm2$imm %in% imm_p_w$ENSEMBL),]
icp2imm2=icp2imm2[which(icp2imm2$imm %in% imm_p_w$ENSEMBL),]
lnc2icp2=lnc2icp2[which(lnc2icp2$lnc %in% lnc2imm2$lnc),]

dir.create("./GSEA后")
setwd("./GSEA后")
write.table(lnc2icp2,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")
write.table(icp2imm2,"icp2imm_result.txt",quote = F,row.names = F,sep="\t")
write.table(lnc2imm2,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
write.table(purity_score,"purity_score.txt",sep = "\t",quote = F)
save(lnc_exp,file="lnc_exp.RData")
save(ICP_exp,file="ICP_exp.RData")
save(IMM_exp,file="IMM_exp.RData")
setwd(path)


