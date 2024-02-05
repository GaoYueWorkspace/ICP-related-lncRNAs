path=getwd()
library(clusterProfiler)
library(org.Hs.eg.db)
library(edgeR)# 差异分析时使用
library(limma)# 可做差异分析
library(ppcor)
library(progress)
library(data.table)
library(tidyr)
cell_list=readRDS("./细胞列表/cell_list.rds")
cell_name=names(cell_list)

#### 文件准备 ###
SYMBOL_ENSEMBL=read.table("./注释文件/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("./注释文件/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("./注释文件/ICP.txt",sep="\t",header=T)
IMM=read.table("./注释文件/imm_gene.txt",sep="\t",header=T)
load("./表达谱/exp.RData")
for(c in 1:length(cell_name)){
  setwd(paste("./",cell_name[c],sep = ""))
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
  setwd(path)
  #### 十、进行GSEA富集分析，对lnc-icp对进行打分 ####
  
  ### 1、读入免疫相关基因的相关通路文件
  gmt<-read.table("./注释文件/imm_genelist.txt",sep = "\t",header = T)
  colnames(gmt)<-c("ont","gene")# ont是生物学通路
  length(table(gmt$ont))
  ### 2、读入IMM基因的表达谱文件
  gmtlast=data.frame()
  name=as.data.frame(rownames(exp))
  colnames(name)="name"
  ### 4、筛选表达谱中含有的IMM基因通路
  for(i in 1:nrow(gmt)){
    if(gmt$gene[i] %in% name$name){
      gmtlast<-rbind(gmtlast,gmt[i,])
    } 
  }
  GO_term=gmtlast[which(gmtlast$ont=="GO"),]
  tmp=SYMBOL_ENSEMBL[which(SYMBOL_ENSEMBL$ENSG%in%gene_300$ENSG),]
  GO_term=GO_term[which(GO_term$gene %in% tmp$SYMBOL),]
  gmtlast=gmtlast[which(!(gmtlast$ont=="GO")),]
  gmtlast=rbind(gmtlast,GO_term)
  gmt=gmtlast # 表达谱中存在于GeneList的基因
  remove(gmtlast,GO_term)
  
  tT=read.table("markers.txt",sep = "\t")
  
  ### 6、制作基因集
  GSEAdata<-data.frame(gene=rownames(tT),logFC=tT$avg_log2FC)
  geneList<-GSEAdata$logFC#第二列可以是folodchange，也可以是logFC
  names(geneList)=GSEAdata$gene 
  geneList=sort(geneList,decreasing = T) #log2FC从高到低排序
  
  ### 7、GSEA富集
  gene_enrichment_w=data.frame()
  
  for(n in 1:length(table(gmt$ont))){
    gmt1<-gmt[which(gmt$ont==unique(gmt$ont)[n]),]# 取出富集到第n个通路上的基因集
    gmts<-gmt1
    if(nrow(gmt1)==1){
      gmts=rbind(gmts,gmt1)
    }else{
      for(m in 1:nrow(gmt1)){
        gmt2<-data.frame(ont=m,gene=gmt1[-m,2])
        gmts<-rbind(gmts,gmt2)
      }
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
    if(nrow(GSEAresult2)==0){
      next
    }
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
  # 提取通路名称
  pathway_names <- gmt[,c(1,2)]
  colnames(pathway_names) <- c("ont", "SYMBOL")
  # 合并通路名称和富集得分
  gene_enrichment_w <- merge(gene_enrichment_w, pathway_names, by = "SYMBOL")
  imm_p_w=merge(imm_gene,gene_enrichment_w,by="SYMBOL")# 其中包括了一个基因可能对应多个通路的情况
  keys = colnames(imm_p_w)[!grepl('w',colnames(imm_p_w))]
  imm_p_w=as.data.table(imm_p_w)
  imm_p_w=imm_p_w[,list(w= mean(w)),keys]# 有多个的取均值
  imm_p_w=as.data.frame(imm_p_w)
  colnames(imm_p_w)[c(1,2,5,6)]=c("SYMBOL","ENSEMBL","ont","w")
  setwd(paste("./",cell_name[c],sep = ""))
  dir.create("./富集分析")
  setwd("./富集分析")
  write.table(imm_p_w,"imm_pagerank_and_enrichment.txt",quote = F,sep = "\t",row.names = F)
  setwd("../")
  remove(allxishu,ES,gene_enrichment_w,geneList,gmt1,gmt2,
         gmts,GSEAdata,GSEAresult2,keys,m,n,i,tT,w,gmt)
  
  ## 8、过滤imm基因
  lnc2imm2=lnc2imm2[which(lnc2imm2$imm %in% imm_p_w$ENSEMBL),]
  icp2imm2=icp2imm2[which(icp2imm2$imm %in% imm_p_w$ENSEMBL),]
  lnc2icp2=lnc2icp2[which(lnc2icp2$lnc %in% lnc2imm2$lnc),]
  
  dir.create("./GSEA后")
  setwd("./GSEA后")
  write.table(lnc2icp2,"lnc2icp_result.txt",quote = F,row.names = F,sep="\t")
  write.table(icp2imm2,"icp2imm_result.txt",quote = F,row.names = F,sep="\t")
  write.table(lnc2imm2,"lnc2imm_result.txt",quote = F,row.names = F,sep="\t")
  save(lnc_exp,file="lnc_exp.RData")
  save(ICP_exp,file="ICP_exp.RData")
  save(IMM_exp,file="IMM_exp.RData")
  setwd(path)
}





















