## GSE72056
exp=read.table("GSE72056.txt")
exp_anno=read.table("GSE72056_anno.txt",header = T,sep = "\t")
pbmc=CreateSeuratObject(counts = exp,project = "GSE72056",min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)# 选择10个主成分
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) # 聚成了7类
pbmc <- RunUMAP(pbmc, dims = 1:10, label = T)
p1 <- DimPlot(pbmc, reduction = "umap")
p1
pbmc <- RunTSNE(pbmc, dims = 1:10)
p2 <- DimPlot(pbmc, reduction = "tsne")
p2
# 假设sce是一个Seurat对象
pbmc <- AddMetaData(pbmc, metadata = exp_anno$cell, col.name = "cell.type")
Idents(pbmc) <- "cell.type"

DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()


path=getwd()

COAD=read.table("./COAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
ESCA=read.table("./ESCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
HNSC=read.table("./HNSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
KIRC=read.table("./KIRC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LIHC=read.table("./LIHC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LUAD=read.table("./LUAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
LUSC=read.table("./LUSC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
PRAD=read.table("./PRAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
READ=read.table("./READ/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
STAD=read.table("./STAD/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
THCA=read.table("./THCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
UCEC=read.table("./UCEC/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
BRCA=read.table("./BRCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
BLCA=read.table("./BLCA/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
SKCM=read.table("./SKCM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
CHOL=read.table("./CHOL/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
GBM=read.table("./GBM/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]
OV=read.table("./OV/包含免疫基因的lnc-icp关系对/LMI_genepair.txt",sep="\t",header = T)[,1:3]

triplet_all=list("BLCA"=BLCA,"BRCA"=BRCA,"CHOL"=CHOL,"COAD"=COAD,"ESCA"=ESCA,
                 "GBM"=GBM,"HNSC"=HNSC,"KIRC"=KIRC,"LIHC"=LIHC,"LUAD"=LUAD,
                 "LUSC"=LUSC,"OV"=OV,"PRAD"=PRAD,"READ"=READ,"SKCM"=SKCM,
                 "STAD"=STAD,"THCA"=THCA,"UCEC"=UCEC)
remove(COAD,BLCA,BRCA,CHOL,ESCA,GBM,HNSC,KIRC,LIHC,
       LUAD,LUSC,OV,PRAD,READ,SKCM,STAD,THCA,UCEC)
setwd(path)
triplet=do.call(rbind,triplet_all)
remove(triplet_all)
triplet=unique(triplet)
name=rownames(exp)
result_lnc=read.table("result_lncRNA.txt",header = T)
result_imm=read.table("result_imm.txt",header = T)
result_icp=read.table("result_icp.txt",header = T)
result_lnc=result_lnc[which(result_lnc$SYMBOL%in%name),]
result_icp=result_icp[which(result_icp$SYMBOL%in%name),]
result_imm=result_imm[which(result_imm$SYMBOL%in%name),]

triplet=triplet[which(triplet$lnc%in%result_lnc$ENSG),]
triplet=triplet[which(triplet$icp%in%result_icp$ENSG),]
triplet=triplet[which(triplet$imm%in%result_imm$ENSG),]

anno=read.table("SYMBOL_ENSEMBL.txt",header = T)
for(i in 1:nrow(triplet)){
  triplet[i,"lnc_anno"]=anno[which(anno$ENSG%in%triplet[i,2]),2]
  triplet[i,"icp_anno"]=anno[which(anno$ENSG%in%triplet[i,1]),2]
  triplet[i,"imm_anno"]=anno[which(anno$ENSG%in%triplet[i,3]),2]
  
}
triplet=triplet[,4:6]
rownames(triplet)=1:nrow(triplet)
triplet[,"id"]=paste(triplet$lnc_anno,triplet$icp_anno,triplet$imm_anno,sep = ";")

change_exp=function(x){
  exp1=as.numeric(exp[x[1],])
  exp2=as.numeric(exp[x[2],])
  exp3=as.numeric(exp[x[3],])
  exp_triplet=(exp1+exp2+exp3)/3
  return(exp_triplet)
}

exp_triplet=apply(triplet,1,change_exp)
exp_triplet=as.data.frame(exp_triplet)
exp_triplet=as.data.frame(t(exp_triplet))
colnames(exp_triplet)=colnames(exp)
rownames(exp_triplet)=triplet$id

triplet_use=read.table("surv_triplet.txt",header = T)
triplet_use[,"id"]=paste(triplet_use$lnc_symbol,triplet_use$icp_symbol,triplet_use$imm_symbol,sep = ";")

genelist=list(triplet_use$id)
ssgsea<- gsva(as.matrix(exp_triplet), genelist,method='ssgsea',kcdf='Poisson',abs.ranking=TRUE)
gsva_matrix <-ssgsea
tsne1<-as.data.frame(pbmc@reductions$tsne@cell.embeddings)
colnames(gsva_matrix)[1]==rownames(tsne1)[1]
tsne1<-cbind(tsne1,t(gsva_matrix))
tsne1$label<-Idents(pbmc)
#标签的坐标
biaox<-aggregate(tsne1$tSNE_1,by=list(tsne1$label),median)
biaoy<-aggregate(tsne1$tSNE_2,by=list(tsne1$label),median)
biao<-data.frame(label=biaox$Group.1,x=biaox$x,y=biaoy[,2])
colnames(tsne1)[3]="ssgsea"
ggplot(tsne1,aes(x=tSNE_1,y=tSNE_2,color=ssgsea))+
  geom_point(size=0.7,pch=16)+
  theme_classic()+
  labs(title = "result lnc")+
  theme(text = element_text(size = 14),plot.title = element_text(hjust = 0.5))+
  scale_colour_gradientn(colours = c("#3498db","#2e86c1","#f1b4ae","#f1948a", "#e74c3c")) +
  annotate("text", x = biao$x , y = biao$y,label = biao$label,colour="black",size=5) 

boxplot(ssgsea ~ label, data = tsne1,outline=F)
