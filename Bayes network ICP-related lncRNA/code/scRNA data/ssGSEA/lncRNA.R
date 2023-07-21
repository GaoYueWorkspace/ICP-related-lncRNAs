# BiocManager::install(version ='3.17')
# devtools::install_github("ZJUFanLab/scCATCH")

library(Seurat)
library(scCATCH)
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(genefilter)
library(ggplot2)
library(SingleR)
library(celldex)
library(tidyverse)

## GSE127471

exp=read.table("GSE127471.txt")
pbmc=CreateSeuratObject(counts = exp,project = "GSE127471",min.cells = 3, min.features = 200)
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
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers=pbmc.markers %>%  group_by(cluster) %>%   slice_max(n = 4, order_by = avg_log2FC)
new.cluster.ids <- c("Megakaryocytes","Monocytes","CD8 T cells","CD4 T cells",
                     "NKT cells","B cells","NK cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
name=row.names(pbmc@meta.data)
name=as.data.frame(name)
name[,"cell"]=pbmc@active.ident
write.table(name,"GSE127471_anno.txt",quote = F,sep = "\t",row.names = F)

count<-as.data.frame(pbmc@assays$RNA@counts)
GENE_LIST<-read.table("result_lncRNA.txt",header = T,sep = "\t")
# genelist=list(GENE_LIST$SYMBOL)
genelist=intersect(GENE_LIST$SYMBOL,rownames(exp))
genelist=list(genelist)
ssgsea<- gsva(as.matrix(count), genelist,method='ssgsea',kcdf='Poisson',abs.ranking=TRUE)
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

################################################################################

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
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%  group_by(cluster) %>%   slice_max(n = 2, order_by = avg_log2FC)
# 假设sce是一个Seurat对象
pbmc <- AddMetaData(pbmc, metadata = exp_anno$cell, col.name = "cell.type")
Idents(pbmc) <- "cell.type"

DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

count<-as.data.frame(pbmc@assays$RNA@counts)
GENE_LIST<-read.table("result_lncRNA.txt",header = T,sep = "\t")
# genelist=list(GENE_LIST$SYMBOL)
genelist=intersect(GENE_LIST$SYMBOL,rownames(exp))
genelist=list(genelist)
ssgsea<- gsva(as.matrix(count), genelist,method='ssgsea',kcdf='Poisson',abs.ranking=TRUE)
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
  geom_point(size=0.8,pch=16)+
  theme_classic()+
  labs(title = "result lnc")+
  theme(text = element_text(size = 14),plot.title = element_text(hjust = 0.5))+
  scale_colour_gradientn(colours = c("#2e86c1","#7fb3d5","#f1948a", "#e74c3c")) +
  annotate("text", x = biao$x , y = biao$y,label = biao$label,colour="black",size=5) 

################################################################################
# GSE117570
exp=read.table("GSE117570.txt")
pbmc=CreateSeuratObject(counts = exp,project = "GSE117570",min.cells = 3, min.features = 200)
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
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
new.cluster.ids <- c("T cell","Macrophage","Stromal cells","Macrophage","T cell",
                     "Stromal cells","B cell","Macrophage","Macrophage","Endothelial cells",
                     "malignant Epithelial cells","malignant Epithelial cells","Macrophage",
                     "malignant Epithelial cells","Endothelial cells","malignant Epithelial cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
name=row.names(pbmc@meta.data)
name=as.data.frame(name)
name[,"cell"]=pbmc@active.ident
write.table(name,"GSE117570_anno.txt",quote = F,sep = "\t",row.names = F)

count<-as.data.frame(pbmc@assays$RNA@counts)
GENE_LIST<-read.table("result_lncRNA.txt",header = T,sep = "\t")
# genelist=list(GENE_LIST$SYMBOL)
genelist=intersect(GENE_LIST$SYMBOL,rownames(exp))
genelist=list(genelist)
ssgsea<- gsva(as.matrix(count), genelist,method='ssgsea',kcdf='Poisson',abs.ranking=TRUE)
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
  geom_point(size=0.8,pch=16)+
  theme_classic()+
  labs(title = "result lnc")+
  theme(text = element_text(size = 14),plot.title = element_text(hjust = 0.5))+
  scale_colour_gradientn(colours = c("#3498db","#2e86c1","#f1b4ae","#f1948a", "#e74c3c")) +
  annotate("text", x = biao$x , y = biao$y,label = biao$label,colour="black",size=5) 
