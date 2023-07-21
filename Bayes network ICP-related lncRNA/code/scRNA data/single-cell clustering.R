path=getwd()
library(readxl)
library(Seurat)
library(tidyverse)
########################################################################################
GSE72056=read.table("单细胞数据/GSE72056_melanoma_single_cell_revised_v2.txt")
GSE72056_anno=GSE72056[1:4,]
GSE72056_anno=as.data.frame(t(GSE72056_anno))
colnames(GSE72056_anno)=GSE72056_anno[1,]
GSE72056_anno=GSE72056_anno[which(GSE72056_anno$`malignant(1=no,2=yes,0=unresolved)`==1),]
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==1),4]="T cell"
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==2),4]="B cell"
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==3),4]="Macro"
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==4),4]="Endo"
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==5),4]="CAF"
GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`==6),4]="NK"
GSE72056_anno=GSE72056_anno[which(GSE72056_anno$`non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)`!=0),]
GSE72056_cell=GSE72056_anno[,c(1,4)]
colnames(GSE72056_cell)=c("sample","cell")
write.table(GSE72056_cell,"GSE72056_anno.txt",quote = F,sep = "\t",row.names = F)

GSE72056=GSE72056[-c(2:4),]
colnames(GSE72056)=GSE72056[1,]
GSE72056=GSE72056[-1,]
GSE72056=GSE72056[!duplicated(GSE72056[,1]),]
rownames(GSE72056)=GSE72056[,1]
GSE72056=GSE72056[,-1]
GSE72056=GSE72056[,GSE72056_cell$sample]
write.table(GSE72056,"GSE72056.txt",sep = "\t",quote = F)

exp=GSE72056
name=colnames(exp)
name=gsub("-", "_", name)
colnames(exp)=name
pbmc=CreateSeuratObject(counts = exp,project = "GSE97681",min.cells = 3, min.features = 200)
# 四、数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 五、细胞聚类
# 1、寻找mark基因,鉴定细胞间表达量高变的基因,选取2000个高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# 2、改变每个基因的表达，使细胞的平均表达为0。缩放每个基因的表达，使细胞间的差异为1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 3、PCA线性降维，默认会返回50个主成分
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 4、PCA可视化
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# 5、主成分数目的选择
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)# 选择10个主成分

# 6、细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) # 聚成了7类


#7、非线性降维——UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10, label = T)
p1 <- DimPlot(pbmc, reduction = "umap")
p1

#8、非线性降维——T-SNE
pbmc <- RunTSNE(pbmc, dims = 1:10)
p2 <- DimPlot(pbmc, reduction = "tsne")
p2

# 六、细胞类型注释
# 1、与所有其他亚群相比，找到每个亚群的标记，仅报告阳性细胞，每类取前两个
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%  group_by(cluster) %>%   slice_max(n = 2, order_by = avg_log2FC)
# 2、细胞注释
library(SingleR)
library(celldex)
# 准备人类初级细胞分类数据
humanRNA=DatabaseImmuneCellExpressionData()
humanRNA2=HumanPrimaryCellAtlasData()
humanRNA3=BlueprintEncodeData()
sce_for_SingleR=GetAssayData(pbmc,slot = "data")
clusters=pbmc@meta.data$seurat_clusters
# 进行细胞注释
pred.humanRNA=SingleR(test = sce_for_SingleR,ref = list(humanRNA,humanRNA2,humanRNA3),
                      labels = list(humanRNA$label.main,humanRNA2$label.main,humanRNA3$label.main),
                      method = "cluster",clusters = clusters,
                      assay.type.test="logcounts",
                      assay.type.ref = "logcounts")
# 提取并保存注释信息
cellType=data.frame(ClusterID=levels(pbmc@meta.data$seurat_clusters),
                    celltype=pred.humanRNA$labels)
pbmc@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),"celltype"]
# 3、注释后的可视化
new.cluster.ids <- cellType$celltype
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
write.table(name,"GSE72056_anno.txt",quote = F,sep = "\t",row.names = F)

##############################################################################


GSE69405=read.table("GSE69405.txt",header = T)
GSE69405=GSE69405[!duplicated(GSE69405$gene_name),]
rownames(GSE69405)=GSE69405$gene_name
exp=GSE69405
exp=exp[,-c(1:3)]
name=as.data.frame(colnames(exp))
name=gsub("\\.", "_", name$`colnames(exp)`)
colnames(exp)=name
write.table(exp,"GSE69405.txt",sep = "\t",quote = F)
pbmc=CreateSeuratObject(counts = exp,project = "GSE97681",min.cells = 3, min.features = 200)
# 四、数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 五、细胞聚类
# 1、寻找mark基因,鉴定细胞间表达量高变的基因,选取2000个高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# 2、改变每个基因的表达，使细胞的平均表达为0。缩放每个基因的表达，使细胞间的差异为1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 3、PCA线性降维，默认会返回50个主成分
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 4、PCA可视化
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# 5、主成分数目的选择
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)# 选择10个主成分

# 6、细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.5) # 聚成了7类


#7、非线性降维——UMAP
pbmc <- RunUMAP(pbmc, dims = 1:5, label = T)
p1 <- DimPlot(pbmc, reduction = "umap")
p1

#8、非线性降维——T-SNE
pbmc <- RunTSNE(pbmc, dims = 1:5)
p2 <- DimPlot(pbmc, reduction = "tsne")
p2

# 六、细胞类型注释
# 1、与所有其他亚群相比，找到每个亚群的标记，仅报告阳性细胞，每类取前两个
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%  group_by(cluster) %>%   slice_max(n = 2, order_by = avg_log2FC)
# 2、细胞注释
library(SingleR)
library(celldex)
# 准备人类初级细胞分类数据
humanRNA=DatabaseImmuneCellExpressionData()
humanRNA2=HumanPrimaryCellAtlasData()
humanRNA3=BlueprintEncodeData()
sce_for_SingleR=GetAssayData(pbmc,slot = "data")
clusters=pbmc@meta.data$seurat_clusters
# 进行细胞注释
pred.humanRNA=SingleR(test = sce_for_SingleR,ref = list(humanRNA),
                      labels = list(humanRNA$label.main),
                      method = "cluster",clusters = clusters,
                      assay.type.test="logcounts",
                      assay.type.ref = "logcounts")
# 提取并保存注释信息
cellType=data.frame(ClusterID=levels(pbmc@meta.data$seurat_clusters),
                    celltype=pred.humanRNA$labels)
pbmc@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),"celltype"]
# 3、注释后的可视化
new.cluster.ids <- cellType$celltype
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
name=row.names(pbmc@meta.data)
name=as.data.frame(name)
name[,"cell"]=pbmc@meta.data$singleR
write.table(name,"GSE69405_anno.txt",quote = F,sep = "\t",row.names = F)

pbmc.data=Read10X("GSE127471")
# 二、创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data,# 未标准化的数据：原始计数或TPMs
                           project = "HCC",# 项目名称
                           min.cells = 3,# 包含至少在这些细胞检测到的features
                           min.features = 200,# 包含至少检测到这些features的细胞
                           names.field = 1,
                           names.delim = "_")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")# 线粒体基因比例
head(pbmc@meta.data, 5)
# 三、数据预处理——QC和细胞筛选
# nFeature_RNA代表每个细胞测到的基因数目，去除掉小于200的，或大于2500的。
# nCount_RNA代表每个细胞测到所有基因的表达量之和。
# percent.mt代表测到的线粒体基因的比例，去除掉大于5%的。
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# 四、数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 五、细胞聚类
# 1、寻找mark基因,鉴定细胞间表达量高变的基因,选取2000个高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)


# 2、改变每个基因的表达，使细胞的平均表达为0。缩放每个基因的表达，使细胞间的差异为1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 3、PCA线性降维，默认会返回50个主成分
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 4、PCA可视化
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# 5、主成分数目的选择
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
ElbowPlot(pbmc)# 选择10个主成分

# 6、细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5) # 聚成了7类


#7、非线性降维——UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10, label = T)
p1 <- DimPlot(pbmc, reduction = "umap")
p1

#8、非线性降维——T-SNE
pbmc <- RunTSNE(pbmc, dims = 1:10)
p2 <- DimPlot(pbmc, reduction = "tsne")
p2

# 六、细胞类型注释
# 1、与所有其他亚群相比，找到每个亚群的标记，仅报告阳性细胞，每类取前两个
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%  group_by(cluster) %>%   slice_max(n = 2, order_by = avg_log2FC)
# 2、细胞注释
library(SingleR)
library(celldex)
# 准备人类初级细胞分类数据
humanRNA=DatabaseImmuneCellExpressionData()
humanRNA2=HumanPrimaryCellAtlasData()
humanRNA3=BlueprintEncodeData()
sce_for_SingleR=GetAssayData(pbmc,slot = "data")
clusters=pbmc@meta.data$seurat_clusters
# 进行细胞注释
pred.humanRNA=SingleR(test = sce_for_SingleR,ref = list(humanRNA,humanRNA2,humanRNA3),
                      labels = list(humanRNA$label.main,humanRNA2$label.main,humanRNA3$label.main),
                      method = "cluster",clusters = clusters,
                      assay.type.test="logcounts",
                      assay.type.ref = "logcounts")
# 提取并保存注释信息
cellType=data.frame(ClusterID=levels(pbmc@meta.data$seurat_clusters),
                    celltype=pred.humanRNA$labels)
pbmc@meta.data$singleR=cellType[match(clusters,cellType$ClusterID),"celltype"]
# 3、注释后的可视化
new.cluster.ids <- cellType$celltype
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()
name=row.names(pbmc@meta.data)
name=as.data.frame(name)
name[,"cell"]=pbmc@meta.data$singleR
write.table(name,"GSE127471_anno.txt",quote = F,sep = "\t",row.names = F)
exp=pbmc@assays$RNA@counts
exp=as.data.frame(exp)
write.table(exp,"GSE127471.txt",sep="\t",quote = F)

GSE117570=read_rds("GSE117570.integrated.rds")
exp=GSE117570@assays$RNA@counts
exp=as.data.frame(exp)
write.table(exp,"GSE117570.txt",sep="\t",quote = F)
name=row.names(GSE117570@meta.data)
name=as.data.frame(name)
name[,"cell"]=GSE117570@meta.data$singleR
write.table(name,"GSE117570_anno.txt",quote = F,sep = "\t",row.names = F)

GSE144945=read_rds("GSE144945_harmony.rds")
exp=GSE144945@assays$RNA@counts
exp=as.data.frame(exp)
write.table(exp,"GSE144945.txt",sep="\t",quote = F)
name=row.names(GSE144945@meta.data)
name=as.data.frame(name)
name[,"cell"]=GSE144945@meta.data$singleR
write.table(name,"GSE144945_anno.txt",quote = F,sep = "\t",row.names = F)

