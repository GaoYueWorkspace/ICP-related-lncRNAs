library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
pbmc.data=Read10X("GSE127471")
pbmc <- CreateSeuratObject(counts = pbmc.data,# 未标准化的数据：原始计数或TPMs
                           project = "HCC",# 项目名称
                           min.cells = 3,# 包含至少在这些细胞检测到的features
                           min.features = 200,# 包含至少检测到这些features的细胞
                           names.field = 1,
                           names.delim = "_")
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
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
remove(all.genes,new.cluster.ids,p1,p2,pbmc.data,markers)

sample_ann=pbmc@meta.data
head(sample_ann)
gene_ann=data.frame(
  gene_short_name=rownames(pbmc@assays$RNA),
  row.names = rownames(pbmc@assays$RNA)
)
head(gene_ann)
ct=as.data.frame(pbmc@assays$RNA@counts)
ct[1:4,1:4]

expr_matrix <-  as.matrix(ct)

#构建CDS对象
pd <- new("AnnotatedDataFrame",data = sample_ann)
fd <- new("AnnotatedDataFrame",data = gene_ann)
#利用Seurat数据创建monocle对象
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())
remove(ct,sample_ann,gene_ann)
HSMM=cds
## 归一化 
HSMM=estimateSizeFactors(HSMM)
HSMM=estimateDispersions(HSMM)
#过滤低质量的细胞,这一操作会在fData(cds)中添加一列num_cells_expressed
HSMM <- detectGenes(HSMM, min_expr = 0.1 )
print(head(fData(HSMM)))
#过滤掉小于在10个细胞中表达的基因，还剩943个基因 
# expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
# print(head(pData(HSMM)))

HSMM@phenoData@data$celltype=pbmc@active.ident
#Run ordering algorithm
# var_genes <- pbmc[["RNA"]]@var.features
# ordering_genes <- var_genes

disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
# 

HSMM <- setOrderingFilter(HSMM, disp.genes)
print(dim(exprs(HSMM)))
plot_ordering_genes(HSMM)
plot_pc_variance_explained(HSMM,return_all = F)


## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM <- reduceDimension(HSMM, 
                        reduction_method="tSNE",
                        max_components=2,
                        num_dim=7,
                        verbose=TRUE)
# HSMM=clusterCells(HSMM,num_clusters = 7)# 无法运行
# 使用seurat的分群
pData(HSMM)$Cluster=pData(HSMM)$celltype

table(pData(HSMM)$Cluster)
# First decide what you want to color your cells by
print(head(pData(HSMM)))

# 寻找差异基因
Sys.time()
diff_test_res <- differentialGeneTest(HSMM,
                                      fullModelFormulaStr = "~Cluster")
Sys.time()
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes=sig_genes[order(sig_genes$pval),]
head(sig_genes[,c("gene_short_name", "pval", "qval")] ) 
cg=as.character(head(sig_genes$gene_short_name)) 

# 挑选差异最显著的基因可视化
plot_genes_jitter(HSMM[cg,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
cg2=as.character(tail(sig_genes$gene_short_name)) 
plot_genes_jitter(HSMM[cg2,],
                  grouping = "Cluster",
                  color_by = "Cluster",
                  nrow= 3,
                  ncol = NULL )
## 拟时序分析
# 第一步: 挑选合适的基因. 有多个方法，例如提供已知的基因集，这里选取统计学显著的差异基因列表
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
head(ordering_genes)
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)


# 第二步: 降维。降维的目的是为了更好的展示数据。函数里提供了很多种方法,
# 不同方法的最后展示的图都不太一样, 其中“DDRTree”是Monocle2使用的默认方法
HSMM <- reduceDimension(HSMM, max_components = 2,
                       method = 'DDRTree')

# 第三步: 对细胞进行排序
HSMM <- orderCells(HSMM)
p1=plot_cell_trajectory(HSMM, 
                     color_by = "Pseudotime",
                     theta = -15,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 4) + theme(legend.position = "right")
p2=plot_cell_trajectory(HSMM, 
                     color_by = "celltype",
                     theta = -15,
                     show_branch_points = FALSE,
                     show_tree = TRUE,
                     cell_size = 4) + theme(legend.position = "right")
p1|p2
dev.off()
# now onwards to pseudotemporal plots
sig_gene_names <- read.table("result_lncRNA.txt",header = T)[,2]
sig_gene_names=sig_gene_names[which(sig_gene_names%in%rownames(HSMM))]
head(sig_gene_names)
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 6, 
                        cores = 4,
                        hmcols = NULL,
                        show_rownames = T,
                        return_heatmap = T)
dev.off()


