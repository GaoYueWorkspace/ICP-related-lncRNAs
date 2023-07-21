# library (utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages ("estimate", repos = rforge, dependencies = TRUE)

load("./表达谱/after_group_exp.RData")
library(AnnotationDbi)
library(org.Hs.eg.db)
symbol2ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(exp), columns = c("SYMBOL", "ENSEMBL"), keytype = "ENSEMBL")
symbol2ensembl <- symbol2ensembl[!duplicated(symbol2ensembl$SYMBOL), ]
symbol2ensembl <- symbol2ensembl[!duplicated(symbol2ensembl$ENSEMBL), ]
symbol2ensembl=na.omit(symbol2ensembl)
exp=exp[symbol2ensembl$ENSEMBL,]
rownames(exp) <- symbol2ensembl$SYMBOL[match(rownames(exp), symbol2ensembl$ENSEMBL)]
id=rownames(exp)
exp=cbind(id,exp)
write.table (exp, file = "PRAD-FR_exp.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 计算肿瘤纯度得分
library (estimate)
filterCommonGenes (input.f = "PRAD-FR_exp.tsv", output.f = "PRAD-FR_10412genes.gct", id = "GeneSymbol")
estimateScore (input.ds = "PRAD-FR_10412genes.gct", output.ds = "PRAD-FR_estimate_score.gct", platform = "illumina")

# 查看结果
scores <- read.table ("PRAD-FR_estimate_score.gct", skip = 2, header = TRUE)
rownames (scores) <- scores[, 1]
scores <- t (scores[, 3:ncol (scores)])
TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])
scores=as.data.frame(scores)
scores[,"TumorPurity"]=TumorPurity

write.table(scores,file = "inte_PRAD.txt",quote = F,sep = "\t")