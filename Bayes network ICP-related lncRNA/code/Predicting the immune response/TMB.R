# ### 安装缺失的R包-CRAN
# if(!requireNamespace("dplyr",quietly = TRUE)) install.packages("dplyr")
# if(!requireNamespace("DT",quietly = TRUE)) install.packages("DT")
# if(!requireNamespace("ggplot2",quietly = TRUE)) install.packages("ggplot2")
# if(!requireNamespace("reshape2",quietly = TRUE)) install.packages("reshape2")
# if(!requireNamespace("plyr",quietly = TRUE)) install.packages("plyr")
# if(!requireNamespace("survival",quietly = TRUE)) install.packages("survival")
# if(!requireNamespace("circlize",quietly = TRUE)) install.packages("circlize")
# if(!requireNamespace("lsmeans",quietly = TRUE)) install.packages("lsmeans")
# if(!requireNamespace("spatstat",quietly = TRUE)) install.packages("spatstat")
# if(!requireNamespace("corrplot",quietly = TRUE)) install.packages("corrplot")
# ### 安装缺失的R包-Bioconductor
# if(!requireNamespace("ComplexHeatmap",quietly = TRUE)) BiocManager::install("ComplexHeatmap")
# if(!requireNamespace("biomaRt",quietly = TRUE)) BiocManager::install("biomaRt")
# if(!requireNamespace("DESeq2",quietly = TRUE)) BiocManager::install("DESeq2")
# if(!requireNamespace("edgeR",quietly = TRUE)) BiocManager::install("edgeR")
# if(!requireNamespace("limma",quietly = TRUE)) BiocManager::install("limma")
# 
# install.packages("DESeq_1.34.1.zip")
# install.packages("./IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL, type = "source")
# load("IMvigor210CoreBiologies/data/cds.RData")
# load("IMvigor210CoreBiologies/data/dat19.RData")
# library(IMvigor210CoreBiologies)
# data("cds")
# 
# 

# if(!requireNamespace("easierData",quietly = TRUE)) BiocManager::install("easierData")
library(easierData)
library(SummarizedExperiment)
data <- get_Mariathasan2018_PDL1_treatment()
names(assays(data))
count <- assay(data,"counts")
phenoData <- as.data.frame(colData(data))


library(tidyverse)
anno=read.table("SYMBOL_ENSEMBL.txt",header = T)
load("genelength.RData")
gene=gene[which(gene$gene_id%in%anno$ENSG),]
for(i in 1:nrow(gene)){
  gene[i,"gene_name"]=anno[which(anno$ENSG%in%gene[i,1]),2]
}

count[,"gene_name"]=rownames(count)
tmp <- inner_join(count,gene,by = 'gene_name') # 匹配gene_count与gene_length
temple=tmp[,1:192]
temple%>% mutate(across(where(is.character), as.numeric))  -> temple
tmp=cbind(tmp[,193],temple,tmp[,195])
colnames(tmp)[1]="gene_name"
colnames(tmp)[194]="length"
tmp=tmp[!duplicated(tmp$gene_name),]
fpkm <- data.frame(row.names = tmp$gene_name)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}
exp=log2(1+fpkm)

write.table(exp,"IMvigor210.txt",quote = F,sep = "\t")
write.table(phenoData,"IMvigor210_pheno.txt",row.names = F,quote = F,sep = "\t")

phenoData=na.omit(phenoData)
exp=exp[,phenoData$pat_id]


TMB=phenoData[,c(2,3)]
TMB[which(TMB$BOR%in%"NR"),"group"]=0
TMB[which(TMB$BOR%in%"R"),"group"]=1


#### 划分训练集和检验集
set.seed(111)
train_index=sample(1:150, size = 0.7 * 150)
test_index=setdiff(1:150, train_index)


train=TMB[train_index,]
test=TMB[test_index,]



library(e1071)
library(pROC)
library(ROCR)

svm_model = svm(group~TMB,data=TMB,knernel = "radial")


summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)

roc.obj=performance(pred.obj,measure = "tpr",x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_TMB=roc.obj
auc.value_TMB=auc.value
##############################################################################
surv_triplet=read.table("surv_triplet.txt",header = T)
gene=union(union(surv_triplet$icp_symbol,surv_triplet$lnc_symbol),surv_triplet$imm_symbol)
exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)
surv_triplet=surv_triplet[which(surv_triplet$icp_symbol%in%gene &
                                  surv_triplet$lnc_symbol%in%gene &
                                  surv_triplet$imm_symbol%in%gene),]
surv_triplet=surv_triplet[,5:7]
change_exp=function(x){
  exp1=as.numeric(exp[x[1],])
  exp2=as.numeric(exp[x[2],])
  exp3=as.numeric(exp[x[3],])
  exp_triplet=(exp1+exp2+exp3)/3
  return(exp_triplet)
}

exp_triplet=apply(surv_triplet,1,change_exp)
exp_triplet=as.data.frame(exp_triplet)

surv_triplet[,"id"]=paste(surv_triplet$lnc_symbol,surv_triplet$icp_symbol,surv_triplet$imm_symbol,sep = ";")
colnames(exp_triplet)=surv_triplet$id

exp_triplet=cbind(TMB,exp_triplet)
exp_triplet=exp_triplet[,-1]

set.seed(111)
train_index=sample(1:150, size = 0.7 * 150)
test_index=setdiff(1:150, train_index)


train=exp_triplet[train_index,]
test=exp_triplet[test_index,]

svm_model = svm(group~.,data=exp_triplet,knernel = "radial")


summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)

roc.obj=performance(pred.obj,measure = "tpr",x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_TMBandtriplet=roc.obj
auc.value_TMBandtriplet=auc.value
#################################################################################
exp_triplet2=exp_triplet[,-1]

set.seed(111)
train_index=sample(1:150, size = 0.7 * 150)
test_index=setdiff(1:150, train_index)


train=exp_triplet2[train_index,]
test=exp_triplet2[test_index,]

svm_model = svm(group~.,data=exp_triplet2,knernel = "radial")


summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred

#绘制ROC曲线
ran_roc <- roc(test$group,as.numeric(svm_pred))
plot(ran_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="skyblue", print.thres=TRUE,main='SVM模型ROC曲线')

pred.obj=prediction(as.numeric(test$svm_pred),test$group)
auc.obj=performance(pred.obj,measure = "auc")
auc.value=as.numeric(auc.obj@y.values)

roc.obj=performance(pred.obj,measure = "tpr",x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_triplet=roc.obj
auc.value_triplet=auc.value


plot(roc.obj_TMB, main = "ROC curve", col = "#509ee0",lwd=2) 
plot(roc.obj_triplet, main = "ROC curve", col = "#e27823",add=T,lwd=2) 
plot(roc.obj_TMBandtriplet, main = "ROC curve", col = "#96bba3",add=T,lwd=2) 


abline(a = 0, b = 1, lty = 2, col = "gray",lwd=1) 

barplot(c(auc.value_triplet,auc.value_TMB,auc.value_TMBandtriplet),ylim = c(0,0.8))
