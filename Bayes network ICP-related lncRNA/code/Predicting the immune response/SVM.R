library(e1071)
library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
library(pROC)

#### 读取五个数据集
Gide_PD1=read.table("Gide_PD1.txt",header = T)
Gide_PD1CTLA4=read.table("Gide_PD1+CTLA4.txt",header = T)
Naive=read.table("Naive.txt",header = T)
Prog=read.table("Prog.txt",header = T)
Vanllen=read.table("Vanllen.txt",header = T)

#### 去除批次效应
tmp=intersect(rownames(Gide_PD1),rownames(Gide_PD1CTLA4))
tmp=intersect(tmp,rownames(Naive))
tmp=intersect(tmp,rownames(Prog))
tmp=intersect(tmp,rownames(Vanllen))


Gide_PD1=Gide_PD1[tmp,]
Gide_PD1CTLA4=Gide_PD1CTLA4[tmp,]
Naive=Naive[tmp,]
Prog=Prog[tmp,]
Vanllen=Vanllen[tmp,]


exp_all=cbind(Gide_PD1,Gide_PD1CTLA4,Naive,Prog,Vanllen)
batch=paste0("batch",rep(c(1,2,3,4,5),c(41,32,25,26,35)))

exp_all=as.matrix(exp_all)
exp <-  ComBat(dat = exp_all, batch = batch, mod = NULL)
exp=as.data.frame(exp)


Gide_PD1_resp=read.table("Gide_PD1_clinical.txt",header = T)
Gide_PD1CTLA4_resp=read.table("Gide_PD1+CTLA4_clinical.txt",header = T)
Naive_resp=read.table("Naive_clinical.txt",header = T)
Prog_resp=read.table("Prog_clinical.txt",header = T)
Vanllen_resp=read.table("Vanllen_clinical.txt",header = T)

surv_triplet=read.table("surv_triplet.txt",header = T)
gene=union(union(surv_triplet$icp_symbol,surv_triplet$lnc_symbol),surv_triplet$imm_symbol)

exp=exp[gene,]
exp=na.omit(exp)
gene=rownames(exp)

surv_triplet=surv_triplet[which(surv_triplet$icp_symbol%in%gene &
                                  surv_triplet$lnc_symbol%in%gene &
                                  surv_triplet$imm_symbol%in%gene),]
surv_triplet=surv_triplet[,5:7]

exp=as.data.frame(t(exp))

os=rbind(Gide_PD1_resp[,2:3],Gide_PD1CTLA4_resp[,2:3],Naive_resp[,2:3],
         Prog_resp[,2:3],Vanllen_resp[,2:3])
exp=cbind(os,exp)
colnames(exp)=gsub("-","\\.",colnames(exp))

surv_triplet2=apply(surv_triplet,1,function(x){gsub("-","\\.",x)})
surv_triplet2=as.data.frame(t(surv_triplet2))

risk=function(x){
  formula <- as.formula(paste("Surv(OS.time, OS) ~", paste(x, collapse = "+")))
  fit <- coxph(formula, data = exp)
  coef <- as.numeric(fit$coef) # 回归系数
  score1=coef[1]*(as.numeric(exp[,x[1]]))
  score2=coef[2]*(as.numeric(exp[,x[2]]))
  score3=coef[3]*(as.numeric(exp[,x[3]]))
  score=score1+score2+score3
  return(score)
}

result=apply(surv_triplet2,1,risk)
name=paste(surv_triplet2$lnc_symbol,surv_triplet2$icp_symbol,surv_triplet2$imm_symbol,sep = ";")
colnames(result)=name

exp_resp=c(Gide_PD1_resp$Response,
           Gide_PD1CTLA4_resp$Response,
           Naive_resp$Response,
           Prog_resp$Response,
           Vanllen_resp$response)

data=cbind(result,exp_resp)
data=as.data.frame(data)
colnames(data)[4]="group"
data$group=factor(data$group)

#### 划分训练集和检验集
set.seed(121)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)


train=data[train_index,]
test=data[test_index,]

# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`,data=train,knernel = "radial")
# svm_model = svm(group~`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
svm_model = svm(group~`DBH.AS1;TIGIT;LTA`+`FAM30A;NCR3;IL16`,data=train,knernel = "radial")
# svm_model = svm(group~`FAM30A;NCR3;IL16`+`ZEB2.AS1;BTK;TNFSF13B`,data=train,knernel = "radial")
# svm_model = svm(group~.,data=train,knernel = "radial")


summary(svm_model)

svm_pred=predict(svm_model,test,decision.values = TRUE)
test$svm_pred = svm_pred
head(test)
table(test$group,test$svm_pred)

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

roc.obj_vector=roc.obj
auc.value_vector=auc.value
