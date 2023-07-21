library(glmnet)
library(sva)
library(data.table)
library(tidyverse)
library(survival)
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

#### 划分训练集和检验集
set.seed(121)
train_index=sample(1:159, size = 0.7 * 159)
test_index=setdiff(1:159, train_index)

result=as.data.frame(t(result))
train=result[,train_index]
test=result[,test_index]
train=as.matrix(train)
test=as.matrix(test)
train=t(train)
test=t(test)

y_train=exp_resp[train_index]
y_test=exp_resp[test_index]

model = glmnet(train,y_train,family = "binomial",alpha = 0.5)
plot(model,xvar="lambda",label=T)
fit_cv <- cv.glmnet(train, y_train, alpha=1, family = 'binomial', type.measure='auc')
plot(fit_cv)

# 交叉验证选择最优lambda
cv_model=cv.glmnet(train,y_train,family="binomial",alpha=1)
best_lambda=cv_model$lambda.min
coef <- as.matrix (coef (model, s = best_lambda))
genes=rownames(coef)[coef!=0]

library(ROCR)
pred=predict(model,newx=test,s=cv_model$lambda.min,type = 'response')
get_confusion_stat <- function(pred,y_test,threshold=0.5){
  # auc
  tmp <- prediction(as.vector(pred),y_test)
  auc <- unlist(slot(performance(tmp,'auc'),'y.values'))
  # statistic
  pred_new <- as.integer(pred>threshold) 
  tab <- table(pred_new,y_test)
  if(nrow(tab)==1){
    print('preds all zero !')
    return(0)
  }
  TP <- tab[2,2]
  TN <- tab[1,1]
  FP <- tab[2,1]
  FN <- tab[1,2]
  accuracy <- round((TP+TN)/(TP+FN+FP+TN),4)
  recall_sensitivity <- round(TP/(TP+FN),4)
  precision <- round(TP/(TP+FP),4)
  specificity <- round(TN/(TN+FP),4)
  # 添加，预测的负例占比（业务解释：去除多少的样本，达到多少的recall）
  neg_rate <- round((TN+FN)/(TP+TN+FP+FN),4)
  re <- list('AUC' = auc,
             'Confusion_Matrix'=tab,
             'Statistics'=data.frame(value=c('accuracy'=accuracy,
                                             'recall_sensitivity'=recall_sensitivity,
                                             'precision'=precision,
                                             'specificity'=specificity,
                                             'neg_rate'=neg_rate)))
  return(re)
}


print(get_confusion_stat(pred,y_test))

pred.obj = prediction(pred, y_test)
auc.obj = performance(pred.obj, measure = "auc")
auc.value = as.numeric(auc.obj@y.values) 
print(auc.value)
roc.obj = performance(pred.obj, measure = "tpr", x.measure = "fpr")
plot(roc.obj, main = "ROC curve", col = "blue") 
abline(a = 0, b = 1, lty = 2, col = "gray") 
text(0.8, 0.2, paste("AUC =", round(auc.value, 3)))

roc.obj_plastic=roc.obj
auc.value_plastic=auc.value
