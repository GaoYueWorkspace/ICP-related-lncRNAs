path=getwd()
library(clusterProfiler)
#library(org.Hs.eg.db)
#library(edgeR)# 差异分析时使用
#library(limma)# 可做差异分析
library(ppcor)
library(progress)
library(data.table)
library(tidyr)
library(dplyr)
SYMBOL_ENSEMBL=read.table("./anno_file/SYMBOL_ENSEMBL.txt",sep="\t",header=T)
lncRNA=read.table("./anno_file/lncRNA_ENSG.txt")
colnames(lncRNA)=c("ENSG","SYMBOL")
ICP=read.table("./anno_file/ICP.txt",sep="\t",header=T)
IMM=read.table("./anno_file/imm_gene.txt",sep="\t",header=T)

cell_list=readRDS("./cell_list/cell_list.rds")
cell_name=names(cell_list)

for(c in 1:length(cell_name)){
  setwd(paste("./",cell_name[c],sep = ""))
  lnc2icp2=read.table("./GSEA/lnc2icp_result.txt",header=T)
  lnc2imm2=read.table("./GSEA/lnc2imm_result.txt",header=T)
  icp2imm2=read.table("./GSEA/icp2imm_result.txt",header=T)
  load("./GSEA/ICP_exp.RData")
  load("./GSEA/lnc_exp.RData")
  load("./GSEA/IMM_exp.RData")
  ICP_exp=as.data.frame(t(ICP_exp))
  lnc_exp=as.data.frame(t(lnc_exp))
  IMM_exp=as.data.frame(t(IMM_exp))
  imm_p_w=read.table("./enrich/imm_pagerank_and_enrichment.txt",sep = "\t",header = T)
  #### 十一、对lncRNA-ICP关系对进行打分,S得分 ####
  gene_pair=lnc2icp2[,c(1,2)]
  
  ### 2、计算每个lnc-icp基因对有多少个imm基因
  imm_table=function(x){
    imm=lnc2imm2[which(lnc2imm2$lnc %in% x[1]),2]
    return(length(imm))
  }
  system.time({
    lnc2icp_imm=apply(gene_pair,1,imm_table)
  })
  lnc2icp_imm=as.data.frame(lnc2icp_imm)
  lnc2icp_imm=cbind(gene_pair,lnc2icp_imm)
  colnames(lnc2icp_imm)=c("lnc","icp","imm_num")
  remove(imm_table)
  gc()
  
  S_calculate=function(x){# x是lnc基因
    # 1、首先要找到每个lnc基因与多少imm基因存在关联
    imm=lnc2imm2[which(lnc2imm2$lnc %in% x[1]),2]
    weight=imm_p_w[which(imm_p_w$ENSEMBL %in% imm),"w"]
    sums=0
    for(i in 1:length(imm)){
      cor1=lnc2imm2[which(lnc2imm2$lnc %in% x[1] & lnc2imm2$imm %in% imm[i]),3]
      if(length(cor1)==0){next}
      p1=lnc2imm2[which(lnc2imm2$lnc %in% x[1] & lnc2imm2$imm %in% imm[i]),4]
      # if(length(p1)==0){next}
      cor2=icp2imm2[which(icp2imm2$imm %in% imm[i] & icp2imm2$icp %in% x[2]),3]
      if(length(cor2)==0){next}
      p2=icp2imm2[which(icp2imm2$imm %in% imm[i] & icp2imm2$icp %in% x[2]),4]
      # if(length(p2)==0){next}
      # 公式：S=w1*(-log10(p1)*sign(cor1))+w2*(-log10(p2)*sign(cor2))+……wi*(-log10(pi)*sign(cori))
      sums=sums+weight[i]*((-log10(p1))*sign(cor1)+(-log10(p2))*sign(cor2))
    }
    return(sums)
  }
  system.time({
    S=apply(gene_pair,1,S_calculate)
  })
  S=as.data.frame(S)# lnc-icp关系对的打分
  S=cbind(gene_pair,S)
  gc()
  #### 十二、进行permutation ####
  permutation_LI=lnc2icp2[,c(1,2)]
  permutation_LI[,"num"]=0
  
  ################################################################################
  ### 3、函数3
  S_calculate=function(x){# x是lnc基因
    # 1、首先要找到每个lnc基因与多少imm基因存在关联
    imm=lnc2imm2[which(lnc2imm2$lnc %in% x[1]),2]
    weight=imm_p_w[which(imm_p_w$ENSEMBL %in% imm),"w"]
    sums=0
    for(i in 1:length(imm)){
      cor1=lnc2imm2[which(lnc2imm2$lnc %in% x[1] & lnc2imm2$imm %in% imm[i]),3]
      if(length(cor1)==0){next}
      p1=lnc2imm2[which(lnc2imm2$lnc %in% x[1] & lnc2imm2$imm %in% imm[i]),4]
      # if(length(p1)==0){next}
      cor2=icp2imm2[which(icp2imm2$imm %in% imm[i] & icp2imm2$icp %in% x[2]),3]
      if(length(cor2)==0){next}
      p2=icp2imm2[which(icp2imm2$imm %in% imm[i] & icp2imm2$icp %in% x[2]),4]
      # if(length(p2)==0){next}
      sums=sums+weight[i]*((-log10(p1))*sign(cor1)+(-log10(p2))*sign(cor2))
    }
    return(sums)
  }
  ################################################################################
  pb <- progress_bar$new(total=1000)
  system.time({
    for(i in 1:1000){# 1000次
      pb$tick()
      # 1、随机扰动IMM基因的样本标签
      col=colnames(IMM_exp)
      IMM_exp2=IMM_exp[,sample(col,length(col))]
      colnames(IMM_exp2)=col
      
      S2=apply(gene_pair,1,S_calculate)
      S2=as.data.frame(S2)# lnc-icp关系对的打分
      S2=cbind(gene_pair,S2)
      
      # 3、比较大小，用于计算P值
      for(j in 1:length(S$S)){
        if(S[j,3]<S2[j,3]){
          permutation_LI[j,3]=permutation_LI[j,3]+1
        }
        else{permutation_LI[j,3]=permutation_LI[j,3]+0}
      }
      ### 随机扰动的结果储存在permuta_num 中
    }
  })# 计时结束
  remove(i,j,S2)
  permutation_LI[,"P"]=permutation_LI$num/1000
  
  # 进行FDR校正
  permutation_LI[,"FDR"]=p.adjust(permutation_LI$P,method = "BH")
  permutation_LI[,"imm_num"]=lnc2icp_imm$imm_num
  dir.create("./permutation_result")
  setwd("./permutation_result")
  write.table(permutation_LI,"permutation_LI.txt",quote = F,row.names = F,sep = "\t")
  setwd("../")
  gc()
  
  #### 十三、提取包含免疫基因的lnc-icp关系对的合格的贝叶斯推断结果 ####
  
  LMI=permutation_LI[which(permutation_LI$FDR<0.01),c(1,2)]# 包含的
  LI=permutation_LI[which(permutation_LI$FDR>=0.01),c(1,2)]# 不包含的
  dir.create("./if_contain_imm")
  setwd("./if_contain_imm")
  write.table(LMI,"LMI.txt",quote = F,row.names = F,sep = "\t")
  write.table(LI,"LI.txt",quote = F,row.names = F,sep = "\t")
  setwd("../")
  
  # 更新lnc2icp基因对
  # 未改动，直接使用了LMI
  lnc2icp=LMI
  
  #### 十三、拼接三元组 ####
  lnc2icp=lnc2icp2[,1:2]
  lnc2imm=lnc2imm2[,1:2]
  icp2imm=icp2imm2[,1:2]
  lnc_icp_imm=merge(lnc2icp,icp2imm,by="icp") 
  lnc_imm_icp=merge(lnc2imm,icp2imm,by="imm") 
  icp_lnc_imm=merge(lnc2icp,lnc2imm,by='lnc') 
  a=c("lnc","imm","icp")
  a=sort(a)
  lnc_icp_imm=lnc_icp_imm[,a]
  lnc_imm_icp=lnc_imm_icp[,a]
  icp_lnc_imm=icp_lnc_imm[,a]
  triplet=rbind(lnc_icp_imm,lnc_imm_icp,icp_lnc_imm) 
  triplet=unique(triplet)# icp,imm,lnc
  
  #### 十四、贝叶斯推断 ####
  bugg <- as.numeric(-Inf)
  library(bnlearn)
  modeSelection=function(lnc,icp,imm){
    learning.test=data.frame(M=lnc,E=icp,G=imm)# 输入准备构建贝叶斯网络的数据
    # 模式1
    LGI = empty.graph(names(learning.test))# 创建空的贝叶斯网络
    modelstring(LGI) = "[M][G|M][E|G]"
    # 模式2
    LIG = empty.graph(names(learning.test))
    modelstring(LIG) = "[M][E|M][G|E]"
    # 模式3
    INDEP = empty.graph(names(learning.test))
    modelstring(INDEP) = "[M][G|M][E|M]"
    # 模式4
    COO = empty.graph(names(learning.test))
    modelstring(COO) = "[M][G][E|M:G]"
    
    #### 计算AIC值
    AICS=c(AIC_LGI=score(LGI,learning.test,type="aic-g"),
           AIC_LIG=score(LIG,learning.test,type="aic-g"),
           AIC_INDEP=score(INDEP,learning.test,type="aic-g"),
           AIC_COO=score(COO,learning.test,type="aic-g")
    )
    # ifelse((bugg %in% AICS),AICS=rep(1,4),AICS)
    if(bugg %in% AICS){
      AICS=rep(1,4)
    }else{AICS=AICS}
    AICS=AICS[order(AICS,decreasing = F)]# 四个AIC值升序排列
    AIC.min=min(AICS)
    RL_AIC=exp(0.5*(min(AICS)-AICS[2]))# 最小值减次小值
    dAIC=AICS-AIC.min
    names(dAIC)=gsub("AIC_","dAIC_",names(AICS))
    index=-0.5*dAIC
    w=exp(index)/sum(exp(index))# 计算四个模型的权重
    names(w)=gsub("AIC_","weight_",names(AICS))
    model=gsub("weight_","",names(which.max(w)))
    r <- c(model=model,AICS,RL_AIC,dAIC,w,
           ME_mi=ci.test(lnc,icp,test = "mi-g")$statistic,# lnc与ICP
           ME_p=ci.test(lnc,icp,test = "mi-g")$p.value,
           MG_mi=ci.test(lnc,imm,test = "mi-g")$statistic,# lnc与Gene
           MG_p=ci.test(lnc,imm,test = "mi-g")$p.value,
           EG_mi=ci.test(icp,imm,test = "mi-g")$statistic,# ICP与Gene
           EG_p=ci.test(icp,imm,test = "mi-g")$p.value,
           CIT_ME_mi=ci.test(lnc,icp,imm,test = "mi-g")$statistic,#lnc,icp|imm
           CIT_ME_p=ci.test(lnc,icp,imm,test = "mi-g")$p.value,
           CIT_MG_mi=ci.test(lnc,imm,icp,test = "mi-g")$statistic,#lnc,imm|icp
           CIT_MG_p=ci.test(lnc,imm,icp,test = "mi-g")$p.value,
           CIT_EG_mi=ci.test(icp,imm,lnc,test = "mi-g")$statistic,#icp,imm|lnc
           CIT_EG_p=ci.test(icp,imm,lnc,test = "mi-g")$p.value)
    if(length(r) <26 ){# 长度小于26赋默认值
      r <- rep(0,26)
    }else{
      r <- c(model=model,AICS,RL_AIC,dAIC,w,
             ME_mi=ci.test(lnc,icp,test = "mi-g")$statistic,
             ME_p=ci.test(lnc,icp,test = "mi-g")$p.value,
             MG_mi=ci.test(lnc,imm,test = "mi-g")$statistic,
             MG_p=ci.test(lnc,imm,test = "mi-g")$p.value,
             EG_mi=ci.test(icp,imm,test = "mi-g")$statistic,
             EG_p=ci.test(icp,imm,test = "mi-g")$p.value,
             CIT_ME_mi=ci.test(lnc,icp,imm,test = "mi-g")$statistic,
             CIT_ME_p=ci.test(lnc,icp,imm,test = "mi-g")$p.value,
             CIT_MG_mi=ci.test(lnc,imm,icp,test = "mi-g")$statistic,
             CIT_MG_p=ci.test(lnc,imm,icp,test = "mi-g")$p.value,
             CIT_EG_mi=ci.test(icp,imm,lnc,test = "mi-g")$statistic,
             CIT_EG_p=ci.test(icp,imm,lnc,test = "mi-g")$p.value)
    }
    return(r)
  }
  # 模式选择 
  pb <- progress_bar$new(total=length(triplet$icp))
  system.time({
    for(n in 1:length(triplet$icp)){
      pb$tick()
      lnc1=lnc_exp[triplet[n,3],]
      icp1=ICP_exp[triplet[n,1],]
      imm1=IMM_exp[triplet[n,2],]
      lnc1=as.double(lnc1)
      icp1=as.double(icp1)
      imm1=as.double(imm1)
      aaa=modeSelection(lnc = lnc1,icp = icp1,imm = imm1)
      triplet[n,4:29]=aaa
    }
  })
  dir.create("./bys_result")
  setwd("./bys_result")
  write.table(triplet,"model_result.txt",quote = F,row.names = F,sep = "\t")
  # 十五、结果统计
  bys_names=c("model","AIC1","AIC2","AIC3","AIC4","RL_AIC","dAIC1","dAIC2",
              "dAIC3","dAIC4","w1","w2","w3","w4",
              "MI_lnc2icp","p_lnc2icp",
              "MI_lnc2imm","p_lnc2imm",
              "MI_icp2imm","p_icp2imm",
              "MI_lncicp|imm","p_lncicp|imm",
              "MI_lncimm|icp","p_lncimm|icp",
              "MI_icpimm|lnc","p_icpimm|lnc")
  # INDEP模式
  INDEP_triplet=triplet[which(triplet[,4] == "INDEP"),] # 
  a1=nrow(INDEP_triplet)
  colnames(INDEP_triplet)[4:29]=bys_names
  INDEP_triplet=INDEP_triplet[which(INDEP_triplet$`p_icpimm|lnc` > 0.05),] # 
  a2=nrow(INDEP_triplet)
  
  # LGI模式
  LGI_triplet=triplet[which(triplet[,4] == "LGI"),] #
  b1=nrow(LGI_triplet)
  colnames(LGI_triplet)[4:29]=bys_names
  LGI_triplet=LGI_triplet[which(LGI_triplet$`p_lncicp|imm` > 0.05),] # 
  b2=nrow(LGI_triplet)
  
  # LIG模式
  LIG_triplet=triplet[which(triplet[,4] == "LIG"),]# 
  c1=nrow(LIG_triplet)
  colnames(LIG_triplet)[4:29]=bys_names
  LIG_triplet=LIG_triplet[which(LIG_triplet$`p_lncimm|icp` > 0.05),] # 
  c2=nrow(LIG_triplet)
  
  # COO模式
  COO_triplet=triplet[which(triplet[,4] == "COO"),] # 
  d1=nrow(COO_triplet)
  colnames(COO_triplet)[4:29]=bys_names
  COO_triplet=COO_triplet[which(COO_triplet$`p_lncimm|icp` > 0.05),] # 
  d2=nrow(COO_triplet)
  
  static_result=data.frame("triple_number"=c(a1,b1,c1,d1),"number"=c(a2,b2,c2,d2))
  static_result[,"rate"]=static_result[,2]/static_result[,1]
  rownames(static_result)=c("INDEP","LGI","LIG","COO")
  write.table(static_result,"static_result.txt",sep="\t",quote = F)
  
  # 合格的基因对
  qualified_genepair=rbind(COO_triplet,INDEP_triplet)
  qualified_genepair=rbind(qualified_genepair,LGI_triplet)
  qualified_genepair=rbind(qualified_genepair,LIG_triplet) # 18402
  
  write.table(qualified_genepair,"qualified_genepair.txt",quote = F,row.names = F,sep = "\t")
  setwd("../")
  
  # 提取贝叶斯、随机扰动合格的三元组
  LMI_genepair=merge(qualified_genepair,LMI,all.y = T)
  LMI_genepair=drop_na(LMI_genepair)
  
  # 统计每个lnc-icp关系对对应几个imm基因
  LMI_gene=LMI_genepair[,c(2,1)]
  LMI_gene=LMI_gene[!duplicated(LMI_gene),]
  LMI_gene=merge(permutation_LI,LMI_gene,all.y = T)
  LMI_gene=drop_na(LMI_gene)
  
  dir.create("./contain_imm_lnc2icp")
  setwd("./contain_imm_lnc2icp")
  write.table(LMI_gene,"LMI_gene.txt",sep = "\t",row.names = F,quote = F)
  write.table(LMI_genepair,"LMI_genepair.txt",sep = "\t",row.names = F,quote = F)
  setwd(path)
}








