#Top gene extraction
  dat <- read.table(myFiles[i],header = T,stringsAsFactors = F)
  res <-data.frame(apply(dat,1,mean))  
  result[which(res[,j] >= quantile(res[,j])[3]),j] <- "top50%"
  result[which(res[,j] >= quantile(res[,j])[4]),j] <- "top25%"
  result[which(res[,j] < quantile(res[,j])[3]),j] <- "bottom"

#Consensus clustering

  out <- read.table(expFiles[i],header = T,stringsAsFactors = F,sep = "\t")
  out <- t(out)
  ddist = dist(out)
  ConsensusClusterObject <- ConsensusClusterPlus(ddist,                                                                
                                                 maxK = 6,                                                                              # set K
                                                 pItem = 0.8,
                                                 reps=5000,                                                                             # set repeats
                                                 clusterAlg = "hc",                                                                     # clustering Algorithm : Hierarchiocal clustering
                                                 innerLinkage = "ward.D2",                                                              # for color coding the clusters use tmyPal = ...
                                                 finalLinkage = "complete",
                                                 plot = 'pdf',                                                                          # write resut to pdf (Alt:png)
                                                 writeTable = TRUE,
                                                 verbose = TRUE)

#DE
    t_test<-t.test(exp_new[j,ICR_high],exp_new[j,ICR_low])$p.value
    ICR_low_ave<-mean(as.numeric(exp_new[j,ICR_low]))
    ICR_high_ave<-mean(as.numeric(exp_new[j,ICR_high]))
    log2foldchange<-ICR_high_ave - ICR_low_ave


#pcc
    PCC <- function(exp,pairss){
      result<-apply(pairss,1,function(x,exp){
        xx<-as.numeric(as.matrix(exp[x[1],]))
        yy<-as.numeric(as.matrix(exp[x[2],]))
        test_can<-cor.test(xx,yy,method="pearson")
        r<-c(x,cor=test_can$estimate,p=test_can[3][[1]])
        return(r)
      },exp)
      result<-as.data.frame(t(result))
    }
	
#survival
    res1 <- coxph(Surv(OS.time,OS) ~.,data=dat)
    r<-c(x,res1$coefficients[(length(res1$coefficients)-1) : length(res1$coefficients)])
    survival <- survdiff(Surv(surv[sett,2],surv[sett,1]) ~ x)

#SVM	
	tObj = tune.svm(scale(ssSKCM[,-c(1,2)]),factor(ssSKCM[,2]),
                probability = TRUE,
                cost=c(0.1,1,2,5,10,100,1000,10000,50000),
                gamma = 10^(-4:4),scale=F)


	
	
	
	
	
	
