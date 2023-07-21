rm(list=ls()) ###清空工作空间
library(readr)
exp_seq <- read_tsv("exp_seq.OV-AU.tsv.gz", "\t", trim_ws = TRUE, guess_max = 50000)
exp_seq=as.data.frame(exp_seq)
colnames(exp_seq)=exp_seq[1,]
exp_seq=exp_seq[-1,]
#readr可以直接读取压缩包内的文件，无需解压。
output1 <- unique(exp_seq[,c("icgc_donor_id","icgc_specimen_id","gene_id","raw_read_count")])
remove(exp_seq)
gc()
library(data.table)
library(dplyr)
output2 <- dcast(data.table(output1),gene_id~icgc_specimen_id,
                 value.var="raw_read_count",fun.aggregate = max)
output2=as.data.frame(output2)
rownames(output2)=output2[,1]
output2=output2[,-1]
#由于有重复的基因ID，因此要对重复的基因进行合并，这里合并的方式是保留最大的值
name=colnames(output2)
name=name[order(name)]
specimen_file <- read_delim(file="specimen.OV-AU.tsv.gz",delim="\t",col_names = TRUE)
sepcimen_simplify_file <-   specimen_file %>% mutate(specimen_type=ifelse(specimen_type=="Primary tumour - solid tissue","C","N")) %>%
  dplyr::select(1,5,7)

sepcimen_simplify_file=as.data.frame(sepcimen_simplify_file)
sepcimen_simplify_file=sepcimen_simplify_file[which(sepcimen_simplify_file$icgc_specimen_id%in%name),]

output2=output2[,name]
sepcimen_simplify_file=sepcimen_simplify_file[order(sepcimen_simplify_file$icgc_specimen_id),]
colnames(output2)=paste(sepcimen_simplify_file$icgc_specimen_id,sepcimen_simplify_file$icgc_donor_id,sepcimen_simplify_file$specimen_type,sep = "_")

# library(AnnotationDbi)
# library(org.Hs.eg.db)
# symbol2ensembl <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(output2), columns = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")
# symbol2ensembl <- symbol2ensembl[!duplicated(symbol2ensembl$SYMBOL), ]
# symbol2ensembl <- symbol2ensembl[!duplicated(symbol2ensembl$ENSEMBL), ]
# symbol2ensembl=na.omit(symbol2ensembl)
# output2=output2[symbol2ensembl$SYMBOL,]
# rownames(output2) <- symbol2ensembl$ENSEMBL[match(rownames(output2), symbol2ensembl$SYMBOL)]

load("genelength.RData")
dat=output2
dat[,"gene_id"]=rownames(dat)
tmp <- inner_join(dat,gene,by = 'gene_id') # 匹配gene_count与gene_length
temple=tmp[,1:111]
temple%>% mutate(across(where(is.character), as.numeric))  -> temple
tmp=cbind(tmp[,112],temple,tmp[,113])
colnames(tmp)[1]="gene_id"
colnames(tmp)[113]="length"
fpkm <- data.frame(row.names = tmp$gene_id)
for (i in 2:(dim(tmp)[2]-1)){
  col <- tmp[[i]]
  N <- sum(col) # 计算每个样本的mapped reads数
  FPKMi <- (col*1e9)/(N*tmp[[dim(tmp)[2]]]) # 计算FPKM值
  FPKMi <- pmax(FPKMi,0) %>% as.data.frame() # 去掉矫正带来的负值
  colnames(FPKMi) <- colnames(tmp)[i]
  fpkm <- cbind(fpkm,FPKMi)
}

data=log2(fpkm+1)

name=colnames(data)
cancer=which(substr(name,start = 18,stop = 18)=="C")
normal=which(substr(name,start = 18,stop = 18)=="N")
name=name[c(normal,cancer)]
data=data[,name]
save(data,file="exp_OV-AU.RData")# 30个正常，81个癌症

remove(col,dat,exp_seq,fpkm,FPKMi,gene,i,N,name,output1,output2,sepcimen_simplify_file,specimen_file,tmp,temple)


