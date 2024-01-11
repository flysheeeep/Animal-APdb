
library(data.table)
library(tidyr)

species<-'your Species'
rel_path = '/path/to/promoter_usage'

tissue_file<-read.table('/path/to/sample/tissue',
                        header = T,
                        sep = ',',
                        quote = '',
                        stringsAsFactors = F)

TF<-fread("/path/to/species/TF_listFromAnimalTFDB", 
          sep = "\t", 
          quote = '', 
          header = T, 
          stringsAsFactors = F,  
          data.table = F)

TF<-TF[c(3,2,4)]
gene_exp <- as.data.frame(fread('/path/to/gene/expression'))
row.names(gene_exp) <-  gene_exp[,1]
#gene_exp <- gene_exp[,-1]
gene_exp <- gene_exp[which(row.names(gene_exp) %in% TF$Ensembl),]
#colnames(TF)<-tolower(colnames(TF))
for (i in unique(tissue_file$tissue)) {
  
  # i="Embryo"
 
  
  tissue = tissue_file[which(tissue_file$tissue == i),]
  if (length(tissue$run) < 5) {
    print(paste0('small tissue:',unique(tissue$tissue),';',length(tissue$run)))
    next
  }else{
    tissue_tf = cbind(gene_exp[,1],gene_exp[tissue$run])
    colnames(tissue_tf)[1] = 'TF'
    tissue_tf = tissue_tf[which(rowMeans(tissue_tf[,2:ncol(tissue_tf)]) >= 1),]
    tissue_tf = tissue_tf[apply(tissue_tf[,2:ncol(tissue_tf)],1,var)!= 0,]
    #relative activity
    tissue_rel = as.data.frame(fread(paste0(rel_path,i,'_modi.txt')))
    tissue_rel$promoterId = paste0('promo_',tissue_rel$promoterId)
    tissue_rel = tissue_rel[,c(3,8:ncol(tissue_rel))]
    tissue_rel = tissue_rel[,-ncol(tissue_rel)]
    row.names(tissue_rel) = tissue_rel$promoterId
    tissue_rel = tissue_rel[,-1]
    tissue_rel = tissue_rel[tissue$run]
    
    rel2tf = data.frame(rownames(tissue_rel))
    colnames(rel2tf) = 'Pos'
    
    
    fun1<-function(tf,id,rel){
      # x=unlist(tf[1,])
      # y=id
      # z=rel
      id$TFID<-tf[1]
      tf<-as.numeric(tf[-1])
      if(length(which(tf>0))>4){
        result<-apply(rel,1,fun2,y=tf)
        result<-cbind(id,t(result))
        rownames(result)<-NULL
        colnames(result)<-c('rel','tf','Rho','FDR')
        return(result)
      }
    }
    fun2<-function(x,y){
      if(length(which(x>0))>4){
        result<-cor.test(x,y,method = 'spearm')
        result<-matrix(c(result$estimate,result$p.value),nrow = 1)
        return(result)
      }else{
        result<-matrix(c(NA,NA),nrow=1)
        return(result)
      }
    }
    result<-apply(tissue_tf,1,fun1,id=rel2tf,rel=tissue_rel)
    result<-do.call(rbind,result)
    rownames(result)<-NULL
    result<-na.omit(result)
    result$FDR<-p.adjust(result$FDR,method="fdr",nrow(result))
    # filename<-paste0('m3-TF/o_',gsub(' ','_',i))
    # write.table(result,filename,quote = F,sep = '\t',col.names = T,row.names = F)
    result[(abs(result$Rho)<0.3 |result$FDR>=0.05) & !is.na(result$Rho),3:4]<-NA
    result<-na.omit(result)
    filename<-paste0('/path/to/save',i,'_tf.txt')
    # write.table(result,filename,quote = F,sep = '\t',col.names = T,row.names = F)
    
  }
  
}

