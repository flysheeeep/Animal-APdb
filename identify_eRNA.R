###eRNA-correlation
##package_needed
library(GenomicRanges)
library(dplyr)
##load path
#path <- as.character(args[1])
path <- '/path/to/save'


eRNA_path <- list.files(path = path, pattern = 'eRNA.csv',full.names = T) 'eRNA expresion file'
project_path <- list.files(path = path, pattern = 'proj.csv',full.names = T)

project <- read.csv(file = project_path,header = T)
#project <- read.csv(file = '/home/xuefy/00aniATSS/05eRNA/00proj/mouse_proj.csv')

#' @load_env
#load(loadpath)

bam <- project$run


#' @load_eRNA_atlas：all_eRNA_with_its_middle_expression
eRNA_atlas <- read.csv(eRNA_path)
tissue_type <- unique(eRNA_atlas$tissue)
eRNA_atlas$start <- eRNA_atlas$middle-3000
eRNA_atlas$end <- eRNA_atlas$middle+3000
eRNA_atlas <- eRNA_atlas[c("bio_project","tissue","erna_id","chr","start","end")]
eRNA_ranges <- GRanges(seqnames = eRNA_atlas$chr,
                       ranges = IRanges(start = eRNA_atlas$start,end = eRNA_atlas$end),
                       strand = '*',
                       erna_id=eRNA_atlas$erna_id)

#' @merge_function

multimerge<-function(dat=list()){
  if(length(dat)<2)return(as.data.frame(dat[[1]]))
  mergedat<-dat[[1]]
  dat[[1]]<-NULL
  for(i in dat){
    mergedat<-merge(mergedat,i,all = T)
  }
  return(mergedat)
}



for (i in tissue_type) {
  print(Sys.time())
  print(paste0('start:',i))
  temp_eRNA <- eRNA_ranges
  tissue_atss <- read.table(file = paste0('/path/to/tissue_info/file/',i,'_modi.txt'),
                            sep = '\t',header = T)
  
  atss_ranges <- GRanges(seqnames = tissue_atss$seqnames,
                         ranges = IRanges(start = tissue_atss$start - 1000000,end = tissue_atss$start + 100000),
                         strand = '*',
                         atss_id=tissue_atss$promoterId)
  
  tissue_run <- project[which(project$tissue==i),]$run
  tissue_proj <- unique(project[which(project$run %in% tissue_run),]$bio_project)
  
  proj_temp <- list()
  for (h in tissue_proj) {
    proj_eRNA <- read.csv(paste0(path,h,'/eRNASampleRPM'),sep = '\t')
    # if (ncol(proj_eRNA)==0) {
    #   print(paste0(h,'no eRNA data'))
    # }
    proj_run <- colnames(proj_eRNA)[(colnames(proj_eRNA) %in% tissue_run)]
    proj_eRNA <- proj_eRNA[c("eRNAID",proj_run)]
    proj_temp[[h]] <- proj_eRNA
    #assign(h,proj_eRNA)
  }
  
  
  tissue_eRNA <- multimerge(dat = proj_temp)
  tissue_eRNA[is.na(tissue_eRNA)] <- 0
  
  
  
  #tissue_eRNA: eRNA expression across samples
  temp_eRNA <- temp_eRNA[which(temp_eRNA$erna_id %in% tissue_eRNA$eRNAID)]
  overlap <- findOverlaps(temp_eRNA,atss_ranges)
  eRNA_overlap <- temp_eRNA$erna_id[queryHits(overlap)]
  atss_overlap <- atss_ranges$atss_id[subjectHits(overlap)]
  eRNA2atss <- data.frame(eRNA_id=eRNA_overlap,atss_id=atss_overlap)
  eRNA2atss <- eRNA2atss[!duplicated(eRNA2atss),]
  eRNA2atss <- group_by(eRNA2atss,atss_id)
  result <- eRNA2atss
  result$Rho <- NA
  result$p <- NA
  result$FDR <- NA
  for (j in 1:nrow(result)) {
    atss_id <- result$atss_id[j]
    eRNA_id <- result$eRNA_id[j]
    atss2exp <- t(tissue_atss[which(tissue_atss$promoterId==atss_id),c(8:(ncol(tissue_atss)-1))])
    eRNA2exp <- t(tissue_eRNA[which(tissue_eRNA$eRNAID == eRNA_id),c(2:ncol(tissue_eRNA))])
    if (length(which(atss2exp !=0)) < 5 | length(which(eRNA2exp !=0))<5 ) {
      print(paste0(i,' :0'))
      next
    }else{
      relation <- merge(atss2exp,eRNA2exp,by = 'row.names')
      relation <- relation[,c(2,3)]
      colnames(relation) <- c('atss','eRNA')
      test <- cor.test(x = relation$atss,y =relation$eRNA,alternative = "two.sided",method = "spearman",)
      result$Rho[j] <- test[["estimate"]][["rho"]]
      result$p[j] <- test[["p.value"]]
    }
  }
  for (k in unique(result$atss_id)) {
    result[result$atss_id==k,]$FDR <- p.adjust(result[result$atss_id==k,]$p,method = 'fdr', n = length(result[result$atss_id==k,]$p))
  }
  result <- result[which(result$FDR<0.05 & abs(result$Rho) >=0.3),]
  tissue_result <- result
  tissue_result$tissue <- i
  tissue_result$seqnames <- '.'
  tissue_result$geneId <- '.'
  tissue_result$start <- '.'
  for (l in 1:nrow(tissue_result)) {
    n <- which(tissue_atss$promoterId==tissue_result[l,]$atss_id)
    tissue_result$seqnames[l] <- tissue_atss$seqnames[n]
    tissue_result$geneId[l] <- tissue_atss$geneId[n]
    tissue_result$start[l] <- tissue_atss$start[n]
  }
  
  tissue_result <- tissue_result[c('tissue','atss_id','seqnames','geneId','start','eRNA_id','Rho','FDR')]
  write.table(tissue_result,file = paste0(path,'/erna/',i,'_erna.txt'),sep = '\t',
              col.names = T,row.names = F,quote = F)
  
  print(Sys.time())
}




