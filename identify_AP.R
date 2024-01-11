

##A complete workflow to identify alternative promoter usage
library(proActiv)

#running time :start
t1=Sys.time()

##gtf file
gtf.file <- list.files("/data/xuefy/00aniATSS/00proActiv/DM/01reference",full.names = TRUE)

promoterAnnotation.gencode <- preparePromoterAnnotation(file = gtf.file,
                                                        species = 'Drosophila_melanogaster')
promoterAnnotation <- promoterAnnotation.gencode

##bam file input

files <- list.files("/data/xuefy/00aniATSS/00proActiv/DM/00bam",full.names = TRUE)


##proActiv estimates promoter counts and activity in a single command
result <- proActiv(files = files,
                   promoterAnnotation = promoterAnnotation,
                   genome = "dm6",
                   ncores = 6)

##result filter
result <- result[complete.cases(assays(result)$promoterCounts),]



##get absolute activity and usage
absolutePromoterActivity <- assays(result)$absolutePromoterActivity

relativePromoterActivity <- assays(result)$relativePromoterActivity
relativePromoterActivity[is.na(relativePromoterActivity)] <- 0


##compile metadata
promoterId <- result@elementMetadata@listData[["promoterId"]]
geneId <- result@elementMetadata@listData[["geneId"]]
seqnames <- result@elementMetadata@listData[["seqnames"]]
start <- result@elementMetadata@listData[["start"]]
strand <- result@elementMetadata@listData[["strand"]]
internalPromoter <- result@elementMetadata@listData[["internalPromoter"]]
promoterPosition <- result@elementMetadata@listData[["promoterPosition"]]


##format transfering 
promoterId <- as.data.frame(promoterId)
geneId <- as.data.frame(geneId)
seqnames <- as.data.frame(seqnames)
start <- as.data.frame(start)


##result
relActRESULT <- cbind(seqnames,geneId,promoterId,start,strand,internalPromoter,promoterPosition,relativePromoterActivity)
absActRESULT <- cbind(seqnames,geneId,promoterId,start,strand,internalPromoter,promoterPosition,absolutePromoterActivity)  

##output

write.csv(relActRESULT,file="relAct.csv",row.names = F)
write.csv(absActRESULT,file="absAct.csv",row.names = F)  

#runnnin time:end
t2=Sys.time()
print(t2)






