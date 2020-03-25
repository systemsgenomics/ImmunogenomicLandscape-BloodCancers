# hm450 data downloaded from firehose:
get.hm450.data=function(genelist, hm450){
  
  genes=paste0("'", paste(genelist, collapse="\\|"), "'")
  temp=tempfile()
  
  system(paste("head -n 1", hm450, ">", temp))
  system(paste("grep", genes, hm450, ">>", temp))
  
  parsemeth=data.table::fread(temp, data.table=F, header = F)
  
  meth_annot=cbind(parsemeth[-1,1], parsemeth[-1,3], parsemeth[-1,4], parsemeth[-1,5])
  
  vals=seq(2, ncol(parsemeth)-1, 4)
  
  meth=data.matrix(parsemeth[2:dim(parsemeth)[1], vals])
  rownames(meth)=meth_annot[,1]
  colnames(meth)=substr(parsemeth[1,vals], 1, 15)
  
  # remove if half is NA:
  meth_annot=meth_annot[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]
  meth=meth[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]
  unlink(temp)
  ret=list(meth, meth_annot)
  names(ret)=c("meth", "meth_annot")
  return(ret)
}

rownameToFirstColumn <- function(DF,colname){
  DF <- as.data.frame(DF)
  DF[,colname] <- row.names(DF)
  DF <- DF[,c(dim(DF)[2],1:(dim(DF)[2]-1))]
  return(DF)
}

get450KProbeMapping <- function(probeIDs, platform='HM450', genome='hg19'){
  hm450 <- FDb.InfiniumMethylation.hg19::getPlatform(platform = platform, genome = genome)
  probes <- hm450[probeIDs]
  
  TSS = FDb.InfiniumMethylation.hg19::getNearestTSS(probes)
  TSS = rownameToFirstColumn(TSS,'methProbeIDs')
  TSS = dplyr::select(TSS, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  data.table::setnames(TSS,c("distance", "nearestGeneSymbol", "nearestTranscript"),
                       c("distanceToTSS", "nearestTSS", "nearestTSS.ID"))
  
  Tx = FDb.InfiniumMethylation.hg19::getNearestTranscript(probes)
  Tx = rownameToFirstColumn(Tx, 'methProbeIDs')
  Tx = dplyr::select(Tx, methProbeIDs, distance, nearestGeneSymbol, nearestTranscript)
  data.table::setnames(Tx,c("distance", "nearestGeneSymbol", "nearestTranscript"),
                       c("distanceToTx", "nearestTx", "nearestTx.ID"))
  
  Annotation = plyr::join_all(list(TSS,Tx), by = 'methProbeIDs', match = 'all')
  
  return(Annotation)
}