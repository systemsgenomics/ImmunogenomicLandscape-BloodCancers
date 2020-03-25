regulon.feats=function(fm, genelist, cnv_annot=NULL, filter.val=3, filtertypes=NULL){
  fmname=rownames(fm)
  
  if(!is.null(filtertypes)){
    fmname=fmname[!substr(fmname, 1,6)%in%filtertypes]
  }
  feature.gene=suppressWarnings(do.call(rbind, strsplit(fmname, ":"))[,3])
  
  if(any(grepl(":", genelist)))genelist[grepl(":", genelist)]=do.call(rbind, strsplit(genelist[grepl(":", genelist)], ":"))[,3]
  
  if(!is.null(cnv_annot)){
    l.genes=strsplit(cnv_annot[,2], ",")
    names(l.genes)=cnv_annot[,1]
    st=stack(l.genes)
    st=st[st[,1]%in%genelist,]
    st[,2]=as.character(st[,2])
  }
  
  nested.features=function(gene, fmname, feature.gene){
    direct=fmname[feature.gene%in%gene]
    if(!is.null(cnv_annot))direct=c(direct, st[st[,1]%in%gene,2])
    
    # order this to numeric, then to gexp-cnv-meth-gnab
    n=unique(direct)
    types=substr(n, 1, 6)
    ord=c("N:GEXP","B:GEXP","B:GNAB", "N:CNVR", "B:CNVR","N:METH", "N:MIRN", "N:SAMP", "B:SAMP", "N:CLIN", "B:CLIN", "N:DRUG")
    ord=ord[ord%in%types]
    
    r=unlist(lapply(ord, function(o) n[types%in%o]))
    
    return(r)
  }
  feats=lapply(genelist, nested.features, fmname, feature.gene)
  
  # Filter if not enough data:
  m=fm[rownames(fm)%in%unlist(feats),]
  s=is.na(m)
  m=m[!rowSums(s)==dim(m)[1]&rowSums(!s)>=filter.val, !colSums(s)==dim(m)[2]]
  
  feats=lapply(feats, function(n)n[n%in%rownames(m)])
  names(feats)=genelist
  
  feats=feats[!unlist(lapply(feats, length))<1]
  
  return(feats)
}

get.feat.within.region=function(INPUT, FEAT){
  
  temp=tempfile()
  temp2=tempfile()
  
  # write out bed format temp file:
  write.table(INPUT, temp, quote = F, row.names = F, col.names = F, sep="\t")
  
  # intersect with features
  system(paste0("bedtools intersect -a ",FEAT, " -b ",temp, " -wa -wb > ", temp2))  # genes in region
  
  # read in R
  d = read.delim(temp2, header=F, stringsAsFactors = F)
  
  d[,4]=gsub(":.*.", "", d[,4])
  
  d2=d[,c(10, 4)]
  
  unlink(temp)
  unlink(temp2)
  
  return(d2)
}

fun.correlation.regulon=function(i, l.regulon.gene, fm, filter.val=3, method.cor="spearman"){
  feats=l.regulon.gene[[i]]
  genename=names(l.regulon.gene)[i]
  x=data.matrix(t(as.matrix(fm[match(feats,rownames(fm)),,drop=F])))
  class(x) <- "numeric"
  
  if(dim(x)[2]<2)return(NULL) # this one did not have pairs to test
  
  # filter if binary has less than n samples
  if(any(grepl("^B:", colnames(x)))){
    check=colSums(x[,grepl("^B:", colnames(x)),drop=F], na.rm = T)
    x=x[,!colnames(x)%in%names(check[check<filter.val]),drop=F]
  }
  
  # filter if numeric has less than n non-zero samples
  if(any(grepl("^N:", colnames(x)))){
    check=colSums(x[,grepl("^N:", colnames(x)),drop=F]!=0, na.rm = T)
    x=x[,!colnames(x)%in%names(check[check<filter.val]),drop=F]
  }
  
  if(dim(x)[2]<2)return(NULL) # this one did not have pairs to test
  
  cor_rho=cor(x,use = "pairwise.complete.obs", method = method.cor)
  cor_pval=cor.test.p.m(x, method = method.cor, use = "pairwise.complete.obs")
  results=flat_cor_mat(cor_rho, cor_pval)
  
  # respect original order:
  results=do.call(rbind, lapply(colnames(x), function(n)results[results[,1]%in%n,]))
  
  rownames(results)=NULL
  
  # filter pairwise result:
  results=results[!duplicated(t(apply(results[,1:2], 1, sort))),]
  
  # add gene to cnvr:
  # results[grepl(":CNVR:", results[,1]),1]=gsub(":CNVR:", paste0(":CNVR:",genename, "@"), results[grepl(":CNVR:", results[,1]),1])
  # results[grepl(":CNVR:", results[,2]),2]=gsub(":CNVR:", paste0(":CNVR:",genename, "@"), results[grepl(":CNVR:", results[,2]),2])
  
  return(results)
}

fun.correlation.extrafeatures=function(i, l.regulon.gene, fm, extrafeatures, filter.val=3, method.cor="spearman"){
  feats=l.regulon.gene[[i]]
  genename=names(l.regulon.gene)[i]
  x=data.matrix(t(as.matrix(fm[match(feats,rownames(fm)),,drop=F])))
  class(x) <- "numeric"
  
  # filter if binary has less than n samples
  if(any(grepl("^B:", colnames(x)))){
    check=colSums(x[,grepl("^B:", colnames(x)),drop=F], na.rm = T)
    x=x[,!colnames(x)%in%names(check[check<filter.val]),drop=F]
  }
  
  # filter if numeric has less than n non-zero samples
  if(any(grepl("^N:", colnames(x)))){
    check=colSums(x[,grepl("^N:", colnames(x)),drop=F]!=0, na.rm = T)
    x=x[,!colnames(x)%in%names(check[check<filter.val]),drop=F]
  }
  
  # remove if in both lists
  extrafeatures=extrafeatures[!extrafeatures%in%feats]
  
  # Filter if not enough data (less than 3 not NA):
  m2=fm[rownames(fm)%in%extrafeatures,]
  s=is.na(m2)
  m2=m2[!rowSums(s)==dim(m2)[2]&rowSums(!s)>=filter.val, ]
  y=data.matrix(t(as.matrix(m2)))
  class(y) <- "numeric"
  
  y=cbind(x,y)
  
  # filter if binary has less than n 3 samples
  if(any(grepl("^B:", colnames(y)))){
    check=colSums(y[,grepl("^B:", colnames(y)),drop=F], na.rm = T)
    y=y[,!colnames(y)%in%names(check[check<filter.val])]
  }
  
  # filter if numeric has less than n non-zero samples
  if(any(grepl("^N:", colnames(y)))){
    check=colSums(y[,grepl("^N:", colnames(y)),drop=F]!=0, na.rm = T)
    y=y[,!colnames(y)%in%names(check[check<filter.val]),drop=F]
  }
  
  if(dim(x)[2]==0|dim(y)[2]==0)return(NULL) # no pairs to test
  
  cor_rho=cor(x,y,use = "pairwise.complete.obs", method = method.cor)
  
  # compute cor pval
  cor_pval=t(apply(x, 2, function(v1){
    apply(y, 2, function(v2){
      if(sum(complete.cases(cbind(v1, v2)))<filter.val)return(NA)
      cor.test.p(v1, v2, method = method.cor, use = "pairwise.complete.obs")
    })
  }))
  
  results=flat_cor_mat(data.matrix(cor_rho), data.matrix(cor_pval))
  
  results=do.call(rbind, lapply(colnames(x), function(n)results[results[,1]%in%n,]))
  
  rownames(results)=NULL
  
  # remove regulonfeatures, they are coming from another test:
  results=results[!(results[,1]%in%colnames(x)&results[,2]%in%colnames(x)),]
  
  # filter pairwise result:
  results=results[!duplicated(t(apply(results[,1:2], 1, sort))),]
  
  # add gene to cnvr: BUG for wilcox
  # results[grepl(":CNVR:", results[,1]),1]=gsub(":CNVR:", paste0(":CNVR:",genename, "@"), results[grepl(":CNVR:", results[,1]),1])
  
  return(results)
}


p.adjust.datatype=function(results, adjust.method="BH", log10=F, prettyNumbers=T, orderdata=T){
  
  # number of features tested:
  datatypes=suppressWarnings(do.call(rbind, strsplit(results[,2], ":"))[,2]) # sometimes more than 3 columns
  datatypes2=suppressWarnings(do.call(rbind, strsplit(results[,1], ":"))[,2])
  datatypes3=apply(cbind(datatypes2,datatypes),1,paste, collapse=":")
  datatypes=datatypes3
  
  results=results[order(datatypes),]
  datatypes=datatypes[order(datatypes)]
  
  vals=seq(datatypes)
  
  vals=unlist(lapply(unique(datatypes), function(d)p.adjust(results$p[datatypes%in%d], method = adjust.method)))
  
  if(log10){
    vals=abs(log10(vals))
  }
  if(prettyNumbers){
    vals=prettyNum(signif(vals,2))
    results[,3]=prettyNum(signif(results[,3],2))
    results[,4]=prettyNum(signif(results[,4],2))
  }else{
    vals=as.numeric(vals)
    results[,3]=as.numeric(results[,3])
    results[,4]=as.numeric(results[,4])
  }
  
  results$adj.p=vals
  results$signifCode=""
  results$signifCode[as.numeric(vals)<0.1]="."
  results$signifCode[as.numeric(vals)<0.05]="*"
  results$signifCode[as.numeric(vals)<0.01]="**"
  results$signifCode[as.numeric(vals)<0.001]="***"
  results$signifCode[as.numeric(vals)<0.0001]="****"
  
  results$datapairs=datatypes
  
  if(orderdata)results=results[order(results[,1], datatypes, -as.numeric(results$cor)>0, as.numeric(results$adj.p), decreasing = F),]
  
  return(results)
}


ft.alt <- function(a, b, c, d, alternative="greater") {
  as.numeric(fisher.test(matrix(c(a, b, c, d), 2), alternative = alternative)$p.value)
}

fisher.test.pwpw=function(results, find, fm, fisher.alternative, cores){
  results.o=results
  # remove gene anno if CNV type, to match with FM
  results.o[,1]=gsub("CNVR:.*.@", "CNVR:", results.o[,1])
  results.o[,2]=gsub("CNVR:.*.@", "CNVR:", results.o[,2])
  
  featsA=as.list(as.data.frame(t(fm[results.o[find,1],])))
  featsB=as.list(as.data.frame(t(fm[results.o[find,2],])))
  
  grid=data.frame(do.call(rbind, mclapply(seq(find), function(i)c(table(featsA[[i]], featsB[[i]])), mc.cores=cores)))
  colnames(grid)=c("a", "b", "c", "d")
  
  if(fisher.alternative=="two.sided"){
    p.ft <- with(grid, HighSpeedStats::ultrafastfet(a, b, c, d))
    results$p[find]=p.ft

  }else{
    p.ft <- with(grid, mapply(ft.alt, a, b, c, d, fisher.alternative))
    results$p[find]=p.ft
  }
  results$test.method[find]=paste0("Fisher´s exact (", fisher.alternative, ")")
  return(results)
}

# pairwise.correlation=function(l.regulon.gene, fm, extrafeatures=NULL, filter.val=3, cores=5, adjust.method="BH", method.cor="spearman", prettyNumbers=T, use.fisher=T, fisher.alternative="two.sided"){
#   
#   #****************************** Test within regulon ******************************
#   r=parallel::mclapply(seq(l.regulon.gene),fun.correlation.regulon,l.regulon.gene, fm, filter.val=filter.val,method.cor=method.cor, mc.cores=cores)
#   results=do.call(rbind, r)
#   results=results[!is.na(results$p),,drop=F]
#   
#   if(!is.null(results)){
#     results$test.method=paste(method.cor, "correlation")
#     results$test.group="gene.features"
#     
#     # if any pairs are binary-binary, add fisher test p-value
#     fa=substr(results[,1], 0, 1)=="B"
#     fb=substr(results[,2], 0, 1)=="B"
#     find=which(fa&fb)
#     findl=fa&fb
#     
#     if(use.fisher&any(findl)){
#       results=fisher.test.pwpw(results, find, fm, fisher.alternative, cores)
#       results1=p.adjust.datatype(results[findl,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#       results2=p.adjust.datatype(results[!findl,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#       results=rbind(results1, results2)
#     }else{
#       results=p.adjust.datatype(results, orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#     }
#   }
#   #*********************************************************************************
#   
#   # Test against extrafeatures if they exist:
#   if(!is.null(extrafeatures)){
#     r=parallel::mclapply(seq(l.regulon.gene),fun.correlation.extrafeatures, l.regulon.gene, fm, extrafeatures=extrafeatures, filter.val=filter.val,method.cor=method.cor, mc.cores=cores)
#     results.e=do.call(rbind, r)
#     results.e=results.e[!is.na(results.e$p),]
#     results.e$test.method=paste(method.cor, "correlation")
#     results.e$test.group="extra.features"
#     
#     # if any pairs are binary-binary, add fisher test p-value
#     fa=substr(results.e[,1], 0, 1)=="B"
#     fb=substr(results.e[,2], 0, 1)=="B"
#     find=which(fa&fb)
#     findl=fa&fb
#     
#     if(use.fisher&any(findl)){
#       results.e=fisher.test.pwpw(results.e, find,fm, fisher.alternative, cores)
#       results1=p.adjust.datatype(results.e[findl,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#       if(sum(!findl)>0){
#         results2=p.adjust.datatype(results.e[!findl,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#         results.e=rbind(results1, results2)
#       }else{
#         results.e=results1
#       }
#     }else{
#       results.e=p.adjust.datatype(results.e, orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
#     }
#     results=rbind(results, results.e)
#   }
#   
#   if(adjust.method=="BH")adj="FDR"
#   if(adjust.method=="bonferroni")adj="FWER"
#   colnames(results)[7]=adj
#   results[grepl("Fisher", results$test.method),3]=""
#   
#   # order columnts:
#   cols=c("featureA", "featureB", "cor", "p", adj, "signifCode", "test.method", "test.group", "datapairs")
#   results=results[,cols]
#   return(results)
# }

filter.pairwise.res=function(results, filter.table=NULL, FDR=0.1){
  if(is.null(filter.table)){
    filter.table=data.frame("pair"=unique(results$datapairs), "FDR"=rep(FDR, length(unique(results$datapairs))), stringsAsFactors = F)
    type1=gsub(":.*.", "", filter.table[,1])
    type2=gsub(".*.:", "", filter.table[,1])
    
    # usually fewer features
    filter.table$FDR[type1=="GNAB"|type2=="GNAB"]=0.25 # mutations can be interesting, high P-value kept
    filter.table$FDR[type1=="SAMP"&type2=="SAMP"]=0.1 # these can be very different
    filter.table$FDR[type1=="CLIN"&type2=="CLIN"]=0.1 # these can be very different
    
    # exclude, not usually meaningful
    filter.table$FDR[type1=="CNVR"&type2=="CNVR"]=0 # remove, highly correlated
    filter.table$FDR[type1=="METH"&type2=="METH"]=0 # remove, highly correlated
    
    # SAMP special rules:
    filter.table$FDR[type1=="SAMP"&type2=="GEXP"]=1e-30 
    filter.table$FDR[type1=="SAMP"&type2=="MIRN"]=1e-30 
    filter.table$FDR[type1=="SAMP"&type2=="METH"]=1e-8 
    
    print("default table used: ")
    print(filter.table)
  }
  
  pairs=results$datapairs
  p=unique(pairs)
  
  # make sure all pairs are found in filter.table
  if(!all(p%in%filter.table[,1]))stop(paste("Not all pairs were found: ", p[!p%in%filter.table[,1],]))
  
  filter.table=filter.table[match(p, filter.table[,1]),]
  
  # p-value per pairs:
  vals=do.call(rbind, lapply(seq(p), function(i){
    
    r=results[pairs%in%p[i],]
    
    filt=as.numeric(r[,5])<as.numeric(filter.table[i,2])
    
    r=r[filt,]
  }))
  
  return(vals)
}

pairwise.correlation=function(l.regulon.gene, fm, extrafeatures=NULL, filter.val=3, cores=5, adjust.method="BH", method.cor="spearman", prettyNumbers=T, use.fisher=T, use.wilcox=T, fisher.alternative="two.sided", wilcox.alternative="two.sided"){
  
  #****************************** Test within regulon ******************************
  r=parallel::mclapply(seq(l.regulon.gene),fun.correlation.regulon,l.regulon.gene, fm, filter.val=filter.val,method.cor=method.cor, mc.cores=cores)
  results=do.call(rbind, r)
  results=results[!is.na(results$p),,drop=F]
  
  if(!is.null(results)){
    results$test.method=paste(method.cor, "correlation")
    results$test.group="gene.features"
    
    # if any pairs are binary-binary, add fisher test p-value
    fa=substr(results[,1], 0, 1)=="B"
    fb=substr(results[,2], 0, 1)=="B"
    find=which(fa&fb)
    findl=fa&fb
    
    # if any pairs are binary-numeric, add wilcoxon test p-value
    fc=substr(results[,1], 0, 1)=="N"
    fd=substr(results[,2], 0, 1)=="B"
    fe=substr(results[,1], 0, 1)=="B"
    ff=substr(results[,2], 0, 1)=="N"
    
    findw=fc&fd|fe&ff
    
    # compute p-value depending on test type:
    if(use.fisher&any(findl)){
      # update B vs B tests:
      results=fisher.test.pwpw(results, find, fm, fisher.alternative, cores)
    }
    
    # compute p-value depending on test type:
    if(use.wilcox&any(findw)) {
      # update B vs N tests:
      results=test.feats.wilcox(features = results, rows = findw, data = fm, wilcox.alternative = wilcox.alternative, cores = cores)
    }
    
    # go through each test and adjust p-value per data type:
    results=do.call(rbind, lapply(unique(results$test.method), function(type){
      p.adjust.datatype(results[results$test.method%in%type,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
    }))
  }
  #*********************************************************************************
  
  # Test against extrafeatures if they exist:
  if(!is.null(extrafeatures)){
    r=parallel::mclapply(seq(l.regulon.gene),fun.correlation.extrafeatures, l.regulon.gene, fm, extrafeatures=extrafeatures, filter.val=filter.val,method.cor=method.cor, mc.cores=cores)
    results.e=do.call(rbind, r)
    results.e=results.e[!is.na(results.e$p),]
    results.e$test.method=paste(method.cor, "correlation")
    results.e$test.group="extra.features"
    
    # if any pairs are binary-binary, add fisher test p-value
    fa=substr(results.e[,1], 0, 1)=="B"
    fb=substr(results.e[,2], 0, 1)=="B"
    find=which(fa&fb)
    findl=fa&fb
    
    # if any pairs are binary-numeric, add wilcoxon test p-value
    fc=substr(results.e[,1], 0, 1)=="N"
    fd=substr(results.e[,2], 0, 1)=="B"
    fe=substr(results.e[,1], 0, 1)=="B"
    ff=substr(results.e[,2], 0, 1)=="N"
    findw=fc&fd|fe&ff
  
    # compute p-value depending on test type:
    if(use.fisher&any(findl)){
    # update B vs B tests:
    results.e=fisher.test.pwpw(results.e, find, fm, fisher.alternative, cores)
    }
    
    # compute p-value depending on test type:
    if(use.wilcox&any(findw)) {
    # update B vs N tests:
    results.e=test.feats.wilcox(features = results.e, rows = findw, data = fm, wilcox.alternative = wilcox.alternative, cores = cores)
    }
    
    # go through each test and adjust p-value per data type
    results.e=do.call(rbind, lapply(unique(results.e$test.method), function(type){
      p.adjust.datatype(results.e[results.e$test.method%in%type,], orderdata=T, adjust.method=adjust.method, prettyNumbers=prettyNumbers)
    }))
    
    results=rbind(results, results.e)
  }
  
  if(adjust.method=="BH")adj="FDR"
  if(adjust.method=="bonferroni")adj="FWER"
  colnames(results)[7]=adj
  results[grepl("Fisher", results$test.method),3]=""
  
  # order columnts:
  cols=c("featureA", "featureB", "cor", "p", adj, "signifCode", "test.method", "test.group", "datapairs")
  results=results[,cols]
  return(results)
}

filter.pairwise.res=function(results, filter.table=NULL, FDR=0.1){
  if(is.null(filter.table)){
    filter.table=data.frame("pair"=unique(results$datapairs), "FDR"=rep(FDR, length(unique(results$datapairs))), stringsAsFactors = F)
    type1=gsub(":.*.", "", filter.table[,1])
    type2=gsub(".*.:", "", filter.table[,1])
    
    # usually fewer features
    filter.table$FDR[type1=="GNAB"|type2=="GNAB"]=0.25 # mutations can be interesting, high P-value kept
    filter.table$FDR[type1=="SAMP"&type2=="SAMP"]=0.1 # these can be very different
    filter.table$FDR[type1=="CLIN"&type2=="CLIN"]=0.1 # these can be very different
    
    # exclude, not usually meaningful
    filter.table$FDR[type1=="CNVR"&type2=="CNVR"]=0 # remove, highly correlated
    filter.table$FDR[type1=="METH"&type2=="METH"]=0 # remove, highly correlated
    
    # SAMP special rules:
    filter.table$FDR[type1=="SAMP"&type2=="GEXP"]=1e-30
    filter.table$FDR[type1=="SAMP"&type2=="MIRN"]=1e-30
    filter.table$FDR[type1=="SAMP"&type2=="METH"]=1e-8
    
    print("default table used: ")
    print(filter.table)
  }
  
  pairs=results$datapairs
  p=unique(pairs)
  
  # make sure all pairs are found in filter.table
  if(!all(p%in%filter.table[,1]))stop(paste("Not all pairs were found: ", p[!p%in%filter.table[,1],]))
  
  filter.table=filter.table[match(p, filter.table[,1]),]
  
  # p-value per pairs:
  vals=do.call(rbind, lapply(seq(p), function(i){
    
    r=results[pairs%in%p[i],]
    
    filt=as.numeric(r[,5])<as.numeric(filter.table[i,2])&!is.na(r[,5])
    
    r=r[filt,]
  }))
  
  return(vals)
}

test.feats.wilcox <- function(features, rows=NA, data, wilcox.alternative="two.sided", cores=1) {
  if(length(rows)>=1) {
    pairs <- features[rows,]
  } else pairs <- features
  
  p=unlist(parallel::mclapply(1:nrow(pairs), function(i){
    binary=substr(pairs[i,2], 0, 1)=="B"
    
    if(binary) {
      
      x <- as.numeric(data[pairs[i,1],])
      y <- as.logical(as.numeric(data[pairs[i,2],]))
      
    }else {
      
      x <- as.numeric(data[pairs[i,2],])
      y <- as.logical(as.numeric(data[pairs[i,1],]))
      
    }
    
    a <- vector()
    b <- vector()
    
    for(j in 1:length(y)) {
      
      if(is.na(y[j])) next
      if(y[j]) a <- c(a, x[j])
      if(!y[j]) b <- c(b, x[j])
      
    }
    
    if((is.numeric(a) && is.numeric(b))&&(length(a)>0&&length(b>0))&&sum(is.na(a))!=length(a)&&sum(is.na(b))!=length(b)) {
      p <- wilcox.test(a,b, alternative = wilcox.alternative)$p.value
    } else {
      p=NA # not enough data
    }
    return(p)
  }, mc.cores=cores))
  
  pairs$test.method=paste0("Wilcoxon test (", wilcox.alternative, ")")
  pairs$p<- as.numeric(as.character(p))
  
  # return original table but with updated p and test.method
  features[rows,]=pairs
  
  return(features)
}
