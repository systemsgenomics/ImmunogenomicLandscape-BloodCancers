
is.binary=function(x) { all(na.omit(x) %in% 0:1) }
is.categorical=function(x) {data.table::uniqueN(na.omit(x))>=2&data.table::uniqueN(na.omit(x))!=length(na.omit(x))&class(x)=="character"}
# datatypes
datatypes=function(d){
  types=sapply(d, class)
  binary=sapply(d, is.binary)
  categorical=sapply(d, is.categorical)
  types[binary]="binary"
  types[categorical]="categorical"
  return(types)
}

make.binary.feats=function(df, datatype="CLIN", prefix=""){
  d=t(df)
  rownames(d)=paste("B:", datatype, ":", prefix, colnames(df), sep="")
  return(d)
}

make.numeric.feats=function(df, datatype="CLIN", prefix=""){
  d=t(df)
  rownames(d)=paste("N:", datatype, ":", prefix, colnames(df), sep="")
  return(d)
}
# make features from clinical
make.features=function(df, datatype="CLIN", prefix="", make.pairwise=T, max.level=50, cores=10){
  bin=NULL
  cat=NULL
  num=NULL
  if(!prefix=="")prefix=paste0(prefix ,"_")
  if(class(df)!="data.frame")stop("df needs to be data.frame")
  
  # clean names
  colnames(df)=FIX_NAME(colnames(df))
  
  check=datatypes(df)
  
  print(check)
  
  # make all binary feats:
  if(any(check=="binary")){
    bdat=df[,check=="binary", drop=F]
    lvl=apply(bdat, 2, function(d){length(unique(d))})
    bdat=bdat[,!lvl<2,drop=F]
    if(dim(bdat)[2]>0){
      bin=data.matrix(make.binary.feats(bdat,datatype, prefix))
    }
  }
  # make all categorical feats:
  if(any(check=="categorical")){
    cdat=df[,check=="categorical", drop=F]
    
    #compute number of levels
    lvl=apply(cdat, 2, function(d){length(unique(d))})
    cdat=cdat[,!(lvl>max.level|lvl==1),drop=F]
    
    if(dim(cdat)[2]>0){
      # cat=do.call(rbind, parallel::mclapply(seq(dim(cdat)[2]), function(i)FUN_MAKE_ALL(annovector = cdat[,i], prefix = colnames(cdat)[i], datatype = datatype, make.pairwise = make.pairwise), mc.cores=cores))
      cat=do.call(rbind, lapply(seq(dim(cdat)[2]), function(i)FUN_MAKE_ALL(annovector = cdat[,i], prefix = colnames(cdat)[i], datatype = datatype, make.pairwise = make.pairwise)))
      
      # as integer, counts
      cat <- apply (cat, c (1, 2), function (x) {
        (as.integer(x))
      })
      
    }
  }  
  # make all numeric feats:
  if(any(check=="numeric"|check=="integer")){
    num=data.matrix(make.numeric.feats(df[,check=="numeric"|check=="integer", drop=F],datatype, prefix))
  }  
  
  dfres=data.frame(rbind(bin, cat, num), stringsAsFactors = F, check.rows = F)
  
  colnames(dfres)=rownames(df)
  
  return(dfres)
}

FIX_NAME=function(name){
  pwname=gsub(" |\\.", "_", name)
  pwname=gsub("[^[:alnum:][:space:]:-]", "_", pwname)
  pwname=gsub("__", "_", pwname)
  pwname=gsub("___", "_", pwname)
  pwname=gsub("____", "_", pwname)
  pwname=gsub("_$", "", pwname)
  pwname=gsub("$_", "", pwname)
  pwname=gsub(" ", "_", pwname)
  return(pwname)
}

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


FUN_BINARYFEAT=function(type, annovector, prefix, datatype){
  logv=t(as.matrix(annovector%in%type))*1
  rownames(logv)=paste("B:", datatype, ":", prefix, type, sep="")
  return(logv)
}

FUN_MAKEPW=function(annovector){
  write.table(annovector, "pw_annovector2.txt", row.names=F, col.names=F, quote=F)
  system("python /research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/scripts/pairwise.py pw_annovector2.txt pw_annovector2_res.txt")
  pw_annovector2=read.delim("pw_annovector2_res.txt", stringsAsFactors = FALSE, header=F)
}

# makes vs. and comparisons
FUN_COMPARISONS=function(i,pw_annovector2, annovector, annovector2, PERCENTAGE, prefix, datatype){
  if(!prefix=="")prefix=paste0(prefix ,"_")
  
  row=pw_annovector2[i,]
  
  type_annovector2=table(annovector2[annovector%in%row[,2]])
  all_annovector2=table(annovector2)
  all_annovector2=all_annovector2[match(names(type_annovector2), names(all_annovector2))]
  
  test=signif(type_annovector2/all_annovector2,2)>=PERCENTAGE
  annovector2_to_correlate=names(type_annovector2)[test] 
  
  type_annovector2=table(annovector2[annovector%in%row[,3]])
  all_annovector2=table(annovector2)
  all_annovector2=all_annovector2[match(names(type_annovector2), names(all_annovector2))]
  
  test=signif(type_annovector2/all_annovector2,2)>=PERCENTAGE
  annovector2_to_correlate2=names(type_annovector2)[test]
  
  na_vector=rep("NA", length(annovector2))
  na_vector[annovector2%in%annovector2_to_correlate]=1
  na_vector[annovector2%in%annovector2_to_correlate2]=0
  
  logv=t(na_vector)
  name=paste(row[,2],"vs", row[,3], sep="_")
  rownames(logv)=paste("B:", datatype, ":", prefix, name, sep="")
  
  # double annovector2
  logv2=t(as.matrix(annovector2%in%annovector2_to_correlate|annovector2%in%annovector2_to_correlate2))*1
  name=paste(row[,2],row[,3], sep="_and_")
  rownames(logv2)=paste("B:", datatype, ":", prefix, name, sep="")
  
  if(length(unique(as.character(logv)))==3){
    return(rbind(logv, logv2))
  }else{
    NULL
  }
}

FUN_MAKE_ALL=function(annovector, prefix="", annovector2=NULL, PERCENTAGE=0, datatype="SAMP", make.pairwise=T){
  if(is.null(annovector2))annovector2=annovector
  if(!prefix=="")prefix=paste0(prefix ,"_")
  
  annovector[annovector=="NA"]=NA
  
  #checking:
  if(length(annovector)!=length(annovector2))stop(paste("annotation vector length differs from cluster length?"))
  
  # these are the types
  a=unique(annovector)
  a=a[!is.na(a)]
  
  # make binary features
  binaryfeats=do.call(rbind, lapply(a, FUN_BINARYFEAT, annovector, prefix, datatype))
  
  if(make.pairwise){
    # make pairwise comparisons
    pw=FUN_MAKEPW(a)  
    
    # make contrasting features out of pairwise
    comp=do.call(rbind, lapply(seq(dim(pw)[1]), FUN_COMPARISONS, pw, annovector, annovector2, PERCENTAGE, prefix, datatype))
  }else{
    comp=NULL
  }
  
  if(!is.null(comp)){
    data=rbind(data.matrix(binaryfeats), data.matrix(comp))
  }else{
    data=data.matrix(binaryfeats)
  }
  return(data)  
}

FUN_MAKE_ALL_LOGICAL=function(annovector, prefix, annovector2){
  
  #checking:
  if(length(annovector)!=length(annovector2))stop(paste("annotation vector length differs from cluster length?"))
  
  # these are the types
  a=unique(annovector)
  a=a[!is.na(a)]
  
  # make binary features
  binaryfeats=lapply(a, FUN_BINARYFEAT, annovector, prefix)
  names(binaryfeats)=a
  return(binaryfeats)  
}

FUN_MAKE_CATEGORICAL=function(annovector, prefix){
  # these are the types
  a=unique(annovector)
  a=a[!is.na(a)]
  
  # categorical, check if levels are too high
  if(length(a)<36){
    clust_fm=t(as.matrix(annovector))          
    rownames(clust_fm)=paste0("C:SAMP:", prefix)
  }else{
    clust_fm=NULL
  }
  return(clust_fm)
}
  
FIND_LOGICAL=function(name, vector){
  data=data.frame(t(grepl(name, vector, fixed = T)*1))
  rownames(data)=paste0("B:CLIN:cytogenetic_", name)
  return(data)
}  
  
make.gam=function(maf, cnv_arm_org, sv=NULL,del.genes, amp.genes, case.list, NAME="data_gam.Rdata"){
  
  # gistic data
  cnv_anno=cnv_arm_org[!grepl("CN values", cnv_arm_org[,1]),(1:9)]
  cnv_arm=cnv_arm_org[!grepl("CN values", cnv_arm_org[,1]),-(1:9)]
  cnv_arm2=cnv_arm_org[!grepl("CN values", cnv_arm_org[,1]),-(1:9)]
  
  filt1=maf$Variant_Classification%in%c("Silent", "Intron", "3'UTR", "3'Flank", "5'UTR", "5'Flank")
  maf=maf[!filt1,]
  
  maf.m=as.data.frame.matrix(table(maf$Hugo_Symbol, maf$Tumor_Sample_Barcode))
  maf.m=maf.m[,match(case.list, colnames(maf.m))]
  
  cnv_arm=cnv_arm[,match(case.list, colnames(cnv_arm))]
  
  n=toupper(cnv_anno$Descriptor)
  n[grepl("Amp",cnv_anno$`Unique Name`)]=paste0(gsub("\\(.*.", "", cnv_anno$`Wide Peak Limits`[grepl("Amp",cnv_anno$`Unique Name`)]), ":", n[grepl("Amp",cnv_anno$`Unique Name`)], ":AMP")
  n[grepl("Del",cnv_anno$`Unique Name`)]=paste0(gsub("\\(.*.", "", cnv_anno$`Wide Peak Limits`[grepl("Del",cnv_anno$`Unique Name`)]), ":", n[grepl("Del",cnv_anno$`Unique Name`)], ":DEL")
  n=gsub(":|-", ".", n)
  
  rownames(cnv_arm)=n
  
  # make CNV annotations, significant genes:
  del=apply(del.genes[-(1:4),-1], 2, function(v)unique(v[!(is.na(v)|v%in%"")]))
  amp=apply(amp.genes[-(1:4),-1], 2, function(v)unique(v[!(is.na(v)|v%in%"")]))
  
  del.name=gsub(":|-", ".", paste0(del.genes[4,-1], ":", toupper(del.genes[1,-1]), ":DEL"))
  amp.name=gsub(":|-", ".", paste0(amp.genes[4,-1], ":", toupper(amp.genes[1,-1]), ":AMP"))
  
  r.d=do.call(rbind, lapply(seq(del.name), function(i)cbind(del.name[i], unlist(del[i]))))
  r.a=do.call(rbind, lapply(seq(amp.name), function(i)cbind(amp.name[i], unlist(amp[i]))))
  
  cnv_annotations=rbind(r.d, r.a)
  
  # make genie in GAM format:
  AMP=t(cnv_arm[grepl("AMP", rownames(cnv_arm)),]==2)
  GAIN=t(cnv_arm[grepl("AMP", rownames(cnv_arm)),]==1)
  DEL=t(cnv_arm[grepl("DEL", rownames(cnv_arm)),]==2)
  LOSS=t(cnv_arm[grepl("DEL", rownames(cnv_arm)),]==1)
  MUT=t(maf.m>0)
  
  colnames(AMP)=paste0(colnames(AMP), ":AMP")
  colnames(DEL)=paste0(colnames(DEL), ":DEL")
  colnames(GAIN)=paste0(colnames(GAIN), ":GAIN")
  colnames(LOSS)=paste0(colnames(LOSS), ":LOSS")
  colnames(MUT)=paste0(colnames(MUT), ":MUT")
  AMP=AMP[,!colSums(AMP)==0]
  MUT=MUT[,!colSums(MUT)==0]
  DEL=DEL[,!colSums(DEL)==0]
  GAIN=GAIN[,!colSums(GAIN)==0]
  LOSS=LOSS[,!colSums(LOSS)==0]
  
  if(!is.null(sv)){
    SV=t(sv[,match(common, colnames(sv))])
    SV=SV[,!colSums(SV)==0]
    colnames(SV)=paste0(gsub("SV_", "", colnames(SV)), ":SV")
    MUT=cbind(MUT, SV)
  }
  
  gam=data.frame(MUT,DEL,AMP,GAIN,LOSS, check.names = F)
  
  save(list = c("gam", "cnv_annotations"), file=NAME)
}  

make.cnv.gam=function(cnv, cores=5){
  res=do.call(rbind, mclapply(rownames(cnv), function(g){
    n=cnv[g,,drop=F]
    n[n==2]="AMP"
    n[n==-2]="DEL"
    n[n==1]="GAIN"
    n[n==-1]="LOSS"
    n[n==0]="DIPLOID"
    
    go.through=unique(as.character(n))
    go.through=go.through[!go.through%in%c("NA", "DIPLOID")|is.na(go.through)]
    dat=do.call(rbind, lapply(go.through, function(name)n%in%name*1))
    rownames(dat)=paste(g, go.through, sep=":")
    colnames(dat)=colnames(cnv)
    return(dat)
  }, mc.cores=cores))
  return(res)
}

