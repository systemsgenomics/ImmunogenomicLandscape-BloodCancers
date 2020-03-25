##Before running this script create pdata.txt file that has 1st column CEL file names, then other sample source/type description columns

##The results are created using affymetrix annotations and alternive cdf from BrainArray: REFSEQ mapping used!
## alternative cdf files can be obtained from http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp
## installation of alternative probe mapping packages (requires AnnotationDBI Bioconductor package!) from shell: sudo R CMD INSTALL package_name
# you need 3 packages from brainarray, cdf probe db, download from link if you don´t have them

setwd("/research/groups/sysgen/raw_data/petri/GSE108474_RAW")
pdata="pdata.txt"

##---------------------------------------------------------

## set these variables to match the normalization you want
useALT_probemap=TRUE #Use alternative probe or not
useGCRMA=T #Use advanced version of GCRMA or not

##---------------------------------------------------------

normalize.hgu.133plus=function(pdata, outputpath=getwd(), useALT_probemap=T, useGCRMA=T){
  ## record time it takes to run
  ptm <- proc.time()

  ##REQUIRED packages

  require(affy)
  require(gcrma)
  require(affyQCReport)
  require(annaffy)

  require(hgu133plus2hsentrezgcdf)
  require(hgu133plus2hsentrezgprobe)
  require(hgu133plus2hsentrezg.db)

  require(hgu133plus2cdf)
  require(hgu133plus2.db)
  require(hgu133plus2probe)

  if(useGCRMA){
  ##---------------------------------------------------------
  ## PROBE AFFINITY CALCULATION

  ## for GC-RMA background correction we need to compute the probe affinity model
  ## check if the files already exist, if not, create them

  ## Using Affymetrix probe mapping
  if(!useALT_probemap){
    if(!file.exists("affinity_human133Plus2.RData")){
      affinity.info.human133plus2=compute.affinities("hgu133Plus2")
      save(affinity.info.human133plus2,file="affinity_human133Plus2.RData")
    }
    load("affinity_human133Plus2.RData")
  }

  ## Using alternative probe mapping (make sure that these libraries have been installed! cdf, probe and .db packages)
  if(useALT_probemap){
    if(!file.exists("affinity_human133Plus2alt.RData")){
      affinity.info.human133plus2alt=compute.affinities("hgu133Plus2_Hs_ENTREZG")
      save(affinity.info.human133plus2alt,file="affinity_human133Plus2alt.RData")
    }
    load("affinity_human133Plus2alt.RData")
  }

  ##---------------------------------------------------------
}
  ##read in pdata.txt and determine filenames to normalize based on it
  pheno=read.AnnotatedDataFrame(pdata)
  Filenames=rownames(pData(pheno))
  gsms=pheno$Accession
  rownames(pData(pheno))=Filenames
  experiment_ID=pheno$gse[1]

  ## AFFYMETRIX PROBE MAPPING
  ## GC-RMA normalization
  if(!useALT_probemap&useGCRMA){
    eset=just.gcrma(filenames=Filenames, phenoData=pheno, normalize=TRUE, type="fullmodel", fast=FALSE, cdfname="hgu133Plus2",affinity.info=affinity.info.human133plus2)
  }

  ## RMA normalization
  if(!useALT_probemap&!useGCRMA){
    eset=just.rma(filenames=Filenames, phenoData=pheno, normalize=TRUE, cdfname="hgu133Plus2")
  }

  ## ALTERNATIVE PROBE MAPPING
  ## GC-RMA normalization
  if(useALT_probemap&useGCRMA){
    eset=just.gcrma(filenames=Filenames, phenoData=pheno, normalize=TRUE, type="fullmodel", fast=FALSE, cdfname="hgu133Plus2_Hs_ENTREZG",affinity.info=affinity.info.human133plus2alt)
  }

  ## RMA normalization
  if(useALT_probemap&!useGCRMA){
    eset=just.rma(filenames=Filenames, phenoData=pheno, normalize=TRUE, cdfname="hgu133Plus2_Hs_ENTREZG")
  }

  ##---------------------------------------------------------
  ## PUT TOGETHER A TABLE WITH ANNOTATIONS

  ##annotation data
  probeIDs=featureNames(eset)

  ## find and remove ctrl probe data
  ctrls=grep("AFFX", probeIDs)
  discard=1:length(probeIDs)%in%ctrls
  probeIDs=probeIDs[!discard]
  # probeIDs=gsub("_at", "", probeIDs)
  eset=eset[!discard,]

  entrezIDs=aafLocusLink(probeIDs, "hgu133plus2hsentrezg.db")
  entrezIDv=as.vector(entrezIDs, mode="character")

  symbols=aafSymbol(probeIDs, "hgu133plus2.db")
  symbolv=as.vector(symbols, mode="character")

  library(org.Hs.eg.db)
  keys=entrezIDv

  symbolv <- mapIds(org.Hs.eg.db,
                    keys=keys,
                    column="SYMBOL",
                    keytype="ENTREZID",
                    multiVals="first")
  all(keys==names(symbolv))

  rm=is.na(symbolv)

  symbolv[is.na(symbolv)]=keys[is.na(symbolv)]

  ##put together in a dataframe
  DF=data.frame(exprs(eset),stringsAsFactors=FALSE)
  colnames(DF)=paste(gsub("_.*.", "", gsms), pheno$Title, sep="_")

  DF=DF[!(duplicated(symbolv)|symbolv=="integer(0)"),]
  rownames(DF)=symbolv[!(duplicated(symbolv)|symbolv=="integer(0)")]
  ##---------------------------------------------------------

  ## variance filter
  gvar=apply(exprs(eset),1,var)
  var_median=median(gvar)
  select=gvar>var_median
  png(paste(experiment_ID, "variance_histogram.png", sep="_"), bg="transparent")
  hist(gvar, main=paste("distribution of gene variance, median : ",var_median,sep=""), xlab="variance")
  dev.off()

  gexprs=apply(exprs(eset),1,mean) #mean exprs level for each gene
  exprs_median=median(gexprs)
  png(paste(experiment_ID, "exprs_histogram.png", sep="_"), bg="transparent")
  hist(gexprs, main=paste("distribution of gene expression, median : ",exprs_median,sep=""), xlab="variance")
  dev.off()

  png(paste(experiment_ID, "exprs_vs_variance.png", sep="_"), bg="transparent")
  plot(gexprs,gvar, main="exprs level vs variance", xlab="gexprs", ylab="variance")
  dev.off()

  ## use variance as filter
  DFfilt=DF[select,]
  ##---------------------------------------------------------
  ## WRITE OUTPUT FILES

  if(!useALT_probemap&useGCRMA){
    output=paste(experiment_ID,"affyID_gcrma_normalized.txt", sep="_")
  }
  if(!useALT_probemap&!useGCRMA){
    output=paste(experiment_ID,"affyID_rma_normalized.txt", sep="_")
  }
  if(useALT_probemap&useGCRMA){
    output=paste(experiment_ID,"symbol_gcrma_normalized.txt", sep="_")
  }
  if(useALT_probemap&!useGCRMA){
    output=paste(experiment_ID,"symbol_rma_normalized.txt", sep="_")
  }

  save(DF, file=file.path(gsub(".txt", ".Rdata", output)))

  write.table(DF,file=output, sep="\t", row.names=T, quote=FALSE)
  write.table(DFfilt, file=paste("filt", output,sep="_"), sep="\t", row.names=T, quote=FALSE)
  ptmf = proc.time() - ptm
  print(ptmf)
}

# run
normalize.hgu.133plus(pdata, useALT_probemap=useALT_probemap, useGCRMA=useGCRMA)
