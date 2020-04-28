GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "common_scripts/pathway_analysis/functions.GSEA.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

files=list.files(path = ".", "AML_subtypes.Rdata")
names(files)=gsub("_subtypes.Rdata", "", files)

wrapper.de.analysis=function(i, files){
  load(files[i])
  
  if(!is.null(dim(coordinates.subtype))){
    subtype=coordinates.subtype$subtype
  }else{
    subtype=coordinates.subtype
  }
  
  # make lv of the subtype:
  lv=get.logical(list(subtype))
  
  res=wrapper.wilcoxtest(rownames(gexp), data = gexp, logicalVectors = list("MDS-like"=lv$`MDS-like`), ALTERNATIVE = "greater", adj.method = "BH", CORES = 8, prettynum = F)
  
  res$Name=names(files)[i]
  res=res[res$FDR<0.05,]
  
}

res=lapply(seq(files), wrapper.de.analysis, files)

AML.res=do.call(rbind, res)

genes.signif=names(table(AML.res$Gene))[table(AML.res$Gene)>1] # at least 2 data sets support

# find for each subtype genes that are upregulated/downregulated:
find=AML.res$FDR<0.001&abs(AML.res$FC)>1&AML.res$Gene%in%genes.signif
AML.res=AML.res[find,]

AML.res=AML.res[order(AML.res$Name, AML.res$FDR, AML.res$FC, decreasing = F),]

write.table(AML.res, "TableS3_Significant_genes_MDSlike_bulk.txt", quote = F, row.names = F, col.names = T, sep="\t")