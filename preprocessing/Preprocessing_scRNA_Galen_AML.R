library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA")

name="AML_Galen"

# read data
files=list.files("/research/groups/biowhat_share/public_data/scRNAseq/GSE116256_Galen_AML/","dem.txt", full.names = T)

files=files[grepl("D0", files)]
files=files
data=do.call(cbind, parallel::mclapply(files, read.delim, row.names=1, mc.cores=8))

files_anno=list.files("/research/groups/biowhat_share/public_data/scRNAseq/GSE116256_Galen_AML/","_AML.*.anno.txt", full.names = T)

anno.list=lapply(files_anno[grepl("D0", files_anno)], function(f){
  mat <- read.delim(file = f, stringsAsFactors = F)
})

batch=gsub("_.*.", "", colnames(data))

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 30, check.pcs=F, plot.umap = T, nFeature.min = 100, nFeature.max = 3000, percent.mitoDNA = 10)

# Each patient:
ids=gsub(".dem.txt|/research/groups/biowhat_share/public_data/scRNAseq/GSE116256_Galen_AML//GSM35879.._","", files)

for(i in seq(files)){
  
  datas=read.delim(files[i], row.names=1)
  
  if(dim(datas)[2]>1000)test3=sc.data.analysis(scmat = datas, name = ids[i], nr.pcs = 15, check.pcs=F, plot.umap = T, nFeature.min = 100, nFeature.max = 3000, percent.mitoDNA = 10)
  
}
