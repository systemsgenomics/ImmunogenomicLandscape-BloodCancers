library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/")

# read data
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Yang_NK_scRNA", "matrix.mtx", full.names = T, recursive = T)
files=files[grepl("BM", files)]
files=dirname(files)
names(files)=basename(files)

data=Read10X(data.dir = files)

library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

name="Yang_NK_BM"

batch=sapply(strsplit(colnames(data), "_"), function(l)paste(l[2:(length(l)-1)], collapse="_"))

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 2500, percent.mitoDNA = 5, resolution = 0.8)

# read data
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Yang_NK_scRNA", "matrix.mtx", full.names = T, recursive = T)
files=files[grepl("PB", files)]
files=dirname(files)
names(files)=basename(files)

data=Read10X(data.dir = files)

name="Yang_NK_PB"

batch=sapply(strsplit(colnames(data), "_"), function(l)paste(l[2:(length(l)-1)], collapse="_"))

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 25, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 2500, percent.mitoDNA = 5, resolution = 1.5)