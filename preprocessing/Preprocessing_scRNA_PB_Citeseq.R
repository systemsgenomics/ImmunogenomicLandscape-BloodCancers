library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA")

name="CD8_Citeseq"

options(future.globals.maxSize= 4991289600)

# files=list.files(".", "matrix.mtx", full.names = T, recursive = T)
files=list.files("/research/groups/sysgen/PROJECTS/students/citeseq/", "matrix.mtx", full.names = T, recursive = T)
files=gsub("matrix.mtx.gz", "", files)
ids=c("cd8_donor1", "cd8_donor2", "cd8_donor3", "cd8_donor4")
names(files)=ids

# Each patient:
# for(i in seq(ids)[c(1,2,3,4)]){
#   datas=Read10X(data.dir = files[i])
#   rownames(datas[[2]])=gsub("_TotalSeqC","", rownames(datas[[2]]))
#   names(datas)=c("RNA", "PROT")
#   test3=sc.data.analysis(scmat = datas, name = ids[i], nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10, cores = 6)
# }

data=Read10X(data.dir = files[c(1,2,3,4)])
rownames(data[[2]])=gsub("_TotalSeqC","", rownames(data[[2]]))
names(data)=c("RNA", "PROT")

# take data with enough both protein counts:
filt=colSums(data[[2]])>2647&colSums(data[[2]])<6664
filt.genes=rowSums(data[[1]]>0)>20

data[[1]]=data[[1]][filt.genes,filt]

data[[2]]=data[[2]][,filt]

batch=substr(colnames(data[[1]]), 1,10)

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="CD8_Citeseq", nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10, cores = 2)



# sbatch
library(Matrix)
library(Seurat)
source("functions.scRNA.analysis.R")

future::plan("multiprocess", workers = 10)
options(future.globals.maxSize= 8991289600)

name="CD8_Citeseq"

options(future.globals.maxSize= 4991289600)

# files=list.files(".", "matrix.mtx", full.names = T, recursive = T)
files=list.files(".", "matrix.mtx", full.names = T, recursive = T)
files=gsub("matrix.mtx.gz", "", files)
ids=c("cd8_donor1", "cd8_donor2", "cd8_donor3", "cd8_donor4")
names(files)=ids

# Each patient:
# for(i in seq(ids)[c(1,2,3,4)]){
#   datas=Read10X(data.dir = files[i])
#   rownames(datas[[2]])=gsub("_TotalSeqC","", rownames(datas[[2]]))
#   names(datas)=c("RNA", "PROT")
#   test3=sc.data.analysis(scmat = datas, name = ids[i], nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10, cores = 6)
# }

data=Read10X(data.dir = files[c(1,2,3,4)])
rownames(data[[2]])=gsub("_TotalSeqC","", rownames(data[[2]]))
names(data)=c("RNA", "PROT")

# take data with enough protein counts:
filt=colSums(data[[2]])>2647&colSums(data[[2]])<6664
filt.genes=rowSums(data[[1]]>0)>20

data[[1]]=data[[1]][filt.genes,filt]

data[[2]]=data[[2]][,filt]

batch=substr(colnames(data[[1]]), 1,10)

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="CD8_Citeseq", nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10, cores = 2)
