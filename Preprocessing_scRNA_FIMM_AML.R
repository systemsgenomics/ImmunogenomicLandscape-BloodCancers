# all samples together (this is what we used in the end)
library(Matrix)
library(Seurat)
library(DESeq2)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA")

# read data
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/FIMM_scRNA_AML/", "matrix.mtx", full.names = T, recursive = T)
files=files[grepl("filtered_gene_bc", files)&grepl("GRCh38", files)]
files=gsub("matrix.mtx", "", files)
ids=c("5249", "5750", "6187", "6333", "5238")
names(files)=ids

counts1=Read10X(data.dir = files)

# read data
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/FIMM_scRNA_AML/", "matrix.mtx", full.names = T, recursive = T)
files=files[grepl("filtered_gene_bc", files)&!grepl("GRCh38", files)]
files=gsub("matrix.mtx", "", files)
ids=c("3667","5897","6386")
names(files)=ids

counts2=Read10X(data.dir = files)

library("HGNChelper")
update=checkGeneSymbols(rownames(counts1))
update2=checkGeneSymbols(rownames(counts2))

genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(counts1)=genelist1
  
genelist2=update2[,1]
genelist2[!is.na(update2[,3])]=update2[!is.na(update2[,3]),3]
rownames(counts2)=genelist2

common=intersect(genelist1, genelist2)

scmat.combined=cbind(counts1[match(common, rownames(counts1)),], counts2[match(common, rownames(counts2)),])

name="FIMM_AML"

batch=gsub("_.*.", "", colnames(scmat.combined))

test1=sc.data.analysis(scmat = scmat.combined, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10)

