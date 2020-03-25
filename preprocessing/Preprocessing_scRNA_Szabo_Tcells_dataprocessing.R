library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/")

# read data Resting cells from bone marrow of Donor 2 after negative selection to enrich CD3+ cell
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Szabo_data/", "filtered.matrix.txt", full.names = T, recursive = T)
files=files[grepl("GSM3589408|GSM3589409|GSM3589414|GSM3589415", files)] # GSM3589408 GSM3589414 no stim, GSM3589409 and GSM3589415 CD3/CD28 T Cell Activator, 

f1=data.table::fread(files[1],data.table=F)
f2=data.table::fread(files[2],data.table=F)
f3=data.table::fread(files[3],data.table=F)
f4=data.table::fread(files[4],data.table=F)

genes=unique(f1[,2])

data=cbind(f1[match(genes, f1[,2]),3:(dim(f1)[2]-1)], f2[match(genes, f2[,2]),3:(dim(f2)[2]-1)], f3[match(genes, f3[,2]),3:(dim(f3)[2]-1)], f4[match(genes, f4[,2]),3:(dim(f4)[2]-1)])
rownames(data)=genes

library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]

data=data[!duplicated(genelist1),]
rownames(data)=unique(genelist1)

name="Szabo_Tcell_BM"

batch=c(rep("GSM3589408", length(3:(dim(f1)[2]-1))), rep("GSM3589409", length(3:(dim(f2)[2]-1))), rep("GSM3589414", length(3:(dim(f3)[2]-1))), rep("GSM3589415", length(3:(dim(f4)[2]-1))))

annot=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Szabo_data/Umap_coords_Tcells.txt", data.table = F, dec=",")

batch=batch[colnames(data)%in%annot$barcode]
data=data[,colnames(data)%in%annot$barcode]

data=data[!rowSums(data)<0.001*dim(data)[2],]

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, name=name, batch.correction.method = "MNNcorrect", nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4500, percent.mitoDNA = 10, resolution = 0.8, cores = 6)




library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/")

# read data Resting cells from bone marrow of Donor 2 after negative selection to enrich CD3+ cell
files=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Szabo_data/", "filtered.matrix.txt.gz", full.names = T, recursive = T)
files=files[grepl("GSM3589418|GSM3589419|GSM3589420|GSM3589421", files)] # GSM3589421 GSM3589419 no stim, GSM3589418 and GSM3589420 CD3/CD28 T Cell Activator, 

f1=data.table::fread(files[1],data.table=F)
f2=data.table::fread(files[2],data.table=F)
f3=data.table::fread(files[3],data.table=F)
f4=data.table::fread(files[4],data.table=F)

data=cbind(f1[,3:(dim(f1)[2]-1)], f2[,3:(dim(f2)[2]-1)], f3[,3:(dim(f3)[2]-1)], f4[,3:(dim(f4)[2]-1)])
data=data[!duplicated(f1[,2]),]

rownames(data)=f1[!duplicated(f1[,2]),2]

name="Szabo_Tcell_PB"

batch=c(rep("GSM3589418", length(3:(dim(f1)[2]-1))), rep("GSM3589419", length(3:(dim(f2)[2]-1))), rep("GSM3589420", length(3:(dim(f3)[2]-1))), rep("GSM3589421", length(3:(dim(f4)[2]-1))))

data=data[!rowSums(data)<0.02*dim(data)[2],]

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, name=name, batch.correction.method = "MNNcorrect", nr.pcs = 50, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4500, percent.mitoDNA = 10, resolution = 0.8, cores = 6)
