# # HumanCellAtlas data
# file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA.h5ad"
# library(Seurat)
# hca <- ReadH5AD(file)
# Idents(hca) <- "louvain"
# save(hca, file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA.Rdata")

# file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HCA_raw.h5ad"
# 
# counts=t(data.table::fread("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/X_allgenes.csv.gz", data.table = F))
# rownames(counts)=scan("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/var_allgenes.txt","genes")
# colnames(counts)=scan("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/data/raw/obs_allgenes.txt","sample")
# save(counts, file="/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA_rawCounts.Rdata")

#******************************************* all cells *********************************************

library(Matrix)
library(Seurat)
library(dplyr)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

load("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA_rawCounts.Rdata")
load("/research/groups/biowhat_share/public_data/scRNAseq/Human_Cell_Atlas_Preview_Datasets/scanpy/hg19/new_2019/cellranger3/HumanCellAtlas_scRNA.Rdata")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/")

name="HumanCellAtlas"

clusters=Idents(hca)
hca[["rnaClusterID"]] <- clusters

#
library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

# try own processing:
batch=gsub("_HiSeq_.*.", "", colnames(scmat))
name=paste0(colnames(counts), "-", gsub("_HiSeq_.*.", "", colnames(counts)))

test1=sc.data.analysis(scmat = counts[rowSums(counts)>1000,name%in%colnames(scmat)], regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="HCA", nr.pcs = 50, check.pcs=F, plot.umap = T, cores = 10, percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000)


# sbatch
library(Matrix)
library(Seurat)
library(dplyr)
source("functions.scRNA.analysis.R")

load("HumanCellAtlas_scRNA_rawCounts.Rdata")
load("HumanCellAtlas_scRNA.Rdata")

name="HumanCellAtlas"

clusters=Idents(hca)
hca[["rnaClusterID"]] <- clusters

#
library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

# try own processing:
batch=gsub("_HiSeq_.*.", "", colnames(scmat))
name=paste0(colnames(counts), "-", gsub("_HiSeq_.*.", "", colnames(counts)))

test1=sc.data.analysis(scmat = counts[,name%in%colnames(scmat)], regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="HCA", nr.pcs = 50, check.pcs=F, plot.umap = T, cores = 10, percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000)


#******************************************* only TNK cells *********************************************
library(Matrix)
library(Seurat)
library(dplyr)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/")

name="HumanCellAtlas_TNK"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/myrun/HCA_scRNA.Rdata")
scmat=scmat[,grepl("Treg|NK|CD", scmat[["SingleR.label"]][,1])]

data=GetAssayData(scmat, assay = "RNA", slot = "counts")

library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

# try own processing:
batch=gsub("_HiSeq_.*.", "", colnames(scmat))

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 50, check.pcs=F, plot.umap = T, percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000, cores = 4, resolution = 1)
#************************************************************************************************************


#******************************************* only NK cells *********************************************
library(Matrix)
library(Seurat)
library(dplyr)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/")

name="HumanCellAtlas_NK"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/myrun/HCA_scRNA.Rdata")
scmat=scmat[,grepl("NK", scmat[["SingleR.label"]][,1])]

# try own processing:
batch=gsub("_HiSeq_.*.", "", colnames(scmat))

data=GetAssayData(scmat, assay = "RNA", slot = "counts")

library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 50, check.pcs=F, plot.umap = T,  percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000, cores = 4, resolution = 1)
#************************************************************************************************************


#******************************************* only T cells *********************************************
library(Matrix)
library(Seurat)
library(dplyr)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/")

name="HumanCellAtlas_Tcells"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/myrun/HCA_scRNA.Rdata")
scmat=scmat[,grepl("Treg|CD", scmat[["SingleR.label"]][,1])]

data=GetAssayData(scmat, assay = "RNA", slot = "counts")

library("HGNChelper")
update=checkGeneSymbols(rownames(data))
genelist1=update[,1]
genelist1[!is.na(update[,3])]=update[!is.na(update[,3]),3]
rownames(data)=genelist1

# try own processing:
batch=gsub("_HiSeq_.*.", "", colnames(scmat))

test1=sc.data.analysis(scmat = data, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 50, check.pcs=F, plot.umap = T, percent.mitoDNA = 10, nFeature.min = 0, nFeature.max = 6000, cores = 4, resolution = 1)
#************************************************************************************************************



# add names and RNA cluster
clusters=Idents(hca)
hca[["rnaClusterID"]] <- clusters

# take 50% of cells from each cluster, and preserve proportions
filt=unlist(lapply(unique(clusters), function(clust){
  set.seed(1)
  ind=which(clusters%in%clust)
  nr.samples=table(clusters)/2
  sample(colnames(hca)[ind], nr.samples[names(nr.samples)%in%clust])
}))

scmat=subset(hca, cells=filt)

# old:
library("SingleR")
cores=8
counts <- data.matrix(GetAssayData(scmat, assay = "RNA"))

# singler = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(blueprint_encode), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)
# singler2 = CreateSinglerObject(counts, annot = NULL, clusters=Idents(scmat), ref.list = list(hpca), project.name=name, fine.tune = T,do.signatures = T, numCores=cores)

# save(singler, file=paste0(name, "_singler_object_blueprint_encode.Rdata"))
# save(singler2, file=paste0(name, "_singler_object_hpca.Rdata"))

load(paste0(name, "_singler_object_blueprint_encode.Rdata"))
load(paste0(name, "_singler_object_hpca.Rdata"))

# add original identifiers
singler$meta.data$orig.ident = scmat@meta.data$orig.ident # the original identities, if not supplied in 'annot'

## if using Seurat v3.0 and over use:
singler$meta.data$xy = scmat@reductions$umap@cell.embeddings # the tSNE coordinates
singler$meta.data$clusters = scmat@active.ident # the Seurat clusters (if 'clusters' not provided)

# these names are not intuitive, replace them, mentioned that these are actually naive http://comphealth.ucsf.edu/SingleR/SupplementaryInformation2.html:
singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD4+ T-cells"]="CD4+ Tn"
singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%"CD8+ T-cells"]="CD8+ Tn"
singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD4+ T-cells"]="CD4+ Tn"
singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%"CD8+ T-cells"]="CD8+ Tn"

# change B-cells
singler$singler[[1]]$SingleR.single$labels[singler$singler[[1]]$SingleR.single$labels%in%c("B-cells","Class-switched memory B-cells", "naive B-cells", "Memory B-cells")]="B-cells"
singler$singler[[1]]$SingleR.clusters$labels[singler$singler[[1]]$SingleR.clusters$labels%in%c("B-cells","Class-switched memory B-cells", "naive B-cells", "Memory B-cells")]="B-cells"

scmat[["SingleR.label.main"]]=singler$singler[[1]]$SingleR.single.main$labels
scmat[["SingleR.label"]]=singler$singler[[1]]$SingleR.single$labels
scmat[["SingleR.label.cluster"]]=paste(singler$singler[[1]]$SingleR.single$labels, Idents(scmat))

# can be added as cluster id
cluster=singler$singler[[1]]$SingleR.clusters$labels
cluster.main=singler$singler[[1]]$SingleR.clusters.main$labels

# order based on cell type and add celltype to cluster:
lineage=c("HSC", "MPP","GMP","CLP","CMP","MEP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs", "Pre-B_cell", "Pro-B_cell", "naive B-cells", "Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes")
lineage=unique(c(lineage,blueprint_encode$types))
lineage=lineage[lineage%in%cluster]

identityVector.samples=as.character(Idents(scmat))
clusters.samples=Idents(scmat)

for(j in seq(cluster)){
  identityVector.samples[clusters.samples%in%rownames(cluster)[j]]=paste(cluster[j], rownames(cluster)[j])
}

# order and cluster identity;
cluster=cluster[order(match(cluster[,1],lineage)),,drop=F]
cat(paste(levels(Idents(scmat)), collapse=","), paste(cluster, collapse=","), sep="\t\t")   

# order seurat clusters too:
Idents(scmat)=factor(identityVector.samples, levels=paste(cluster, rownames(cluster)))

r4=DimPlot(object = scmat, group.by = "SingleR.label", reduction = 'umap' , label = TRUE,  pt.size = 0.5, repel = T)  + NoLegend() + NoAxes()
ggplot2::ggsave(r4, filename = paste0(name, "_UMAP_singleR.pdf"), width = 6, height = 6)


save(singler, file=paste0(name, "_singler_object.Rdata"))
cat("\nSingleR done", sep="\n\n")

# add immunoscores to FM format data matrix
fm.f <- data.matrix(GetAssayData(scmat, assay = "RNA"))

# add gm based score:
add.scores=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))

gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
  dat3=fm.f[rownames(fm.f)%in%add.scores[[i]],]
  gm=t(apply(dat3, 2, gm_mean)) # done to normalized values
  rownames(gm)=names(add.scores)[i]
  return(gm)
}))

# also add to seurat object:
for(i in seq(add.scores)){
  scmat[[names(add.scores)[i]]] <- gm.objects[i,]
}

# gene expression and scores to fm:
rownames(fm.f)=paste0("N:GEXP:", rownames(fm.f))
rownames(gm.objects)=paste0("N:SAMP:", rownames(gm.objects))
fm=rbind(gm.objects, fm.f)

# DE analysis:
# markers.all=FindAllMarkers(scmat, only.pos = T)#, ...)

save(list=c("fm", "scmat"), file=paste0(name, "_scRNA.Rdata"))
# save(list=c("fm", "scmat", "markers.all"), file=paste0(name, "_scRNA.Rdata"))


r5=DimPlot(object = scmat, reduction = 'umap', label=T)
ggplot2::ggsave(r5, filename = paste0(name, "_UMAP_clusters.pdf"))