source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")
# tools:
library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/visualisation/plotting_functions.R")
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/")

#********************************* merge HCA and FIMM AML, Yang and normalize together **************************************

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/FIMM_AML_TNK_scRNA.Rdata")
TNK_FIMM=scmat
DefaultAssay(TNK_FIMM)="RNA"
TNK_FIMM[["batch"]][,1][!grepl("3667|5897", TNK_FIMM[["batch"]][,1])]="AML other"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/HumanCellAtlas_NK_scRNA.Rdata")
TNK_HCA=scmat
DefaultAssay(TNK_HCA)="RNA"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/Yang_NK_BM_scRNA.Rdata")
TNK_Yang=scmat
DefaultAssay(TNK_Yang)="RNA"

scmat.merged=merge(x = TNK_FIMM, y = list(TNK_HCA, TNK_Yang))
scmat.merged[["Cluster"]]=c(paste("AML", as.character(Idents(TNK_FIMM))), paste("HCA", as.character(Idents(TNK_HCA))), paste("YangNK", as.character(Idents(TNK_Yang))))

scmat.merged=scmat.merged[,grepl("NK", scmat.merged[["SingleR.label"]][,1])]

name="FIMM_AML_HCA_Yang_NK"

batch=scmat.merged[["batch"]][,1]

batch=factor(batch, levels=c(grep("MantonBM", unique(batch), value=T), grep("_BM", unique(batch), value=T), "3667", "5897", "AML other"))

gexp=GetAssayData(scmat.merged, assay="RNA", slot="counts")
common=intersect(intersect(rownames(TNK_FIMM), rownames(TNK_HCA)), rownames(TNK_Yang))
gexp=gexp[rownames(gexp)%in%common, ]
test1=sc.data.analysis(scmat = gexp, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 25, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 2500, percent.mitoDNA = 10, resolution = 0.8)

#********************************** Same to T cells, combine with Szabo data **********************************

source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")
# tools:
library(Matrix)
library(Seurat)
source("/research/users/ppolonen/git_home/common_scripts/visualisation/plotting_functions.R")
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/")

#********************************* merge HCA and FIMM AML and normalize together **************************************
load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/FIMM_AML_TNK_scRNA.Rdata")
TNK_FIMM=scmat
DefaultAssay(TNK_FIMM)="RNA"
find=table(TNK_FIMM[["batch"]][,1])
find=names(find)[find<500]
TNK_FIMM[["batch"]][,1][TNK_FIMM[["batch"]][,1]%in%find]="AML other"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/HCA_scRNA/HumanCellAtlas_Tcells_scRNA.Rdata")
TNK_HCA=scmat
DefaultAssay(TNK_HCA)="RNA"

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/Szabo_Tcell_BM_scRNA.Rdata")
TNK_Yang=scmat
DefaultAssay(TNK_Yang)="RNA"

# scmat.merged=merge(x = TNK_FIMM, y = list(TNK_HCA, TNK_Yang))
scmat.merged=merge(x = TNK_FIMM, y = list(TNK_HCA))

scmat.merged=scmat.merged[,grepl("CD|Treg", scmat.merged[["SingleR.label"]][,1])]

name="FIMM_AML_HCA_T"

batch=scmat.merged[["batch"]][,1]

# sort based on their size
batch=factor(batch, levels=names(table(batch))[order(-table(batch), names(table(batch)))])

gexp=GetAssayData(scmat.merged, assay="RNA", slot="counts")
common=intersect(intersect(rownames(TNK_FIMM), rownames(TNK_HCA)), rownames(TNK_Yang))
gexp=gexp[rownames(gexp)%in%common, ]

test1=sc.data.analysis(scmat = gexp, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name=name, nr.pcs = 25, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 2500, percent.mitoDNA = 10, resolution = 0.5)
