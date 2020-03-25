library(Matrix)
library(Seurat)
library(dplyr)
source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA")

name="CLL_D0"

# DATA:
m=data.matrix(Matrix::readMM("/research/groups/biowhat_share/public_data/scRNAseq/GSE111014_CLL/GSE111014_matrix.mtx"))
colnames(m)=scan("/research/groups/biowhat_share/public_data/scRNAseq/GSE111014_CLL/GSE111014_barcodes.tsv", "a")
genes=data.table::fread("/research/groups/biowhat_share/public_data/scRNAseq/GSE111014_CLL/GSE111014_genes.tsv", data.table = F, header = F)[,2]
rownames(m)=make.unique(genes)

# filter
m=m[,grepl("d0", colnames(m))]

# using the function:
batch=gsub(".*.-seq_", "", colnames(m))

test1=sc.data.analysis(scmat = m, regress.cell.label = batch, batch.correction.method = "MNNcorrect", name="CLL_D0", nr.pcs = 22, check.pcs=F, plot.umap = T, nFeature.min = 200, nFeature.max = 4000, percent.mitoDNA = 10)