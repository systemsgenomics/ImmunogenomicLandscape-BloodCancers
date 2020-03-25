GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))

# tools:
library(Matrix)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

#****************************** plot Yang scores and immunoscores, singleR etc ******************************
geneset1=unique(data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/Yang_2019.txt", data.table = F, header = F))
geneset=unique(data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/Yang_2019_allgenes.txt", data.table = F, header = T)) # 41467_2019_11947_MOESM5_ESM

geneset=geneset[order(geneset$cluster, -geneset$avg_logFC),]
geneset=geneset[!geneset$avg_logFC<0,]
geneset=do.call(rbind, lapply(unique(geneset$cluster), function(clu)head(geneset[geneset$cluster%in%clu,], 10)))

add=c("CD7", "NKG7","KLRC1","NCR1", "NCAM1", "FCER1G","ZBTB16","SH2D1B", "SYK","IL32")

geneset.l=lapply(unique(geneset[,2]), function(g)geneset[geneset[,2]%in%g,1])
names(geneset.l)=unique(geneset[,2])
# add gm based score:
add.scores=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))
add.scores=append(geneset.l, add.scores)

# Integrated dataset:
load("FIMM_AML_HCA_Yang_NK_scRNA.Rdata")

# make sure that using logNormalized data and scaled to 10000:
scmat <- NormalizeData(scmat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
fm.f <- data.matrix(GetAssayData(scmat, assay = "RNA", slot="data"))

sample=gsub("_.*.", "", colnames(scmat))
sample[grepl("Manton", sample)]="HCA BM"
sample[grepl("GSM", sample)]="Yang BM"
sample[!grepl("3667|5249|5897", sample)&!grepl("HCA|Yang", sample)]="AML other"
sample[grepl("3667|5249|5897", sample)&!grepl("HCA|Yang", sample)]="MDS-like"
scmat[["sample"]]=sample

scmat[["sample.cluster"]]=paste(sample, Idents(scmat))

# take cluster assignment for AML cells:
scmat2=FindNeighbors(scmat[,grepl("MDS-like", sample)], dims = 1:25, reduction = "mnnCorrect")
scmat2=FindClusters(scmat2, resolution = 0.5, algorithm = 1)

AML_cluster_cell=Idents(scmat2)
save(AML_cluster_cell, file="Yang_HCA_AML_cluster_cell.Rdata")

pdf("FigureS3N_NK_FIMM_AML_reclustered.pdf", width = 6, height = 5)
DimPlot(scmat, pt.size = 0.25, cells.highlight = colnames(scmat2)[Idents(scmat2)%in%0], cols.highlight = "#006f00")
DimPlot(scmat, pt.size = 0.25, cells.highlight = colnames(scmat2)[Idents(scmat2)%in%1], cols.highlight = "#006f00")
DimPlot(scmat, pt.size = 0.25, cells.highlight = colnames(scmat2)[Idents(scmat2)%in%2], cols.highlight = "#006f00")
dev.off()

# add gm based score:
gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
  dat3=fm.f[rownames(fm.f)%in%add.scores[[i]],]
  gm=log2(t(apply(dat3, 2, gm_mean))) # done to normalized values
  rownames(gm)=names(add.scores)[i]
  return(gm)
}))

# also add to seurat object:
for(i in seq(add.scores)){
  scmat[[names(add.scores)[i]]] <- gm.objects[i,]
}

# pdf("Yang_markers_FIMM_AML_HCA_Yang_NKcells.pdf", width=10, height=10)
# DimPlot(scmat, label = T)
# DimPlot(scmat, group.by = "SingleR.label", label = T)
# DimPlot(scmat, group.by = "batch", label = T)
# FeaturePlot(scmat, features = names(add.scores), cols=c("grey75", "red"))
# FeaturePlot(scmat, features = add[add%in%rownames(scmat)], cols=c("grey75", "red"))
# VlnPlot(scmat, names(add.scores), pt.size = 0)
# VlnPlot(scmat, add, pt.size = 0)
# dev.off()

pdf("FigureS3N_Yang_markers_FIMM_AML_HCA_Yang_NKcells_Scores_scaled.pdf", width=16, height=10)
p1=FeaturePlot(scmat, features = names(add.scores), combine = F, min.cutoff = 0.1, max.cutoff = 1.5)
fix.sc <- scale_color_gradientn( colours = c('grey75', 'red'),  limits = c(0.1, 1.5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

pdf("Figure3J_Yang_markers_FIMM_AML_HCA_Yang_NKcells_Scores_blend.pdf", width=20, height=5)
FeaturePlot(scmat, features = names(add.scores)[c(1,5)], min.cutoff = 0.2, blend = T,cols = c("grey74", "blue", "red"))
FeaturePlot(scmat, features = names(add.scores)[c(4,5)], min.cutoff = 0.2, blend = T,cols = c("grey74", "blue", "red"))
FeaturePlot(scmat, features = names(add.scores)[c(4,10)], min.cutoff = 0.2, blend = T,cols = c("grey74", "blue", "red"))
dev.off()

# plot MDS-like samples to map:
scmat[["MDSlike"]]=ifelse(grepl("5897|3667|5249", scmat[["batch"]][,1]), "MDS-like", "other")
scmat[["MDSlike"]][grepl("BM", scmat[["batch"]][,1]),1]="normal BM"

lv = scmat[["MDSlike"]][,1]%in%c("MDS-like", "other")
# plot.scatter.seurat(scmat = scmat, colors.group = data.frame("V1"=factor(c("normal BM", "other", "MDS-like")), "V2"=c("grey75", "orange", "#21a366"), stringsAsFactors = F), seurat.feature = "MDSlike", name="NK cells (AML, HCA, Yang)",lv = lv, cores=1, add.density = T, add.proportions = T, SIZE = 1, rasterize=F, width = 64*4, height = 74*4, text.size=10*4)

xy=data.frame(Embeddings(scmat, "umap"))

pdf("Figure3J_Density_MDSlike_NKcell.pdf", width = 5, height = 4)
ggplot(xy[lv,], aes(x=UMAP_1, y=UMAP_2) ) +
  stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +
  viridis::scale_fill_viridis()
ggplot(xy[scmat[["MDSlike"]]=="normal BM",], aes(x=UMAP_1, y=UMAP_2) ) +
  stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +
  viridis::scale_fill_viridis()
dev.off()