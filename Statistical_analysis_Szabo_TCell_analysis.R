GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))

library(Matrix)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

#****************************** plot Szabo scores and immunoscores, singleR etc ******************************

a=fread("Szabo_Tcell_markers.txt", data.table = F)
rownames(a)=make.unique(a[,1])
dat=data.matrix(a[,-1])
dat[is.na(dat)]=0
# plot markers in

# compare NK cells
geneset=data.frame("genes"=a[,1], "name"=NA, stringsAsFactors = F)
geneset$name[1:70]="Treg"
geneset$name[71:140]="CD4 Tn / Tcm resting"
geneset$name[141:210]="CD4 / CD8 resting"
geneset$name[211:280]="IFN response"
geneset$name[281:350]="Proliferation"
geneset$name[351:420]="CD8 Cytotoxic"
geneset$name[421:490]="CD8 Cytokine"

geneset2=data.frame("genes"=a[,1], "name"=NA, stringsAsFactors = F)
geneset2$name[1:70]="Resting"
geneset2$name[71:140]="Resting"
geneset2$name[141:210]="Resting"
geneset2$name[211:280]="Activated"
geneset2$name[281:350]="Activated"
geneset2$name[351:420]="Activated"
geneset2$name[421:490]="Activated"

geneset.l=lapply(unique(geneset[,2]), function(g)geneset[geneset[,2]%in%g,1])
names(geneset.l)=unique(geneset[,2])

geneset.l2=lapply(unique(geneset2[,2]), function(g)geneset2[geneset2[,2]%in%g,1])
names(geneset.l2)=unique(geneset2[,2])

geneset.exhaustion=c("PDCD1", "CTLA4", "LAYN", "LAG3", "HAVCR2", "CD244", "CD160")
geneset.cytokine=c("CCL3", "CCL4", "XCL1", "XCL2", "IFNG", "CSF2", "IL10","TNFRSF9")
geneset.cytotoxic=c("CCL5", "GZMK", "GNLY", "EOMES", "ZEB2", "ZNF683", "KLRG1", "NKG7")

geneset.l3=list("Exhaustion"=geneset.exhaustion, "Cytokine"=geneset.cytokine, "Cytotoxic"=geneset.cytotoxic)

# add gm based score:
add.scores=list(HLAIScore=c("B2M", "HLA-A", "HLA-B","HLA-C"), HLAIIScore=c("HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1","HLA-DRA","HLA-DRB1"), CytolyticScore=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"))
add.scores=append(append(geneset.l, append(geneset.l2, add.scores)), geneset.l3)



load("Szabo_Tcell_BM_scRNA.Rdata")
szabo=scmat
szabo <- NormalizeData(szabo, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
fm.szabo <- data.matrix(GetAssayData(szabo, assay = "RNA", slot="data"))

# same in FIMM data:
load("FIMM_AML_HCA_T_scRNA.Rdata")

scmat=FindClusters(scmat, resolution = 0.5)

# make sure that using logNormalized data and scaled to 10000:
scmat <- NormalizeData(scmat, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
fm.f <- data.matrix(GetAssayData(scmat, assay = "RNA", slot="data"))

sample=gsub("_.*.", "", colnames(scmat))
sample[grepl("Manton", sample)]="HCA BM"
sample[!grepl("3667|5249|5897", sample)&!grepl("HCA", sample)]="AML other"
scmat[["sample"]]=sample

# take cluster assignment for AML cells:
scmat2=FindNeighbors(scmat[,grepl("3667|5897|5249", sample)], dims = 1:25, reduction = "mnnCorrect")
scmat2=FindClusters(scmat2, resolution = 0.5, algorithm = 1)

AML_cluster_Tcell=Idents(scmat2)
save(AML_cluster_Tcell, file="Szabo_HCA_AML_cluster_cell.Rdata")

pdf("Figure_S3M_FIMM_AML_reclustered_Tcell.pdf", width = 6, height = 5)
DimPlot(scmat, pt.size = 0.1, cells.highlight = colnames(scmat2)[Idents(scmat2)%in%0], cols.highlight = "#006f00") # cytotoxic
DimPlot(scmat, pt.size = 0.1, cells.highlight = colnames(scmat2)[Idents(scmat2)%in%3], cols.highlight = "#006f00") # Cytokine
DimPlot(scmat2, label = T) # all clusters
dev.off()

geneset.l=lapply(unique(geneset[,2]), function(g){
  genes=geneset[geneset[,2]%in%g,1]
  
  #take only significant:
  genes=genes[genes%in%markers.all$gene[grepl("CD|Treg", markers.all$cluster)]]
})
names(geneset.l)=paste(unique(geneset[,2]), "filtered")

# add.scores=append(add.scores, geneset.l)

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

# add gm based score:
gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
  dat3=fm.szabo[rownames(fm.szabo)%in%add.scores[[i]],]
  gm=log2(t(apply(dat3, 2, gm_mean))) # done to normalized values
  rownames(gm)=names(add.scores)[i]
  return(gm)
}))

# also add to seurat object:
for(i in seq(add.scores)){
  szabo[[names(add.scores)[i]]] <- gm.objects[i,]
}

# pdf("FigureS3M_Szabo_markers_FIMM_AML_HCA_Tcells_Scores.pdf", width=16, height=10)
# FeaturePlot(scmat, features = names(add.scores), cols=c("grey75", "red"), min.cutoff = 0.1)
# FeaturePlot(szabo, features = names(add.scores), cols=c("grey75", "red"), min.cutoff = 0.1)
# dev.off()

pdf("Figure3M_Szabo_markers_FIMM_AML_HCA_Tcells_Scores_scaled.pdf", width=16, height=14)
p1=FeaturePlot(scmat, features = names(add.scores), combine = F, min.cutoff = 0.2, max.cutoff = 1)
fix.sc <- scale_color_gradientn( colours = c('grey75', 'red'),  limits = c(0.2, 1))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

p1=FeaturePlot(szabo, features = names(add.scores), combine = F, min.cutoff = 0.2, max.cutoff = 1)
fix.sc <- scale_color_gradientn( colours = c('grey75', 'red'),  limits = c(0.2, 1))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
dev.off()

pdf("Figure3I_FigureS3M_Szabo_markers_FIMM_AML_HCA_Tcells_Scores_blend.pdf", width=20, height=5)
FeaturePlot(scmat, features = names(add.scores)[c(14,15)], blend.threshold = 0, blend = T,cols = c("grey74", "darkblue", "red"))
FeaturePlot(scmat, features = names(add.scores)[c(8,9)], blend = T, blend.threshold = 0.6, cols = c("grey74", "darkblue", "red"))
dev.off()

colors.group=data.table::fread("colors_lineage.txt", data.table = F, header = F)
colors.group=colors.group[grepl("CD|Treg", colors.group$V1),]
plot.multi.scatter.scRNA(data.frame("Sample"="T cells HCA AML", "Data"="FIMM_AML_HCA_T_scRNA.Rdata", "Clusters"="", "Type"="", stringsAsFactors = F), colors.group = colors.group, seurat.feature = "SingleR.label", name="FigureS3M_singleR", cores=9, SIZE = 0.1, rasterize=F, width = 77*4, height = 74*4, text.size = 10*4)

# plot MDS-like samples to map:
scmat[["MDSlike"]]=ifelse(grepl("5897|3667|5249", scmat[["batch"]][,1]), "MDS-like", "other")
scmat[["MDSlike"]][grepl("BM", scmat[["batch"]][,1]),1]="normal BM"

lv = scmat[["MDSlike"]][,1]%in%c("MDS-like", "other")
# plot.scatter.seurat(scmat = scmat, colors.group = data.frame("V1"=c("normal BM", "other", "MDS-like"), "V2"=c("grey75", "orange", "#21a366"), stringsAsFactors = F), seurat.feature = "MDSlike", name="T cells (AML, HCA)",lv = lv, cores=1, add.density = T, add.proportions = T, SIZE = 1, rasterize=F, width = 64*4, height = 74*4, text.size=10*4)

xy=data.frame(Embeddings(scmat, "umap"))

pdf("Figure3I_Density_MDSlike_Tcell.pdf", width = 5, height = 4)
ggplot(xy[lv,], aes(x=UMAP_1, y=UMAP_2) ) +
stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +
  viridis::scale_fill_viridis()
ggplot(xy[scmat[["MDSlike"]]=="normal BM",], aes(x=UMAP_1, y=UMAP_2) ) +
  stat_density2d(aes(fill = ..density..), contour = F, geom = 'tile') +
  viridis::scale_fill_viridis()
dev.off()