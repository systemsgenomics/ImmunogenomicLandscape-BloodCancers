source("/research/users/ppolonen/git_home/common_scripts/scRNA/functions.scRNA.analysis.R")
# tools:
library(Matrix)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# Compare MDS in blasts between FIMM AML and then compare to HCA
load("HCA_scRNA.Rdata")

MDS=get(load("MDS_genesets.Rdata"))

Idents(scmat)=scmat[["SingleR.label"]]
markers=FindAllMarkers(scmat, features = MDS$MDS_signature_all_filt[MDS$MDS_signature_all_filt%in%rownames(scmat)], only.pos = T)

for(i in 1:5){
  markers[,i]=prettyNum(signif(markers[,i],2))
}
write.table(markers[markers$avg_logFC>0,c("gene",	"cluster",	"avg_logFC",	"pct.1",	"pct.2",	"p_val",	"p_val_adj")], "TableS3_MDS_signature_markers_celltype_HCA.txt", quote = F, row.names = F, sep="\t")



#********************************** Genes expressed in MDS Blasts vs other blasts **************************************
markers=markers[markers$avg_logFC>0,]
erythroid=markers$gene[markers$cluster%in%c("Erythrocytes")]
HSC=markers$gene[markers$cluster%in%c("HSC", "MPP", "MEP")]
HSC=HSC[!HSC%in%erythroid]

scmat[["SingleR.label2"]]=factor(scmat[["SingleR.label"]][,1], levels=c("HSC","MPP","MEP", "CLP","CMP","GMP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes"))

pdf("HCA_HSC_genes.pdf", width = 5, height = 4)
plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = unique(HSC), cols = c("white", "red"), dot.scale = 5)
dev.off()

pdf("HCA_erythroid_genes.pdf", width = 5, height = 8)
plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = unique(erythroid), cols = c("white", "red"), dot.scale = 5)
dev.off()

pdf("HCA_HSC_erythroid_genes_filt.pdf", width = 4, height = 3)
plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = c("DPPA4","CRHBP","CSF3R","BEX1","EREG", "CA1", "HBD", "HEMGN", "MYL4", "SNCA"), cols = c("white", "red"), dot.scale = 6, scale.max = 50, scale.min = 0)
dev.off()

# same in FIMM MDS samples
load("FIMM_AML_scRNA.Rdata")

Idents(scmat)=scmat[["SingleR.label"]]
markers=FindAllMarkers(scmat, features = MDS$MDS_signature_all_filt[MDS$MDS_signature_all_filt%in%rownames(scmat)], only.pos = T, logfc.threshold = 0.025)
excl=table(scmat[["SingleR.label"]])<0.001*dim(scmat)[2]

for(i in 1:5){
  markers[,i]=prettyNum(signif(markers[,i],2))
}

write.table(markers[markers$avg_logFC>0&!markers$cluster%in%names(excl)[excl],c("gene",	"cluster",	"avg_logFC",	"pct.1",	"pct.2",	"p_val",	"p_val_adj")], "TableS3_MDS_signature_markers_celltype_AML.txt", quote = F, row.names = F, sep="\t")

scmat[["SingleR.label2"]]=factor(scmat[["SingleR.label"]][,1], levels=c("HSC","MPP","MEP", "CLP","CMP","GMP","Monocytes","DC","Macrophages","Macrophages M1","Macrophages M2","CD4+ Tn","CD4+ T-cells", "CD4+ Tcm", "CD4+ Tem","CD8+ Tn","CD8+ T-cells","CD8+ Tcm","CD8+ Tem","NK cells","Tregs","naive B-cells","Memory B-cells","Class-switched memory B-cells","Plasma cells","Endothelial cells","Neutrophils","Eosinophils","Fibroblasts","Smooth muscle","Erythrocytes","Megakaryocytes"))

# pdf("FIMM_AML_HSC_genes.pdf", width = 5, height = 4)
# plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = unique(HSC), cols = c("white", "red"), dot.scale = 5)
# dev.off()
# 
# pdf("FIMM_AML_erythroid_genes.pdf", width = 5, height = 8)
# plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = unique(erythroid), cols = c("white", "red"), dot.scale = 5)
# dev.off()

pdf("FigureS3L.pdf", width = 4, height = 3)
plot.DotPlot(scmat[,grepl("HSC|MPP|MEP|CMP|GMP|Mono|Eryth", scmat[["SingleR.label2"]][,1])], group.by = "SingleR.label2", features = c("DPPA4","CRHBP","CSF3R","BEX1","EREG", "CA1", "HBD", "HEMGN", "MYL4", "SNCA"), cols = c("white", "red"), dot.scale = 6, scale.max = 50, scale.min = 0)
dev.off()
