# created with 3.3.3 Rtsne_0.13
# first make cancer map and validate clusters:
# gexp annot:
source("/research/users/ppolonen/git_home/Hemap_published_scripts/Cancermap.R")
source("/research/users/ppolonen/git_home/Hemap_published_scripts/Plot_color_vector.R")
source("/research/users/ppolonen/git_home/Hemap_published_scripts/Plot_cancermap_clusters.R")
source("/research/users/ppolonen/git_home/Hemap_published_scripts/Example_useCases/useCase3/Find_correlated_genes.R")
source("/research/users/ppolonen/git_home/common_scripts/statistics/functions_statistics.R")
source("/research/users/ppolonen/git_home/common_scripts/statistics/statistics_wrappers.R")

# cancermap coordinates:
library(reshape2)
library(gridExtra)
library(GSVA)
library(RColorBrewer)
library(ggplot2)
library(Rtsne)
library(LPCM)
library(parallel)
library(ComplexHeatmap)
library(circlize)


setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")

# gexp data
data=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data9544_with_gene_symbols.RData"))

# annotations
annot = get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Hemap_immunology_Annotations.Rdata"))
gexp.hemap=t(data[rownames(data)%in%annot[annot$colorClass=="AML"&annot$Sample.type=="Cancer",1],]) # GSE16716,GSE24080 and GSE19784
annot.sub=annot[annot$colorClass=="AML"&annot$Sample.type=="Cancer",]

coord_hemap=read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/webresource/data9544_AML_cancermap_cluster_15pct_coordinates.txt", header=T, stringsAsFactors=F)
peaks_hemap=read.table("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/webresource/data9544_AML_cancermap_cluster_15pct_cluster_centroids.txt", sep="\t", header=T, stringsAsFactors=F)

coord_hemap=coord_hemap[coord_hemap$ID%in%annot.sub$GSM.identifier..sample.,]

a1=Plot_cancermap_clusters(X = coord_hemap, peaks = peaks_hemap)

#*********************** make genesets from TCGA AML: ********************

matrix=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/DUFVA_TCGA_AML_FM_meth.Rdata")) # NAME: matrix
gex=data.matrix(matrix[grep("GEXP", rownames(matrix)),])
rownames(gex)=make.unique(gsub("N:GEXP:|:chr.*.|::.*.","", rownames(gex)))
coord=read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/webresource/TCGA_AML_cancermap_cluster_15pct_coordinates.txt", header=T, stringsAsFactors=F)
peaks=read.table("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/webresource/TCGA_AML_cancermap_cluster_15pct_cluster_centroids.txt", sep="\t", header=T, stringsAsFactors=F)

matrix=matrix[,match(gsub("\\.", "-", substr(coord$ID, 1, 12)), substr(colnames(matrix), 1, 12))]
gex=gex[,match(gsub("\\.", "-", substr(coord$ID, 1, 12)), substr(colnames(gex), 1, 12))]

coord$X1.5..cluster=as.numeric(coord$X1.5..cluster)
clusters=coord$X1.5..cluster

clusters[coord$X1.5..cluster%in%c(1)]="PML-RARA"
clusters[coord$X1.5..cluster%in%c(2)]="CBFB-MYH11"
clusters[coord$X1.5..cluster%in%c(3)]="CMP-like"
clusters[coord$X1.5..cluster%in%c(4)]="MDS-like"
clusters[coord$X1.5..cluster%in%c(5)]="CEBPA"
clusters[coord$X1.5..cluster%in%c(6)]="Monocyte-like"
clusters[coord$X1.5..cluster%in%c(7)]="RUNX1-RUNX1T1"
clusters[grepl("MLL", as.character(matrix["C:SAMP:GENETICS:::::",]))]="Monocyte-like-MLL"

coord$subtype=clusters

coord$ID=gsub("\\.", "-", substr(coord$ID, 1, 12))

# make gene sets
genesets.cor=unlist(mclapply(unique(clusters), Find_correlated_genes_new, gex, clusters, mc.cores=7), recursive=F)

gexp=gex
coordinates.subtype=coord
  
# AML data:
save(list = c("gexp", "coordinates.subtype"), file="TCGA_AML_subtypes.Rdata")

#*******************************************************************************

# run GSVA
gsva_es <- gsva(as.matrix(gexp.hemap), method="gsva", genesets.cor, mx.diff=F, tau=0.25, verbose=T, min.sz=5, max.sz=500, parallel.sz=8)

#********************************************* Step 3: Plot GSVA scores to hemap cancermap ************************************************************
feats=rownames(gsva_es)

gsva_es2=t(scale(t(gsva_es)))

Plot_GSVA_scores=function(feat, data_plot, VALUE, SIZE, CLUSTER_CENTRE, coord, peaks){
  
  # Find specified feature(feat) from data_plot
  data=data_plot[rownames(data_plot)%in%feat,]
  
  # Transform matrix to numeric. 
  data=as.numeric(data)
  
  # Color vector for gradient colors from blue to red.
  rbPal <- colorRampPalette(c('blue','red'))
  
  # Adjust data for gradient colors.
  data=c(data, 2, -2) # adjust range
  datCol <- rbPal(10)[as.numeric(cut(data,breaks = 10))]
  datCol=datCol[-c(length(datCol)-1, length(datCol))]
  data=data[-c(length(data)-1, length(data))]
  
  # Samples below cutoff colored grey.
  datCol[abs(data)<VALUE]="grey75"
  front=abs(data)>VALUE
  
  # Prepare coordinate data for plotting.
  dat2show <- cbind(coord$x, coord$y)
  df=as.data.frame(dat2show)
  colnames(df) = c("X1","X2")
  
  # Generate plot title. 
  cutoff=paste("GSVA score >", VALUE)
  plotname=paste(feat, cutoff, sep="\n")
  
  # Call actual plotting function.
  drawFig(df, CLUSTER_CENTRE, datCol, front, plotname, SIZE, peaks) 
}

# Plot GSVA scores to Hemap cancermap, obtain plotting data
p.all=lapply(feats, Plot_GSVA_scores, gsva_es2, VALUE=0.5, SIZE=0.4, CLUSTER_CENTRE=T, coord_hemap, peaks_hemap)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0("Hemap_AML_TCGA_multipage.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)

coordinates.subtype=coord_hemap
# checked the genetics and combined the subtypes as larger subtypes:
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(14)]="PML-RARA"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(6)]="CBFB-MYH11"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(23,8,9,26,13,1)]="CMP-like"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(24,2,22,15,12,3,19)]="MDS-like"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(4)]="CEBPA"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(16,29,5,25,28)]="Monocyte-like"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(11,17)]="RUNX1-RUNX1T1"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(7)]="Monocyte-like-MLL"
coordinates.subtype$subtype[coord_hemap$X1.5..cluster%in%c(18,20,21,30, 10)]=NA

gexp=as.matrix(gexp.hemap)

# AML data:
save(list = c("gexp", "coordinates.subtype"), file="Hemap_AML_subtypes.Rdata")

pl=gsva_es[,order(match(coordinates.subtype$subtype, c("PML-RARA", "CBFB-MYH11","CMP-like", "MDS-like", "CEBPA", "Monocyte-like", "Monocyte-like-MLL", "RUNX1-RUNX1T1")))]
pl=t(scale(t(pl)))
pl[pl>2]=2

ha=HeatmapAnnotation("TSNE"=as.character(coordinates.subtype$subtype)[order(coordinates.subtype$subtype)])

pdf("TCGA_cluster_hemap_AML.pdf", height = 8, width = 8)
Heatmap(pl,top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
Heatmap(pl[grepl("UP", rownames(pl)),], top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
dev.off()


# Also subtype the beatAML:
load("/research/groups/sysgen/PROJECTS/LEUKEMIA/petri_wrk/BeatAML/BeatAML_gexp.Rdata")
load("/research/groups/sysgen/PROJECTS/LEUKEMIA/petri_wrk/BeatAML/BeatAML_fm_annot.Rdata")
coord.beatAML=read.delim("/research/groups/sysgen/PROJECTS/LEUKEMIA/petri_wrk/BeatAML/cancermap_BeatAML_15pct_genes_BH-SNE_mean-shift_BW1.5.txt", header=T, stringsAsFactors = F)
peaks.beatAML=read.delim("/research/groups/sysgen/PROJECTS/LEUKEMIA/petri_wrk/BeatAML/cancermap_BeatAML_15pct_genes_BH-SNE_mean-shift_BW1.5_cluster_centroids.txt", header=T, stringsAsFactors = F)

annot=annot[!is.na(annot$TCGA_coord),]

# run GSVA
gsva_es <- gsva(as.matrix(gexp), method="gsva", genesets.cor, mx.diff=F, tau=0.25, verbose=T, min.sz=5, max.sz=500, parallel.sz=8)

#********************************************* Step 3: Plot GSVA scores to hemap cancermap ************************************************************
feats=rownames(gsva_es)

gsva_es2=t(scale(t(gsva_es)))

# Plot GSVA scores to Hemap cancermap, obtain plotting data
p.all=lapply(feats, Plot_GSVA_scores, gsva_es2, VALUE=0.5, SIZE=0.4, CLUSTER_CENTRE=T, coord.beatAML, peaks.beatAML)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0("BeatAML_TCGA_multipage.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)

coordinates.subtype=coord.beatAML
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(13)]="PML-RARA"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(11)]="CBFB-MYH11"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(2,7,6)]="CMP-like"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(9, 14,3, 12, 5)]="MDS-like"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(8)]="CEBPA"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(4,1)]="Monocyte-like"
coordinates.subtype$subtype[coord.beatAML$cluster%in%c(10)]="RUNX1-RUNX1T1"
coordinates.subtype$subtype[annot$finalFusion%in%"MLLT3-KMT2A"]="Monocyte-like-MLL"

coordinates.subtype$cluster[coordinates.subtype$cluster%in%6]="CMP-like-HLAII.low"

# AML data:
save(list = c("gexp", "coordinates.subtype"), file="BeatAML_subtypes.Rdata")

pl=gsva_es[,order(match(coordinates.subtype$subtype, c("PML-RARA", "CBFB-MYH11","CMP-like", "MDS-like", "CEBPA", "Monocyte-like", "Monocyte-like-MLL", "RUNX1-RUNX1T1")))]
pl=t(scale(t(pl)))
pl[pl>2]=2

ha=HeatmapAnnotation("TSNE"=as.character(coordinates.subtype$subtype)[order(coordinates.subtype$subtype)])

pdf("TCGA_cluster_beatAML.pdf", height = 8, width = 8)
Heatmap(pl,top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
Heatmap(pl[grepl("UP", rownames(pl)),], top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
dev.off()


