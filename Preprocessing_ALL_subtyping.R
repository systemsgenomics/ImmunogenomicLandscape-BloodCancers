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

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")

# gexp data
data=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data9544_with_gene_symbols.RData"))

# annotations
annot = get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Hemap_immunology_Annotations.Rdata"))
gexp.hemap=t(data[rownames(data)%in%annot[annot$colorClass=="pre-B-ALL"&annot$Sample.type=="Cancer",1],])
annot.sub=annot[annot$colorClass=="pre-B-ALL"&annot$Sample.type=="Cancer",]

load("Hemap_ALL_subtypes.Rdata")

# cancermap from pre-B-ALL, gu used about 1000 genes ~ 5% var.
coord_hemap=CancerMap(data = t(as.matrix(gexp.hemap)), name = "Hemap_pre-B-ALL", VAR = 5, BW = 1.25, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")

plot.scatter(x=coord_hemap$x, y = coord_hemap$y, group =  coordinates.subtype$subtype, namev = coordinates.subtype$subtype, main = "Hemap ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)

coord_hemap=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/cancermap_Hemap_pre-B-ALL_5pct_genes_BH-SNE_mean-shift_BW1.75.txt", header=T, stringsAsFactors = F)
peaks_hemap=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/cancermap_Hemap_pre-B-ALL_5pct_genes_BH-SNE_mean-shift_BW1.75_cluster_centroids.txt", header=T, stringsAsFactors = F)

a1=Plot_cancermap_clusters(X = coord_hemap, peaks = peaks_hemap)

#*********************** make genesets from St judes ALL: ********************
load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/PecanALL_subtypes.Rdata")

# clusters=coordinates.subtype$subtype
# # make gene sets
# genesets.cor=unlist(mclapply(unique(clusters), Find_correlated_genes_new, gexp, clusters, mc.cores=7), recursive=F)
# 
# clusters2=annot$`2nd subtype`
# 
# # make gene sets
# genesets.cor.2=unlist(mclapply(unique(clusters2), Find_correlated_genes_new, gexp, clusters2, mc.cores=7), recursive=F)
# 
# genesets.cor=append(genesets.cor, genesets.cor.2)


# use predefined gene sets from Gu et al:

genesets.up=lapply(3:16, function(i){
  genes=openxlsx::read.xlsx("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL/B-ALL-subtype-signature.xlsx", sheet = i, colNames = TRUE)
  genes=genes$geneName[genes$log2FoldChange>1&genes$padj<0.001]
  genes=genes[!is.na(genes)]
})

names(genesets.up)=c("PAX5 P80R", "PAX5alt", "IKZF1 N159Y", "Hyperdiploid", "ETV6-RUNX1", "KMT2A", "Ph", "TCF3-PBX1", "Low hypodiploid", "DUX4", "ZNF384", "MEF2D", "iAMP21", "BCL2 MYC")
  

genesets.down=lapply(3:16, function(i){
  genes=openxlsx::read.xlsx("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL/B-ALL-subtype-signature.xlsx", sheet = i, colNames = TRUE)
  genes=genes$geneName[genes$log2FoldChange>1&genes$padj<0.001]
  genes=genes[!is.na(genes)]
  genes=head(genes, 20)
})

names(genesets.down)=c("PAX5 P80R", "PAX5alt", "IKZF1 N159Y", "Hyperdiploid", "ETV6-RUNX1", "KMT2A", "Ph", "TCF3-PBX1", "Low hypodiploid", "DUX4", "ZNF384", "MEF2D", "iAMP21", "BCL2 MYC")

# cancermap from pre-B-ALL, gu used about 1000 genes ~ 5% var.
coord_pecan=CancerMap(data = t(as.matrix(gexp)), name = "Pecan_pre-B-ALL", VAR = 10, BW = 1.75, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")
plot.scatter(x=coord_pecan$x, y = coord_pecan$y, group =  coordinates.subtype$`primary subtype`, namev = coordinates.subtype$`primary subtype`, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)

# gsva_es <- gsva(as.matrix(gexp), method="gsva", genesets.up, mx.diff=F, tau=0.25, verbose=T, min.sz=5, max.sz=10500, parallel.sz=8)
# save(gsva_es, file="gsva_es_Pecan_ALL_pecan_theirsets.Rdata")

pecan_gsva=get(load("gsva_es_Pecan_ALL_pecan_theirsets.Rdata"))

# plot clusters
coord_pecan=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/cancermap_Pecan_pre-B-ALL_10pct_genes_BH-SNE_mean-shift_BW1.75.txt", header=T, stringsAsFactors = F)
peaks_pecan=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/cancermap_Pecan_pre-B-ALL_10pct_genes_BH-SNE_mean-shift_BW1.75_cluster_centroids.txt", header=T, stringsAsFactors = F)

a1=Plot_cancermap_clusters(X = coord_pecan, peaks = peaks_pecan)

feats=rownames(pecan_gsva)

gsva_es2=t(scale(t(pecan_gsva)))

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
p.all=lapply(feats, Plot_GSVA_scores, gsva_es2, VALUE=0.5, SIZE=0.2, CLUSTER_CENTRE=T, coord_hemap, peaks_hemap)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0("PECAN_ALL_theirsubtypes_multipage.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)

#************************ change the name to match with their gene sets: ************************

# cohort=openxlsx::read.xlsx("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Pecan_ALL/B-ALL-subtype-signature.xlsx", sheet = 1, colNames = TRUE)
# 
# gexp=gexp[,match(cohort$patient, colnames(gexp))]
# 
# # cancermap from pre-B-ALL, gu used about 1000 genes ~ 5% var.
# coordinates.subtype=CancerMap(data = t(as.matrix(gexp)), name = "Pecan_pre-B-ALL_835", VAR = 10, BW = 1.75, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")
# 
# coordinates.subtype$subtype=cohort$primary.subtype
# plot.scatter(x=coordinates.subtype$x, y = coordinates.subtype$y, group =  coordinates.subtype$subtype, namev = coordinates.subtype$subtype, main = "Pecan ALL", rasterize = F, width = 70*2, height = 74*2, SIZE = 0.5)
# 
# save(list = c("gexp", "coordinates.subtype"), file = "Pecan_ALL_835_subtypes.Rdata")
# 

#*******************************************************************************

# run GSVA
gsva_es <- gsva(as.matrix(gexp.hemap), method="gsva", genesets.up, mx.diff=F, tau=0.25, verbose=T, min.sz=5, max.sz=10500, parallel.sz=8)
save(gsva_es, file="gsva_es_hemap_ALL_pecan_theirsets.Rdata")
load("gsva_es_hemap_ALL_pecan_theirsets.Rdata")

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
p.all=lapply(feats, Plot_GSVA_scores, gsva_es2, VALUE=0.5, SIZE=0.2, CLUSTER_CENTRE=T, coord_hemap, peaks_hemap)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0("Hemap_ALL_PECAN_ALL_theirsubtypes_multipage.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)

coordinates.subtype=coord_hemap
# checked the genetics and combined the subtypes as larger subtypes:
coordinates.subtype$subtype[coord_hemap$cluster%in%c(23)]="PAX5 P80R"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(12)]="PAX5alt"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(17,14,31,20,6)]="ETV6-RUNX1"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(9)]="DUX4"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(21,4,22, 27)]="KMT2A"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(1,29)]="ZNF384"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(13,10,25,24,33)]="Ph"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(19,16,32)]="Hyperdiploid"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(2,15)]="TCF3-PBX1"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(26)]="MEF2D"
coordinates.subtype$subtype[coord_hemap$cluster%in%c(3, 30, 5,11,8,18,28,7)]="Other"


Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_BCR.ABL1==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_hyperdiploid==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_hypodiploid==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_TEL.AML1==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_TCF3_PBX1==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_IgH_cMYC==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_MLL==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_other==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_TELdel==1, "darkgreen", "grey75"))
Plot_color_vector(coord_hemap, color = ifelse(annot.sub$GENETICS_preBALL_trisomy_chr4_chr10==1, "darkgreen", "grey75"))

aml.gexp=as.matrix(gexp)

ord=c("PAX5 P80R", "PAX5alt","ETV6-RUNX1","DUX4", "KMT2A", "ZNF384", "Ph", "Hyperdiploid", "TCF3-PBX1", "MEF2D", "other")
pl=gsva_es[,order(match(coordinates.subtype$subtype, ord))]
pl=t(scale(t(pl)))
pl[pl>2]=2

ha=HeatmapAnnotation("TSNE"=as.character(coordinates.subtype)[order(match(coordinates.subtype, ord))])

pdf("PECAN_ALL_subtype_hemap_ALL.pdf", height = 8, width = 8)
Heatmap(pl[match(ord, rownames(pl)),],top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
dev.off()

ha=HeatmapAnnotation("TSNE"=as.character(coordinates.subtype$subtype)[order(match(coordinates.subtype$subtype, ord))], df = annot.sub[order(match(coordinates.subtype$subtype, ord)),grepl("GENETICS_preBALL", colnames(annot.sub))])

pdf("PECAN_ALL_subtype_hemap_ALL_validation.pdf", height = 8, width = 16)
Heatmap(pl[match(ord, rownames(pl)),],top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
dev.off()


gexp=gexp.hemap

save(list = c("gexp", "coordinates.subtype"), file="Hemap_ALL_subtypes.Rdata")
