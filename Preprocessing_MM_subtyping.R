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

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/")

# gexp data
data=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data9544_with_gene_symbols.RData"))

# annotations
annot = get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/Hemap_immunology_Annotations.Rdata"))
gexp.hemap=t(data[rownames(data)%in%annot[annot$colorClass=="MM"&annot$Sample.type=="Cancer"&annot$GSE.identifier..experiment.%in%c("GSE16716,GSE24080"),1],]) # GSE16716,GSE24080 and GSE19784
annot.sub=annot[annot$colorClass=="MM"&annot$Sample.type=="Cancer"&annot$GSE.identifier..experiment.%in%c("GSE16716,GSE24080"),]
  
# remove gender spesific genes
res2=wrapper.wilcoxtest(genelist = rownames(gexp.hemap), data = gexp.hemap, logicalVectors = list("gender"=annot.sub$GENDER=="male"), ALTERNATIVE = "two.sided", adj.method = "BH", CORES = 8, prettynum = F)
rm2=res2[res2$FDR<0.05,1]

coord=CancerMap(data = t(as.matrix(gexp.hemap[!rownames(gexp.hemap)%in%rm2,])), name = "Hemap_MM", VAR = 12.5, BW = 1.5, perplexity = 30, PATH_OUTPUT = "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/HEMAP_MM/")

coord.hemap=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/HEMAP_MM/cancermap_Hemap_MM_12.5pct_genes_BH-SNE_mean-shift_BW1.5.txt", header=T, stringsAsFactors = F)
peaks.hemap=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/HEMAP_MM/cancermap_Hemap_MM_12.5pct_genes_BH-SNE_mean-shift_BW1.5_cluster_centroids.txt", header=T, stringsAsFactors = F)

a1=Plot_cancermap_clusters(X = coord, peaks = peaks)

#*********************** make genesets from COMPASS MM: ********************
load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/MM_COMPASS_FM.Rdata")
load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/MM_COMPASS_ANNOT.Rdata")
annot=annot[!is.na(annot$subtype.cluster),]

gex=fm[grepl("N:GEXP", rownames(fm)),match(rownames(annot), colnames(fm))]
rownames(gex)=gsub("N:GEXP:", "", rownames(gex))
gex=gex[rownames(gex)%in%rownames(gexp),]

#clusters
clusters=annot$subtype.cluster[rownames(annot)%in%colnames(gex)]

# make gene sets
genesets.cor=unlist(mclapply(sort(unique(clusters)), Find_correlated_genes_new, gex, clusters, mc.cores=7), recursive=F)

# save data:
coordinates.subtype=annot[,c(109:112, 108)]
colnames(coordinates.subtype)=c("ID", "x", "y", "subtype", "cluster")

coordinates.subtype$subtype[coordinates.subtype$subtype%in%"Hyperdiploid"]="Hyperdiploid_gain11q"
coordinates.subtype$subtype[coordinates.subtype$subtype%in%"Hyperdiploid_amp1q"]="Hyperdiploid_gain1q"

gexp=as.matrix(gex)

save(list = c("gexp", "coordinates.subtype"), file="CoMMpass_MM_subtypes.Rdata")

#*******************************************************************************
gs=c(genesets.cor)

append.list=function(vector.lists){
  unlist(lapply(vector.lists, list), recursive=F)
}

gs=append.list(gs)

# run GSVA
gsva_es <- gsva(as.matrix(gexp.hemap), method="gsva", gs, mx.diff=F, tau=0.25, verbose=T, min.sz=5, max.sz=500, parallel.sz=8)

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
p.all=lapply(feats, Plot_GSVA_scores, gsva_es2, VALUE=0.5, SIZE=0.7, CLUSTER_CENTRE=T, coord, peaks)
p.all=p.all[!sapply(p.all, is.null)] # remove empty just in case

# Save PDF figure (A4) with multiple panels
ggsave(paste0("Hemap_MM_COMPASS_multipage_GSE16716,GSE24080.pdf"),
       do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)


coord2=coord.hemap
# checked the genetics and combined the subtypes as larger subtypes:
coord2$subtype[coord$cluster%in%c(10)]="MAF_Ig"
coord2$subtype[coord$cluster%in%c(4)]="WHSC1_FGFR3_Ig"
coord2$subtype[coord$cluster%in%c(1,3,9)]="CCND1_Ig"
coord2$subtype[coord$cluster%in%c(2,7,5)]="Hyperdiploid_gain11q"
coord2$subtype[coord$cluster%in%c(8)]="TRAF3_Aberrated"
coord2$subtype[coord$cluster%in%c(6)]="Hyperdiploid_gain1q"


pl=gsva_es[,order(coord2$subtype)]
pl=t(scale(t(pl)))
pl[pl>2]=2

# add CGA number:
t.df = read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/t.antigen_df.txt", stringsAsFactors=F, header=T)
CGA=c(unique(t.df[,1]))

profile=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/revision_2018/petri/mixtureM_profile.Rdata"))
profile[profile==-1] = 0
profile2=profile[,colnames(profile)%in%annot.sub$GSM.identifier..sample.]

# take only high expressed into account
profile2[data.matrix(t(data[rownames(data)%in%annot.sub$GSM.identifier..sample.,]))<5]=0

expressed_testis_num=colSums(profile2[rownames(profile2)%in%unique(t.df$gene),])

ha=HeatmapAnnotation(n.CGA=anno_barplot(expressed_testis_num), "TSNE"=as.character(coord2$cluster)[order(coord2$cluster)])

pdf("CoMMpass_cluster_hemapMM.pdf", height = 8, width = 8)
Heatmap(pl,top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)
Heatmap(pl[grepl("UP", rownames(pl)),],top_annotation = ha, cluster_rows = F, cluster_columns=F, show_column_names = F, cluster_column_slices = T)

pl=lapply(unique(coord2$cluster), function(clu)expressed_testis_num[coord2$cluster%in%clu])
names(pl)=unique(coord2$cluster)
boxplot(pl)
dev.off()

coordinates.subtype=coord2
gexp=as.matrix(gexp.hemap)

save(list = c("gexp", "coordinates.subtype"), file="Hemap_MM_subtypes.Rdata")

