library(GSVA)
GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "pathway_analysis/functions.GSEA.R"))
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "statistics/useful_functions.R"))

# cancermap coordinates:
library(reshape2)
library(gridExtra)
library(GSVA)
library(RColorBrewer)
library(ggplot2)
library(Rtsne)
library(LPCM)

# plotting function
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

load("MDS_genesets.Rdata")
geneset=list("MDS_signature"=geneset$MDS_signature_all_filt)

# go through each cohort and plot GSVA score:
setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

files=list.files(path = ".", "subtypes.Rdata")
names(files)=gsub("_subtypes.Rdata", "", files)

files=files[grepl("AML", names(files))]

p.all=lapply(files, function(f){
  load(f)
  
  score=gsva(data.matrix(gexp), geneset, mx.diff=F, tau=0.25, parallel.sz=4)
  
  plot.val=t(scale(t(score)))
  
  a=Plot_GSVA_scores(feat = "MDS_signature", data_plot = t(scale(t(score))), VALUE = 0.5,SIZE = 1, coord = coordinates.subtype, peaks = NULL, CLUSTER_CENTRE = F)
  
})

# Save PDF figure (A4) with multiple panels
ggsave("Figure3E_FigureS3D.pdf", do.call(marrangeGrob, list(grobs=p.all, nrow=4, ncol=3)), width = 210, height = 297, units = "mm", dpi=150)
