GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "pathway_analysis/functions.GSEA.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)

annot = get(load("Hemap_immunology_Annotations.Rdata"))

# gexp data
data=t(get(load("data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]

t.df=t.df[order(t.df[,3]),]
t.df=t.df[t.df$name%in%"Cancer_Myeloma",]
# rank patients by number of testis antigens expressed
# profiles
profile=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/revision_2018/petri/mixtureM_profile.Rdata"))
profile[profile==-1] = 0
profile2=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]

# take only high expressed into account
profile2[data.matrix(data)<5]=0

expressed_testis_num=colSums(profile2[rownames(profile2)%in%unique(t.df$gene),])
feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0_Antigens"
feat_class[expressed_testis_num>=1&expressed_testis_num<=4]="1to4_Antigens"
feat_class[expressed_testis_num>=5&expressed_testis_num<=6]="5to6_Antigens"
feat_class[expressed_testis_num>=7]="over7_Antigens"

# plot CGA to hemap:

colorv=rep("grey75", length(annot$y))

colorf=colorRamp2(breaks = c(0, 2, 4, 6, 7),colors = c("grey75", "#f4acac", "#f77777", "#e23030", "#800000"))

Plot_color_vector=function(X, color, NAME=NULL, SIZE=0.8, TITLE="tSNE plot", PATH_OUTPUT=getwd(), peaks=NULL) {
  
  # if cluster centers are not defined as inputs, do not plot them.
  CLUSTER_CENTRE=ifelse(is.null(peaks), F, T)
  
  # Color vector
  datCol=color
  
  # Prepare input matrix for plotting.
  dat2show <- cbind(X$x, X$y)
  df=as.data.frame(dat2show)
  colnames(df) = c("X1","X2")
  
  # If for some reason we use a color vector which includes blanks,
  # create a color vector without them to be plotted at the front.
  # (Increases visibility of the "important" colored samples.)
  front=!datCol==""
  
  # Call actual plotting function.
  p=drawFig(df, CLUSTER_CENTRE, datCol, front, TITLE, SIZE, peaks)
  
  if(!is.null(NAME)){
    # Write out print quality figure as PDF.
    ggsave(paste0(PATH_OUTPUT, NAME, "_singlepage.pdf"), p, width = 210, height = 210, units = "mm", dpi=150)
  }
  return(p)
}

pdf("FigS6B_TSNE_CGA_hemap.pdf", width = 6.5, height = 6)
Plot_color_vector(X = data.frame("x"=annot$x, "y"=annot$y), color = colorf(expressed_testis_num), NAME = "FigS6B_TSNE_CGA_hemap.pdf")
dev.off()

# Only hemap MM
load("Hemap_MM_subtypes.Rdata")

# plot CGA to hemap:
expressed_testis_num=colSums(profile2[rownames(profile2)%in%unique(t.df$gene),])
feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0_Antigens"
feat_class[expressed_testis_num>=1&expressed_testis_num<=2]="1to2_Antigens"
feat_class[expressed_testis_num>=3&expressed_testis_num<=4]="3to4_Antigens"
feat_class[expressed_testis_num>=5&expressed_testis_num<=6]="5to6_Antigens"
feat_class[expressed_testis_num>=7]="over7_Antigens"


colorf=colorRamp2(breaks = c(0, 2, 4, 6, 7),colors = c("grey75", "#f4acac", "#f77777", "#e23030", "#800000"))

pdf("FigS6F_TSNE_CGA_hemap_MM_subtypes.pdf", width = 6.5, height = 6)
Plot_color_vector(X = data.frame("x"=coordinates.subtype$x, "y"=coordinates.subtype$y), color = colorf(expressed_testis_num[match(coordinates.subtype$ID, names(expressed_testis_num))]), NAME = "TSNE_CGA_hemap.pdf", SIZE = 4)
dev.off()