GIT_HOME="/research/users/ppolonen/git_home/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/pathway_analysis/functions.GSEA.R"))

# FigureS2 A-B related analysis
setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# annotations
annot=get(load("Hemap_immunology_Annotations.Rdata"))

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/data9544_with_gene_symbols.RData")
data = t(data[rownames(data)%in%annot$GSM.identifier..sample.,])

subclass2=get.logical(annovector = list(annot$subclasses), filterv = annot$CELLS_SORTED==0&annot$Sample.type%in%c("Prolif", "Cancer"))
tbly2=get.logical(annovector = list(annot$tbLY), filterv = annot$CELLS_SORTED==0&annot$Sample.type%in%c("Prolif", "Cancer"))
tbly2=tbly2[!names(tbly2)%in%"Lymphoma_BCL"]

list_cancers=lapply(c(subclass2, tbly2), list)
list_cancers=unlist(list_cancers, recursive=F)

l=c("AML", "BCL_DLBCL", "CLL", "CML", "Lymphoma_BCL_CHL", "Lymphoma_BCL_MCL", "Lymphoma_BCL_FL", "pre-B-ALL", "Prolif_Myeloproliferative_MDS", "T-ALL")
logicalvectors=list_cancers[names(list_cancers)%in%l]

#*********************** run GSEA using cytolyticScore as ranking metric on each of the disease: ***********************

GENESETS=file.path(getwd(), "HALLMARKS.gmt")
WD=file.path(getwd(), "GSEA")
OUTDIR=file.path(getwd(), "GSEA")

names(logicalvectors)=gsub("-", "_", names(logicalvectors))

res=parallel::mclapply(seq(logicalvectors), function(i){
  filt=logicalvectors[[i]]
  command=run.GSEA(data=data[,filt], cls.vector = as.numeric(annot$CytolyticScore[filt]), datatype = "N", GENESETS = GENESETS, dataname = paste0(names(logicalvectors)[i]), clsname = "CytolyticScore", WD=WD, OUTDIR=OUTDIR)
  try(system(command))
}, mc.cores=5)

#***********************************************************************************************************************

#************** run GSEA using cytolyticScore as ranking metric on TNK cells vs other normal ******************
cls=annot$immunoNormals%in%c("CD8+Tcell", "NKCell")
filt=!annot$immunoNormals%in%c("")

command=run.GSEA(data=data[,filt], cls.vector = as.numeric(cls[filt]), datatype = "B", GENESETS = GENESETS, dataname = "CD8_NK_pathways", clsname = "CD8_NK", WD=WD, OUTDIR=OUTDIR)
try(system(command))

#***********************************************************************************************************************

# combine results to same matrix:
datp=list.files(path = OUTDIR, pattern = "gsea_report_for_1_*.*.xls|gsea_report_for_feat_pos.*..xls", recursive = T, full.names = T)
datn=list.files(path = OUTDIR, pattern = "gsea_report_for_0_*.*.xls|gsea_report_for_feat_neg.*..xls", recursive = T, full.names = T)

library(parallel)

# filter to contain only significant
pwnames=scan(pipe("cut -f1 HALLMARKS.gmt"), "genesets")

FUN=function(f){
  
  dat=try(read.delim(f, stringsAsFactors=F, header=T), silent = T)
  dat$FDR.q.val[dat$FDR.q.val==0]=0.0001
  # variable.1, variable.2, features (y name), id (x name)
  d=data.frame("features"=gsub("-MSIGDB_HALLMARKS", "", dat$NAME), "variable.1"=dat$NES, "variable.2"=-log10(as.numeric(dat$FDR.q.val)), "id"=gsub("_CytolyticScore_.*.|_pathways_.*.|Prolif_Myeloproliferative_|Lymphoma_BCL_|^/BCL_|/","", gsub(eval(OUTDIR), "", f)))
  
}

pwp=mclapply(datp, FUN, mc.cores=10)
pwp2=mclapply(datn, FUN, mc.cores=10)

dat=do.call(rbind, pwp)
dat2=do.call(rbind, pwp2)

dat_comb=rbind(dat, dat2)

dat_comb$id=factor(gsub("_", "-", dat_comb$id), levels=c("DLBCL","CHL", "FL", "MCL", "CLL", "CML", "AML", "MDS", "pre-B-ALL", "T-ALL", "CD8-NK"))

pdf("FigureS2A.pdf", width = 5.5, height = 6.5)
plot.DotPlot.df(data.plot = dat_comb, cols = c("blue", "white", "red"), fontsize = 8, name.variable.1 = "NES", name.variable.2 = "-log10 FDR", dot.scale = 3)
dev.off()

# plot correlation plot between protein and Pathway
fm=get(load("TCGA_DLBCL_FM_DUFVA.Rdata"))

# run GSVA for certain pathways and add to matrix:
gexp=data.matrix(fm[grepl("N:GEXP", rownames(fm)),])
rownames(gexp)=gsub(":.*.", "", gsub("N:GEXP:", "",rownames(gexp)))

# gene sets:
GENESETS="Combined_pathway_signatures_2017_filtered_robust.gmt"

# get GSVA visualization for the pathways
library(GSVA)
library(parallel)

# Geneset list
Onc.pathways=read.delim(GENESETS, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(GENESETS, sep = '\t'), na.rm = T)), fill = TRUE)

# Make list
listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
B=A[!A==""&!A=="NA"]}, mc.cores=6)
names(listA) <- toupper(Onc.pathways[,1])
viz_scores=gsva(expr = gexp, gset.idx.list = listA, parallel.sz=8, tau=0.25)
rownames(viz_scores)=paste0("N:GSVA:", rownames(viz_scores))
fm.add=rbind(fm, viz_scores)

fm.add=fm.add[!(!grepl("^B:GNAB:*.*chr*.*y_n_somatic", rownames(fm.add))&grepl("^B:GNAB:", rownames(fm.add))),!is.na(fm.add["N:RPPA:PARP1:chr1:226548392:226595780:-:PARP_cleaved",])]

annotdf=data.frame("CytolyticScore"=as.numeric(fm["N:SAMP:CytolyticScore",]), "CASP7"=as.numeric(fm["N:RPPA:CASP7:chr10:115438942:115490662:+:Caspase-7_cleavedD198",]), "IRF1"=as.numeric(fm["N:RPPA:IRF1:chr5:131817301:131826490:-:IRF-1",]))
rownames(annotdf)=colnames(fm)

annotdf=annotdf[!is.na(annotdf$CASP7),]

pdf("FigureS2B.pdf", height = 2.5, width = 2.5)
ggscatter(annotdf, x = "CytolyticScore", y = "CASP7",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "pearson", label.sep = "\n") +
  xlab("Cytolytic Score") +
  ylab("CASP7")

ggscatter(annotdf, x = "CytolyticScore", y = "IRF1",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "pearson", label.sep = "\n") +
  xlab("Cytolytic Score") +
  ylab("IRF1")
dev.off()
