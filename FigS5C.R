GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "pathway_analysis/functions.GSEA.R"))
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "statistics/useful_functions.R"))
library(Matrix)
library(Seurat)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(parallel)
library(ggplot2)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")


# analyze FIMM and Galen AML and find costim associations to certain clusters.
co.stim=data.table::fread("costim_ligands_final.txt", data.table = F)[,c(1,3,5)]
co.stimR=data.table::fread("costim_ligands_final.txt", data.table = F)
co.stimR.f=c(unlist(strsplit(co.stimR$`Receptor gene`, ", ")), "LAG3")
co.stim=rbind(co.stim, c("FGL1", "LAG3", "Inhibitory"))

load("AML_Galen_scRNA.Rdata")
galen=scmat

load("HCA_scRNA.Rdata")
HCA=scmat

load("FIMM_AML_scRNA.Rdata")

# plot from each dataset all significant genes:
subtype.order=c("MDS-like", "Progenitor-like", "Monocyte-like", "Monocyte-like-MLL", "CEBPA", "RUNX1-RUNX1T1", "CBFB-MYH11", "PML-RARA")

fimm=scmat[,scmat[["SingleR.label"]][,1]%in%c("HSC", "MPP", "GMP", "CMP", "MEP", "Monocytes", "Erythrocytes")]
galen.sub=galen[,galen[["SingleR.label"]][,1]%in%c("HSC", "MPP", "GMP", "CMP", "MEP", "Monocytes", "Erythrocytes")]
hca=HCA[,HCA[["SingleR.label"]][,1]%in%c("HSC", "MPP", "GMP", "CMP", "MEP", "Monocytes", "Erythrocytes")]

Idents(fimm)=fimm[["SingleR.label"]][,1]
Idents(galen.sub)=factor(galen.sub[["SingleR.label"]][,1], levels=c("HSC", "MPP", "GMP", "CMP", "MEP", "Monocytes", "Erythrocytes"))
Idents(hca)=hca[["SingleR.label"]][,1]

# do DE gene tests:
DE.FIMM=FindAllMarkers(object = fimm, features = genelist[genelist%in%rownames(scmat)], only.pos = T, logfc.threshold = 0.15)
DE.galen=FindAllMarkers(object = galen.sub, features = genelist[genelist%in%rownames(galen)], only.pos = T, logfc.threshold = 0.15)
DE.HCA=FindAllMarkers(object = hca, features = genelist[genelist%in%rownames(HCA)], only.pos = T, logfc.threshold = 0.15)

all=rbind(DE.FIMM, DE.galen, DE.HCA)
all=all[order(match(all$cluster, c("HSC", "MPP", "GMP", "CMP", "MEP", "Monocytes", "Erythrocytes"))),]

a=table(all$cluster, all$gene)
genelist.filt=c(colnames(a[,apply(a, 2, function(v)any(v>1))]), "C10orf54")

# these were selected for the figure:
genelist.filt=c("CD34", "CLEC2B", "TNFRSF14", "CD84", "VSIR", "C10orf54", "CD68", "CD48","CD86","ENTPD1") 

pdf("FigS5C_VISTA_FIMM_dotplot.pdf", width = 4.5, height = 2.5)
plot.DotPlot(fimm, features = genelist.filt[genelist.filt%in%rownames(scmat)], cols=c("grey75", "red"), scale.max = 50, scale.min = 10, dot.scale = 4)
plot.DotPlot(galen.sub, features = genelist.filt[genelist.filt%in%rownames(galen)], cols=c("grey75", "red"), scale.max = 50, scale.min = 10, dot.scale = 4)
plot.DotPlot(hca, features = genelist.filt[genelist.filt%in%rownames(HCA)], cols=c("grey75", "red"), scale.max = 50, scale.min = 10, dot.scale = 4)
dev.off()
