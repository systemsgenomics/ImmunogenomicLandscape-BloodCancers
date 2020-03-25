# AML correlations below
GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
library(data.table)
library(parallel)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/DUFVA_TCGA_AML_FM_meth.Rdata"))
fm=fm[!(!grepl("^B:GNAB:*.*chr*.*y_n_somatic", rownames(fm))&grepl("^B:GNAB:", rownames(fm))),]

# add methylation global in heatmap:
meth_global=fread("450array_methylation_GMM_profile_3m.txt", data.table=F)
meth_global[,1]=gsub("\\.", "-", meth_global[,1])

fm=fm[,match(meth_global$V1, colnames(fm))]

meth=t(meth_global[,-1])

hlaII=as.numeric(fm["N:SAMP:HLAIIScore",])

clusters=as.character(fm["C:SAMP:cancermap_cluster",])
clusters2=as.character(fm["C:SAMP:cancermap_cluster",])
clusters2[clusters%in%3&hlaII<8]="3 CIITA methylated"

clu=get.logical(list(clusters2))

a=plot.boxplot("high_pct", logicalVectors = clu, names.lv = names(clu), data = meth, spread = F, ylab = "% High Methylated", color.v = c("#e8c5e4", "#ab3724", "#93a891","#93a891", "#917e99", "#7cb7d8", "#f18aad", "#95dbb2"))

ggsave(plot = a, filename = "FigS4D_Global_hypermethylation_HLAIIlow_TCGA_clusters.pdf", width = 4, height = 3)

wilcox.test(x = meth["high_pct",clu[[8]]], y = meth["high_pct",!clu[[8]]], alternative = "greater")
0.00173
