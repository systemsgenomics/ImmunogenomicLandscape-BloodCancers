GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))

library(data.table)
library(parallel)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")


fm=get(load("REDDY_DLBCL_fm.Rdata"))
name=""

# fm=fm[,fm["B:CLIN:ABC_GCB_RNAseq_GCB",]==1]
# name="_GCB"
#  
# fm=fm[,fm["B:CLIN:ABC_GCB_RNAseq_ABC",]==1]
# name="_ABC"

genelist=gsub(".:SAMP:", "", grep("Score", rownames(fm), value=T))

fm.rm=fm[grepl("CNVR", rownames(fm)),]
filt=rowSums(fm.rm>0.7|fm.rm<(-0.7), na.rm = T)<dim(fm.rm)[2]*0.025

extrafeatures=c(grep("CLIN|GNAB|SAMP:", rownames(fm), value=T), rownames(fm.rm)[filt])
extrafeatures=extrafeatures[!(grepl("Score",extrafeatures)&grepl("B:", extrafeatures))]
extrafeatures=extrafeatures[!grepl("CGA|Cytolytic|pos_Cells|neg_Cells",extrafeatures)]


# HLA Score associations:
l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene[grep("HLA", names(l.regulon.gene))], fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS4_Reddy_DLBCL_HLA_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Cytolytic Score associations:
results=pairwise.correlation(l.regulon.gene[grep("CytolyticScore", names(l.regulon.gene))], fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater", use.wilcox = T)
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS3_Reddy_DLBCL_CytScore_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Immunomodulatory associations:
d=fread("costim_ligands_final.txt", data.table = F)
genelist=unique(d[,1])
l.regulon.gene=regulon.feats(fm, genelist, filtertypes = "B:GEXP")

results=pairwise.correlation(l.regulon.gene = l.regulon.gene, fm, extrafeatures, filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS5_Reddy_DLBCL_costim_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************
