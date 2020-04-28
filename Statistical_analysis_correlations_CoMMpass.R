GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))

library(data.table)
library(parallel)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("MM_COMPASS_FM.Rdata"))

genelist=gsub(".:SAMP:", "", grep("Score", rownames(fm), value=T))

name=""

# nested correlations:
extrafeatures=c(grep(":CLIN:|:GNAB:|:CNVR:|:SAMP:", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl("vs|and", extrafeatures)]


# HLA Score associations:
l.regulon.gene=regulon.feats(fm, genelist)
results=pairwise.correlation(l.regulon.gene[grep("HLA", names(l.regulon.gene))], fm, extrafeatures[!(grepl("Score",extrafeatures))],filter.val = 5, cores=6, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS4_CoMMpass_MM_HLA_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Cytolytic Score associations:
d=fread("costim_ligands_final.txt", data.table = F)
genelist=unique(d[,1])
l.regulon.gene=regulon.feats(fm, genelist, filtertypes = "B:GEXP")

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores=6, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS5_CoMMpass_MM_costim_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# CGA associations:
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
genelist=c(unique(t.df[,1]))

l.regulon.gene=regulon.feats(fm, genelist)

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures, filter.val = 5, cores=6, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS6_CoMMpass_MM_antigen_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# n.CGA associations:
genelist=gsub(".:SAMP:", "", grep("SAMP.*.nCGA|SAMP.*._CGA", rownames(fm), value=T))
l.regulon.gene=regulon.feats(fm, genelist)

results=pairwise.correlation(l.regulon.gene, fm, c(extrafeatures[!(grepl("catCGA|cancermap",extrafeatures))]),filter.val = 5, cores=5, adjust.method = "BH")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS6_CoMMpass_MM_nCGA_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Subtype assocuated mutations
genelist=gsub(".:SAMP:", "", grep("SAMP.*cancermap", rownames(fm), value=T))
genelist=genelist[!grepl("vs|and", genelist)]

l.regulon.gene=regulon.feats(fm, genelist)

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures[!grepl("SAMP.*cancermap|_vs_|_and_", extrafeatures)],filter.val = 5, cores=6, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("Subtyping_CoMMpass_MM_cancermap_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************