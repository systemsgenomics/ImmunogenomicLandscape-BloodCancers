GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))

library(data.table)
library(parallel)


setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("GSE98588_fm.Rdata"))
name=""

# fm=fm[,fm["B:SAMP:COO_byGEP_GCB",]==1]
# name="_GCB"
# fm=fm[,fm["B:SAMP:COO_byGEP_ABC",]==1]
# name="_ABC"

genelist=gsub(".:SAMP:", "", grep(":SAMP:.*.Score", rownames(fm), value=T))
genelist=genelist[!grepl("and|vs", genelist)]

# nested correlations:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:CNVR:|^B:SAMP:|^N:SAMP:", rownames(fm), value=T))
extrafeatures=extrafeatures[!(grepl("Score",extrafeatures)&grepl("B:", extrafeatures))]
extrafeatures=extrafeatures[!grepl("CGA|Cytolytic",extrafeatures)]

cnv_annot=fread("41591_2018_16_MOESM8_ESM_CNV_ANNOT.txt", data.table=F)


# HLA Score associations:
l.regulon.gene=regulon.feats(fm, genelist, cnv_annot)
results=pairwise.correlation(l.regulon.gene[grep("HLA", names(l.regulon.gene))], fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS4_GSE98588_DLBCL_HLA_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Cytolytic Score associations:
results=pairwise.correlation(l.regulon.gene[grep("CytolyticScore", names(l.regulon.gene))], fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS3_GSE98588_DLBCL_CytScore_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# Immunomodulatory associations:
d=fread("costim_ligands_final.txt", data.table = F)
genelist=unique(d[,1])
l.regulon.gene=regulon.feats(fm, genelist, cnv_annot, filtertypes = "B:GEXP")

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS5_GSE98588_DLBCL_costim_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************


# CGA associations:
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
genelist=c(gsub("N:SAMP:|B:SAMP:", "", grep("nCGA|catCGA", rownames(fm), value=T)), unique(t.df[,1]))
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|B:GNAB|^B:CNVR:|^B:SAMP:|^N:SAMP:", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl("nCGA|catCGA|cat2CGA|vs|and", extrafeatures)]

l.regulon.gene=regulon.feats(fm, genelist, cnv_annot, filtertypes="N:GEXP")

fm=fm[,!colnames(fm)%in%"DLBCL_LS2208"] # testicular DLBCL removed

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, paste0("TableS6_GSE98588_DLBCL_antigen_correlations", name, ".tsv"), sep ="\t")
#***********************************************************************