GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))
library(data.table)
library(parallel)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

load("BeatAML_fm.Rdata")
load("BeatAML_fm_annot.Rdata")

filtv=annot$specimenType%in%"Bone Marrow Aspirate"
fm=fm[,filtv]

# correlations against these features:
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|^B:GNAB|^B:SAMP:|^N:SAMP:", rownames(fm), value=T))
extrafeatures=extrafeatures[!(grepl("_and_|_vs_|ic50|cancermap_cluster", extrafeatures)&!grepl("TCGA", extrafeatures))]
extrafeatures=extrafeatures[!(grepl("_and_|_vs_|ic50", extrafeatures))]

genelist=gsub("N:SAMP:|B:SAMP:|:::::|N:GEXP:", "", grep("Score", rownames(fm), value=T))
genelist=genelist[!(grepl("_and_|_vs_", genelist))]

l.regulon.gene=regulon.feats(fm, genelist)


# Cytolytic Score associations:
results=pairwise.correlation(l.regulon.gene[grepl("CytolyticScore", l.regulon.gene)], fm, extrafeatures[!grepl("Cytolytic",extrafeatures)], filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS3_BeatAML_cytscore_correlations.tsv", sep ="\t")
#***********************************************************************


# HLA Score and CIITA associations:
l.regulon.gene=regulon.feats(fm, c(genelist, "CIITA"))
extrafeatures=c(grep("^B:CLIN:|^N:CLIN:|^B:GNAB|^B:SAMP:|^N:SAMP:|N:GEXP", rownames(fm), value=T))
extrafeatures=extrafeatures[!(grepl("__", extrafeatures))]

results=pairwise.correlation(l.regulon.gene[grepl("HLA|CIITA", l.regulon.gene)], fm, c(extrafeatures[!grepl("SAMP:HLA|CIITA",extrafeatures)]), filter.val = 5, cores=10, adjust.method = "BH", fisher.alternative = "greater")
results2=filter.pairwise.res(results)

# filter GEXP-GEXP pairs more:
results2=results2[!(results2$datapairs=="GEXP:GEXP"&as.numeric(results2$FDR)>1e-20),]

#***********************************************************************
fwrite(results2, "TableS4_BeatAML_HLAscore_correlations.tsv", sep ="\t")
#***********************************************************************


# Immunomodulatory associations:
d=fread("costim_ligands_final.txt", data.table = F)
genelist=unique(d[,1])

l.regulon.gene=regulon.feats(fm, genelist)
l.regulon.gene[41]="N:GEXP:ENTPD1"

results=pairwise.correlation(l.regulon.gene, fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS5_BeatAML_coStim_correlations.tsv", sep ="\t")
#***********************************************************************