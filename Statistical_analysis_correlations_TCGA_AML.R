GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))

library(data.table)
library(parallel)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("DUFVA_TCGA_AML_FM_meth.Rdata"))
fm=fm[!(!grepl("^B:GNAB:*.*chr*.*y_n_somatic", rownames(fm))&grepl("^B:GNAB:", rownames(fm))|grepl("OS.months..3.31.12|TCGA.Patient.ID|X.PB.Blast|X.BM.Blast", rownames(fm))),]

genelist=grep(":SAMP:CytolyticScore|:SAMP:HLAIScore|:SAMP:HLAIIScore", rownames(fm), value=T)

# nested correlations:
extrafeatures=c(grep("^B:GNAB|GENETICS|^B:CLIN:|^N:CLIN:|CNVR|MLL\\|FISH_test_component|MUTATION_RATE|GENOME_FRAGMENTATION_RATE|B:SAMP:cancermap_cluster_.$", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl(",|GSVA", extrafeatures)]

FEAT="refGene_coding_transcripts.bed"
bed=grep("N:CNVR", rownames(fm), value=T)
bed=cbind(do.call(rbind, strsplit(bed, ":"))[,c(4,5,6)], bed)
cnv_annot=get.feat.within.region(bed, FEAT)

FEAT="refGene_non_coding_transcripts.bed"
bed=grep("N:CNVR", rownames(fm), value=T)
bed=cbind(do.call(rbind, strsplit(bed, ":"))[,c(4,5,6)], bed)
cnv_annot2=get.feat.within.region(bed, FEAT)

cnv_annot=rbind(cnv_annot, cnv_annot2)
cnv_annot=cnv_annot[!grepl("Gistic", cnv_annot[,1]),]


# Cytolytic Score associations:
l.regulon.gene=regulon.feats(fm, genelist, cnv_annot)
results=pairwise.correlation(l.regulon.gene[grep("CytolyticScore", names(l.regulon.gene))], fm, extrafeatures, filter.val = 5, cores=10, adjust.method = "BH", use.fisher = T, fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS3_TCGA_AML_cytScore_correlations.tsv", sep ="\t")
#***********************************************************************


# HLA Score associations:
extrafeatures=c(grep("^B:GNAB|GENETICS|^B:CLIN:|^N:CLIN:|CNVR|MLL\\|FISH_test_component|MUTATION_RATE|GENOME_FRAGMENTATION_RATE|B:SAMP:cancermap_cluster_.$|^N:METH:.*.chr.*.|N:GEXP:", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl(",|GSVA", extrafeatures)]

l.regulon.gene=regulon.feats(fm, genelist, cnv_annot)
results=pairwise.correlation(l.regulon.gene[grep("HLA", names(l.regulon.gene))], fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH", use.fisher = T, fisher.alternative = "greater")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS4_TCGA_AML_HLA_correlations.tsv", sep ="\t")
#***********************************************************************


# Immunomodulatory associations:
fm=fm[!(grepl("METH", rownames(fm))&!grepl("mean", rownames(fm))),]

d=fread("costim_ligands_final.txt", data.table = F)
genelist=unique(d[,1])

extrafeatures=c(grep("^B:GNAB|GENETICS|^B:CLIN:|^N:CLIN:|CNVR|MLL\\|FISH_test_component|MUTATION_RATE|GENOME_FRAGMENTATION_RATE|B:SAMP:cancermap_cluster_.$", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl(",|GSVA", extrafeatures)]

l.regulon.gene=regulon.feats(fm, genelist, cnv_annot)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS5_TCGA_AML_costim_LR_correlations.tsv", sep ="\t")
#***********************************************************************


# CIITA associations:
fm=fm[!(grepl("METH", rownames(fm))&!grepl("mean", rownames(fm))),]

genelist="CIITA"

extrafeatures=c(grep("^B:GNAB|GENETICS|^B:CLIN:|^N:CLIN:|CNVR|MLL\\|FISH_test_component|MUTATION_RATE|GENOME_FRAGMENTATION_RATE|B:SAMP:cancermap_cluster_.$", rownames(fm), value=T))
extrafeatures=extrafeatures[!grepl(",|GSVA", extrafeatures)]

l.regulon.gene=regulon.feats(fm, genelist, cnv_annot)
results=pairwise.correlation(l.regulon.gene, fm, extrafeatures,filter.val = 5, cores=10, adjust.method = "BH")
results2=filter.pairwise.res(results)

#***********************************************************************
fwrite(results2, "TableS4_TCGA_AML_CIITA_correlations.tsv", sep ="\t")
#***********************************************************************