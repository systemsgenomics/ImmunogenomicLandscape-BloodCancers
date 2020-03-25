GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "statistics/useful_functions.R"))
library(parallel)
source("/research/users/ppolonen/git_home/common_scripts/featurematrix/compute.pairwise.R")
source("/research/users/ppolonen/git_home/common_scripts/statistics/functions_statistics.R")

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/")

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/TCGA_AML_meth_probes_genelist.Rdata")

fm=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/TCGA_AML_FM_DUFVA.Rdata"))
fm=fm[!(!grepl("^B:GNAB:*.*chr*.*y_n_somatic", rownames(fm))&grepl("^B:GNAB:", rownames(fm))),!is.na(fm[1,])]

head(meth_annot_extra)

rownames(meth)=apply(cbind(meth_annot_extra[,3],meth_annot_extra[,1], meth_annot_extra[,2]), 1, paste, collapse=":")

methfm=data.matrix(make.features(data.frame(t(meth), check.names = F), datatype="METH"))
colnames(methfm)=colnames(meth)

methfm=methfm[,match(colnames(fm), colnames(methfm))]
colnames(methfm)=colnames(fm)

# exclude probes with low variance:
f=apply(methfm, 1, sd, na.rm=T)<0.1 # low SD
methfm=methfm[!f,]
meth_annot_extra=meth_annot_extra[!f,]

fm2=rbind(methfm, fm[grepl("GEXP", rownames(fm)),])

l.regulon.gene=regulon.feats(fm2, unique(meth_annot_extra[,3]))
l.regulon.gene=l.regulon.gene[sapply(l.regulon.gene, function(l)any(grepl("GEXP", l))&any(grepl("METH", l)))]

res=pairwise.correlation(l.regulon.gene=l.regulon.gene, fm = fm2, prettyNumbers = F)
res_neg=res[res$datapairs=="GEXP:METH"&res$p<0.05&res$cor<=0,]
res_pos=res[res$datapairs=="GEXP:METH"&res$p<0.05&res$cor>0,]

ind=rownames(methfm)%in%res_pos$featureB
methmean_pos=aggregate_by_mean(t(methfm[ind,]), xs = meth_annot_extra[ind,3])
colnames(methmean_pos)=paste(colnames(methmean_pos), "mean_positive_cor", sep=":")
ind=rownames(methfm)%in%res_neg$featureB
methmean_neg=aggregate_by_mean(t(methfm[ind,]), xs = meth_annot_extra[ind,3])
colnames(methmean_neg)=paste(colnames(methmean_neg), "mean_negative_cor", sep=":")

methmean=cbind(methmean_pos, methmean_neg)

methmeanfm=data.matrix(make.features(data.frame(methmean, check.names = F), datatype="METH"))
colnames(methmeanfm)=colnames(fm)

fm=rbind(fm, methmeanfm)

save(fm, file="DUFVA_TCGA_AML_FM_meth.Rdata")


setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/")

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/TCGA_DLBCL_meth_probes_genelist.Rdata")

fm=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/TCGA_DLBCL_FM_DUFVA.Rdata"))

head(meth_annot_extra)

rownames(meth)=apply(cbind(meth_annot_extra[,3],meth_annot_extra[,1], meth_annot_extra[,2]), 1, paste, collapse=":")

methfm=data.matrix(make.features(data.frame(t(meth), check.names = F), datatype="METH"))
colnames(methfm)=colnames(meth)

methfm=methfm[,match(colnames(fm), colnames(methfm))]
colnames(methfm)=colnames(fm)

# exclude probes with low variance:
f=apply(methfm, 1, sd, na.rm=T)<0.1 # low SD
methfm=methfm[!f,]
meth_annot_extra=meth_annot_extra[!f,]

fm2=rbind(methfm, fm[grepl("GEXP", rownames(fm)),])

l.regulon.gene=regulon.feats(fm2, unique(meth_annot_extra[,3]))
l.regulon.gene=l.regulon.gene[sapply(l.regulon.gene, function(l)any(grepl("GEXP", l))&any(grepl("METH", l)))]

library(parallel)
source("/research/users/ppolonen/git_home/common_scripts/featurematrix/compute.pairwise.R")
source("/research/users/ppolonen/git_home/common_scripts/statistics/functions_statistics.R")

res=pairwise.correlation(l.regulon.gene=l.regulon.gene, fm = fm2, prettyNumbers = F)
res_neg=res[res$datapairs=="GEXP:METH"&res$p<0.05&res$cor<0,]
res_pos=res[res$datapairs=="GEXP:METH"&res$p<0.05&res$cor>0,]

ind=rownames(methfm)%in%res_pos$featureB
methmean_pos=aggregate_by_mean(t(methfm[ind,,drop=F]), xs = meth_annot_extra[ind,3])
colnames(methmean_pos)=paste(colnames(methmean_pos), "mean_positive_cor", sep=":")
ind=rownames(methfm)%in%res_neg$featureB
methmean_neg=aggregate_by_mean(t(methfm[ind,,drop=F]), xs = meth_annot_extra[ind,3])
colnames(methmean_neg)=paste(colnames(methmean_neg), "mean_negative_cor", sep=":")

methmean=cbind(methmean_pos, methmean_neg)

methmeanfm=data.matrix(make.features(data.frame(methmean, check.names = F), datatype="METH"))
colnames(methmeanfm)=colnames(fm)

fm=rbind(fm, methmeanfm)

save(fm, file="DUFVA_TCGA_DLBCL_FM_meth.Rdata")
