library(FDb.InfiniumMethylation.hg19)

# TCGA AML
hm450="/research/groups/biowhat_share/public_data/TCGA/firehose/Meth_data/gdac.broadinstitute.org_LAML.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/LAML.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"

tnf=c("TNFRSF1A", "TNFRSF1B", "TNFRSF10A", "TNFRSF10B","TNFRSF10C","TNFRSF10D","RIPK1","FADD")
hlagenes=c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1", "CIITA")
co.stim=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/costim_ligands_final.txt", data.table = F)
t.df = read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/t.antigen_df.txt", stringsAsFactors=F, header=T)

genelist=c(t.df[,1], co.stim[,1], "CIITA", hlagenes,tnf)

# first extract the genes from the meth data:
source("~/git_home/common_scripts/raw_data_processing/hm450_meth_analysis.R")

# get the objects back


# better annotations
meth_annot_extra=get450KProbeMapping(meth_annot[,1])

save(list = c("meth", "meth_annot_extra"), file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/TCGA_AML_meth_probes_genelist.Rdata")


# TCGA DLBCL
hm450="/research/groups/biowhat_share/public_data/TCGA/firehose/Meth_data/gdac.broadinstitute.org_DLBC.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/DLBC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"

# first extract the genes from the meth data:
source("~/git_home/common_scripts/raw_data_processing/hm450_meth_analysis.R")

# get the objects back
x=get.hm450.data(genelist, hm450)
meth=x[[1]]
meth_annot=x[[2]]

# better annotations
meth_annot_extra=get450KProbeMapping(meth_annot[,1])

save(list = c("meth", "meth_annot_extra"), file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/TCGA_DLBCL_meth_probes_genelist.Rdata")







library(FDb.InfiniumMethylation.hg19)

# TCGA AML
hm450="/research/groups/biowhat_share/public_data/TCGA/firehose/Meth_data/gdac.broadinstitute.org_LAML.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/LAML.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"

hlagenes=c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1", "CIITA")
co.stim=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/costim_ligands_final.txt", data.table = F)
t.df = read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/t.antigen_df.txt", stringsAsFactors=F, header=T)

parsemeth=data.table::fread(hm450, data.table=F)
meth_annot=cbind(parsemeth[-1,1], parsemeth[-1,3], parsemeth[-1,4], parsemeth[-1,5])

vals=seq(2, ncol(parsemeth)-1, 4)

meth=data.matrix(parsemeth[2:dim(parsemeth)[1], vals])
rownames(meth)=meth_annot[,1]
colnames(meth)=substr(parsemeth[1,vals], 1, 15)

# remove if half is NA:
meth_annot=meth_annot[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]
meth=meth[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]

# better annotations
meth_annot_extra=get450KProbeMapping(meth_annot[,1])

save(list = c("meth", "meth_annot_extra"), file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/TCGA_AML_meth_probes_genelist_all.Rdata")


# TCGA DLBCL
hm450="/research/groups/biowhat_share/public_data/TCGA/firehose/Meth_data/gdac.broadinstitute.org_DLBC.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015110100.0.0/DLBC.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"

# first extract the genes from the meth data:
source("~/git_home/common_scripts/raw_data_processing/hm450_meth_analysis.R")

# get the objects back
parsemeth=data.table::fread(hm450, data.table=F)
meth_annot=cbind(parsemeth[-1,1], parsemeth[-1,3], parsemeth[-1,4], parsemeth[-1,5])

vals=seq(2, ncol(parsemeth)-1, 4)

meth=data.matrix(parsemeth[2:dim(parsemeth)[1], vals])
rownames(meth)=meth_annot[,1]
colnames(meth)=substr(parsemeth[1,vals], 1, 15)

# remove if half is NA:
meth_annot=meth_annot[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]
meth=meth[!rowSums(is.na(meth))>dim(meth)[2]*0.5,]


# better annotations
meth_annot_extra=get450KProbeMapping(meth_annot[,1])

save(list = c("meth", "meth_annot_extra"), file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/TCGA_DLBCL_meth_probes_genelist_all.Rdata")
