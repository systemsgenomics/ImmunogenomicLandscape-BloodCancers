# data downloaded 2019 from:
# https://gdc.cancer.gov/about-data/publications/panimmune
# https://gdc.cancer.gov/node/905/
# TCGA pancancer FM

GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "statistics/useful_functions.R"))
library(parallel)

# annotations
clin=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/CoMMpass_IA13_FlatFiles/MMRF_CoMMpass_IA13_PER_PATIENT.csv", data.table=F)
surv=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/CoMMpass_IA13_FlatFiles/MMRF_CoMMpass_IA13_STAND_ALONE_SURVIVAL.csv", data.table=F)

# OMICS data
rna=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_E74GTF_HtSeq_Gene_Counts.txt", data.table=F)
rna2=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_E74GTF_Cufflinks_Gene_FPKM.txt", data.table=F)
rna2=rna2[,-2]

# Mutation data
mut=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_All_Canonical_Variants_ENSG_Mutation_Counts.txt", data.table=F)
maf=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_All_Canonical_NS_Variants.txt", data.table=F)
# maf=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_All_Canonical_Variants.txt", data.table=F)

# CNV
cnv.long=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_CNA_LongInsert_FISH_CN_All_Specimens.txt", data.table=F)

# fusion
fusion=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_TophatFusion_Results.txt", data.table=F)
ig=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_LongInsert_Canonical_Ig_Translocations.txt", data.table=F)
rnaig=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/MMRF_CoMMpass_IA13a_RNAseq_Canonical_Ig_Translocations.txt", data.table=F)

#*************************************** filter and process each data type ******************************************

clin.filt=clin[,!colnames(clin)%in%c("PUBLIC_ID")]
rownames(clin.filt)=clin$PUBLIC_ID
clin.filt=data.frame(clin.filt, stringsAsFactors = F)

clin.filt$R_ISS=as.character(clin.filt$R_ISS)
clin.filt$IMWG_Risk_Class=as.character(clin.filt$IMWG_Risk_Class)

surv.filt=data.frame("STATUS"=0, "OS"=surv$lstalive)
surv.filt$STATUS[!is.na(surv$deathdy)]=1
surv.filt$OS[!is.na(surv$deathdy)]=surv$deathdy[!is.na(surv$deathdy)]
rownames(surv.filt)=surv$public_id

#******************************************************************************************************************

mut.m=data.matrix((table(maf$`ANN[*].GENE`, maf$Sample)>0)*1)
colnames(mut.m)=gsub("_._BM|_._PB", "", colnames(mut.m))
mut.m=mut.m[,!duplicated(colnames(mut.m))]

# # all mutations, is mutated or not:
# mut.m=data.matrix(mut[,-1]>0)*1
# match.gene=intersect(genes$gene_id, mut[,1])
# genes2=genes[match(match.gene, genes$gene_id),]
# mut.m=mut.m[match(match.gene, mut[,1]),]
# mut=mut[match(match.gene, mut[,1]),]
# rownames(mut.m)=genes2$gene_name[match(match.gene, mut[,1])]

# convert to symbol:
library(data.table)
# http://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz same file was used in the analysis, so genes will match
genes <- fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/Homo_sapiens.GRCh37.74.gtf")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )

# the problem is the attributes column that tends to be a collection
# of the bits of information you're actually interested in
# in order to pull out just the information I want based on the 
# tag name, e.g. "gene_id", I have the following function:
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

# this is how to, for example, extract the values for the attributes of interest (here: "gene_id")
genes$gene_id <- unlist(mclapply(genes$attributes, extract_attributes, "gene_id", mc.cores=10))
genes$gene_name <- unlist(mclapply(genes$attributes, extract_attributes, "gene_name", mc.cores=10))
genes=unique(genes[,10:11])
genes=genes[!is.na(genes$gene_name),]

#******************************************************************************************************************

# fusions:
rnaig_filt=rnaig[,grepl("Call", colnames(rnaig))]
colnames(rnaig_filt)=gsub("_Call", "_Ig_translocation", colnames(rnaig_filt))

name=gsub("_._BM|_._PB", "", rnaig$Specimen_ID)

rnaig_m=do.call(rbind, lapply(unique(name), function(n){
  (colSums(rnaig_filt[name%in%n,])>0)*1
}))

rownames(rnaig_m)=unique(name)


ig.filt=ig[,grepl("CALL", colnames(ig))]
colnames(ig.filt)=gsub("_CALL", "_Ig_translocation", colnames(ig.filt))
name=gsub("_._BM|_._PB", "", ig$Study_Visit_iD)
           
ig_m=do.call(rbind, lapply(unique(name), function(n){
  (colSums(ig.filt[name%in%n,])>0)*1
}))

rownames(ig_m)=unique(name)

a=fusion[gsub("_._BM|_._PB", "", fusion$ID)%in%coord$ID[coord$cluster=="5"],]
b=(table(paste(a$left.Gene, a$right.Gene), gsub("_._BM|_._PB", "", a$ID))>0)*1

sort(rowSums(b))

b2=(table(paste(fusion$left.Gene, fusion$right.Gene), gsub("_._BM|_._PB", "", fusion$ID))>0)*1
b3=b2==1

test=do.call(rbind, apply(b3, 1, function(lv2)fisher.2x2(lv1 = colnames(b2)%in%coord$ID[coord$cluster=="5"], lv2, alternative = "greater")))

#******************************************************************************************************************

# GEXP
rna.filt=rna[,-1]
rna.filt2=rna2[,-1]
rna.filt=rna.filt[,grepl("1_BM", colnames(rna.filt))]
rna.filt2=rna.filt2[,grepl("1_BM", colnames(rna.filt2))]

match.gene=intersect(genes$gene_id, rna[,1])
genes2=genes[match(match.gene, genes$gene_id),]
rna.filt=rna.filt[match(match.gene, rna[,1]),]
rna=rna[match(match.gene, rna[,1]),]
rownames(rna.filt)=make.unique(genes2$gene_name[match(match.gene, rna[,1])])

match.gene=intersect(genes$gene_id, rna2[,1])
genes2=genes[match(match.gene, genes$gene_id),]
rna.filt2=rna.filt2[match(match.gene, rna2[,1]),]
rna2=rna2[match(match.gene, rna2[,1]),]
rownames(rna.filt2)=make.unique(genes2$gene_name[match(match.gene, rna2[,1])])
rna.filt2=rna.filt2[,match(colnames(rna.filt), colnames(rna.filt2))]

# filter:
filt=rowSums(edgeR::cpm(rna.filt)>1)>dim(rna.filt)[2]*0.025

# normalize library size
DGE <- edgeR::DGEList(rna.filt[filt,])
DGE <- edgeR::calcNormFactors(DGE)

# voom transform
gexp=limma::voom(DGE, plot=T)$E
colnames(gexp)=gsub("_._BM|_._PB", "", colnames(gexp))

#******************************************************************************************************************
cnv.long.filt=cnv.long[,-1]
rownames(cnv.long.filt)=cnv.long[,1]
cnv.long.filt=cnv.long.filt[,!grepl("percent", colnames(cnv.long.filt))]
colnames(cnv.long.filt)=gsub("SeqWGS_Cp_", "", colnames(cnv.long.filt))
cnv.long.filt=cnv.long.filt[grepl("1_BM", rownames(cnv.long.filt)),]
rownames(cnv.long.filt)=gsub("_._BM|_._PB", "", rownames(cnv.long.filt))

#******************************************************************************************************************
# immunoscores:
library(circlize)

dat_a3=gexp[rownames(gexp)%in%c("B2M",
                                "HLA-A",
                                "HLA-B",
                                "HLA-C"),]

dat3=2^dat_a3+0.01
gm1=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm1)="HLAIScore"

dat_a3=gexp[rownames(gexp)%in%c("HLA-DMA",
                                "HLA-DMB",
                                "HLA-DPA1",
                                "HLA-DPB1",
                                "HLA-DRA",
                                "HLA-DRB1"),]

dat3=2^dat_a3+0.01
gm2=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm2)="HLAIIScore"

dat_a3=gexp[rownames(gexp)%in%c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"),]

dat3=2^dat_a3+0.01
gm3=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm3)="CytolyticScore"

classification1=data.frame(t(rep("medium", length(gm3))))
zscore=as.numeric(scale(t(gm3)))
classification1[zscore>=1]="high"
classification1[zscore<=(-1)]="low"
rownames(classification1)="CytolyticScore" 
colnames(classification1)=colnames(gexp)

classification2=data.frame(t(rep("medium", length(gm1))))
zscore=as.numeric(scale(t(gm1)))
classification2[zscore>=1]="high"
classification2[zscore<=(-1)]="low"
rownames(classification2)="HLAIScore" 
colnames(classification2)=colnames(gexp)

classification3=data.frame(t(rep("medium", length(gm2))))
zscore=as.numeric(scale(t(gm2)))
classification3[zscore>=1]="high"
classification3[zscore<=(-1)]="low"
rownames(classification3)="HLAIIScore" 
colnames(classification3)=colnames(gexp)

classification=data.frame(t(rbind(classification1,classification2,classification3)), stringsAsFactors = F)

# CGA:
t.df = read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/t.antigen_df.txt", stringsAsFactors=F, header=T)
genelist=c(unique(t.df[,1]))

# rank patients by number of testis antigens expressed
expressed_testis_num=data.frame(t(colSums(gexp[rownames(gexp)%in%genelist,]>3))) # good cutoff point in fpkm and cpm, a lot of noise otherwise
# expressed_testis_num=data.frame(t(colSums(log2(rna.filt2[rownames(rna.filt2)%in%genelist,])>3)))

colnames(expressed_testis_num)=colnames(gexp)
rownames(expressed_testis_num)="nCGA"

feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0_Antigens"
feat_class[expressed_testis_num>=1&expressed_testis_num<=2]="1to2_Antigens"
feat_class[expressed_testis_num>=3&expressed_testis_num<=4]="3to4_Antigens"
feat_class[expressed_testis_num>=5&expressed_testis_num<=6]="5to6_Antigens"
feat_class[expressed_testis_num>=7]="over7_Antigens"

rownames(feat_class)="catCGA"

immunoscores=as.data.frame(t(rbind(gm1, gm2,expressed_testis_num)),stringsAsFactors=F)
immunoscores$catCGA=as.character(feat_class)

#****************************************************************************************************************************
#*************************************** Make a feature matrix from each data type ******************************************
#****************************************************************************************************************************
# clustering:
coord=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/cancermap_CoMMpass_12.5pct_genes_BH-SNE_mean-shift_BW1.txt", header=T, stringsAsFactors = F)
peaks=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/cancermap_CoMMpass_12.5pct_genes_BH-SNE_mean-shift_BW1_cluster_centroids.txt", header=T, stringsAsFactors = F)

coord$cluster[coord$cluster%in%c(21)]="CGA_Prolif"

clust=make.features(data.frame("cancermap_cluster"=as.character(coord$cluster), stringsAsFactors = F), datatype="SAMP", make.pairwise = F)
colnames(clust)=coord$ID

# clustering larger:
coord.2=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/cancermap_CoMMpass_12.5pct_genes_BH-SNE_mean-shift_BW2.txt", header=T, stringsAsFactors = F)
peaks.2=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/MM_compass/cancermap_CoMMpass_12.5pct_genes_BH-SNE_mean-shift_BW2_cluster_centroids.txt", header=T, stringsAsFactors = F)

coord2=coord.2

# checked the genetics and combined the subtypes as larger subtypes:
coord2$cluster[coord.2$cluster%in%c(6)]="MAF_Ig"
coord2$cluster[coord.2$cluster%in%c(3)]="WHSC1_FGFR3_Ig"
coord2$cluster[coord.2$cluster%in%c(1,2)]="CCND1_Ig"
coord2$cluster[coord.2$cluster%in%c(4)]="Hyperdiploid"
coord2$cluster[coord.2$cluster%in%c(5)]="Hyperdiploid_amp1q"
coord2$cluster[coord.2$cluster%in%c(7)]="TRAF3_Aberrated"

clust2=make.features(data.frame("cancermap_subtypes_"=as.character(coord2$cluster), stringsAsFactors = F), datatype="SAMP")
colnames(clust2)=coord$ID


# clinical and surv
clindatfm=make.features(df = clin.filt, datatype="CLIN", prefix="", make.pairwise = F)
survdatfm=make.features(df = surv.filt, datatype="CLIN", prefix="", make.pairwise = F)

# OMICS
gexpfm=make.features(df = data.frame(t(gexp), check.names = F), datatype="GEXP", prefix="", make.pairwise = F)
colnames(gexpfm)=gsub("_._BM|_._PB", "", colnames(gexp))

# Genetics
cnvfm=make.features(df = data.frame(cnv.long.filt, check.names = F), datatype="CNVR", prefix="", make.pairwise = F)
colnames(cnvfm)=gsub("_._BM|_._PB", "", rownames(cnv.long.filt))

mutfm=make.features(as.data.frame(t(mut.m), check.names = F), datatype="GNAB", prefix="")
colnames(mutfm)=gsub("_._BM|_._PB", "", colnames(mut.m))

# immunology
immunoscoresfm=make.features(immunoscores, datatype="SAMP", prefix="")
colnames(immunoscoresfm)=gsub("_._BM|_._PB", "", colnames(gexp))

immunoscoresfm2=make.features(classification, datatype="SAMP", prefix="")
colnames(immunoscoresfm2)=gsub("_._BM|_._PB", "", colnames(gexp))

rnaig_fm=make.features(data.frame(rnaig_m), datatype="CNVR", prefix="")
ig_fm=make.features(data.frame(ig_m), datatype="CNVR", prefix="")


annot=cbind(clin.filt,surv.filt[match(rownames(clin.filt), rownames(surv.filt)),], "cluster"=coord$cluster[match(rownames(clin.filt),coord$ID)], "subtype"=coord2[match(rownames(clin.filt),coord2$ID),])
l.fm=list(clust,clust2, clindatfm, survdatfm, gexpfm, cnvfm, mutfm, immunoscoresfm, immunoscoresfm2, rnaig_fm, ig_fm)

library(data.table)
fm=rbindlist(l.fm, use.names=T, fill=T)

fm=data.frame(fm, stringsAsFactors=F, check.names = F)
rownames(fm)=unlist(lapply(l.fm, rownames))

save(gexp, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/MM_COMPASS_GEXP.Rdata")
save(fm, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/MM_COMPASS_FM.Rdata")
save(annot, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_COMPASS/MM_COMPASS_ANNOT.Rdata")