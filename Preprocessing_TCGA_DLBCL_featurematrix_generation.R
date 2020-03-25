setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_DLBCL/")

source("~/git_home/common_scripts/featurematrix/functions_generate_fm.R")

# fm data and clusters
fm = read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/dlbc.newMerge.all.31012017.tsv", stringsAsFactors=F, row.names = 1, header=T)
annot = read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/TCGA_DLBCL_CLUSTERS.tsv", stringsAsFactors=F, header=T)
fm[grepl("N:GEXP:", rownames(fm)),]=log2(data.matrix(fm[grepl("N:GEXP:", rownames(fm)),])+0.01)
colnames(fm)=substr(gsub("\\.", "-", colnames(fm)), 1, 15)
fm=fm[,1:48] # excluding normals

# GSVA data
# dat=log2(data.matrix(fm[grepl("N:GEXP:", rownames(fm)),])+0.01)
# dat=dat[,!apply(dat, 2, function(v)any(is.na(v)))]
# rownames(dat)=gsub("N:GEXP:|:chr*.*", "", rownames(dat))
# colnames(dat)=gsub("\\.", "-", colnames(dat))
# save(dat, file="/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/TCGA_DLBCL/TCGA_DLBCL_RNASEQ.Rdata")

load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/TCGA_DLBCL/GSVA/TCGA_DLBCL_all_samples_Combined_pathway_drug_signatures_2017_GSVA.Rdata")
rownames(gsva_es)=gsub(" ", "_", rownames(gsva_es))
gsva_es=gsva_es[,match(colnames(fm), colnames(gsva_es))]
colnames(gsva_es)=colnames(fm)

setwd("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/")

# cell fractions
load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/MCP_counter_results.Rdata")
load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_DLBCL/CIBERSORT_results.Rdata")

ExampleEstimates=ExampleEstimates[,match(colnames(fm), colnames(ExampleEstimates))]
colnames(ExampleEstimates)=colnames(fm)

results_fm=results_fm[,match(colnames(fm), colnames(results_fm))]
colnames(results_fm)=colnames(fm)

# Clusters
clusters_cancermap=annot$cluster

# annotated and cancermap clusters
cluster_cancermap=FUN_MAKE_ALL(clusters_cancermap, "cancermap_cluster", clusters_cancermap, 0)
colnames(cluster_cancermap)=substr(gsub("\\.", "-", annot$ID), 1, 15)

# cancermap
cluster_cancermap=cluster_cancermap[,match(colnames(fm), colnames(cluster_cancermap))]
colnames(cluster_cancermap)=colnames(fm)

# categorical
class_cancermap=FUN_MAKE_CATEGORICAL(annot$cluster, "cancermap_cluster")
colnames(class_cancermap)=substr(gsub("\\.", "-", annot$ID), 1, 15)
class_cancermap=class_cancermap[,match(colnames(fm), colnames(class_cancermap)), drop=F]
colnames(class_cancermap)=colnames(fm)

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#*********************************** compute geometric mean ********************************
gexp=data.matrix(fm[grepl("N:GEXP:", rownames(fm)),])
gexp=gexp[,!apply(gexp, 2, function(v)any(is.na(v)))]

dat_a=gexp[grep("GZMA|PRF1|GNLY|GZMH|GZMM", rownames(gexp)),]
dat=2^dat_a+0.01
rownames(dat)
gm=log2(t(apply(dat, 2, gm_mean, zero.propagate = F)))
rownames(gm)="N:SAMP:CytolyticScore"

classification1=data.frame(t(rep("medium", length(gm))))
zscore=as.numeric(scale(t(gm)))
classification1[zscore>=1]="high"
classification1[zscore<=(-1)]="low"
rownames(classification1)="CytolyticScore" 
colnames(classification1)=colnames(gexp)

# also HLA
dat_a2=gexp[grep("N:GEXP:B2M:|N:GEXP:HLA-A:|N:GEXP:HLA-B:|N:GEXP:HLA-C:", rownames(gexp)),]

dat2=2^dat_a2+0.01
rownames(dat2)
gm2=log2(t(apply(dat2, 2, gm_mean, zero.propagate = F)))
rownames(gm2)="N:SAMP:HLAIScore"

classification2=data.frame(t(rep("medium", length(gm))))
zscore=as.numeric(scale(t(gm2)))
classification2[zscore>=1]="high"
classification2[zscore<=(-1)]="low"
rownames(classification2)="HLAIScore" 
colnames(classification2)=colnames(gexp)

# also HLAII
dat_a3=gexp[grep("N:GEXP:HLA-DMA:|N:GEXP:HLA-DMB:|N:GEXP:HLA-DPA1:|N:GEXP:HLA-DPB1:|N:GEXP:HLA-DRA:|N:GEXP:HLA-DRB1:", rownames(gexp)),]

dat3=2^dat_a3+0.01
rownames(dat3)
gm3=log2(t(apply(dat3, 2, gm_mean, zero.propagate = F)))
rownames(gm3)="N:SAMP:HLAIIScore"

classification3=data.frame(t(rep("medium", length(gm))))
zscore=as.numeric(scale(t(gm3)))
classification3[zscore>=1]="high"
classification3[zscore<=(-1)]="low"
rownames(classification3)="HLAIIScore" 
colnames(classification3)=colnames(gexp)

classification=data.frame(t(classification1), t(classification2), t(classification3), check.names = F, stringsAsFactors = F)

cat_immunoscores=make.features(df = classification, datatype = "SAMP", make.pairwise = F)

cat_feats=data.frame(t(fm[grep("^C:CLIN|^C:SAMP", rownames(fm)),]) ,stringsAsFactors = F, check.names = F)
cat_feats=make.features(df = cat_feats, datatype = "SAMP", make.pairwise = F)
rownames(cat_feats)=gsub("B:SAMP:|:::::","", rownames(cat_feats))
rownames(cat_feats)=gsub("^C","B", rownames(cat_feats))

cat_clin=make.features(df = classification, datatype = "SAMP", make.pairwise = F)


# excluding categorical here, they slow everything down!
l.data_list=list(gm, gm2,gm3,cat_immunoscores)
data_list=data.frame(do.call(rbind, l.data_list))
colnames(data_list)=substr(gsub("\\.", "-", colnames(data_list)), 1, 15)

#*******************************************************
sum_mutations=t(data.frame(colSums(data.matrix(fm[grep("code_potential_somatic", rownames(fm)),]), na.rm = T)))
rownames(sum_mutations)="N:SAMP:MUTATION_RATE"
sum_mutations=sum_mutations[,match(colnames(fm), colnames(sum_mutations)), drop=F]
colnames(sum_mutations)=colnames(fm)

# sum certain type mutation
mafs_org <- read.delim("/research/groups/biowhat_share/public_data/TCGA/firehose/mafs/gdac.broadinstitute.org_DLBC-TP.Mutation_Assessor.Level_4.2016012800.0.0/DLBC-TP.maf.annotated", header=T, stringsAsFactors=F)
mafs=mafs_org

# filter to coding:
table(mafs$type)
mafs=mafs[!mafs$is_silent==1,]
mafs=mafs[mafs$type,]

gt=unique(mafs$type)

# sapply(gt, function(g){
#   a=mafs[mafs$type%in%g,]
#   
#   a$Tumor_Sample_Barcode=substr(a$Tumor_Sample_Barcode, 1, 15)
#   mafs_table=colSums(table(a$Hugo_Symbol, a$Tumor_Sample_Barcode))
#   mafs_table=mafs_table[match(colnames(fm), names(mafs_table))]
#   
#   if(sum(is.na(mafs_table))>2)return(NULL)
#   
#   cor.test(as.numeric(mafs_table), as.numeric(gm2), use="complete.obs", method="spearman")
# })

mafs$Tumor_Sample_Barcode=substr(mafs$Tumor_Sample_Barcode, 1, 15)
mafs_table=colSums(table(mafs$Hugo_Symbol, mafs$Tumor_Sample_Barcode))
mafs_table=mafs_table[match(colnames(fm), names(mafs_table))]

FILE="/research/groups/biowhat_share/public_data/TCGA/firehose/cnvr_segments/gdac.broadinstitute.org_DLBC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/DLBC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"
segments=read.delim(FILE, header=T, stringsAsFactors=F, sep="\t")
segments=segments[!substr(segments$Sample,14,15)=="11",]

segments$Sample=substr(segments$Sample,1,15)
samples=unique(segments$Sample)

library(parallel)
fragments=unlist(mclapply(samples, function(sample){
  A=segments[segments$Sample%in%sample,]
  segments_fragmented=A[abs(A$Segment_Mean)>=0.3,]
  segments_not_fragmented=A[abs(A$Segment_Mean)<0.3,]
  segments_fragmented=sum(as.numeric(abs(segments_fragmented[,4]-segments_fragmented[,3])))
  segments_not_fragmented=sum(as.numeric(abs(segments_not_fragmented[,4]-segments_not_fragmented[,3])))
  segments_fragmented/(segments_not_fragmented+segments_fragmented)}, mc.cores=10))

fragments=data.frame(t(fragments[match(colnames(fm), samples)]))
colnames(fragments)=colnames(fm)
rownames(fragments)="N:SAMP:GENOME_FRAGMENTATION_RATE"

cor.test(as.numeric(fragments), as.numeric(gm2), use="complete.obs", method="spearman")
cor.test(as.numeric(mafs_table), as.numeric(gm), use="complete.obs", method="spearman")

cor.test(as.numeric(sum_mutations), as.numeric(gm), use="complete.obs", method="spearman")
cor.test(as.numeric(gm), as.numeric(gm2), use="complete.obs", method="spearman")

plot(as.numeric(gm), as.numeric(fragments))
plot(as.numeric(gm), as.numeric(mafs_table))

plot(as.numeric(gm2), as.numeric(fragments))
plot(as.numeric(gm), as.numeric(gm3))

# *****************************************************************************************
l.data_list=list(fm, data_list, gsva_es, ExampleEstimates, results_fm, sum_mutations, fragments, cat_feats)
fm_dat=data.frame(do.call(rbind, l.data_list), stringsAsFactors = F)

matrix=fm_dat
colnames(matrix)=substr(gsub("\\.", "-", colnames(matrix)), 1, 15)

A=apply(cluster_cancermap, 1, unique)

B=unlist(lapply(A, function(d)sum(d%in%c(1,0))>=2))

if(!all(B))stop("Check comparisons, impossible comparisons made")

save(matrix, file="TCGA_DLBCL_FM_DUFVA.Rdata")

write.table(t(c("N:SAMP", as.character(colnames(matrix)))), file="TCGA_DLBCL_FM_DUFVA.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(matrix, file="TCGA_DLBCL_FM_DUFVA.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)
