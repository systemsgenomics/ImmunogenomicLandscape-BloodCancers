setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/")

source("~/git_home/common_scripts/featurematrix/functions_generate_fm.R")

# fm data and clusters
load("TCGA_AML_FM.Rdata")

# annotations
annot = read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/TCGA_AML/TCGA_AML_cancermap_cluster_15pct_coordinates.txt", stringsAsFactors=F, header=T)
colnames(fm)=substr(gsub("\\.", "-", colnames(fm)), 1, 15)

# GSVA data
# dat=data.matrix(fm[grepl("N:GEXP:", rownames(fm)),])
# dat=dat[,!apply(dat, 2, function(v)any(is.na(v)))]
# rownames(dat)=gsub("N:GEXP:|:chr*.*", "", rownames(dat))
# colnames(dat)=gsub("\\.", "-", colnames(dat))
# save(dat, file="TCGA_AML_RNASEQ.Rdata")

load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/TCGA_AML/GSVA/TCGA_AML_all_samples_Combined_pathway_drug_signatures_2017_GSVA.Rdata")
rownames(gsva_es)=gsub(" ", "_", rownames(gsva_es))
gsva_es=gsva_es[,match(colnames(fm), colnames(gsva_es))]
colnames(gsva_es)=colnames(fm)

# TCGA survival update
surv <- read.delim("TCGA_clinical_annotations_updated_survival.txt", header=T, stringsAsFactors=F)
A=gsub("\\.", "-", substr(colnames(fm), 9, 12))
surv=surv[match(A, surv$TCGA.Patient.ID),]

D=t(surv)
rownames(D)=paste(ifelse(sapply(surv, class)=="character", "C", "N"), ":CLIN:", colnames(surv), ":::::", sep="")
colnames(D)=colnames(fm)

# vital status
living <- read.delim("TCGA_clinical_annotations_status.txt", header=T, stringsAsFactors=F)
A=gsub("\\.", "-", substr(colnames(fm), 1, 12))

living=living[match(A, living[,1]),]

E=t(living[,2])
rownames(E)="C:CLIN:vital_status_TCGA_paper:::::"
colnames(E)=colnames(fm)

# cell fractions
ciber <- t(read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_AML/CIBERSORT-Results.txt", row.names = 1, header=T, stringsAsFactors=F))
rownames(ciber)=paste0("N:SAMP:DUFVA_CIBERSORT_", gsub(" |-|\\.", "_", rownames(ciber)), ":::::")
colnames(ciber)=gsub("\\.", "-", substr(colnames(ciber), 1, 15))
ciber=ciber[,match(colnames(fm), colnames(ciber))]
colnames(ciber)=colnames(fm)

load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/TCGA_AML/MCP_counter_results.Rdata")
rownames(ExampleEstimates)=paste0("N:SAMP:DUFVA_", gsub(" |-|\\.", "_", rownames(ExampleEstimates)), ":::::")
ExampleEstimates=ExampleEstimates[,match(colnames(fm), substr(gsub("\\.", "-", colnames(ExampleEstimates)), 1, 15))]
colnames(ExampleEstimates)=colnames(fm)

# Clusters
clusters_cancermap=annot$X1.5..cluster

# annotated and cancermap clusters
cluster_cancermap=FUN_MAKE_ALL(clusters_cancermap, "cancermap_cluster", clusters_cancermap, 0)
colnames(cluster_cancermap)=substr(gsub("\\.", "-", annot$ID), 1, 15)

# cancermap
cluster_cancermap=cluster_cancermap[,match(colnames(fm), colnames(cluster_cancermap))]
colnames(cluster_cancermap)=colnames(fm)

# categorical
class_cancermap=FUN_MAKE_CATEGORICAL(clusters_cancermap, "cancermap_cluster")
colnames(class_cancermap)=substr(gsub("\\.", "-", annot$ID), 1, 15)
class_cancermap=class_cancermap[,match(colnames(fm), colnames(class_cancermap)), drop=F]
colnames(class_cancermap)=colnames(fm)

classes=as.character(class_cancermap)[1:173]
size=c(length(unique(classes)))
col1 <- colorRampPalette(brewer.pal(min(size, 12), "Set1"))(size)
feat_col1=unlist(lapply(classes, function(c){ind=unique(classes)%in%c
col1[ind]}))

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NA)
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

# make indicator features:
cat_immunoscores=make.features(df = classification, datatype = "SAMP", make.pairwise = F)

# excluding categorical here, they slow everything down!
l.data_list=list(gm, gm2, gm3, cat_immunoscores)
data_list=data.matrix(do.call(rbind, l.data_list))
colnames(data_list)=substr(gsub("\\.", "-", colnames(data_list)), 1, 15)

# 
data_list=data_list[,match(colnames(fm), colnames(data_list))]
colnames(data_list)=colnames(fm)

# mutation rate
sum_mutations=t(data.frame(colSums(data.matrix(fm[grep("y_n", rownames(fm)),]), na.rm = T)))
find=substr(colnames(sum_mutations),14,15)=="11"|substr(colnames(sum_mutations),14,15)=="01"
sum_mutations[find]=NA
rownames(sum_mutations)="N:SAMP:MUTATION_RATE"
colnames(sum_mutations)=colnames(fm)

mafs <- read.delim("/research/groups/biowhat_share/public_data/TCGA/firehose/mafs/gdac.broadinstitute.org_LAML-TB.Mutation_Assessor.Level_4.2015082100.0.0/LAML-TB.maf.annotated", header=T, stringsAsFactors=F)
mafs$Tumor_Sample_Barcode=substr(mafs$Tumor_Sample_Barcode, 1, 15)

mafs_table=colSums(table(mafs$Hugo_Symbol, mafs$Tumor_Sample_Barcode))
mafs_table=mafs_table[match(colnames(fm), names(mafs_table))]

FILE="/research/groups/biowhat_share/public_data/TCGA/firehose/cnvr_segments/gdac.broadinstitute.org_LAML.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/LAML.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"
segments=read.delim(FILE, header=T, stringsAsFactors=F, sep="\t")
segments=segments[!substr(segments$Sample,14,15)=="11",]

segments$Sample=substr(segments$Sample,1,15)
samples=unique(segments$Sample)

library(parallel)
fragments=unlist(mclapply(samples, function(sample){A=segments[segments$Sample%in%sample,]
segments_fragmented=A[abs(A$Segment_Mean)>0.3,]
segments_not_fragmented=A[abs(A$Segment_Mean)<0.3,]
segments_fragmented=sum(as.numeric(abs(segments_fragmented[,4]-segments_fragmented[,3])))
segments_not_fragmented=sum(as.numeric(abs(segments_not_fragmented[,4]-segments_not_fragmented[,3])))
segments_fragmented/(segments_not_fragmented+segments_fragmented)}, mc.cores=10))
options(scipen=999)

fragments=data.frame(t(fragments[match(colnames(fm), samples)]))
colnames(fragments)=colnames(fm)
rownames(fragments)="N:SAMP:GENOME_FRAGMENTATION_RATE"

cor(as.numeric(fragments), as.numeric(data_list[1,]), use="complete.obs",method="spearman")

cor.test(as.numeric(fragments), data_list[1,], use="complete.obs", method="spearman")
cor.test(as.numeric(sum_mutations)[1:173], data_list[1,1:173], use="complete.obs", method="spearman")
cor.test(as.numeric(gm), as.numeric(gm2), use="complete.obs", method="spearman")
cor.test(data_list[5,], data_list[1,], use="complete.obs", method="spearman")
boxplot(lapply(seq(7), function(n)c(na.omit(as.numeric(fragments[,clusters_cancermap==n])))))
boxplot(lapply(seq(7), function(n)c(na.omit(as.numeric(data_list[1,clusters_cancermap==n])))))

# *****************************************************************************************
l.data_list=list(fm, data_list, gsva_es, ExampleEstimates, class_cancermap, cluster_cancermap, sum_mutations, fragments, D, E)
fm_dat=data.frame(do.call(rbind, l.data_list), stringsAsFactors = F)

matrix=fm_dat
colnames(matrix)=substr(gsub("\\.", "-", colnames(matrix)), 1, 15)

save(matrix, file="TCGA_AML_FM_DUFVA.Rdata")

write.table(t(c("N:SAMP", as.character(colnames(matrix)))), file="TCGA_AML_FM_DUFVA.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(matrix, file="TCGA_AML_FM_DUFVA.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)
