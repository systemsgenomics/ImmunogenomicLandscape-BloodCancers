library(methylSig)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Instructions to use the package (v 0.1):
# http://sartorlab.ccmb.med.umich.edu/sites/default/files/methylSig.pdf
# https://qcb.ucla.edu/wp-content/uploads/sites/14/2017/02/Workshop-6-WGBS-D1.pdf # check also 2 and 3
# note, v 0.1 of the package is not available in https://github.com/sartorlab/methylSig. Current version of the package is very different.

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# annotations:
annot=read.delim("GSE86952_gsm2patient.txt", header = T, stringsAsFactors = F)
annot$X.Sample_title=gsub("AML blast | AML blast ", "", annot$X.Sample_title)

# use same samples as cancermap
coords=read.delim("cancermap_GSE86952_AML_15pct_genes_BH-SNE_mean-shift_BW1.5.txt", header=T, stringsAsFactors = F)
coords=coords[match(annot$X.Sample_title, as.character(coords$ID)),]

annot=annot[!is.na(coords[,1]),]

# ********************* differential methylation analysis: **************************
fileList=list.files("/research/groups/sysgen/raw_data/petri/GSE86952/", pattern = ".txt", full.names = T)
match_fileL=gsub("*.*Sample_|.mincov10.txt|_repeat", "", fileList)
fileList=fileList[match_fileL%in%annot$X.Sample_title]
match_fileL=match_fileL[match_fileL%in%annot$X.Sample_title]

library(devtools)
install_github('sartorlab/methylSig')

hlalow=scan("GSE86952_AML_RRBS/HLAII_low_samples_Zscore_minus1.txt", "a")

# # gexp:
sample.id = match_fileL
treatment = (annot$X.Sample_characteristics_ch1.5=="cytogenetics: t(15;17)"|annot$X.Sample_title%in%hlalow)*1

# Done already
meth <- methylSigReadData(fileList, sample.ids = sample.id, assembly = "hg19",treatment = treatment, context = "CpG", destranded=TRUE)
save(meth, file="meth_hlalow_pml.Rdata") # meth with pml-rara

# start here
meth=get(load("meth_hlalow_pml.Rdata"))

myDiffSigboth = methylSigCalc(meth, groups=c(1,0), min.per.group=5, num.cores=12)
save(myDiffSigboth, file="myDiffSigboth_hlalow_pml.Rdata") # meth with pml-rara