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

meth=get(load("meth_hlalow_pml.Rdata"))
myDiffSigboth=get(load("myDiffSigboth_hlalow_pml.Rdata"))

#******************************* visualize location **********************************

# wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
# gunzip *.gz
refGeneInfo = getRefgeneInfo(system.file("annotation", "refGene.txt",
                                         package = "methylSig"))


refGeneAnn = refGeneAnnotation(refGeneInfo, myDiffSigboth)

cpgInfo = getCpGInfo("cpgIslandExt.txt")

# DE methylated cites:
myDiffDMCs = myDiffSigboth[myDiffSigboth[,"qvalue"] < 0.05
                           & abs(myDiffSigboth[,"meth.diff"]) > 25,]

# Annotate
cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)

# annotate refgene
refGeneInfo = getRefgeneInfo("refGene.txt")

refGeneAnn = refGeneAnnotation(refGeneInfo, myDiffSigboth)

# what kin are DMCs
refGeneAnnotationPlot(refGeneAnn,main="ALL",
                      priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))

refGeneAnnotationPlot(refGeneAnnDmc, main="DMC",
                      priority=c("promoter","cds", "noncoding", "5'utr", "3'utr"))

# TF cites:
tfbsInfo = getTFBSInfo("ENCODE_AwgTfbs.hg19.txt")

DMCIndex = (myDiffSigboth[,"qvalue"] < 0.05
            & abs(myDiffSigboth[,"meth.diff"]) > 25)

pvalue = methylSig.tfbsEnrichTest(myDiffSigboth, DMCIndex, tfbsInfo, plot = F)

highest=sort(pvalue)
plot_h=data.frame("n"=names(highest)[highest<1E-5], "v"=-1*log10(highest)[highest<1E-5]) # wrote this out

write.table(rownames(plot_h), "ENCODE_AwgTfbs.hg19_inDMC.txt", row.names = F, col.names = F, quote = F, sep="\t")

system("grep -f Binomial_TFs_enriched_DMC.txt ENCODE_AwgTfbs.hg19.txt >ENCODE_AwgTfbs.hg19_inDMC.txt")

tfbsInfo2 = getTFBSInfo("ENCODE_AwgTfbs.hg19_inDMC.txt")

# CIITA and enhancers: chr16:10960573-11019050
# CIITA only chr16:10970055-11018840

pdf("FigureSE_CIITA_methylation_RRBS.pdf", width = 8, height = 12)

p <-ggplot(plot_h, aes(n, v))
p=p +geom_bar(stat = "identity", aes(fill = n), position = "dodge") +
  xlab("TFbinding_site") + ylab("-log10 P-value") +
  ggtitle("Binomial test, DMC vs TFbinding ENCODE") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))
p

methylSigPlot(meth, "chr16", c(10960573, 11019050), groups=c(1,0),
              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, myDiff=myDiffSigboth,
              tfbsInfo=tfbsInfo, tfbsDense=T, sigQ=0.01, noGroupEst = F, noDataLine=T, cex=0.1)


methylSigPlot(meth, "chr16", c(10960055, 10990055), groups=c(1,0),
              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, myDiff=myDiffSigboth,
              tfbsInfo=tfbsInfo, tfbsDense=T, sigQ=0.01, noGroupEst = F, noDataLine=T, cex=0.1)

# only significant binding cites
methylSigPlot(meth, "chr16", c(10960055, 10990055), groups=c(1,0),
              cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, myDiff=myDiffSigboth,
              tfbsInfo=tfbsInfo2, tfbsDense=F, sigQ=0.01, noGroupEst = F, noDataLine=T, cex=0.1)

dev.off()

coord=data.frame(refGeneInfo[[2]], stringsAsFactors = F)

coord=coord[grepl("^HLA", coord$gene)&!grepl("_", coord$chr)&!grepl("HLA-F", coord$gene),]
coord$chr=as.character(coord$chr)

for(i in 1:dim(coord)[1]){
  pdf(paste0(getwd(), "/", coord$gene[i], "_methylation_RRBS.pdf"), width = 8, height = 12)
  
  
  methylSigPlot(meth, chr = coord[i,1], loc.range = c(coord[i,3]-2000, coord[i,4]+2000), groups=c(1,0),
                cpgInfo=cpgInfo, refGeneInfo=refGeneInfo, myDiff=myDiffSigboth,
                tfbsInfo=tfbsInfo, tfbsDense=T, sigQ=0.01, noGroupEst = F, noDataLine=T, cex=0.1)
  dev.off()
  
}


#******************************* visualize as heatmap **********************************

# summarize per group:

# take DE meth regions and extract from data
x = meth[,"numCs"]/meth[, "coverage"]
colnames(x) = meth@sample.ids
rownames(x) = meth@data.ids

listInMeth = match(myDiffDMCs@data.ids, meth@data.ids)
y = x[listInMeth,]
colnames(y) = meth@sample.ids

hlalow_samples=colnames(y)%in%hlalow
HLAlow_mean=rowMeans(y[,hlalow_samples], na.rm=T)
HLArest_mean=rowMeans(y[,!hlalow_samples], na.rm=T)

y2=data.frame("chr"=as.character(myDiffDMCs@data.chr), "start"=myDiffDMCs@data.start, "end"=myDiffDMCs@data.end,"ID"=myDiffDMCs@data.ids, stringsAsFactors=F)

refgene=data.frame(refGeneInfo[[2]], stringsAsFactors = F)
i <- sapply(refgene, is.factor)
refgene[i] <- lapply(refgene[i], as.character)

refgene$start[refgene$strand=="+"]=refgene$start[refgene$strand=="+"]-1000
refgene$end[refgene$strand=="-"]=refgene$end[refgene$strand=="-"]+1000
refgene$start[refgene$start<0]=0

# intersect this with methylation regions, return name of the region and gene, transcript name:
bed=cbind(refgene[,1], refgene[,3], refgene[,4], refgene[,8], refgene[,9])

temp1=tempfile()
temp2=tempfile()
temp3=tempfile()

write.table(bed, temp1, row.names = F, col.names = F, quote = F, sep="\t")
write.table(y2, temp2, row.names = F, col.names = F, quote = F, sep="\t")

command=paste0("bedtools intersect -a ",temp1, " -b ",temp2, " -wa -wb |awk -v OFS=\"\t\" -F\"\t\" '{print $4, $5, $9}' |sort -u > ", temp3)
cat(command,"\n")
try(system(command))

matched_data=read.delim(temp3, header = F, stringsAsFactors = F)
unlink(c(temp1, temp2, temp3))

y4=y[y2$ID%in%matched_data[grepl("CIITA|HLA-D", matched_data[,2]),3],]

dat=y[,order(colnames(y)%in%hlalow, decreasing = T)]
dat_sub=dat[rownames(dat)%in%rownames(y4),]

fab=annot$X.Sample_characteristics_ch1.4
cyto=annot$X.Sample_characteristics_ch1.5
npm=annot$X.Sample_characteristics_ch1.18=="npm1: Pos"
pml=annot$X.Sample_characteristics_ch1.9=="t(15;17): Pos"
cbfb=annot$X.Sample_characteristics_ch1.11=="inv(16): Pos"
runx1_t1=annot$X.Sample_characteristics_ch1.10=="t(8;21): Pos"
mll=annot$X.Sample_characteristics_ch1.6=="t(v;11q23): Pos"
del5q=grepl("cytogenetics: -5", annot$X.Sample_characteristics_ch1.5)
idh1=annot$X.Sample_characteristics_ch1.23=="idh1: Pos"
idh2=annot$X.Sample_characteristics_ch1.24=="idh2: Pos"
dnmt3a=annot$X.Sample_characteristics_ch1.25=="dnmt3a: Pos"
dnmt3a_type=annot$X.Sample_characteristics_ch1.26
evi1=annot$X.Sample_characteristics_ch1.27=="evi1 overexpression: Pos"
kit=annot$X.Sample_characteristics_ch1.22=="kit: Pos"
cebpa=annot$X.Sample_characteristics_ch1.21 =="cebpa: Sil"
cebpaDM=annot$X.Sample_characteristics_ch1.21 =="cebpa: DM" 
kras=annot$X.Sample_characteristics_ch1.20=="kras: Pos"
nras=annot$X.Sample_characteristics_ch1.19=="nras: Pos"
flt1=annot$X.Sample_characteristics_ch1.17=="flt3 tkd: Pos"
flt2=annot$X.Sample_characteristics_ch1.16=="flt3 itd: Pos"
normal=annot$X.Sample_characteristics_ch1.15=="normal karyotype: Pos"
complex=grepl("complex|Complex", annot$X.Sample_characteristics_ch1.5)

load("GSE6891_immunoscores.Rdata")
immunoscores=immunoscores[,match(annot$X.Sample_title, colnames(immunoscores))]
colnames(immunoscores)=annot$X.Sample_title
coords=read.delim("cancermap_GSE86952_AML_15pct_genes_BH-SNE_mean-shift_BW1.5.txt", header=T, stringsAsFactors = F)
coords=coords[match(annot$X.Sample_title, as.character(coords$ID)),]

hlalow=scan("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE86952_AML_RRBS/HLAII_low_samples_Zscore_minus1.txt", "a")
hlalow.s=(pml|annot$X.Sample_title%in%hlalow)*1
ha=data.frame("HLAlow"=hlalow.s,pml,cebpa,cebpaDM,idh1,idh2,dnmt3a,mll, cbfb,runx1_t1,del5q, fab,cyto,npm,evi1,kit,kras,nras,flt1,flt2,normal,complex, t(immunoscores), "cluster"=as.character(coords$cluster), "CIITA Meth Mean"=colMeans(dat_sub, na.rm = T), stringsAsFactors = F)

# remove MDS samples:
dat=dat[,!is.na(coords$cluster)]
dat_sub=dat_sub[,!is.na(coords$cluster)]

ha=ha[!is.na(ha$cluster),]

# order based on mutual exclusivity and HLA groups
m=cbind(ha$HLAlow, ha$pml, ha$cebpa, ha$cebpaDM, ha$idh1, ha$idh2, ha$dnmt3a, ha$mll, ha$cbfb, ha$runx1_t1, ha$del5q,  ha$npm, ha$kras, ha$nras, ha$flt1, ha$flt2)
rownames(m)=rownames(ha)
colnames(m)=c("HLAIIlow", "PML-RARA", "CEBPA", "CEBPA_DM", "IDH1", "IDH2", "DNMT3A","MLL", "CBFB-MYH11", "RUNX1_RUNX1T1", "-5/7(q)", "NPM1", "KRAS", "NRAS", "FLT3tkd", "FLT3itd")

m=m[do.call(order, -as.data.frame(m)),,drop=F]
sample_ord=rownames(m)

# order samples
ha=ha[match(sample_ord, rownames(ha)),]
dat=dat[,match(sample_ord, colnames(dat))]
dat_sub=dat_sub[,match(sample_ord, colnames(dat_sub))]

# add row annotations:
cpgInfo = getCpGInfo("cpgIslandExt.txt")
cpgAnn = cpgAnnotation(cpgInfo,myDiffSigboth)

rowannot=cpgAnn[myDiffSigboth@data.ids%in%rownames(dat_sub)]
rowannot=as.character(rowannot)
rownames(dat_sub)=apply(cbind(rowannot,matched_data[match(rownames(dat_sub), matched_data[,3]),2]), 1, paste, collapse="_")

# simplify annotations for heatmap:
ha2=data.frame(ha$HLAlow,ha$fab,ha$DUFVA_HLAI_SCORE, ha$DUFVA_HLAII_SCORE, ha$DUFVA_CYTOLYTIC_SCORE,ha$CIITA, ha$normal,ha$complex, ha$cluster, ha$CIITA.Meth.Mean)

# Final plot:
pdf("FigureS4F_ComplexHeatmap_DMC_methylation.pdf", height = 15, width = 15)
Heatmap(dat_sub, top_annotation = HeatmapAnnotation(df = ha2), name = "DMC", cluster_columns = F, cluster_rows = F, col = colorRamp2(c(0, 0.5, 1), c("brown","lightgray", "orange")), use_raster = F)+
  Heatmap(t(scale(data.frame(ha$DUFVA_HLAI_SCORE, ha$DUFVA_HLAII_SCORE, ha$DUFVA_CYTOLYTIC_SCORE, ha$CIITA))),name = "DMC", cluster_columns = F, cluster_rows = F, use_raster = F)
dev.off()

source("/research/users/ppolonen/git_home/common_scripts/visualisation/oncoprint_memo.R")
pdf("FigureS4F_ComplexHeatmap_DMC_methylation_genetics.pdf", height = 15, width = 15)
oncoPrint(t(m), sort=F)
dev.off()