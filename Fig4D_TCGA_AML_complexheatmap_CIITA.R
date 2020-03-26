GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))

library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("TCGA_AML_FM_DUFVA.Rdata"))

gexp=t(scale(t(data.matrix(fm[grepl("CIITA|HLA-D|Score", rownames(fm))&grepl("GEXP|N:SAMP:.*.Score", rownames(fm)),]))))
gexp=gexp[!grepl("HLA-DO", rownames(gexp)),]

parsemeth=read.delim("data_meth_parse.txt", header = F, stringsAsFactors = F)
parsemeth=parsemeth[grep("CIITA|HLA-D|CD274|CD40|RAET1E|PDCD1LG2|TCGA-", parsemeth[,3]),]

meth_annot=cbind(parsemeth[-1,1], parsemeth[-1,3], parsemeth[-1,4], parsemeth[-1,5])
meth=data.matrix(parsemeth[2:dim(parsemeth)[1],grep("\\.", parsemeth[4,])])
rownames(meth)=paste(meth_annot[,2],meth_annot[,1], sep="_")
colnames(meth)=substr(parsemeth[1,grep("\\.", parsemeth[4,])], 1, 15)

meth=meth[!rowSums(is.na(meth))>100,]
gexp=gexp[,colnames(gexp)%in%colnames(meth)]
fm=fm[,colnames(fm)%in%colnames(meth)]

meth=meth[,match(colnames(gexp), colnames(meth))]
colnames(meth)=colnames(gexp)

meth=meth[!grepl("HLA-DO|LG", rownames(meth)),]

meth=meth[order(rownames(meth)),]

ciita_gexp=data.matrix(fm[grepl("CIITA", rownames(fm))&grepl("GEXP", rownames(fm)),!is.na(gexp[1,])])
HLAII_gexp=data.matrix(fm[grepl("HLAIIScore", rownames(fm))&grepl("N:SAMP", rownames(fm)),!is.na(gexp[1,])])

dat=data.frame("CIITA"=as.numeric(ciita_gexp),"HLAIIscore"=as.numeric(HLAII_gexp))
dat$grp=factor("group")
myColors <- c("#E41A1C")

names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(name = "grp",values = myColors)

main="CIITA - HLAII score"

p=ggplot(dat, aes(CIITA, HLAIIscore, colour = grp)) +
  geom_point(size=1) + colScale +
  theme_classic(base_size = 12) +
  scale_y_continuous(breaks = pretty(dat$HLAIIscore, n = 6))+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        text=element_text()
  )+
  theme(axis.text.x = element_text(colour="grey20",size=10,face="plain"),
        axis.text.y = element_text(colour="grey20",size=10,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=10,face="plain"),
        axis.title.y = element_text(colour="grey20",size=10,face="plain"))  

ggsave(p, filename = "CIITA_HLAII_correlation.pdf", width = 2 , height = 2.3)

# get CEBPA annotations:
mut=read.delim("CEBPA_mutations.txt", header = T, stringsAsFactors = F)
cebpa=table(substr(mut$Tumor_Sample_Barcode, 1, 15))
cebpa_sil=colnames(gexp)%in%names(cebpa)[cebpa==1]
cebpa_dm=colnames(gexp)%in%names(cebpa)[cebpa>1]

# make complex heatmap:
# make logical vectors with fab and genetic subtype
is.na.s=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(Complex_Cytogenetics|GENETICS):::::",]))

complex=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(Complex_Cytogenetics|GENETICS):::::",!(is.na(is.na.s))]))
normal=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(Normal_Karyotype|GENETICS):::::",!(is.na(is.na.s))]))
mll=as.logical(as.numeric(colSums(data.matrix(fm[rownames(fm)%in%c("B:SAMP:I(MLL_translocation,_poor_risk|GENETICS):::::", "B:SAMP:I(MLL_translocation,_t(9;11)|GENETICS):::::"),!(is.na(is.na.s))]))>0))
runx=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(RUNX1-RUNX1T1|GENETICS):::::",!(is.na(is.na.s))]))
cbfb=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(CBFB-MYH11|GENETICS):::::",!(is.na(is.na.s))]))
pml=as.logical(as.numeric(fm[rownames(fm)%in%"B:SAMP:I(PML-RARA|GENETICS):::::",!(is.na(is.na.s))]))

# mut
FLT3=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:FLT3:chr13:28577411:28674729:-:y_n_somatic",!(is.na(is.na.s))]))
cebpa_sil=cebpa_sil[!(is.na(is.na.s))]
cebpa_dm=cebpa_dm[!(is.na(is.na.s))]
RUNX=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:RUNX1:chr21:36160098:37357047:-:y_n_somatic",!(is.na(is.na.s))]))
TP53=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:TP53:chr17:7565097:7590863:-:y_n_somatic",!(is.na(is.na.s))]))
NPM1=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:NPM1:chr5:170814708:170837888:+:y_n_somatic",!(is.na(is.na.s))]))
IDH1=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:IDH1:chr2:209100953:209119806:-:y_n_somatic",!(is.na(is.na.s))]))
IDH2=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:IDH2:chr15:90627212:90645708:-:y_n_somatic",!(is.na(is.na.s))]))
TET=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:TET1:chr10:70320117:70454239:+:y_n_somatic",!(is.na(is.na.s))]))
DNMT3A=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:DNMT3A:chr2:25455845:25565459:-:y_n_somatic",!(is.na(is.na.s))]))
FLT3=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:FLT3:chr13:28577411:28674729:-:y_n_somatic" ,!(is.na(is.na.s))]))
del5q7q=as.logical(as.numeric(fm[rownames(fm)%in%"B:CLIN:I(-5_or_del(5q)|FISH_test_component):::::" ,!(is.na(is.na.s))]))|as.logical(as.numeric(fm[rownames(fm)%in%"B:CLIN:I(-7_or_del(7q)|FISH_test_component):::::" ,!(is.na(is.na.s))]))
NRAS=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:NRAS:chr1:115247085:115259515:-:y_n_somatic" ,!(is.na(is.na.s))]))
KRAS=as.logical(as.numeric(fm[rownames(fm)%in%"B:GNAB:KRAS:chr12:25358180:25403863:-:y_n_somatic" ,!(is.na(is.na.s))]))

d.l=list(pml,cebpa_sil,cebpa_dm,IDH1,IDH2, TET, DNMT3A, mll,cbfb, runx, del5q7q, NPM1,NRAS,KRAS, FLT3, RUNX, TP53, complex, normal)
names(d.l)=c("PML_RARA","CEBPA","CEBPA_DM", "IDH1", "IDH2", "TET3", "DNMT3A", "MLL", "CBFB_MYH11", "RUNX1_RUNX1T1", "-5/7q", "NPM1", "NRAS", "KRAS", "FLT3", "RUNX1","TP53", "Complex_Karyotype", "Normal_Karyotype")
clin=do.call(cbind, d.l)*1

fab=t(fm[rownames(fm)%in%"C:CLIN:leukemia_french_american_british_morphology_code:::::",!(is.na(is.na.s)),drop=F])

fab=data.frame(fab, stringsAsFactors = F)
colnames(fab)=c("FAB")

# sort everything based on CIITA exp:
clin=clin[order(gexp[1,!(is.na(is.na.s))], decreasing = T),]
fab=fab[order(gexp[1,!(is.na(is.na.s))], decreasing = T),,drop=F]
meth=meth[,order(gexp[1,!(is.na(is.na.s))], decreasing = T)]
gexp=gexp[,order(gexp[1,!(is.na(is.na.s))], decreasing = T)]

ha2 = HeatmapAnnotation(df = data.frame(clin))
ha = HeatmapAnnotation(df = data.frame(fab[match(colnames(gexp), rownames(fab)),]))

gexp[gexp>2]=2
gexp[gexp<(-2)]=-2

rownames(gexp)=do.call(rbind, strsplit(rownames(gexp), ":"))[,3]

pdf("Figure4D_TCGA_AML_ComplexHeatmap_CIITA_gexp.pdf", height = 15, width = 30)
Heatmap(data.matrix(gexp[grep("CIITA|HLA-D", rownames(gexp)),]), bottom_annotation = ha, name = "CIITA", cluster_columns = F, cluster_rows = F, show_column_names = F)
dev.off()

corres=cor(t(meth), gexp[1,], use = "complete.obs")

res=p.adjust(apply(meth, 1, cor.test.p, y = gexp[1,], method="pearson"), method="bonferroni")

filt=names(res)[res<0.0001]

meth2=meth[match(filt, rownames(meth)),]
meth_annot2=meth_annot[match(filt, paste(meth_annot[,2],meth_annot[,1], sep="_")),]

meth2=meth2[order(meth_annot2[,2]),]
meth_annot2=meth_annot2[order(meth_annot2[,2]),]

# take only highest anti-correlated from significant:
probes=sapply(unique(meth_annot2[,2]), function(g){
  a=rownames(meth2[meth_annot2[,2]%in%g,])
  names(which.min(corres[rownames(corres)%in%a,]))
}
)

meth3=t(scale(t(meth2[rownames(meth2)%in%probes,])))
rownames(meth3)=unique(meth_annot2[,2])

pdf("Figure4D_TCGA_AML_ComplexHeatmap_CIITA_meth.pdf", height = 15, width = 30)
Heatmap(meth3, name = "CIITA", bottom_annotation = ha, cluster_columns = F, cluster_rows = F, col = colorRamp2(c(-0.5, 0, 1), c("brown","grey", "orange")), show_column_names = F)
dev.off()

# order based on mutual exclusivity and HLA groups
# fab[fab[,1]%in%"Not_Classified",1]=NA
# m=cbind((gexp[17,]<(-1))*1, clin,(fab[,1]%in%c("M3"))*1, (fab[,1]%in%c("M0_Undifferentiated", "M1", "M2"))*1, (fab[,1]%in%c("M4", "M5"))*1)
# colnames(m)[1]="HLAIIlow"
# colnames(m)[21]="FAB_M3"
# colnames(m)[22]="FAB_M0_M1_M2"
# colnames(m)[23]="FAB_M4_M5"
# 
# m=m[do.call(order, -as.data.frame(m)),,drop=F]
# sample_ord=rownames(m)
# 
# # order samples
# m=t(m[match(sample_ord, rownames(m)),])
# gexp2=gexp[,match(sample_ord, colnames(gexp))]
# meth4=meth3[,match(sample_ord, colnames(meth3))]
# 
# pdf("TCGA_AML_ComplexHeatmap_gexp.pdf", height = 15, width = 30)
# Heatmap(data.matrix(gexp2), top_annotation = HeatmapAnnotation(df = data.frame(fab[match(sample_ord, rownames(fab)),])), name = "CIITA", cluster_columns = F, cluster_rows = F, show_column_names = F)
# dev.off()
# 
# pdf("TCGA_AML_ComplexHeatmap_meth.pdf", height = 15, width = 30)
# Heatmap(meth4, name = "CIITA", cluster_columns = F, cluster_rows = F, col = colorRamp2(c(-0.5, 0, 1), c("brown","grey", "orange")), show_column_names = F)
# dev.off()
#
# source("/research/users/ppolonen/git_home/common_scripts/visualisation/oncoprint_memo.R")
# pdf("TCGA_AML_ComplexHeatmap_genetics_HLA.pdf", height = 15, width = 15)
# oncoPrint(m, sort=F)
# dev.off()