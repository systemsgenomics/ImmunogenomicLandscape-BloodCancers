# Tools
GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))
library(RColorBrewer)
library(survival)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

files=list.files(pattern = "tableS7")
files=files[!grepl("_CHOP", files)] # CHOP treated DLBCL not shown, as it has very different survival to RCHOP

univariate.results=lapply(files[grepl("signif", files)], fread, data.table=F)

univariate.results.all=lapply(files[!grepl("signif", files)], fread, data.table=F)

names(univariate.results)=files[grepl("signif", files)]
names(univariate.results.all)=files[!grepl("signif", files)]

AML=do.call(rbind, univariate.results[grep("AML", names(univariate.results))])
AML$Name=gsub("Cancer_Myeloma", "Hemap_MM", AML$Name)
AML$Name=gsub("_RCHOP", "", AML$Name)

DLBCL=do.call(rbind, univariate.results[grep("DLBCL", names(univariate.results))])
DLBCL$Name=gsub("Cancer_Myeloma", "Hemap_MM", DLBCL$Name)
DLBCL$Name=gsub("_RCHOP", "", DLBCL$Name)

MM=do.call(rbind, univariate.results[grep("MM", names(univariate.results))])
MM$Name=gsub("Cancer_Myeloma", "Hemap_MM", MM$Name)
MM$Name=gsub("_RCHOP", "", MM$Name)

# all
pancan=do.call(rbind, univariate.results.all)
pancan$Name=gsub("Cancer_Myeloma", "Hemap_MM", pancan$Name)
pancan$Name=gsub("_RCHOP", "", pancan$Name)

fun_forestplot=function(data, NAME="data", BOX=0.1,cex=2, colorv="black"){
  library(forestplot)
  
  txt.df=data.frame("Feature"=data$Feature, "Cohort"=gsub("_", " ", data$Name), "FDR"=signif(data$Adj.P, 2), stringsAsFactors = F)
  txt.df$Feature[duplicated(txt.df$Feature)]=""
  txt.df=rbind(c("Feature", "Cohort", "FDR"), txt.df)
  
  coef.df=data.frame("HR"=data$`exp(coef)`, "lower .95"=data$`lower .95`, "upper .95"=data$`upper .95`)
  coef.df=rbind(c(NA, NA, NA), coef.df)
  
  xticks=seq(ifelse(min(data$`lower .95`)<0.5, 0, 0.5), min(c(5, max(data$`upper .95`))), by = 0.1)
  
  attr(xticks, "labels") = xticks%in%seq(-4, 4, by=0.5)
  
  forestplot(txt.df,coef.df, new_page = T, zero = c(0.98, 1.02), 
             clip =c(-1, 2), is.summary = c(T, rep(F, length(data$Name))),
             xticks=xticks,
             boxsize=BOX,
             xlab="HR",
             col=fpColors(box=colorv),
             txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Helvetica", cex=cex*1.25),
                                           gpar(fontfamily = "Helvetica", col = "black", cex=cex)),
                              summary = gpar(fontfamily = "Helvetica", cex=cex*1.5),
                              ticks = gpar(fontfamily = "Helvetica", cex=cex),
                              xlab  = gpar(fontfamily = "Helvetica", cex = cex*1.5)))
  
}
cat=c("Subtype","ImmunoScores","Inhibitory ligand", "Stimulatory ligand", "Stromal/cancer gene (Rho > 0)", "Stromal/cancer gene (Rho < 0)","CTL/NK gene","Clinical", "CGA")
cols=data.frame("name"=c("Subtype","ImmunoScores","Inhibitory ligand", "Stimulatory ligand", "Stromal/cancer gene (Rho > 0)", "Stromal/cancer gene (Rho < 0)","CTL/NK gene","Clinical", "CGA", "MDS-signature gene"),
                "color"=c("#acb839","#5e2883","#1f78b4","#b2df8a","#377eb8","grey50","#e41a1c","brown", "#d7a85b", "indianred"), stringsAsFactors = F)

count.signif=table(AML$Feature)
signif.feats=names(count.signif)[count.signif>1]
AML.filt=AML[AML$Feature%in%signif.feats,]
AML.filt=AML.filt[order(match(AML.filt$Type, cat), AML.filt$Feature, AML.filt$Adj.P),]

count.signif=table(MM$Feature)
signif.feats=names(count.signif)[count.signif>1]
MM.filt=MM[MM$Feature%in%signif.feats,]

count.signif=table(DLBCL$Feature)
signif.feats=names(count.signif)[count.signif>1]
DLBCL.filt=DLBCL[DLBCL$Feature%in%signif.feats,]

Fig7_AML=AML[AML$Feature%in%c("MDS-like", "CLEC2B", "Monocyte-like","C10orf54"),]
Fig7_MM=MM[MM$Feature%in%c("Freq-CGA", "MAGEA1", "MAGEB2"),]
Fig7_DLBCL=DLBCL[DLBCL$Feature%in%c("ABC","C1QA","C1QB","C1QC","CCL18", "CD163", "LYZ"),]
Fig7_pancan=pancan[pancan$Feature%in%c("CD274", "C10orf54", "CD58", "CD84"),]
Fig7_pancan=Fig7_pancan[Fig7_pancan$Adj.P<0.2,]

all=rbind(Fig7_pancan)

all=all[order(match(all$Feature, c("CD274", "C10orf54", "CD58", "CD84", "HLAII")), -all$`exp(coef)`),]

col=cols[match(unique(all$Type), cols$name),2]
names(col)=unique(all$Type)

pdf("FigS7B_univariate_pancan_forestPlot.pdf", height = 4, width = 5)
fun_forestplot(all, "AML, MM, DLBCL", BOX=0.5, cex=0.5, colorv ="black")
dev.off()

# show even if not significant:
Fig7_AML=pancan[pancan$Feature%in%c("MDS-like", "CLEC2B", "Monocyte-like","C10orf54")&grepl("AML", pancan$Name),]
Fig7_MM=pancan[pancan$Feature%in%c("Freq-CGA", "MAGEA1", "MAGEB2")&grepl("MM", pancan$Name),]
Fig7_DLBCL=pancan[pancan$Feature%in%c("C1QA","C1QB","C1QC", "CCL18","CD163", "LYZ")&grepl("DLBCL", pancan$Name),]

all=rbind(Fig7_AML, Fig7_MM, Fig7_DLBCL)

all=all[order(match(all$Feature, c(c("MDS-like", "CLEC2B", "Monocyte-like","C10orf54"), c("Freq-CGA", "MAGEA1", "MAGEB2"), c("C1QA","C1QB","C1QC", "CCL18", "CD163", "LYZ"))), -all$`exp(coef)`),]

col=cols[match(unique(all$Type), cols$name),2]
names(col)=unique(all$Type)

pdf("Fig7A_univariate_forestPlot.pdf", height = 5, width = 5)
fun_forestplot(all, "AML, MM, DLBCL", BOX=0.5, cex=0.5, colorv ="black")
dev.off()

# draw KM curve from a few examples:
load("Hemap_AML_survival_data.Rdata")

pdf("Fig7B_KM_VISTA.pdf", width = 2, height = 2.5)
fun.kapplanMeier(TIME = TIME[logicalv[[1]]], STATUS = STATUS[logicalv[[1]]], CONTINUOUS  = as.numeric(gexp[logicalv[[1]],"C10orf54"]), conf.int = F, MONTHS=F, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "VISTA")
dev.off()

# draw KM curve from a few examples:
load("CoMMpass_survival_data.Rdata")

pdf("FigS7A_CGA_CoMMpass.pdf", width = 2, height = 2.5)
fun.kapplanMeier(TIME = TIME[logicalv[[1]]], STATUS = STATUS[logicalv[[1]]], GROUPS  = as.numeric(immunoscore[logicalv[[1]],"Freq.CGA"]>=7), conf.int = F, MONTHS=T, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "MAGEA1")
dev.off()

# draw KM curve from a few examples:
load("Hemap_MM_survival_data.Rdata")

pdf("FigS7A_CGA_Hemap.pdf", width = 2, height = 2.5)
fun.kapplanMeier(TIME = TIME[logicalv[[1]]], STATUS = STATUS[logicalv[[1]]], GROUPS  = as.numeric(immunoscore[logicalv[[1]],"Freq.CGA"]>=7), conf.int = F, MONTHS=T, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "MAGEA1")
dev.off()

# draw KM curve from a few examples:
load("Hemap_DLBCL_survival_data.Rdata")

pdf("Fig7D_KM_CD163_LYZ.pdf", width = 2, height = 2.5)
fun.kapplanMeier(TIME = TIME[logicalv[[5]]], STATUS = STATUS[logicalv[[5]]], CONTINUOUS  = as.numeric(gexp[logicalv[[5]],"CD163"]), conf.int = F, MONTHS=F, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "CD163")
fun.kapplanMeier(TIME = TIME[logicalv[[5]]], STATUS = STATUS[logicalv[[5]]], CONTINUOUS  = as.numeric(gexp[logicalv[[5]],"LYZ"]), conf.int = F, MONTHS=F, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "LYZ")
dev.off()

# load("GSE98588_DLBCL_survival_data.Rdata")
# 
# pdf("Fig7_KM_CGA_DLBCL.pdf", width = 2, height = 2.5)
# fun.kapplanMeier(TIME = TIME[logicalv[[1]]], STATUS = STATUS[logicalv[[1]]], GROUPS  = as.numeric(immunoscore[logicalv[[1]],"Freq.CGA"])>=4*1, conf.int = F, MONTHS=F, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = "75th_25th_percentile", INDIVIDUAL_GROUPS=F, NAME = "#CGA", COLORV =c("grey75", "red"))
# dev.off()