# survival analysis for all scores:
# many coefficients in survival analysis, using regularization and elastic net to select features for cox model.

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

validate.model=function(dataset,coef, ind=1, months=F, summary="3quantiles", filt=NULL){
  load(dataset)
  
  if(!is.null(filt)){
    logicalv=list(filt&logicalv[[ind]])
    ind=1
  }
  
  time=TIME[logicalv[[ind]]]
  time[time==0]=0.1
  status=STATUS[logicalv[[ind]]]
  
  genelist=rownames(coef)
  
  test.data=cbind(gexp[logicalv[[ind]], colnames(gexp)%in%genelist,drop=F], immunoscore[logicalv[[ind]], colnames(immunoscore)%in%genelist,drop=F], samp[logicalv[[ind]], colnames(samp)%in%genelist,drop=F])
  
  if(!(all(genelist%in%colnames(test.data)))){
    warning(paste("Features not found from test data:", paste(genelist[!genelist%in%colnames(test.data)], collapse = ",")))
    coef=coef[rownames(coef)%in%colnames(test.data),,drop=F]
    genelist=rownames(coef)
  }
  
  m=test.data[,match(genelist, colnames(test.data))]
  riskPI=as.numeric(coef) %*% data.matrix(t(m))
  
  print(summary(coxph(Surv(time, status) ~ PI.test, data.frame("PI.test"=as.numeric(riskPI)))))
  
  # plot validation data:
  print(fun.kapplanMeier(time, status, CONTINUOUS = as.numeric(riskPI), conf.int = F, MONTHS=months, PVAL=1,LWD = 0.5, CONTINUOUS_SUMMARY = summary, INDIVIDUAL_GROUPS=F, NAME = ""))
  
  return(riskPI)
}

validate.model.cox=function(dataset,coef, ind=1, months=F, summary="3quantiles", filt=NUL, NAME=""){
  load(dataset)
  
  if(!is.null(filt)){
    logicalv=list(filt&logicalv[[ind]])
    ind=1
  }
  
  time=TIME[logicalv[[ind]]]
  time[time==0]=0.1
  status=STATUS[logicalv[[ind]]]
  
  genelist=rownames(coef)
  
  test.data=cbind(gexp[logicalv[[ind]], colnames(gexp)%in%genelist,drop=F], immunoscore[logicalv[[ind]], colnames(immunoscore)%in%genelist,drop=F], samp[logicalv[[ind]], colnames(samp)%in%genelist,drop=F])
  
  if(!(all(genelist%in%colnames(test.data)))){
    warning(paste("Features not found from test data:", paste(genelist[!genelist%in%colnames(test.data)], collapse = ",")))
    coef=coef[rownames(coef)%in%colnames(test.data),,drop=F]
    genelist=rownames(coef)
  }
  
  m=test.data[,match(genelist, colnames(test.data))]
  riskPI=as.numeric(coef) %*% data.matrix(t(m))
  
  fit=coxph(Surv(time, status) ~ PI.test, data.frame("PI.test"=as.numeric(riskPI)))

  a=summary(fit)
  
  pval=a$coefficients[,5]
  coef=a$conf.int[,c(1,3,4)]
  univ=data.frame(rownames(a$coefficients),t(coef), pval, a$concordance[1],NAME, stringsAsFactors = F)
  rownames(univ)=NULL
  colnames(univ)=c("Feature", "exp(coef)", "lower .95", "upper .95", "P", "concordance", "Name")
  return(univ)
}


cols=data.frame("name"=c("Subtype","ImmunoScores","Inhibitory ligand", "Stimulatory ligand", "Stromal/cancer gene (Rho > 0)", "Stromal/cancer gene (Rho < 0)","CTL/NK gene","Clinical", "CGA", "MDS-signature gene"),
                "color"=c("#acb839","#5e2883","#1f78b4","#b2df8a","#377eb8","grey50","#e41a1c","brown", "#d7a85b", "indianred"), stringsAsFactors = F)


Plot.model=function(dataset, coef, type.feat, ind=1, NAME=NULL){
  load(dataset)
  
  time=TIME[logicalv[[ind]]]
  time[time==0]=0.1
  status=STATUS[logicalv[[ind]]]
  
  genelist=rownames(coef)
  
  gexp=gexp[logicalv[[ind]], colnames(gexp)%in%genelist,drop=F]
  immunoscore=immunoscore[logicalv[[ind]], colnames(immunoscore)%in%genelist,drop=F]
  samp=samp[logicalv[[ind]], colnames(samp)%in%genelist,drop=F]
  
  test.data=cbind(gexp, immunoscore, samp)
  
  if(!(all(genelist%in%colnames(test.data)))){
    warning(paste("Features not found from test data:", paste(genelist[!genelist%in%colnames(test.data)], collapse = ",")))
    coef=coef[rownames(coef)%in%colnames(test.data),,drop=F]
    genelist=rownames(coef)
  }
  
  m=test.data[,match(genelist, colnames(test.data))]
  riskPI=as.numeric(coef) %*% data.matrix(t(m))
  
  gene.annot=data.frame(type.feat[match(rownames(coef), type.feat[,1]),], "HR"=as.numeric(exp(coef)), stringsAsFactors = F)
  gene.annot=gene.annot[order(gene.annot$Type, -gene.annot$HR),]
  fm.m=t(m)
  fm.m=fm.m[match(gene.annot$Feature, rownames(fm.m)),]
  
  rownames(fm.m)[rownames(fm.m)%in%colnames(gexp)]=paste0("N:GEXP:", rownames(fm.m)[rownames(fm.m)%in%colnames(gexp)])
  rownames(fm.m)[rownames(fm.m)%in%colnames(immunoscore)]=paste0("N:GSVA:", rownames(fm.m)[rownames(fm.m)%in%colnames(immunoscore)])
  rownames(fm.m)[rownames(fm.m)%in%gene.annot$Feature[gene.annot$Type%in%c("Subtype")]]=paste0("B:SAMP:", rownames(fm.m)[rownames(fm.m)%in%gene.annot$Feature[gene.annot$Type%in%c("Subtype")]])
  rownames(fm.m)[rownames(fm.m)%in%gene.annot$Feature[gene.annot$Type%in%c("Clinical")]]=paste0("N:CLIN:", rownames(fm.m)[rownames(fm.m)%in%gene.annot$Feature[gene.annot$Type%in%c("Clinical")]])
  
  rownames(gene.annot)=rownames(fm.m)
  HR=prettyNum(signif(gene.annot[,3], 3))
  names(HR)=rownames(gene.annot)
  
  annotdf=data.frame("OS"=status, "RiskPI"=t(riskPI), stringsAsFactors = F)
  rownames(annotdf)=colnames(fm.m)
  
  plot.complexHM.fm(feats = rownames(gene.annot), text_annot = HR, feats.barplot = c("RiskPI"), split.columns = F, annotdf = annotdf, fm.f = fm.m, order_columns = colnames(fm.m)[order(as.numeric(riskPI))], use_raster = F, order_rows = F, NAME = NAME)

}

fun_forestplot=function(data, NAME="data", BOX=0.1,cex=2, colorv="black"){
  library(forestplot)
  
  txt.df=data.frame("Feature"=data$Feature, "Cohort"=gsub("_", " ", data$Name), "P"=signif(data$P, 2), stringsAsFactors = F)
  txt.df$Feature[duplicated(txt.df$Feature)]=""
  txt.df=rbind(c("Feature", "Cohort", "P"), txt.df)
  
  coef.df=data.frame("HR"=data$`exp(coef)`, "lower .95"=data$`lower .95`, "upper .95"=data$`upper .95`)
  coef.df=rbind(c(NA, NA, NA), coef.df)
  
  xticks=seq(ifelse(min(data$`lower .95`)<0.5, 0, 0.5), min(c(5, max(data$`upper .95`))), by = 0.5)
  
  attr(xticks, "labels") = xticks%in%seq(-4, 4, by=1)
  
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

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

#**************************************** training data DLBCL **************************************** 
load("Hemap_DLBCL_survival_data.Rdata")

# significant genes:
files=list.files(pattern = "tableS7")

# significant genes:
univariate.results=lapply(files[grepl("signif", files)], fread, data.table=F)
names(univariate.results)=files[grepl("signif", files)]

DLBCL=do.call(rbind, univariate.results[grep("DLBCL", names(univariate.results))])
DLBCL$Name=gsub("_RCHOP", "", DLBCL$Name)
DLBCL$Feature=gsub("-", ".", DLBCL$Feature)

# genes with FDR<0.2, beta to same direction and observed in 2 cohorts:
genes=table(DLBCL$Feature, DLBCL$`exp(coef)`>1)>1
genelist=rownames(genes)[rowSums(genes)==1]

# combine data:
time=TIME[logicalv[[5]]]
status=STATUS[logicalv[[5]]]

regression.data=cbind(gexp[logicalv[[5]], colnames(gexp)%in%genelist,], immunoscore[logicalv[[5]], colnames(immunoscore)%in%genelist,], samp[logicalv[[5]], colnames(samp)%in%genelist,])
all(genelist%in%colnames(regression.data))

# results_dlbcl=fun.cox.elasticnet(DATA_ORG = regression.data, time, status, summary.km = "3quantile", cores = 10, REPEATS = 100, percentage = 0, nfold = 10, min.elnet = 0, max.elnet = 0.1)
# save(results_dlbcl, file="Hemap_DLBCL_cox_datasets_filt_adj20_revision.Rdata")
load("Hemap_DLBCL_cox_datasets_filt_adj20_revision.Rdata")

# revision
# [1] 0.1
# [1] 0.0745611
# [1] 11.60345

annot=get(load("GSE98588_annot.Rdata"))

pdf("FigS7D_DLBCL_model_cox.pdf", height = 2.5, width = 2)
riskInd_Chapuy=validate.model("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients)
riskInd_Hemap=validate.model("Hemap_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, 5)
dev.off()

DLBCL_risk1=validate.model.cox("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, filt=annot$IPI%in%c(0, 1, 2), NAME = "DLBCL IPI 0-2")
DLBCL_risk2=validate.model.cox("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, filt=annot$IPI%in%c(3), NAME = "DLBCL IPI 3")
DLBCL_risk3=validate.model.cox("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, filt=annot$IPI%in%c(4,5), NAME = "DLBCL IPI 4-5")
DLBCL_risk4=validate.model.cox("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, filt=annot$COO_byGEP=="ABC", NAME = "DLBCL ABC")
DLBCL_risk5=validate.model.cox("GSE98588_DLBCL_survival_data.Rdata", results_dlbcl$coefficients, filt=annot$COO_byGEP=="GCB", NAME = "DLBCL GCB")

type.feat=unique(DLBCL[,c("Feature", "Type")])
type.feat=type.feat[order(type.feat[,2]),]

Plot.model(dataset = "GSE98588_DLBCL_survival_data.Rdata", coef = results_dlbcl$coefficients, type.feat = type.feat, ind = 1,NAME = "Fig7F_Chapuy_model")
# Plot.model(dataset = "Reddy_DLBCL_survival_data.Rdata", coef = results_dlbcl$coefficients, type.feat = type.feat, ind = 1,NAME = "Reddy_model")

#**************************************** training data MM **************************************** 
load("Hemap_MM_survival_data.Rdata")

# significant genes:
files=list.files(pattern = "tableS7")

# significant genes:
univariate.results=lapply(files[grepl("signif", files)], fread, data.table=F)
names(univariate.results)=files[grepl("signif", files)]

MM=do.call(rbind, univariate.results[grep("MM", names(univariate.results))])
MM$Feature=gsub("-", ".", MM$Feature)

# genes with FDR<0.2, beta to same direction and observed in 2 cohorts:
# use GSE19784 Cancer_Myeloma as training and GSE16716,GSE24080 Cancer_Myeloma as test set
genes=table(MM$Feature, MM$`exp(coef)`>1)>1
genelist=rownames(genes)[rowSums(genes)==1]
genelist=genelist[!genelist%in%"WHSC1_FGFR3_Ig"]

# combine data:
ind=1
time=TIME[logicalv[[ind]]]
time[time==0]=0.1
status=STATUS[logicalv[[ind]]]

regression.data=cbind(gexp[logicalv[[ind]], colnames(gexp)%in%genelist,drop=F], immunoscore[logicalv[[ind]], colnames(immunoscore)%in%genelist,drop=F], samp[logicalv[[ind]], colnames(samp)%in%genelist,drop=F])
all(genelist%in%colnames(regression.data))

# results_MM=fun.cox.elasticnet(DATA_ORG = regression.data, time, status, summary.km = "3quantile", cores = 10, REPEATS = 100, percentage = 0, nfold = 10)
# save(results_MM, file="Hemap_MM_cox_datasets_filt_adj20_revision.Rdata")
# [1] 0.05
# [1] 0.1089963
# [1] 14.10799

load("Hemap_MM_cox_datasets_filt_adj20_revision.Rdata")

pdf("FigS7C_MM_model_cox.pdf", height = 2.5, width = 2)
riskInd_CoMMpass=validate.model(dataset = "Hemap_MM_survival_data.Rdata", coef = results_MM$coefficients, ind=1, summary="75th_25th_percentile")
riskInd_CoMMpass=validate.model("CoMMpass_survival_data.Rdata", results_MM$coefficients, summary="75th_25th_percentile", months = T)
dev.off()

type.feat=unique(MM[,c("Feature", "Type")])
type.feat=type.feat[order(type.feat[,2]),]

Plot.model(dataset = "CoMMpass_survival_data.Rdata", coef = results_MM$coefficients, type.feat = type.feat, ind = 1,NAME = "Fig7E_CoMMpass_model")

load("CoMMpass_MM_subtypes.Rdata")

subtype=coordinates.subtype[match(colnames(riskInd_CoMMpass)[order(-riskInd_CoMMpass)],coordinates.subtype$ID),]
subtype$subtype[subtype$cluster=="CGA_Prolif"]="CGA_Prolif"

col=data.frame("subtype"=c("CCND1_Ig", "WHSC1_FGFR3_Ig", "Hyperdiploid_gain11q", "Hyperdiploid_gain1q", "MAF_Ig", "TRAF3_Aberrated", "CGA_Prolif"),
"color"=c("#e41a1b", "#357eb8", "#5eb45b", "#9b53a4", "#ff7d00", "#f8f875", "darkred"), stringsAsFactors = F)

load("CoMMpass_survival_data.Rdata")
MM_risk1=validate.model.cox("CoMMpass_survival_data.Rdata", results_MM$coefficients, filt=samp$ISS1==1, NAME = "CoMMpass ISS 1")
MM_risk2=validate.model.cox("CoMMpass_survival_data.Rdata", results_MM$coefficients, filt=samp$ISS2==1, NAME = "CoMMpass ISS 2")
MM_risk3=validate.model.cox("CoMMpass_survival_data.Rdata", results_MM$coefficients, filt=samp$ISS3==1, NAME = "CoMMpass ISS 3")

#**************************************** training data AML **************************************** 
load("Hemap_AML_survival_data.Rdata")

# significant genes:
files=list.files(pattern = "tableS7")

# significant genes:
univariate.results=lapply(files[grepl("signif", files)], fread, data.table=F)
names(univariate.results)=files[grepl("signif", files)]

AML=do.call(rbind, univariate.results[grep("AML", names(univariate.results))])
AML$Feature=gsub("-", ".", AML$Feature)

# genes with FDR<0.2, beta to same direction and observed in 2 cohorts:
genes=table(AML$Feature, AML$`exp(coef)`>1)>1
genelist=rownames(genes)[rowSums(genes)==1]
genelist=gsub("-", ".", genelist)

# combine data:
time=TIME[logicalv[[1]]]
time[time==0]=0.1
status=STATUS[logicalv[[1]]]

regression.data=cbind(gexp[logicalv[[1]], colnames(gexp)%in%genelist,drop=F], immunoscore[logicalv[[1]], colnames(immunoscore)%in%genelist,drop=F], samp[logicalv[[1]], colnames(samp)%in%genelist,drop=F])
all(genelist%in%colnames(regression.data))

# results_AML=fun.cox.elasticnet(DATA_ORG = regression.data, time, status, summary.km = "3quantile", cores = 10, REPEATS = 100, percentage = 0, nfold = 5)

# model revision
# [1] 0.06
# [1] 0.1982066
# [1] 11.88285

# save(results_AML, file="Hemap_AML_cox_datasets_filt_adj20_revision.Rdata")
load("Hemap_AML_cox_datasets_filt_adj20_revision.Rdata")
exp(results_AML$coefficients)

pdf("FigS7E_Risk_AML_validation.pdf", height = 2.5, width = 2)
riskInd_Hemap_AML=validate.model("Hemap_AML_survival_data.Rdata", results_AML$coefficients, summary = "75th_25th_percentile")
riskInd_BeatAML=validate.model("BeatAML_survival_data.Rdata", results_AML$coefficients, months=T, summary = "75th_25th_percentile")
riskInd_TCGA_AML=validate.model("TCGA_AML_survival_data.Rdata", results_AML$coefficients, summary = "75th_25th_percentile")
dev.off()

type.feat=unique(AML[,c("Feature", "Type")])
type.feat=type.feat[order(type.feat[,2]),]

# Plot.model(dataset = "Hemap_AML_survival_data.Rdata", coef = results_AML$coefficients, type.feat = type.feat, ind = 1,NAME = "Hemap_AML_model")
Plot.model(dataset = "BeatAML_survival_data.Rdata", coef = results_AML$coefficients, type.feat = type.feat, ind = 1,NAME = "FigS7G_BeatAML_model")
# Plot.model(dataset = "TCGA_AML_survival_data.Rdata", coef = results_AML$coefficients, type.feat = type.feat, ind = 1,NAME = "TCGA_AML_model")

load("BeatAML_survival_data.Rdata")
annot=get(load("BeatAML_fm_annot.Rdata"))

AML_risk1=validate.model.cox("BeatAML_survival_data.Rdata", results_AML$coefficients, filt=grepl("Adverse", annot$ELN2017), NAME = "BeatAML Adverse")
AML_risk2=validate.model.cox("BeatAML_survival_data.Rdata", results_AML$coefficients, filt=annot$ELN2017=="Intermediate", NAME = "BeatAML Intermediate")
AML_risk3=validate.model.cox("BeatAML_survival_data.Rdata", results_AML$coefficients, filt=grepl("Favorable", annot$ELN2017), NAME = "BeatAML Favorable")

dat=rbind(DLBCL_risk1, DLBCL_risk2, DLBCL_risk3, DLBCL_risk4, DLBCL_risk5, MM_risk1, MM_risk2, MM_risk3, AML_risk1, AML_risk2, AML_risk3)

pdf("FigS7F_Subset_data_model.pdf", width = 4, height = 3)
fun_forestplot(dat, "AML, MM, DLBCL", BOX=0.5, cex=0.5, colorv ="black")
dev.off()