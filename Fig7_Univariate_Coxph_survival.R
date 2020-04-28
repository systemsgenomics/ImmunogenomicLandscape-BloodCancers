GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/statistics_wrappers.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/compute.pairwise.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))

library(RColorBrewer)
library(survival)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(survMisc)
library(survminer)
library(plyr)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

#******************************************** Hemap *********************************************

annot = get(load("Hemap_immunology_Annotations.Rdata"))
annot=annot[!is.na(annot$OS_Time),]

# survival time and status
TIME=annot$OS_Time
STATUS=as.numeric(annot$OS_Status)
TIME2=annot$PFS_Time
STATUS2=as.numeric(annot$PFS_Status)

#************************************** Process gene expression data ************************************
profile=data.matrix(get(load("mixtureM_profile.Rdata")))
profile[profile==-1] = 0
profile2=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]
data=t(get(load("data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]

# take only high expressed into account, for CGA score
profile2[data<5]=0

#********************************************************************************************************


#********************************** Make necessary gene lists ***************************************

# other lists from each analysis
cga = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)[,1]
co.stim=fread("costim_ligands_final.txt", data.table = F)

hlagenes=c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1", "CIITA")
immunoscore=c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")

# for AML:
MDS=get(load("MDS_genesets.Rdata"))

# significant per disease
microenv=get(load("Hemap_cytolytic_correlated_genes_TableS2_onlysignif.Rdata"))

#****************************************************************************************************
# data for survival
df=data.frame("time"=annot$OS_Time, "status"=annot$OS_Status, "HLAI"=annot$HLAIScore, "HLAII"=annot$HLAIIScore, "CytolyticScore"=annot$CytolyticScore, "Freq.CGA"=colSums(profile2[rownames(profile2)%in%cga,]), stringsAsFactors = F)

#********************************* just AML *********************************
filterv = annot$subclasses%in%"AML"&annot$CELLS_SORTED==0
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

filterv = annot$subclasses%in%"AML"&annot$CELLS_SORTED==0
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="Hemap_AML"
logicalv=append(logicalv2, logicalv)
logicalv=logicalv[lapply(logicalv, sum)>2]

# HLA, cytolytic
DATA=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "CytolyticScore"=scale(annot$CytolyticScore))
aml_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
aml_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = F, pretty=F))

DATA_clin=data.frame("Age"=annot$AGE, "Gender"=annot$GENDER)
aml_res_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))
aml_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = F, pretty=F))

DATA_hla=data.frame(t(data[hlagenes,]))
aml_res_HLA=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_hla,TIME,STATUS, univariate = T, pretty=F))

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
aml_res_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# MDS:
DATA=data.frame(scale(t(data[rownames(data)%in%c(MDS$MDS_signature_all_filt),])))
aml_res_MDS=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# microenvironment
ME=microenv[microenv$disease%in%"AML",]
DATA=data.frame(scale(t(data[rownames(data)%in%ME$gene,])))
aml_res_microenvironment=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
aml_res_microenvironment[aml_res_microenvironment$Adj.P<0.1,]

# subtype:
load("Hemap_AML_subtypes.Rdata")
coordinates.subtype=coordinates.subtype[!is.na(coordinates.subtype$subtype),]

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
DATA=data.frame(do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1)
colnames(DATA)=gsub("-", ".", unique(coordinates.subtype$subtype))

aml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(aml_res, aml_res_subtype, aml_res_costim, aml_res_MDS, aml_res_microenvironment, aml_res_clin)
tableS7=tableS7[tableS7$Name%in%c("Hemap_AML"),]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%unique(coordinates.subtype$subtype)]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%aml_res_MDS[,1]]="MDS-signature gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"
tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_AML.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_AML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "CytolyticScore"=scale(annot$CytolyticScore))

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
subtypes=do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1 
colnames(subtypes)=gsub("-", ".", unique(coordinates.subtype$subtype))
samp=data.frame(subtypes, "Age"=annot$AGE, "Gender"=annot$GENDER)

save(list = c("gexp","immunoscore", "TIME", "STATUS", "logicalv", "samp"), file="Hemap_AML_survival_data.Rdata")

#************************************* just MM *************************************
filterv = annot$subclasses%in%"Cancer_Myeloma"&!is.na(annot$MM_ISS)
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

# just MM
filterv = annot$subclasses%in%"Cancer_Myeloma"&!is.na(annot$MM_ISS)
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="MM_all"
logicalv=append(logicalv2, logicalv)
names(logicalv)=c("Hemap_MM", "GSE19784_Hemap_MM", "GSE16716,GSE24080_Hemap_MM")

data.test=data.frame("time"=TIME[logicalv[[1]]], "status"=STATUS[logicalv[[1]]])

ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
           xlab = "months", 
           ylab = "Overall survival probability")


DATA=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "ISS3"=(annot$MM_ISS==3)*1, "ISS1"=(annot$MM_ISS==1)*1, "Freq.CGA"=as.numeric(df$Freq.CGA),"Age"=annot$AGE, "Gender"=annot$GENDER, check.names = F)

mm_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
# mm_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA[,3:5],TIME,STATUS, univariate = F, pretty=F))

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
mm_res_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

DATA=data.frame(scale(t(data[rownames(data)%in%cga,])))
mm_res_antigen=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# subtype:
load("Hemap_MM_subtypes.Rdata")

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
DATA=data.frame(do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1)
colnames(DATA)=gsub("-", "_", unique(coordinates.subtype$subtype))

mm_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(mm_res, mm_res_costim, mm_res_antigen, mm_res_subtype)
tableS7=tableS7[tableS7$Name%in%c("GSE16716,GSE24080_Hemap_MM"),]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")


tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%gsub("-", "_", unique(coordinates.subtype$subtype))]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"
tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_MM_GSE24080.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_MM_GSE24080_signif.tsv", sep="\t")


# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(mm_res, mm_res_costim, mm_res_antigen)
tableS7=tableS7[tableS7$Name%in%c("GSE19784_Hemap_MM"),]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%gsub("-", "_", unique(coordinates.subtype$subtype))]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"
tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_MM_GSE19784.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_MM_GSE19784_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "Freq.CGA"=as.numeric(df$Freq.CGA), check.names = F)

samples=lapply(unique(coordinates.subtype$subtype), function(type)coordinates.subtype$ID[coordinates.subtype$subtype%in%type])
subtype=do.call(cbind, lapply(samples, function(id)annot$GSM.identifier..sample.%in%id))*1
colnames(subtype)=c(gsub("-", ".", unique(coordinates.subtype$subtype)), colnames(subtype))
samp=data.frame(subtype, "Age"=annot$AGE, "Gender"=annot$GENDER, "ISS3"=(annot$MM_ISS==3)*1, "ISS1"=(annot$MM_ISS==1)*1)

save(list = c("gexp","immunoscore", "TIME", "STATUS", "logicalv", "samp"), file="Hemap_MM_survival_data.Rdata")

#***************************** just DLBCL ************************************
# remove genes not in chapyu dataset:
gexp=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_symbol_rma_normalized.Rdata"))

ME=microenv[microenv$disease%in%"DLBCL"&microenv$category=="Stromal/cancer gene (Rho > 0)"&microenv$Rho>0.4,]
microenv_dlbcl=ME[,1]
microenv_dlbcl=microenv_dlbcl[microenv_dlbcl%in%rownames(gexp)]

filterv = annot$subclasses%in%"BCL_DLBCL"&!is.na(annot$dlbcl_ipi)
logicalv=get.logical(annovector = list(annot$GSE.identifier..experiment.), filterv = filterv)
names(logicalv)=unique(paste(names(logicalv), annot$subclasses[filterv]))

filterv = annot$subclasses%in%"BCL_DLBCL"
logicalv2=get.logical(annovector = list(annot$subclasses), filterv = filterv)
names(logicalv2)="BCL_DLBCL_all"

logicalv=append(logicalv2, logicalv)
logicalv=logicalv[lapply(logicalv, sum)>2]

RCHOP=list(annot$Chemotherapy_RCHOP==1&!is.na(STATUS)&!TIME==0)
CHOP=list(annot$Chemotherapy_CHOP==1&!is.na(STATUS)&!TIME==0)
names(RCHOP)="Hemap_DLBCL_RCHOP"
names(CHOP)="Hemap_DLBCL_CHOP"

logicalv <- append(logicalv, append(RCHOP, CHOP))

# HLA, cytolytic
DATA=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "CytolyticScore"=scale(annot$CytolyticScore), "IPI_0to1"=(annot$dlbcl_ipi%in%c(0,1))*1, "IPI_4to5"=(annot$dlbcl_ipi%in%c(4,5))*1, "Freq.CGA"=as.numeric(df$Freq.CGA),  "ABC"=(grepl("ABC", annot$tbLY))*1, "GCB"=(grepl("GCB", annot$tbLY))*1,check.names = F)

dlbcl_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
# dlbcl_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA[,4:6],TIME,STATUS, univariate = F, pretty=F))

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
dlbcl_res_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# microenv
DATA=data.frame(scale(t(data[rownames(data)%in%c(microenv_dlbcl),])))
dlbcl_res_microenv=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

DATA=data.frame(scale(t(data[rownames(data)%in%cga,])))
dlbcl_res_antigen=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(dlbcl_res, dlbcl_res_costim, dlbcl_res_microenv, dlbcl_res_antigen)
tableS7=tableS7[tableS7$Name%in%c("Hemap_DLBCL_RCHOP"),]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[tableS7$Feature%in%c("ABC", "GCB")]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_DLBCL_RCHOP.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_DLBCL_RCHOP_signif.tsv", sep="\t")


# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(dlbcl_res, dlbcl_res_costim, dlbcl_res_microenv, dlbcl_res_antigen) 
tableS7=tableS7[tableS7$Name%in%c("Hemap_DLBCL_CHOP"),]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[tableS7$Feature%in%c("ABC", "GCB")]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_DLBCL_CHOP.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Hemap_DLBCL_CHOP_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(annot$HLAIScore), "HLAII"=scale(annot$HLAIIScore), "Freq.CGA"=as.numeric(df$Freq.CGA), check.names = F)

samp=data.frame("Age"=annot$AGE, "Gender"=annot$GENDER, "IPI_0to1"=(annot$dlbcl_ipi%in%c(0,1))*1, "IPI_4to5"=(annot$dlbcl_ipi%in%c(4,5))*1,"ABC"=(grepl("ABC", annot$tbLY))*1, "GCB"=(grepl("GCB", annot$tbLY))*1)

save(list = c("gexp","immunoscore","samp", "TIME", "STATUS", "logicalv"), file="Hemap_DLBCL_survival_data.Rdata")


#******************************************** TCGA AML *****************************************************
fm_org=get(load("TCGA_AML_FM_DUFVA.Rdata"))
fm=fm_org[,!is.na(fm_org["N:SAMP:CytolyticScore",])]

risks=c("N:CLIN:Age:::::",
        "C:CLIN:acute_myeloid_leukemia_calgb_cytogenetics_risk_category:::::" ,
        "C:CLIN:FISH_test_component:::::",
        "B:GNAB:NPM1:chr5:170814708:170837888:+:y_n_somatic",
        "B:GNAB:FLT3:chr13:28577411:28674729:-:y_n_somatic",
        "B:GNAB:CEBPA:chr19:33790840:33793430:-:y_n_somatic",
        "B:GNAB:TP53:chr17:7565097:7590863:-:y_n_somatic")

df=t(fm[rownames(fm)%in%risks,])
colnames(df)=do.call(rbind, strsplit(colnames(df), ":"))[,3]

data=data.matrix(fm[grepl("GEXP", rownames(fm)),])
rownames(data)=make.unique(do.call(rbind, strsplit(rownames(data), ":"))[,3])


OS=as.numeric(fm["N:CLIN:OS.months..3.31.12:::::",])
TIME=OS
STATUS=as.numeric(fm["C:CLIN:vital_status_TCGA_paper:::::",]=="DECEASED")
PFS=as.numeric(fm["N:CLIN:EFS.months....4.30.13:::::",])

DATA=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])))

logicalv=list("TCGA_AML"=rep(T, dim(DATA)[1]))

TCGA_AML_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
# TCGA_AML_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = F, pretty=F))

# clinical
DATA_clin=data.frame("Age"=scale(as.numeric(fm["N:CLIN:Age:::::",])), "Blast.percentage"=scale(as.numeric(fm["N:CLIN:X.BM.Blast:::::",])), "Gender"=as.character(fm["C:CLIN:Sex:::::",]))
TCGA_AML_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
TCGA_AML_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

DATA_hla=data.frame(scale(t(data[rownames(data)%in%hlagenes,])))
TCGA_AML_HLA=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_hla,TIME,STATUS, univariate = T, pretty=F))

# MDS
DATA=data.frame(scale(t(data[rownames(data)%in%c(MDS$MDS_signature_all_filt),])))
TCGA_AML_MDS=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# mikroenvironment
ME=microenv[microenv$disease%in%"AML",]
DATA=data.frame(scale(t(data[rownames(data)%in%ME$gene,])))
TCGA_res_microenvironment=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
TCGA_res_microenvironment[TCGA_res_microenvironment$Adj.P<0.1,]

# subtype:
load("TCGA_AML_subtypes.Rdata")

DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)

TCGA_aml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(TCGA_AML_res, TCGA_AML_costim, TCGA_AML_MDS, TCGA_res_microenvironment, TCGA_aml_res_subtype, TCGA_AML_clin)
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%unique(coordinates.subtype$subtype)]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%aml_res_MDS[,1]]="MDS-signature gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_TCGA_AML.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_TCGA_AML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])))

samp=cbind(DATA_clin, DATA_subtypes)

save(list = c("gexp","immunoscore", "samp", "TIME", "STATUS", "logicalv"), file="TCGA_AML_survival_data.Rdata")

#********************************************************* Compass MM *********************************************************
fm=get(load("MM_COMPASS_FM.Rdata"))
annot=get(load("MM_COMPASS_ANNOT.Rdata"))

annot=annot[match(colnames(fm), rownames(annot)),]

data=fm[grepl("N:GEXP:", rownames(fm)),]
rownames(data)=gsub("N:GEXP:", "", rownames(data))

data.mut=fm[grepl("B:GNAB:", rownames(fm)),]
rownames(data.mut)=gsub("B:GNAB:", "", rownames(data.mut))
data.mut.filt=data.mut[!rowSums(data.mut, na.rm = T)<10,]

TIME=as.numeric(fm["N:CLIN:OS",])
STATUS=as.numeric(fm["B:CLIN:STATUS",])

STATUS[TIME>1825&!is.na(STATUS)&STATUS==1]=0 # change to 4year survival, 5 year sharp drop
TIME[TIME>1825]=1825

# compute antigen scores
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
t.df=t.df[order(t.df[,3]),]
genelist=t.df[,1]

data.test=data.frame("time"=TIME[logicalv[[1]]]*0.0328767, "status"=STATUS[logicalv[[1]]])

ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
           xlab = "months", 
           ylab = "Overall survival probability")

expressed_testis_num=as.numeric(fm["N:SAMP:nCGA",])

l.regulon.gene=regulon.feats(fm, co.stim[,1])

hlagenes=c("B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1", "CIITA")

DATAcostim=scale(data.frame(t(gexp[rownames(gexp)%in%co.stim[,1],])))
DATAhla=scale(data.frame(t(gexp[rownames(gexp)%in%hlagenes,])))
DATAcga=scale(data.frame(t(gexp[rownames(gexp)%in%t.df[,1],])))

DATA=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "ISS1"=as.numeric(fm["B:CLIN:R_ISS_1",]),"ISS2"=as.numeric(fm["B:CLIN:R_ISS_2",]),"ISS3"=as.numeric(fm["B:CLIN:R_ISS_3",]), "Freq.CGA"=as.numeric(expressed_testis_num), stringsAsFactors = F)

logicalv=list(!is.na(expressed_testis_num))

names(logicalv)="CoMMPass"

# plot overall survival
r=lapply(seq(logicalv), function(i){
  data.test=data.frame("time"=TIME[logicalv[[i]]], "status"=STATUS[logicalv[[i]]])
  
  ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
             xlab = "months", 
             ylab = "Overall survival probability",title=(names(logicalv)[i])
  )
})

names(r)=names(logicalv)

pdf("CoMMpass_MM_cohorts.pdf", width =5, height = ceiling(length(r)/2)*2.75)
plots.together=arrange_ggsurvplots(r, print = TRUE, ncol = 2, nrow = ceiling(length(r)/2))
dev.off()

commpass_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
commpass_res_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = F, pretty=F))

# costim, no mut
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
commpass_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# mutations
mut=t(data.mut.filt[rownames(data.mut.filt)%in%c(t.df[,1]),])
colnames(mut)=paste0("MUT:", colnames(mut))
DATA=data.frame(scale(t(data[rownames(data)%in%c(t.df[,1]),])), mut)

commpass_antig=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# subtype:
load("CoMMpass_MM_subtypes.Rdata")
coordinates.subtype=coordinates.subtype[match(colnames(fm), coordinates.subtype$ID),]

coordinates.subtype$subtype[coordinates.subtype$cluster=="CGA_Prolif"&!is.na(coordinates.subtype$cluster)]="CGA_Prolif"
DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)
commpass_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(commpass_res, commpass_costim, commpass_antig, commpass_subtype)
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

# annotate these genes, needed later:
tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[genelist_signif[,1]%in%commpass_subtype$Feature]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_CoMMpass.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_CoMMpass_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "Freq.CGA"=as.numeric(expressed_testis_num), stringsAsFactors = F)

samp=cbind(data.frame("ISS1"=as.numeric(fm["B:CLIN:R_ISS_1",]),"ISS2"=as.numeric(fm["B:CLIN:R_ISS_2",]),"ISS3"=as.numeric(fm["B:CLIN:R_ISS_3",]), stringsAsFactors = F), DATA_subtypes)

save(list = c("gexp","immunoscore", "samp", "TIME", "STATUS", "logicalv"), file="CoMMpass_survival_data.Rdata")

#********************************************************* Chapyu DLBCL *********************************************************

fm=get(load("GSE98588_fm.Rdata"))
annot=get(load("GSE98588_annot.Rdata"))

TIME=as.numeric(fm["N:CLIN:OS",])
PFS=as.numeric(fm["N:CLIN:PFS",])
STATUS=as.numeric(fm["B:CLIN:OS_STAT",])
STATUS2=as.numeric(fm["B:CLIN:PFS_STAT",])

# compute antigen scores
genelist=t.df[,1]

data=get(load("GSE98588_symbol_rma_normalized.Rdata"))
colnames(data)=gsub("_DLBCL", "", colnames(data))
gexp=data

load("GSE98588_DLBCL_mixtureM_profile.Rdata")
profile[profile==-1] = 0
profile[data.matrix(data)<5]=0

expressed_testis_num=as.numeric(colSums(profile[rownames(profile)%in%genelist,]))

cnv_annot=fread("41591_2018_16_MOESM8_ESM_CNV_ANNOT.txt", data.table=F)

l.regulon.gene=regulon.feats(fm, c(co.stim[,1], hlagenes), cnv_annot)

ME=microenv[microenv$disease%in%"DLBCL"&microenv$category=="Stromal/cancer gene (Rho > 0)"&microenv$Rho>0.4,]

# all costim-cytolytic-CGA correlated mutations:
costim_feats=data.matrix(fm[rownames(fm)%in%unlist(l.regulon.gene),])
rownames(costim_feats)=sapply(rownames(costim_feats), function(n){
  if(!grepl("CNVR", n))return(n)
  a=names(l.regulon.gene)[sapply(l.regulon.gene, function(a)any(a%in%n))]
  if(length(a))paste0(n,"@", paste(a, collapse=","))
})

DATAmut=data.frame(t(costim_feats[grepl("GNAB|CNVR", rownames(costim_feats)),]))

DATA=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])), "Freq.CGA"=as.numeric(expressed_testis_num), stringsAsFactors = F)

logicalv=list(rep(T, dim(DATA)[1]))
names(logicalv)="GSE98588_DLBCL"

GSE98588_DLBCL_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
# GSE98588_DLBCL_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA[,4:6],TIME,STATUS, univariate = F, pretty=F))

# Clinical
DATA_clin=data.frame("IPI_0to1"=annot$IPI%in%c(0,1)*1,"IPI_4to5"=annot$IPI%in%c(4,5)*1, "ABC"=(annot$COO_byGEP=="ABC")*1 ,"GCB"=(annot$COO_byGEP=="GCB")*1,  stringsAsFactors = F)
GSE98588_DLBCL_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# immuno-editing:
GSE98588_DLBCL_immunoediting=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATAmut,TIME,STATUS, univariate = T, pretty=F))
GSE98588_DLBCL_immunoediting$Feature=gsub("B.GNAB.", "MUT:", GSE98588_DLBCL_immunoediting$Feature)
GSE98588_DLBCL_immunoediting$Feature=gsub("HLA.", "HLA-", GSE98588_DLBCL_immunoediting$Feature)
GSE98588_DLBCL_immunoediting$Feature=gsub("B.CNVR.", "", GSE98588_DLBCL_immunoediting$Feature)

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])), DATAmut, "costim.mut"=(rowSums(DATAmut)))
GSE98588_DLBCL_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# monocyte
DATA=data.frame(scale(t(data[rownames(data)%in%c(ME[,1]),])))
GSE98588_DLBCL_microenvironment=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

DATA=data.frame(scale(t(data[rownames(data)%in%c(t.df[,1]),])))
GSE98588_DLBCL_antigen=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(GSE98588_DLBCL_res,GSE98588_DLBCL_clin, GSE98588_DLBCL_costim, GSE98588_DLBCL_microenvironment, GSE98588_DLBCL_antigen, GSE98588_DLBCL_immunoediting) 
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")


tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%colnames(DATAmut)]="Immune Editing mutation"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[tableS7$Feature%in%c("ABC", "GCB")]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_GSE98588_DLBCL.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_GSE98588_DLBCL_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAI"=scale(as.numeric(fm["N:SAMP:HLAIScore",])), "HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])), "Freq.CGA"=as.numeric(expressed_testis_num), stringsAsFactors = F)

samp=DATA_clin

save(list = c("gexp","immunoscore", "samp", "TIME", "STATUS", "logicalv"), file="GSE98588_DLBCL_survival_data.Rdata")

#********************************************************* Reddy DLBCL *********************************************************
fm=get(load("REDDY_DLBCL_fm.Rdata"))
annot=get(load("REDDY_DLBCL_annot.Rdata"))

data=data.matrix(fm[grepl("GEXP", rownames(fm)),])
rownames(data)=gsub("N:GEXP:", "", rownames(data))

TIME=as.numeric(annot$Overall.Survival.years)
STATUS=(as.numeric(annot$Censored)==0)*1 # 0 indicates no censoring, meaning that the death was observed, whereas a 1 indicates that the patient was alive

ME=microenv[microenv$disease%in%"DLBCL"&microenv$category=="Stromal/cancer gene (Rho > 0)"&microenv$Rho>0.4,]

l.regulon.gene=regulon.feats(fm, c(co.stim[,1], hlagenes))

# all costim-cytolytic-CGA correlated mutations:
costim_feats=data.matrix(fm[rownames(fm)%in%unlist(l.regulon.gene),])
rownames(costim_feats)=sapply(rownames(costim_feats), function(n){
  if(!grepl("CNVR", n))return(n)
  a=names(l.regulon.gene)[sapply(l.regulon.gene, function(a)any(a%in%n))]
  # if(length(a))paste0(n,"@", paste(a, collapse=","))
})

DATAmut=data.frame(t(costim_feats[grepl("GNAB|CNVR", rownames(costim_feats)),]))

DATA=data.frame("HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])), stringsAsFactors = F)

logicalv=list("Reddy_DLBCL"=!is.na(data[1,]), "Reddy_DLBCL_ABC"=!is.na(data[1,])&annot$ABC.GCB..RNAseq.=="ABC", "Reddy_DLBCL_GCB"=!is.na(data[1,])&annot$ABC.GCB..RNAseq.=="GCB")

Reddy_DLBCL_res=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
Reddy_DLBCL_multivar=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = F, pretty=F))

# Clinical
DATA_clin=data.frame("IPI_0to1"=annot$IPI%in%c(0,1)*1,"IPI_4to5"=annot$IPI%in%c(4,5)*1, "ABC"=(annot$ABC.GCB..RNAseq.=="ABC")*1 ,"GCB"=(annot$ABC.GCB..RNAseq.=="GCB")*1,  stringsAsFactors = F)
Reddy_DLBCL_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# immuno-editing:
Reddy_DLBCL_immunoediting=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATAmut,TIME,STATUS, univariate = T, pretty=F))
Reddy_DLBCL_immunoediting$Feature=gsub("B.GNAB.", "MUT:", Reddy_DLBCL_immunoediting$Feature)
Reddy_DLBCL_immunoediting$Feature=gsub("HLA.", "HLA-", Reddy_DLBCL_immunoediting$Feature)
Reddy_DLBCL_immunoediting$Feature=gsub("N.CNVR.", "", Reddy_DLBCL_immunoediting$Feature)

# costim
DATA=data.frame(scale(t(data[rownames(data)%in%c(co.stim[,1]),])))
Reddy_DLBCL_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# microenvironment
DATA=data.frame(scale(t(data[rownames(data)%in%c(ME[,1]),])))
Reddy_DLBCL_microenvironment=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(Reddy_DLBCL_res,Reddy_DLBCL_clin[!is.na(Reddy_DLBCL_clin$`exp(coef)`),], Reddy_DLBCL_immunoediting, Reddy_DLBCL_costim, Reddy_DLBCL_microenvironment) 
tableS7=tableS7[tableS7$Name%in%"Reddy_DLBCL",]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")


tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="Clinical", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%Reddy_DLBCL_immunoediting$Feature]="Immune Editing mutation"
genelist_signif$type[genelist_signif[,1]%in%cga]="CGA"
genelist_signif$type[tableS7$Feature%in%c("ABC", "GCB")]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Reddy_DLBCL.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_Reddy_DLBCL_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(data)))
immunoscore=data.frame("HLAII"=scale(as.numeric(fm["N:SAMP:HLAIIScore",])), "CytolyticScore"=scale(as.numeric(fm["N:SAMP:CytolyticScore",])), stringsAsFactors = F)

samp=DATA_clin

save(list = c("gexp","immunoscore", "samp", "TIME", "STATUS", "logicalv"), file="Reddy_DLBCL_survival_data.Rdata")


#************************************ beatAML:
load("BeatAML_fm.Rdata")
annot=get(load("BeatAML_fm_annot.Rdata"))

gexp=fm[grepl("GEXP", rownames(fm)),]
rownames(gexp)=gsub("N:GEXP:", "", rownames(gexp))

annot$vitalStatus[annot$vitalStatus=="Unknown"]=NA
TIME=as.numeric(annot$overallSurvival)
STATUS=as.numeric(annot$vitalStatus=="Dead")

STATUS[TIME>1825&!is.na(STATUS)&STATUS==1]=0 # change to 5year survival
TIME[TIME>1825]=1825
TIME[TIME==0]=0.1

filtv=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(TIME)|is.na(STATUS)|is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv2=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("CMP-like","MDS-like","Monocyte-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv3=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("MDS-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv4=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("CMP-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")
filtv5=annot$specimenType%in%"Bone Marrow Aspirate"&!(is.na(annot$TCGA_coord))&annot$TCGA_coord%in%c("Monocyte-like")&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown")

logicalv=list("BeatAML"=!(is.na(TIME)|is.na(STATUS)|is.na(annot$TCGA_coord))&annot$causeOfDeath%in%c("Alive", "Dead-Disease", "Dead-Unknown"), "beatAML_BMonly"=filtv, "normal_karyotype"=filtv2, "MDS-like"=filtv3,  "CMP-like"=filtv4,  "Monocyte-like"=filtv5)

# plot overall survival
r=lapply(seq(logicalv), function(i){
  ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.frame("time"=TIME[logicalv[[i]]], "status"=STATUS[logicalv[[i]]])), 
             xlab = "months", 
             ylab = "Overall survival probability",title=(names(logicalv)[i])
  )
})

names(r)=names(logicalv)

pdf("BeatAML_cohorts.pdf", width =5, height = ceiling(length(r)/2)*2.75)
plots.together=arrange_ggsurvplots(r, print = TRUE, ncol = 2, nrow = ceiling(length(r)/2))
dev.off()

data.test=data.frame("time"=TIME[logicalv[[1]]]*0.0328767, "status"=STATUS[logicalv[[1]]])

ggsurvplot(survfit(Surv(time, status) ~ 1, data = data.test), 
           xlab = "months", 
           ylab = "Overall survival probability")

# HLA, cytolytic
DATA=data.frame(scale(t(fm[grepl("N:SAMP:.*.Score", rownames(fm)),])), check.names = F)
colnames(DATA)=c("HLAI", "HLAII", "CytolyticScore")
beatAML=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))

# Clinical:
DATA_clin=data.frame(t(fm[grepl("ELN2017", rownames(fm)),]),"Age"=as.numeric(fm[grepl("ageAtDiagnosis", rownames(fm)),,drop=F]), t(fm[grepl("BM_Transplant", rownames(fm)),,drop=F]),t(fm[grepl("B:CLIN:priorMDS_TRUE", rownames(fm)),,drop=F]), t(fm[grepl("in_PB|in_BM|B:CLIN:is|finalFusion_Complex", rownames(fm)),,drop=F]))
DATA_clin$Age[is.na(DATA_clin$Age)]=median(DATA_clin$Age, na.rm = T) #one observation --> set to median

colnames(DATA_clin)=gsub("..CLIN.", "", colnames(DATA_clin))
beatAML_clin=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_clin,TIME,STATUS, univariate = T, pretty=F))

# stroma
ME=microenv[microenv$disease=="AML",]
DATA=data.frame(t(gexp[rownames(gexp)%in%ME$gene,]))
beatAML_microenv=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
beatAML_microenv[beatAML_microenv$P<0.01,]

# MDS
DATA=data.frame(t(gexp[rownames(gexp)%in%MDS$MDS_signature_all_filt,]))
beatAML_MDS=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
beatAML_MDS[beatAML_MDS$P<0.05,]

# costim
DATA=data.frame(t(gexp[rownames(gexp)%in%co.stim[,1],]))
beatAML_costim=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA,TIME,STATUS, univariate = T, pretty=F))
beatAML_costim[beatAML_costim$P<0.05,]

# subtype:
load("BeatAML_subtypes.Rdata")
coordinates.subtype$subtype[coordinates.subtype$subtype%in%"Progenitor-like"]="CMP-like"
coordinates.subtype=coordinates.subtype[match(colnames(fm), coordinates.subtype$ID),]

DATA_subtypes=data.frame(do.call(cbind, get.logical(list(coordinates.subtype$subtype)))*1)

beataml_res_subtype=do.call(rbind, lapply(seq(logicalv), fun.get.cox, logicalv, DATA_subtypes,TIME,STATUS, univariate = T, pretty=F))

fun.kapplanMeier(TIME[logicalv[[1]]], STATUS[logicalv[[1]]],GROUPS=coordinates.subtype$subtype[logicalv[[1]]], MONTHS=T, PVAL=1, INDIVIDUAL_GROUPS=F,LWD = 1, NAME = "Prognostic Index - validation")

# make table S6, adjusted p-value set here to correct for number of comparisons in total:
tableS7=rbind(beatAML, beataml_res_subtype, beatAML_costim, beatAML_MDS, beatAML_microenv)
tableS7=tableS7[tableS7$Name%in%"BeatAML",]
tableS7$Adj.P=p.adjust(tableS7$P, method="BH")

tableS7[,2]=prettyNum(tableS7[,2])
tableS7[,3]=prettyNum(tableS7[,3])
tableS7[,4]=prettyNum(tableS7[,4])
tableS7[,5]=prettyNum(tableS7[,5])
tableS7[,6]=prettyNum(tableS7[,6])
tableS7[,8]=prettyNum(tableS7[,8])

# annotate these genes, needed later:
genelist_signif=data.frame(tableS7[,1], type="", stringsAsFactors = F)
genelist_signif$type[genelist_signif[,1]%in%c("HLAI", "HLAII", "CytolyticScore", "Freq.CGA")]="ImmunoScores"
genelist_signif$type[genelist_signif[,1]%in%beataml_res_subtype$Feature]="Subtype"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho > 0)",1]]="Stromal/cancer gene (Rho > 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"Stromal/cancer gene (Rho < 0)",1]]="Stromal/cancer gene (Rho < 0)"
genelist_signif$type[genelist_signif[,1]%in%ME[ME$category%in%"CTL/NK gene",1]]="CTL/NK gene"
genelist_signif$type[genelist_signif[,1]%in%MDS$MDS_signature_all_filt]="MDS-signature gene"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Inhibitory", co.stim$`Immune checkpoint function`),1]]="Inhibitory ligand"
genelist_signif$type[genelist_signif[,1]%in%co.stim[grepl("Stimulatory", co.stim$`Immune checkpoint function`),1]]="Stimulatory ligand"

tableS7$Type=genelist_signif$type
tableS7=tableS7[order(tableS7$Type),]
tableS7[,1]=gsub("\\.", "-", tableS7[,1])

data.table::fwrite(tableS7[,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_beatAML.tsv", sep="\t")
data.table::fwrite(tableS7[tableS7$Adj.P<0.2,c(1,11,10,2,3,4,5,6,8,9)], "tableS7_beatAML_signif.tsv", sep="\t")

# save cohort and survival
gexp=data.frame(scale(t(gexp)))
immunoscore=data.frame(scale(t(fm[grepl("N:SAMP:.*.Score", rownames(fm)),])), check.names = F)
colnames(immunoscore)=c("HLAI", "HLAII", "CytolyticScore")

samp=cbind(DATA_clin, DATA_subtypes)

save(list = c("gexp","immunoscore", "samp", "TIME", "STATUS", "logicalv"), file="BeatAML_survival_data.Rdata")

#******************************** make excel table
files=list.files(pattern = "tableS7")
files=files[!grepl("_CHOP", files)]

univariate.results=lapply(files[grepl("signif", files)], fread, data.table=F)
names(univariate.results)=gsub("tableS7_|_signif.tsv", "", files[grepl("signif", files)])
univariate.results=univariate.results[c(grep("DLBCL", names(univariate.results)), grep("MM", names(univariate.results)), grep("AML", names(univariate.results)))]

require(openxlsx)

write.xlsx(univariate.results, file = "TableS7_cox.xlsx", )
