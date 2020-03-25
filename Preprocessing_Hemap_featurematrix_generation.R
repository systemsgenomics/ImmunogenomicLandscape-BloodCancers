#***************************************************************************************************
#******************************* Make immunology FM **********************************************
#***************************************************************************************************
library(mclust)
library(data.table)
library(parallel)

source("/research/users/ppolonen/git_home/common_scripts/featurematrix/functions_generate_fm.R")

# WD
setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/")

# GEXP
gexp=t(get(load("data9544_with_gene_symbols.RData")))

# annotation table
check=get(load("Hemap_immunology_Annotations_8304.Rdata"))
annot=read.delim("anno_coord_data9544_15pct_bw2.5_updated.txt", stringsAsFactors=F, header=T)

# listing files, not needed anymore
# f=list.files("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/Annotations_immunology/", pattern = ".csv", full.names = T)
# f=f[-1]
# annot_normal=do.call(rbind, lapply(f, read.csv, header=F, stringsAsFactors=F, skip=1))

#******************************** filtering *************************************
exvivotreatments_allowed <- c("none", "na", "control", "activation", "differentiation", "differentiation followed by activation", " differentiation followed by LPS+IFNg", "differentiation followed by IFN", "differentiation followed by inflammatory cytokines", "differentiation followed by LPS", "differentiation followed by CD40L", "differentiation followed by Poly(I:C)", "differentiation (IL-4)", "differentiation with EPO", "anti-IgM", "IL-2", "IL-3", "IL-3+CpG")

rm_healthy = annot$Sample.type%in%c("NonCancerHealthy")&!(annot$In.vivo.treatment%in%c("na", "none", "no")&annot$Ex.vivo.treatment%in%exvivotreatments_allowed)
rm_cancer_prolif = annot$Sample.type%in%c("Cancer", "Prolif")&!annot$Ex.vivo.treatment=="none"
rm_celline = annot$Sample.type%in%c("CellLine")&!annot$Ex.vivo.treatment%in%c("none", "na", "control")
rm_treatments = annot$GSE.identifier..experiment.%in%c("GSE26661")|annot$GSM.identifier..sample.%in%c("GSM425497", "GSM425499")
rm_noncancer = !annot$Sample.type%in%c("Cancer", "Prolif", "CellLine", "NonCancerHealthy")
rm_mm_pbmc = annot$colorClass=="MM"&!annot$Sample.isolation=="CD138+ plasma cells"
# cluster exclusion
# AML,GSE7538,24
# CLL,GSE18866,GSE9250,32
# LP,GSE12453,GSE7345,6
# MM,GSE24147,GSE24522,18
# MP,GSE12079,GSE15811,13
# TCL,GSE14879,16 microdissected, do not remove

rm_cluster_differ=annot$GSE.identifier..experiment.%in%c("GSE7538", "GSE9250", "GSE12453", "GSE18866","GSE24522", "GSE7345", "GSE12079", "GSE15811", "GSE24147")&annot$Sample.type=="Cancer"

annot_left_out = annot

annot_left_out$reason_removed[rm_healthy]="Normal Cell sample, ex-vivo or in-vivo treated"
annot_left_out$reason_removed[rm_cluster_differ]="Cancer sample, outlier cluster"
annot_left_out$reason_removed[rm_celline]="Cell line sample, ex-vivo treated"
annot_left_out$reason_removed[rm_cancer_prolif|rm_treatments]="Cancer sample, ex-vivo treated"
annot_left_out$reason_removed[rm_noncancer]="NonCancer sample, not healthy"
annot_left_out$reason_removed[rm_mm_pbmc]="NonCancer sample, not healthy"

annot_left_out <- annot_left_out[(rm_noncancer|rm_healthy|rm_cancer_prolif|rm_celline|rm_treatments|rm_cluster_differ|rm_mm_pbmc),]
fwrite(annot_left_out, "hemap_1072_leftout_reasonremoved.txt", sep="\t")

annot <- annot[!(rm_healthy|rm_cancer_prolif|rm_celline|rm_treatments|rm_cluster_differ|rm_noncancer|rm_mm_pbmc),]
gexp <- gexp[,!(rm_healthy|rm_cancer_prolif|rm_celline|rm_treatments|rm_cluster_differ|rm_noncancer|rm_mm_pbmc)]

newSamples=annot[!annot[,1]%in%check[,1],]
fwrite(newSamples, "new_samples_in_hemap_8472.txt", sep="\t")

dim(check[check[,1]%in%annot[,1],])


#****************************************************************************************
# annot[rm_healthy&annot$GSM.identifier..sample.%in%annot_normal[,1],] # excluded from previous normals

# new
annot$CLASS[annot$Sample.type%in%c("NonCancerHealthy")]=gsub("NonCancerHealthy_|NonCancerHealthy|na_","", gsub("_StemCell_|Myeloid_|_na_na|1_na_na|2_na_na|Lymphoid_|_G.*.|_na_G0|_na_G1.*.|T-|B-|M1-Differentiating|Erythroid_","", annot$CLASS[annot$Sample.type%in%c("NonCancerHealthy")]))
annot$CLASS[annot$Sample.type%in%c("NonCancerHealthy")][annot$CLASS[annot$Sample.type%in%c("NonCancerHealthy")]=="na"]="LymphNode_GerminalCentre"
annot$CLASS[annot$CLASS=="CD4+Tcell"]="RestingCD4+Tcell"
annot$CLASS[annot$CLASS=="CD8+Tcell"]="RestingCD8+Tcell"

annot$MAINCLASS[annot$Sample.type%in%c("NonCancerHealthy")]="NonCancerHealthy"
annot$colorClass[annot$Sample.type%in%c("NonCancerHealthy")][annot$CLASS[annot$Sample.type%in%c("NonCancerHealthy")]=="na"]="Lymphoid"

# old
# annot$CLASS[match(annot_normal[,1], annot$GSM.identifier..sample.)]=annot_normal[,5]
# annot$MAINCLASS[match(annot_normal[,1], annot$GSM.identifier..sample.)]=annot_normal[,3]
# annot$colorClass[match(annot_normal[,1], annot$GSM.identifier..sample.)]=annot_normal[,4]

#*********************************** compute geometric mean ********************************
# get certain genes GEXP
rownames(gexp)=paste("N:GEXP:", rownames(gexp), sep="")

dat_a=gexp[grep("GZMA|PRF1|GNLY|GZMH|GZMM", rownames(gexp)),]
# dat_a=matrix[grep("GZMA|PRF1", rownames(matrix)),]
dat=2^dat_a+0.01
rownames(dat)
gm1=log2(t(apply(dat, 2, gm_mean)))
rownames(gm1)="CytolyticScore"

# also HLA
dat_a2=gexp[rownames(gexp)%in%c("N:GEXP:B2M", "N:GEXP:HLA-A", "N:GEXP:HLA-B", "N:GEXP:HLA-C"),]

dat2=2^dat_a2+0.01
rownames(dat2)
gm2=log2(t(apply(dat2, 2, gm_mean)))
rownames(gm2)="HLAIScore"

# also HLAII
dat_a3=gexp[rownames(gexp)%in%c("N:GEXP:HLA-DMA",
                "N:GEXP:HLA-DMB",
                "N:GEXP:HLA-DPA1",
                "N:GEXP:HLA-DPB1",
                "N:GEXP:HLA-DRA",
                "N:GEXP:HLA-DRB1"),]

dat3=2^dat_a3+0.01
rownames(dat3)
gm3=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm3)="HLAIIScore"

classification1=data.frame(t(rep("medium", length(gm1))), stringsAsFactors = F)
zscore=as.numeric(scale(t(gm1)))
classification1[zscore>=1]="high"
classification1[zscore<=(-1)]="low"
rownames(classification1)="CytolyticScore" 
colnames(classification1)=colnames(gexp)

classification2=data.frame(t(rep("medium", length(gm2))), stringsAsFactors = F)
zscore=as.numeric(scale(t(gm2)))
classification2[zscore>=1]="high"
classification2[zscore<=(-1)]="low"
rownames(classification2)="HLAIScore" 
colnames(classification2)=colnames(gexp)

classification3=data.frame(t(rep("medium", length(gm3))), stringsAsFactors = F)
zscore=as.numeric(scale(t(gm3)))
classification3[zscore>=1]="high"
classification3[zscore<=(-1)]="low"
rownames(classification3)="HLAIIScore" 
colnames(classification3)=colnames(gexp)

classification=data.frame(t(rbind(classification1,classification2,classification3)), stringsAsFactors = F)


immunoscores=as.data.frame(t(rbind(gm1, gm2, gm3)))
immunoscoresfm=make.features(immunoscores, datatype="SAMP", prefix="")
colnames(immunoscoresfm)=colnames(gexp)

immunoscores_class_fm=make.features(classification, datatype="SAMP", prefix="")
colnames(immunoscores_class_fm)=colnames(gexp)

# excluding categorical here, they slow everything down!
l.data_list=list(gexp, immunoscoresfm, immunoscores_class_fm)
data_list=data.frame(do.call(rbind, l.data_list))

# ******************************** Infer cell fractions ***********************************

# cibersort
results=read.delim("CIBERSORT-Results.txt", row.names = 1, header=T, stringsAsFactors = F)
colnames(results)=paste0("N:SAMP:CIBERSORT_", gsub(" |-|\\.", "_", colnames(results)), "")
cibersort=t(results)
cibersort=cibersort[,colnames(cibersort)%in%colnames(gexp)]

MCP=get(load("MCP_counter_data.Rdata"))
rownames(MCP)=paste0("N:SAMP:", rownames(MCP), "")
MCP=MCP[,colnames(MCP)%in%colnames(cibersort)]

l.fractions=list(cibersort, MCP)
fractions=data.frame(do.call(rbind, l.fractions), stringsAsFactors = F)

#********************************** Clinical data ******************************

files=list.files(".", pattern=".info.tsv", full.names = T)

surv_data=do.call(rbind, lapply(files, read.delim, header=T, stringsAsFactors=F))
surv_data=surv_data[surv_data[,1]%in%colnames(gexp),]
surv_data=surv_data[!is.na(surv_data[,3]),]

surv_d=t(surv_data[,2:3])
colnames(surv_d)=surv_data[,1]
rownames(surv_d)=c("N:CLIN:OS_Time", "B:CLIN:OS_Status")

surv_d=surv_d[,match(colnames(gexp), colnames(surv_d))]
colnames(surv_d)=colnames(gexp)

#******************************* tumor percentage ********************************
Sys.setlocale(locale="C")
sorteds=read.delim("sorted_samples.txt", stringsAsFactors=F, header=F)

gsm=unlist(lapply(1:dim(sorteds)[1], function(i){
  annot$GSM.identifier..sample.[annot$Sample.isolation%in%sorteds[i,1]&annot$colorClass%in%sorteds[i,2]]
}))
add=annot$GSM.identifier..sample.[grepl("CD303", annot$Sample.isolation)]

sorted=t(annot$GSM.identifier..sample.%in%c(gsm, add))

rownames(sorted)="B:CLIN:CELLS_SORTED"
colnames(sorted)=colnames(gexp)

sorted[,grepl("Padiatr", annot$Additional.notes)]=1
# sorted[,annot$Additional.notes=="The leukemic blasts were sorted based on CD41, CD7, CD117, CD33, and CD34 antibodies as previously described (Klin. Padiatr. 217, 126-134)."] = 1

tumor_per=gsub("blast%: |>=|%|blast cell percentage: |t_cell_purity: |>|;|blast count, % of sample, -1=unavailable : ","", annot$Purity.Tumor.cell.content)
tumor_per[tumor_per=="high"]=80
tumor_per[as.numeric(tumor_per)>100]=100
tumor_per[tumor_per=="-1"]=0
tumor_per[tumor_per%in%c("n./a.", "na")]=NA

malt=read.delim("clinical_annotations_MALT.txt", stringsAsFactors=F, header=T)
replace=malt[match(colnames(gexp), malt$GSMID),]
tumor_per[grepl("Percentage of tumor", tumor_per)]=replace$X..Tumor[grepl("Percentage of tumor", tumor_per)]

tumor_percentage=data.matrix(t(as.numeric(tumor_per)))
rownames(tumor_percentage)="N:SAMP:BLAST_TUMOR_PERCENTAGE"
colnames(tumor_percentage)=colnames(gexp)

# T-cell percentages
T_per=data.matrix(t(as.numeric(replace$X..T.cells)))
rownames(T_per)="N:SAMP:TCELL_PERCENTAGE"
colnames(T_per)=colnames(gexp)

tissue_per=data.matrix(t(as.numeric(replace$X..Lung)))
rownames(tissue_per)="N:SAMP:TISSUE_PERCENTAGE"
colnames(tissue_per)=colnames(gexp)


#*************************** adding some lymphoma annotation *************************************
anno = read.delim("GSE10846_series_matrix_info_ipi_clean.txt", stringsAsFactors=F, header=T)
anno=anno[anno[,1]%in%annot[,1],]

anno=anno[match(annot[,1], anno[,1]),]

annot$In.vivo.treatment[annot[,1]%in%anno[,1]]=gsub("*.*: |;", "", anno$chemotherapy[annot[,1]%in%anno[,1]])
annot$In.vivo.treatment[annot$In.vivo.treatment%in%"NA"]=NA

annot$dlbcl_ipi=NA
annot$dlbcl_ipi[annot[,1]%in%anno[,1]]=anno$ipi[annot[,1]%in%anno[,1]]

CHOP=t(annot$In.vivo.treatment%in%"CHOP-Like Regimen"*1)
RCHOP=t(annot$In.vivo.treatment%in%"R-CHOP-Like Regimen"*1)
CHOP[is.na(annot$In.vivo.treatment)]=NA
RCHOP[is.na(annot$In.vivo.treatment)]=NA
rownames(CHOP)="B:CLIN:Chemotherapy_CHOP"
rownames(RCHOP)="B:CLIN:Chemotherapy_RCHOP"
colnames(CHOP)=colnames(gexp)
colnames(RCHOP)=colnames(gexp)

# add cytogenetic information
genetics_org=read.delim("AML_preBALL_cytogenetics_vectors.txt", stringsAsFactors=F, header=T, row.names=1)
genetics=genetics_org[rownames(genetics_org)%in%colnames(gexp),]
genetics=t(genetics)*1
rownames(genetics)=paste("B:CLIN:", "GENETICS_", rownames(genetics), "", sep="")

cytogenetic=annot$Cytogenetics
cytogenetic[cytogenetic%in%c("na", "n/a", "", " ")]=NA
cytogenetic[grepl("without|unknown|remainingcytogenetics|no del13q|crlf2 fish: Normal|crlf2 fish: n/a", cytogenetic)]=NA
cytogenetic[grepl("ormal", cytogenetic)]="normal_karyotype"
cytogenetic[grepl("MLL", cytogenetic)]="MLL"
cytogenetic=gsub("complex aberrant karyotype", "complex karyotype", cytogenetic)
cytogenetic=gsub("hyperdiploid karyotype", "hyperdiploid", cytogenetic)
cytogenetic=gsub("TAL$", "TAL1", cytogenetic)
cytogenetic=gsub("remaining cytogenetics|other abNormalities", "other", cytogenetic)
cytogenetic=gsub("fish:|trisomy 8 |;deletion|; complex karyotype|, complex karyotype|: positive| chromosomal aberrations|/API2-MALT1|/API2-MALT1 negative|deletion *.*: negative|/IGH-MALT1|, plus other| plus other", "", cytogenetic)
cytogenetic=gsub("trisomy ", "trisomy", cytogenetic)
cytogenetic=gsub("TEL deleted", "TEL_deleted", cytogenetic)
cytogenetic=gsub("p13.1", "p13", cytogenetic)
cytogenetic=gsub(";$", "", cytogenetic)
cytogenetic=gsub("\\+ ", "+", cytogenetic)
cytogenetic=gsub("i\\(", "inv(", cytogenetic)
cytogenetic=gsub("complex karyotype", "complex_karyotype", cytogenetic)
cytogenetic=gsub("^ ", "", cytogenetic)
annot$cytogenetic_clean=cytogenetic

cytogenetic_terms=sort(unique(unlist(strsplit(cytogenetic, " "))))

cytogenetic_terms=unlist(strsplit(cytogenetic, "; |; |  |, | |/"))
cytogenetic_terms=gsub(" ", "", cytogenetic_terms)

terms=table(cytogenetic_terms)
shared=names(terms)[terms>5]

cytogenetics=do.call(rbind, mclapply(shared,FIND_LOGICAL, cytogenetic, mc.cores=8))
colnames(cytogenetics)=colnames(gexp)
#***********************

# age annotations
age=gsub("*.*: |^ |;| yr| age| years|-.*.$| Years|d 32.8|d 54.8", "", annot$Age)
age[age%in%c("na", "n/a", "", "not available")]=NA
age[grepl("month|Month", age)]=signif(as.numeric(gsub(" months.*.| months| Months", "", age[grepl("month|Month", age)]))/12, 2)
age[age=="Adult"]=30
age[age=="Children"]=5
age[age=="pediatric"]=1
age=t(as.numeric(age))
age[age>100]=NA
rownames(age)="N:CLIN:AGE"
colnames(age)=colnames(gexp)

# gender annotations
gender=toupper(annot$Gender)
gender[gender%in%c("GENDER: NOT AVAILABLE", "GENDER: NA;", "SEX: UNKNOWN;")]=NA
gender=gsub(" |;", "", gender)
gender=gsub("GENDER:|SEX/AGE:|/.*.|SEX:", "", gender)
gender[gender%in%c("F", "FEMALE", "WOMAN")]="female"
gender[gender%in%c("M", "MALE", "MAN")]="male"
gender=t(gender)
rownames(gender)="C:CLIN:GENDER"
colnames(gender)=colnames(gexp)

# race annotations
race=toupper(annot$Race)
race[grep("AGE", race)]=NA
race[grep("OTHER", race)]="OTHER"
race[grep("AFRICAN|RACE: AA;|RACE: B;", race)]="AFRICAN"
race[grep("HISPANIC|RACE: H; ", race)]="HISPANIC"
race[grepl("EUROPEAN|CAUCASIAN|ANGLO-AMERICAN|WHITE|RACE: W;|RACE: C;", race)]="EUROPEAN"
race[grep("ASIAN", race)]="ASIAN"
race[!grepl("ASIAN|AFRICAN|HISPANIC|EUROPEAN|OTHER", race)]=NA

race=t(race)
rownames(race)="C:CLIN:RACE"
colnames(race)=colnames(gexp)

#*************************** adding some myeloma annotation *************************************
annomm = read.delim2("GSE24080_MM_clininfo_GSMid_clean.txt", stringsAsFactors=F, header=T)
annomm2 = read.delim2("GSE19784_MM_survival_GSMid_iss_clean.txt", stringsAsFactors=F, header=T)

annomm$b2m=gsub("<0.5", "0.5", annomm$b2m)
annomm$b2m=as.numeric(gsub(",", ".", annomm$b2m))

annomm$aspc=as.numeric(gsub(",", ".", annomm$aspc))
annomm$bmpc=as.numeric(gsub(",", ".", annomm$bmpc))

# first combine the two:
library(data.table)
combmm=data.frame(rbindlist(list(annomm, annomm2), fill = TRUE), stringsAsFactors = F)

combmm=combmm[match(colnames(gexp), combmm$accession),]

# now add these vectors to annot table
age[colnames(gexp)%in%combmm$accession]=signif(combmm$age[colnames(gexp)%in%combmm$accession], 3)
gender[colnames(gexp)%in%combmm$accession]=combmm$sex[colnames(gexp)%in%combmm$accession]

combmm$race[combmm$race%in%"other"]="OTHER"
combmm$race[combmm$race%in%"white"]="EUROPEAN"
race[colnames(gexp)%in%combmm$accession]=combmm$race[colnames(gexp)%in%combmm$accession]

# time and status:
surv_d[1,colnames(gexp)%in%combmm$accession]=signif(combmm$os_time[colnames(gexp)%in%combmm$accession], 3)
surv_d[2,colnames(gexp)%in%combmm$accession]=signif(combmm$os_censor[colnames(gexp)%in%combmm$accession], 3)

# pfs time and status
pfs=cbind(combmm$pfs_time, combmm$pfs_censor)
pfs=t(pfs)
colnames(pfs)=colnames(gexp)
rownames(pfs)=c("N:CLIN:PFS_Time", "B:CLIN:PFS_Status")

# other myeloma annotations:
otherMM_C=t(combmm[,c(4,8)])
rownames(otherMM_C)=paste0("C:CLIN:MM_", toupper(rownames(otherMM_C)), "")
colnames(otherMM_C)=colnames(gexp)

# numeric myeloma:
otherMM_N=t(combmm[,c(9:20)])
rownames(otherMM_N)=paste0("N:CLIN:MM_", toupper(rownames(otherMM_N)), "")
colnames(otherMM_N)=colnames(gexp)

add_mm=t(combmm$cytogenetic_abnormalities)
colnames(add_mm)=colnames(gexp)
rownames(add_mm)=c("B:CLIN:MM_CYTOGENETIC_ABNORMALITIES")

l.clin=list(sorted, surv_d, pfs, tumor_percentage,T_per,tissue_per, CHOP, RCHOP, gender, age,race, cytogenetics, genetics, add_mm, otherMM_N, otherMM_C)
clin=data.frame(do.call(rbind, l.clin), stringsAsFactors = F)

#****************************** make annotation clusters ******************************************
annot$acute=rep("other", dim(annot)[1])
annot$acute[annot$colorClass=="AML"|annot$colorClass=="pre-B-ALL"|annot$colorClass=="T-ALL"]="acute_leukemias"
annot$acute[annot$colorClass=="CLL"|annot$colorClass=="CML"]="chronic_leukemias"
annot$acute[grepl("NonCancer", annot$MAINCLASS)]="NonCancer"
annot$disease=rep("other", dim(annot)[1])
annot$disease[grepl("NonCancer", annot$Sub.maps.available)]="NonCancer"
annot$CLASS2=annot$CLASS
annot$CLASS2[grepl("CellLine_Myeloma", annot$CLASS2)]="CellLine_Myeloma"
annot$CLASS2[grepl("CellLine_Leukemia", annot$CLASS2)]="CellLine_Leukemia"
annot$CLASS2[grepl("CellLine_Lymphoma", annot$CLASS2)]="CellLine_Lymphoma"
annot$CLASS2[grepl("CellLine_mix", annot$CLASS2)]="CellLine_mix"

findthese=c("NonCancer", "Cancer_Leukemia", "Cancer_Myeloma", "Cancer_Lymphoma","CellLine_Leukemia","CellLine_Lymphoma","CellLine_Myeloma","Prolif_Lymphoproliferative_ALPS","Prolif_Lymphoproliferative_MPN", "Prolif_Myeloproliferative_LCH_LC", "Prolif_Myeloproliferative_MDS")

for(f in findthese){
  annot$disease[grepl(f, annot$MAINCLASS)]=f
}

annot$subclasses=rep("other", dim(annot)[1])
findthese=c("NonCancer", "Cancer_Myeloma", "Cancer_Lymphoma", "CellLine_Leukemia", "CellLine_Lymphoma","CellLine_Myeloma","Prolif_Lymphoproliferative_ALPS","Prolif_Lymphoproliferative_MPN", "Prolif_Myeloproliferative_LCH_LC", "Prolif_Myeloproliferative_MDS")
findthese2=c("T-ALL", "pre-B-ALL", "AML","CML","CLL", "BCL", "TCL", "B-Lymphoid", "T-Lymphoid","Lymphoid", "Myeloid", "Erythroid", "StemCell")

for(f in findthese){
  annot$subclasses[grepl(f, annot$MAINCLASS)]=f
}
for(f in findthese2){
  annot$subclasses[grepl(f, annot$colorClass)]=f
}
DLBCL=c("Cancer_Lymphoma_BCL_DLBCL_ABC", "Cancer_Lymphoma_BCL_DLBCL_GCB", "Cancer_Lymphoma_BCL_DLBCL_na")
annot$subclasses[annot$CLASS2%in%DLBCL]="BCL_DLBCL"


annot$CLASS=gsub("_na|_check|_testicular", "",annot$CLASS)

# lymphoma annotations
bLY=(1:nrow(annot)%in%grep("Lymphoma_BCL",annot$CLASS))&(!1:nrow(annot)%in%grep("CellLine",annot$CLASS))&(!1:nrow(annot)%in%grep("NonCancer",annot$CLASS))
tLY=(1:nrow(annot)%in%grep("Lymphoma_TCL",annot$CLASS))&(!1:nrow(annot)%in%grep("CellLine",annot$CLASS))&(!1:nrow(annot)%in%grep("NonCancer",annot$CLASS))
annot$CLASS=gsub("Cancer_", "", annot$CLASS)

annot$tbLY=NA
annot$tbLY[bLY|tLY]=annot$CLASS[bLY|tLY]

# these are the terms to look for
# table(annot$disease)
# table(annot$colorClass)
# table(annot$acute)
# table(annot$subclasses)
# table(annot$tbLY)

#*******************************************************************************************
#*************************** annotation clusters ********************************************

#****************************** make clusters ******************************************

# make immunological normal annotation vectors
annot$plotNormals = "Other"

lv=!annot$Sample.type%in%"NonCancerHealthy"
annot$plotNormals[lv]=""

lv=grepl("RestingBcell|NaiveBcell|MemoryBcell|BcellActivated", annot$CLASS)
annot$plotNormals[lv]="B cell"

lv=grepl("GerminalCentre", annot$CLASS)
annot$plotNormals[lv]="Germinal centre cell"

lv=grepl("PlasmaBcell", annot$CLASS)
annot$plotNormals[lv]="Plasma cell"

lv=grepl("Tcell|NaturalKillerCell", annot$CLASS)
annot$plotNormals[lv]="T/NK cell"

lv=grepl("DendriticCell", annot$CLASS)
annot$plotNormals[lv]="Dendritic cell"

lv=grepl("Langerhans", annot$CLASS)
annot$plotNormals[lv]="Langerhans cell"

lv=grepl("Eryth|Platelet", annot$CLASS)
annot$plotNormals[lv]="Erythroid"

lv=grepl("Monocyte", annot$CLASS)
annot$plotNormals[lv]="Monocyte"

lv=grepl("Macrophage", annot$CLASS)
annot$plotNormals[lv]="Macrophage"

lv=grepl("Neutrophil", annot$CLASS)
annot$plotNormals[lv]="Neutrophil"

lv=grepl("MyeloidProgenitor", annot$CLASS)
annot$plotNormals[lv]="Myeloid progenitor"

lv=grepl("HematopoieticStemCell", annot$CLASS)
annot$plotNormals[lv]="HSC"

lv=grepl("^Mononuclear", annot$CLASS)
annot$plotNormals[lv]="PBMC"

lv=grepl("LymphNode", annot$CLASS)
annot$plotNormals[lv]="Lymph node"

HLAplot_normals <- c("B cell", "Plasma cell", "T/NK cell", "Dendritic cell", "Erythroid", "Monocyte", "Macrophage", "Neutrophil", "Myeloid progenitor", "HSC")
cytolyticplot_normals <- c("PBMC", "Lymph node")
costimplot_normals <- c(HLAplot_normals, "PBMC", "Lymph node", "Langerhans cell", "Germinal centre cell")

annot$immunoNormals=annot$Category.specifying.lineage.tumor.origin

lv=grepl("CD8|CD8+TcellActivated", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="CD8+Tcell"

lv=grepl("NaturalKiller", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="NKCell"

lv=grepl("M2-Macrophage", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="M2-Macrophage"

lv=grepl("M1-Macrophage", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="M1-Macrophage"

lv=grepl("DendriticCell", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="DendriticCell"

lv=grepl("Monocyte", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="Monocyte"

lv=grepl("CD4+", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="CD4+Tcell"

lv=!annot$Sample.type%in%"NonCancerHealthy"
annot$immunoNormals[lv]=""

lv=grepl("Eryth|Platelet", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="Erythroid"

lv=grepl("CD3", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="Tcell"

lv=grepl("GerminalCentre", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="GerminalCentreCell"

lv=grepl("^Tcell$|^TcellActivated$|^TcellResting$", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="Tcell"

lv=grepl("^ActivatedBcell$|^RestingBcell$|^BcellActivated$", annot$Category.specifying.lineage.tumor.origin)
annot$immunoNormals[lv]="Bcell"

# annotated clusters
tbLY=FUN_MAKE_ALL(annot$tbLY, "annotated_class", annot$tbLY, 0)
subclasses=FUN_MAKE_ALL(annot$subclasses, "annotated_class", annot$subclasses, 0)
acute_chronic=FUN_MAKE_ALL(annot$acute, "annotated_class", annot$acute, 0)
colorClass=FUN_MAKE_ALL(annot$colorClass, "annotated_class", annot$colorClass, 0)
disease=FUN_MAKE_ALL(annot$disease, "annotated_class", annot$disease, 0)
tbLY=FUN_MAKE_ALL(annot$tbLY, "annotated_class", annot$tbLY, 0)
fullclass=FUN_MAKE_ALL(annot$CLASS2, "annotated_class", annot$CLASS2, 0)
immunoclass=FUN_MAKE_ALL(annot$immunoNormals, "annotated_class_immunoNormals", annot$immunoNormals, 0)

l.comparisons=list(disease, colorClass, acute_chronic,subclasses, tbLY, fullclass, immunoclass)
comparisons=do.call(rbind, l.comparisons)
comparisons=data.frame(data.matrix(comparisons[!duplicated(rownames(comparisons)),]), stringsAsFactors = F)
colnames(comparisons)=colnames(gexp)

# test if all rows are fine, should be >1 values
A=apply(comparisons, 1, unique)

B=unlist(lapply(A, function(d)sum(d%in%c(1,0))>=2))

if(!all(B))stop("Check comparisons, impossible comparisons made")

# categorical feats
class1=FUN_MAKE_CATEGORICAL(annot$tbLY, "annotated_class_BCL_TCL")
class2=FUN_MAKE_CATEGORICAL(annot$colorClass, "annotated_class_colorclass")
class3=FUN_MAKE_CATEGORICAL(annot$disease, "annotated_class_disease")
class4=FUN_MAKE_CATEGORICAL(annot$immunoNormals, "annotated_class_immunoNormals")

l.comparisons=list(class1, class2, class3, class4)
comparisons_cat=do.call(rbind, l.comparisons)
comparisons_cat=data.frame(comparisons_cat[!duplicated(rownames(comparisons_cat)),], stringsAsFactors = F)

colnames(comparisons_cat)=colnames(gexp)


#****************************************************************************************************
#******************** Next we start to create features of these individual data types ***************
#****************************************************************************************************

l.fm=list(data_list,clin,fractions, comparisons)

library(data.table)
fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

matrix=fm

# also add clinicaldata to annotations
numclin=t(data.matrix(clin[!grepl("^C:", rownames(clin)),]))
numchr=t(clin[grepl("^C:", rownames(clin)),])
colnames(numclin)=gsub(".:CLIN:|.:GEXP:|", "", colnames(numclin))
colnames(numchr)=gsub(".:CLIN:|.:GEXP:|", "", colnames(numchr))

annot_add=data.frame(numclin, numchr, stringsAsFactors = F)

fractions2=t(fractions)
colnames(fractions2)=gsub("N:SAMP:", "", colnames(fractions2))
annot2=data.frame(annot, annot_add, "CytolyticScore"=as.numeric(gm1) ,"HLAIScore"=as.numeric(gm2), "HLAIIScore"=as.numeric(gm3),classification,fractions2, stringsAsFactors = F)


# annot
clusters=read.delim("AML_15pct_BHSNE_mean-shift.txt", stringsAsFactors=F, header=T)
clusters=clusters[clusters$ID%in%annot$GSM.identifier..sample.,]
clusters_cancermap=clusters$X1.5..cluster

matrix_sub=matrix[,colnames(matrix)%in%clusters$ID]
annot_sub=annot[annot$GSM.identifier..sample.%in%clusters$ID,]

# TCGA clusters
cluster_mapping=read.delim("Table_TCGA_cluster_AML_cluster_assignment.txt", header=T, stringsAsFactors=F, sep="\t")
TCGA_cluster=rep("NA", dim(annot_sub)[1])

TCGA_cluster[clusters_cancermap%in%cluster_mapping[1,2]]="TCGA_AML_cluster_1"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[2,2]]="TCGA_AML_cluster_2"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[3:5,2]]="TCGA_AML_cluster_3"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[6:10,2]]="TCGA_AML_cluster_4"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[11,2]]="TCGA_AML_cluster_5"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[12:15,2]]="TCGA_AML_cluster_6"
TCGA_cluster[clusters_cancermap%in%cluster_mapping[16:17,2]]="TCGA_AML_cluster_7"

n=lapply(unique(TCGA_cluster), function(i)annot_sub$GSM.identifier..sample.[TCGA_cluster%in%i])
names(n)=unique(TCGA_cluster)
n=n[!unique(TCGA_cluster)%in%"NA"]
save(n, file="Hemap_immunology_TCGA_clusters.Rdata")

save(matrix, file="Hemap_immunology_fm.Rdata")
save(annot2, file="Hemap_immunology_Annotations.Rdata")
write.table(annot2,"Hemap_immunology_Annotations.tsv", sep="\t", col.names=T, row.names=F, quote=FALSE)

write.table(t(c("N:SAMP", as.character(colnames(matrix)))), file="Hemap_immunology_fm.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(matrix, file="Hemap_immunology_fm.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)

# make a small fix here to harmonize survival data to months:
load("Hemap_immunology_Annotations.Rdata")

unique(cbind(annot2[!is.na(annot2$OS_Time), c(2,4)]))
annot2$OS_Time[!is.na(annot2$OS_Time)&annot2[,2]%in%c("GSE10846,GSE11318", "GSE10846", "GSE10846,GSE17372", "GSE11877")]=annot2$OS_Time[!is.na(annot2$OS_Time)&annot2[,2]%in%c("GSE10846,GSE11318", "GSE10846", "GSE10846,GSE17372", "GSE11877")]*12

# for myeloma, transform data to 5y survival to compare data sets:
modify=!is.na(annot2$OS_Time)&annot2[,2]%in%c("GSE16716,GSE24080")
find=annot2$OS_Time>60&modify

# change status to alive if dead later
find2=annot2$OS_Status==1&modify
annot2$OS_Status[find&find2]=0
annot2$OS_Time[find]=60

save(annot2, file="Hemap_immunology_Annotations.Rdata")

#****************************************************************************************
# This FM can then be used as a backbone for other FMs. GSVA and clusters must be added
#****************************************************************************************

matrix=get(load("Hemap_immunology_fm.Rdata"))
annot=get(load("Hemap_immunology_Annotations.Rdata"))

#********************************** Full map **********************************
clusters=read.delim("anno_coord_data9544_15pct_bw2.5_updated.txt", stringsAsFactors=F, header=T)
clusters=clusters[clusters$ID%in%annot$GSM.identifier..sample.,]

#*********************************** GSVA input ****************************
gsva=get(load("data8238_dufva_immunological_genes_updated_2016_GSVA_geneperm_lean_eFDR.Rdata"))
bindea=get(load("data8238_all_samples_dufva_bindea_2013_geneset_GSVA.Rdata"))
load("data8238_all_samples_Combined_pathway_signatures_210616_GSVA.Rdata")
gsva_es=rbind(gsva, bindea, gsva_es)
gsva_es=gsva_es[!duplicated(rownames(gsva_es)),]

# match cols gsva
gsva_es=data.frame(gsva_es[,match(colnames(matrix), colnames(gsva_es))])
colnames(gsva_es)=colnames(matrix)

#**************************** Cancermap clusters ****************************
clusters_cancermap=clusters$X2.5..cluster
cluster_cancermap=FUN_MAKE_ALL(clusters_cancermap, "cancermap_cluster", clusters_cancermap, 0)
subclasses=FUN_MAKE_ALL(annot$subclasses, "cancermap_cluster", clusters_cancermap, 0.8)
acute_chronic=FUN_MAKE_ALL(annot$acute, "cancermap_cluster", clusters_cancermap, 0.8)
colorClass=FUN_MAKE_ALL(annot$colorClass, "cancermap_cluster", clusters_cancermap, 0.8)
disease=FUN_MAKE_ALL(annot$disease, "cancermap_cluster", annot$disease, 0.8)
fullclass=FUN_MAKE_ALL(annot$CLASS2, "cancermap_cluster", annot$CLASS2, 0.8)

class_cancermap=FUN_MAKE_CATEGORICAL(clusters_cancermap, "cancermap_cluster")

l.comparisons=list(cluster_cancermap, subclasses, acute_chronic, colorClass, disease,fullclass, class_cancermap)
comparisons_cat=do.call(rbind, l.comparisons)
comparisons_cat=data.frame(comparisons_cat[!duplicated(rownames(comparisons_cat)),], stringsAsFactors = F)
colnames(comparisons_cat)=colnames(matrix)

# combine
l.fm=list(matrix, gsva_es, comparisons_cat)

library(data.table)
fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

# remove rows with few values or NAs
rm=apply(fm, 1, function(v)all(is.na(v)))
fm=fm[!rm,]

save(fm, file="Hemap_immunology_fm_cancermap.Rdata")

write.table(t(c("N:SAMP", as.character(colnames(fm)))), file="Hemap_immunology_fm_cancermap.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(fm, file="Hemap_immunology_fm_cancermap.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)


#********************************** Lymphoma **********************************
matrix=get(load("Hemap_immunology_fm.Rdata"))
annot=get(load("Hemap_immunology_Annotations.Rdata"))

# annot
clusters=read.delim("Hemap_Lymphoma_15pct_genes_BHSNE_mean-shift.txt", stringsAsFactors=F, header=T)
clusters=clusters[clusters$ID%in%annot$GSM.identifier..sample.,]

matrix_sub=matrix[,colnames(matrix)%in%clusters$ID]
annot_sub=annot[annot$GSM.identifier..sample.%in%clusters$ID,]

load("data9544_LYMPHOMA_all_samples_Combined_pathway_drug_signatures_2017_GSVA.Rdata")
rownames(gsva_es)=gsub(" ", "_", rownames(gsva_es))

# match cols gsva
gsva_es=data.frame(gsva_es[,match(colnames(matrix_sub), colnames(gsva_es))])
colnames(gsva_es)=colnames(matrix_sub)

clusters_cancermap=clusters$X1.5..cluster

# comparisons
cluster_subtypes=FUN_MAKE_ALL(annot_sub$CLASS, "cancermap_cluster", clusters_cancermap, 0.8)
cluster_cancermap=FUN_MAKE_ALL(clusters_cancermap, "cancermap_cluster", clusters_cancermap, 0)
cluster_BCL_TCL=FUN_MAKE_ALL(annot_sub$tbLY, "cancermap_cluster", clusters_cancermap, 0.8)
class_cancermap=FUN_MAKE_CATEGORICAL(clusters_cancermap, "cancermap_cluster")

l.comparisons=list(cluster_subtypes, cluster_cancermap, cluster_BCL_TCL, class_cancermap)
comparisons_cat=do.call(rbind, l.comparisons)
comparisons_cat=data.frame(data.matrix(comparisons_cat[!duplicated(rownames(comparisons_cat)),]), stringsAsFactors = F)


# combine
l.fm=list(matrix_sub, data.frame(gsva_es), comparisons_cat)

library(data.table)
fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

# remove rows with few values or NAs
rm=apply(fm, 1, function(v)all(is.na(v)))
fm=fm[!rm,]

save(fm, file="Hemap_LYMPHOMA_immunology_fm.Rdata")

write.table(t(c("N:SAMP", as.character(colnames(fm)))), file="Hemap_LYMPHOMA_immunology_fm.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(fm, file="Hemap_LYMPHOMA_immunology_fm.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)

#********************************** AML **********************************

matrix=get(load("Hemap_immunology_fm.Rdata"))
annot=get(load("Hemap_immunology_Annotations.Rdata"))

# annot
clusters=read.delim("AML_15pct_BHSNE_mean-shift.txt", stringsAsFactors=F, header=T)
clusters=clusters[clusters$ID%in%annot$GSM.identifier..sample.,]

clusters_cancermap=clusters$X1.5..cluster

matrix_sub=matrix[,colnames(matrix)%in%clusters$ID]
annot_sub=annot[annot$GSM.identifier..sample.%in%clusters$ID,]

# GSVA input
gsva=get(load("data9544_AML_all_samples_Combined_pathway_drug_signatures_2017_GSVA.Rdata"))

# match cols gsva
gsva_es=data.frame(gsva_es[,match(colnames(matrix_sub), colnames(gsva_es))])
colnames(gsva_es)=colnames(matrix_sub)

# TCGA clusters
cluster_mapping=read.delim("Table_TCGA_cluster_AML_cluster_assignment.txt", header=T, stringsAsFactors=F, sep="\t")

annot_sub$TCGA_cluster=rep("NA", dim(annot_sub)[1])

annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[1,2]]="TCGA_AML_cluster_1"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[2,2]]="TCGA_AML_cluster_2"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[3:5,2]]="TCGA_AML_cluster_3"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[6:10,2]]="TCGA_AML_cluster_4"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[11,2]]="TCGA_AML_cluster_5"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[12:15,2]]="TCGA_AML_cluster_6"
annot_sub$TCGA_cluster[clusters_cancermap%in%cluster_mapping[16:17,2]]="TCGA_AML_cluster_7"

# comparisons
cluster_TCGA=FUN_MAKE_ALL(annot_sub$TCGA_cluster, "cancermap_cluster", annot_sub$TCGA_cluster, 0.9)
cluster_cancermap=FUN_MAKE_ALL(clusters_cancermap, "cancermap_cluster", clusters_cancermap, 0)
cluster_subtypes=FUN_MAKE_ALL(annot_sub$CLASS, "cancermap_cluster", clusters_cancermap, 0.9)
class_cancermap=FUN_MAKE_CATEGORICAL(clusters_cancermap, "cancermap_cluster")
class_TCGA=FUN_MAKE_CATEGORICAL(annot_sub$TCGA_cluster, "cancermap_cluster")

l.comparisons=list(cluster_TCGA, cluster_cancermap, cluster_subtypes, class_cancermap, class_TCGA)
comparisons_cat=do.call(rbind, l.comparisons)
comparisons_cat=data.frame(data.matrix(comparisons_cat[!duplicated(rownames(comparisons_cat)),]), stringsAsFactors = F)

# combine
l.fm=list(matrix_sub, data.frame(gsva_es), comparisons_cat)

library(data.table)
fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

# remove rows with few values or NAs
rm=apply(fm, 1, function(v)all(is.na(v)))
fm=fm[!rm,]

save(fm, file="Hemap_AML_immunology_fm.Rdata")

write.table(t(c("N:SAMP", as.character(colnames(fm)))), file="Hemap_AML_immunology_fm.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
write.table(fm, file="Hemap_AML_immunology_fm.tsv", sep="\t", col.names=F, row.names=T, quote=FALSE, append=T)

# write.table(all_isolation, file="Hemap_all_isolation.tsv", sep="\t", col.names=F, row.names=F, quote=FALSE, append=F)
