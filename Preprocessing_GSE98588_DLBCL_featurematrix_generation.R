GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/")

gexp=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_symbol_rma_normalized.Rdata"))
colnames(gexp)=gsub("_DLBCL", "", colnames(gexp))

# run mixturemodel:
# source(file.path(GIT_HOME, "statistics/mixtureModel.R"))
# profile=do.call(rbind, parallel::mclapply(rownames(gexp), mixtureM.gene, d.matrix = data.matrix(gexp), FIX = T, mc.cores=15))
# save(profile, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_DLBCL_mixtureM_profile.Rdata")

profile=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_DLBCL_mixtureM_profile.Rdata"))
profile[profile==-1] = 0
profile[data.matrix(gexp)<5]=0

# save(immunoscores, file="GSE98588_immunoscores.Rdata")
library(data.table)
clin=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_clinical.txt", stringsAsFactors = F, na.strings = "na")
mut=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_genetic.txt", row.names = 1, stringsAsFactors = F, na.strings = "na")[,-1]
samp=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_sampleinfo.txt", row.names = 1, stringsAsFactors = F, na.strings = "na")[,-1]
signature=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/GSE98588_mutsig.txt", row.names = 1, stringsAsFactors = F)[,-1]

samp=samp[match(colnames(gexp),samp$sample_geo_accession),]
clin=clin[match(colnames(gexp),clin$sample_geo_accession),]
colnames(gexp)=clin$individual_id
colnames(profile)=clin$individual_id
mut=mut[,match(colnames(gexp), colnames(mut))]
names_mut=rownames(mut)
mut=data.frame(t(mut), check.names = F)
signature=signature[match(colnames(gexp),rownames(signature)),]

samp$Normal.coverage=as.numeric(gsub("\\*", "", samp$Normal.coverage))
samp[samp=="no normal"]=NA
samp$tumor_in_normal_estimate=as.numeric(samp$tumor_in_normal_estimate)
samp$tumor_in_normal_post_calls=as.numeric(samp$tumor_in_normal_post_calls)
samp$tumor_in_normal_pre_calls=as.numeric(samp$tumor_in_normal_pre_calls)
samp$Mutations.Recovered.By.deTiN=as.numeric(samp$Mutations.Recovered.By.deTiN)
samp=samp[,!colnames(samp)%in%c("TissueSite", "MoniChapuy.Cancer.Cell.series.", "Lohr.et.al.series.", "Cohort")]

samp=as.data.frame(samp, stringsAsFactors = F)

cnv=mut[,86:158]
cnv2=mut[,86:158]
mut=mut[,1:85]
mut2=mut[,1:85]

colnames(cnv2)[grepl("SV", colnames(cnv2))]=gsub("SV_", "", gsub("/", "_", paste0(gsub("SV:", "", colnames(cnv2)[grepl("SV", colnames(cnv2))]), ":SV")))
colnames(cnv)[grepl("SV", colnames(cnv))]=gsub("SV_", "", gsub("/", "_", paste0(gsub("SV:", "", colnames(cnv)[grepl("SV", colnames(cnv))]), ":SV")))
colnames(cnv)=gsub("\\.", "_", toupper(colnames(cnv)))
colnames(cnv2)=gsub("\\.", "_", toupper(colnames(cnv2)))

# mut[mut==1|mut==2]=1
mut[mut==1]=0 # take only nonsynonymous
mut[mut==2]=1 # take only nonsynonymous

mut2[mut2==1]="synonymous"
mut2[mut2==2]="nonsynonymous"
mut2[mut2==0]="absent"

amp=cnv[,grep("AMP", colnames(cnv))]
del=cnv[,grep("DEL", colnames(cnv))]

amp[amp==1]="GAIN"
amp[amp==2]="AMP"
del[del==1]="LOSS"
del[del==2]="DEL"

cnv=cbind(amp, del)
cnv[cnv==0]="absent"

cnv2[cnv2==3]=1
cnv2[cnv2==1]=1
cnv2[cnv2==2]=1
cnv2[cnv2==0]=0

colnames(cnv2)=gsub("DEL", "DEL_LOSS", colnames(cnv2))
colnames(cnv2)=gsub("AMP", "AMP_GAIN", colnames(cnv2))
colnames(cnv2)=gsub("SV_", "AMP_GAIN", colnames(cnv2))

colnames(cnv)=gsub(":DEL|\\(BCL2\\)|:AMP", "", colnames(cnv))
colnames(cnv2)=gsub("\\(BCL2\\)", "", colnames(cnv2))

# remove if no both del-loss or amp-gain:
take=apply(cnv, 2, function(v)length(unique(v))==3)
# cnv2=cnv2[,grepl("SV", colnames(cnv2))]

# compute antigen scores
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
t.df=t.df[order(t.df[,3]),]

co.stim=fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/costim_ligands_final.txt", data.table = F)

# data:
genelist=unique(t.df[,1])

# rank patients by number of testis antigens expressed
expressed_testis_num=data.frame(t(colSums(profile[rownames(profile)%in%genelist,])))
colnames(expressed_testis_num)=colnames(gexp)
rownames(expressed_testis_num)="nCGA"

expressed_testis_num2=as.numeric(expressed_testis_num)
expressed_testis_num2[expressed_testis_num2>=4]="over4_antigens"
expressed_testis_num2[expressed_testis_num2>=1&expressed_testis_num2<4]="1to3_antigens"
expressed_testis_num2[expressed_testis_num2==0]="0_antigens"

# expressed_testis_num3=as.numeric(expressed_testis_num)
# expressed_testis_num3[expressed_testis_num3>=2]="2orMore_antigens"
# expressed_testis_num3[expressed_testis_num3<2]="0to1_antigens"

expressed_testis_num_c=data.frame(rbind(expressed_testis_num2), stringsAsFactors = F)
rownames(expressed_testis_num_c)=c("catCGA")


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

# get GSVA visualization for the pathways
library(GSVA)
library(parallel)

# # Geneset list
# GENESETS="/research/work/ppolonen/genesets/Combined_pathway_signatures_2017_filtered_robust.gmt"
# Onc.pathways=read.delim(GENESETS, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(GENESETS, sep = '\t'), na.rm = T)), fill = TRUE)
# 
# # Make list
# listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
# B=A[!A==""&!A=="NA"]}, mc.cores=6)
# 
# names(listA) <- Onc.pathways[,1]
# 
# viz_scores=gsva(data.matrix(gexp), gset.idx.list = listA, parallel.sz=8, tau=0.25)
# 
# pwscores=as.data.frame(t(viz_scores))

# pwfm=make.features(df = pwscores, datatype="SAMP", prefix="")
# colnames(pwfm)=colnames(viz_scores)

immunoscores=as.data.frame(t(rbind(gm1, gm2, gm3,expressed_testis_num)),stringsAsFactors=F)

clindatfm=make.features(df = clin, datatype="CLIN", prefix="")
colnames(clindatfm)=colnames(gexp)

immunoscoresfm=make.features(immunoscores, datatype="SAMP", prefix="")
colnames(immunoscoresfm)=colnames(gexp)

immunoscoresfm2=make.features(as.data.frame(data.frame(t(expressed_testis_num_c), classification, stringsAsFactors = F), stringsAsFactors = F), datatype="SAMP", prefix="")
colnames(immunoscoresfm2)=colnames(gexp)

sampfm=make.features(samp, datatype="SAMP", prefix="")
colnames(sampfm)=colnames(gexp)

mutfm=make.features(mut, datatype="GNAB", prefix="")
colnames(mutfm)=colnames(gexp)

mutfm=mutfm[!grepl("_vs_|_and_|_absent", rownames(mutfm)),]

mutfm2=make.features(mut2, datatype="GNAB", prefix="")
colnames(mutfm2)=colnames(gexp)

mutfm2=mutfm2[!grepl("_vs_|_and_|_absent", rownames(mutfm2)),]

cnvfm=make.features(cnv, datatype="CNVR", prefix="")
colnames(cnvfm)=colnames(gexp)
cnvfm=cnvfm[!grepl("_vs_|_and_|_absent", rownames(cnvfm)),]
rownames(cnvfm)=gsub("_GAIN", ":GAIN", rownames(cnvfm))
rownames(cnvfm)=gsub("_AMP", ":AMP", rownames(cnvfm))
rownames(cnvfm)=gsub("_DEL", ":DEL", rownames(cnvfm))
rownames(cnvfm)=gsub("_LOSS", ":LOSS", rownames(cnvfm))

cnvfm2=make.features(cnv2, datatype="CNVR", prefix="")
colnames(cnvfm2)=colnames(gexp)

signaturefm=make.features(signature, datatype="SAMP", prefix="")
colnames(signaturefm)=colnames(gexp)

gexpfm=data.frame(gexp)
rownames(gexpfm)=paste0("N:GEXP:", rownames(gexp))

gexpfm2=data.frame(profile)
rownames(gexpfm2)=paste0("B:GEXP:", rownames(gexpfm2))

l.fm=list(immunoscoresfm,immunoscoresfm2, clindatfm, sampfm, mutfm, cnvfm,cnvfm2, signaturefm, gexpfm, gexpfm2[!rowSums(gexpfm2)==0,])

library(data.table)
fm=rbindlist(l.fm, use.names=F, fill=F)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

annot=data.frame(clin, immunoscores, samp, stringsAsFactors = F)

save(fm, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/GSE98588_fm.Rdata")
save(annot, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/GSE98588_annot.Rdata")


# Chapuy data
fm=get(load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/GSE98588_fm.Rdata"))

load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/GSE98588_annot.Rdata")

gexp=fm[grepl("N:GEXP", rownames(fm)),]
rownames(gexp)=gsub("N:GEXP:", "", rownames(gexp))

coordinates.subtype=annot$COO_byGEP

save(list = c("gexp", "coordinates.subtype"), file="Chapuy_DLBCL_subtypes.Rdata")

# make annotation table for these cnv
cnv_annot=fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/41591_2018_16_MOESM8_ESM.txt", data.table=F)
cnv_annot=data.frame("Name"=cnv_annot$Name, "hgnc_symbol"=cnv_annot$hgnc_symbol, stringsAsFactors = F)

cnv_annot[cnv_annot[,2]%in%"01-Mar",2]="MARCH1"
cnv_annot[cnv_annot[,2]%in%"02-Mar",2]="MARCH2"

cnv_annot[,1]=gsub(":DEL", ":DEL_LOSS", cnv_annot[,1])
cnv_annot[,1]=gsub(":AMP", ":AMP_GAIN", cnv_annot[,1])
cnv_annot[,1]=gsub("\\(BCL2\\)", "", cnv_annot[,1])
cnv_annot[,1]=paste0("B:CNVR:", FIX_NAME(cnv_annot[,1]))

cnv_annot2=cnv_annot
cnv_annot3=cnv_annot
cnv_annot2[,1]=gsub("AMP_|DEL_", "", cnv_annot2[,1])
cnv_annot3[,1]=gsub("_GAIN|_LOSS", "", cnv_annot3[,1])

cnv_annot=rbind(cnv_annot, cnv_annot2, cnv_annot3)
cnv_annot=cnv_annot[cnv_annot[,1]%in%rownames(fm),]

fwrite(cnv_annot, "/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/GSE98588/41591_2018_16_MOESM8_ESM_CNV_ANNOT.txt")
