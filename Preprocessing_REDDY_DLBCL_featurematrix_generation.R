GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "featurematrix/functions_generate_fm.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/DLBCL_CELL")

mut=data.table::fread("Mutation_data.txt", data.table = F, dec = ",")
cnvr=data.table::fread("cnvr_data.txt", data.table = F, dec = ",")
clin=data.table::fread("clinical_data.txt", data.table = F, dec = ",")
gexp=data.table::fread("GeneExpression_FPKM_QuantileNorm_624Samples.txt", data.table = F, dec = ",")
ihc=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/DLBCL_CELL/DLBCL_IHC_results_220117_Reddy_ID.txt", data.table = F)

#************************************************ filter data: ************************************************
clin.filt=clin[,-1]
rownames(clin.filt)=clin[,1]

clin.filt$`Tumor Purity`[clin.filt$`Tumor Purity`%in%"< 30%"]="under_30"
clin.filt$`Tumor Purity`[clin.filt$`Tumor Purity`%in%"30 to 70%"]="30to70"
clin.filt$`Tumor Purity`[clin.filt$`Tumor Purity`%in%"70% or more"]="over_70"

mut.filt=mut[,45:dim(mut)[2]]
mut.filt=do.call(rbind, lapply(unique(mut$Gene.refGene), function(gene)(colSums(mut.filt[mut$Gene.refGene%in%gene,])>0)*1))
rownames(mut.filt)=unique(mut$Gene.refGene)

colnames(mut.filt)==clin$`Sample  ID`

cnvr.filt=data.matrix(cnvr[,4:dim(cnvr)[2]])
rownames(cnvr.filt)=cnvr[,1]

gexp.filt=data.matrix(gexp[-1,3:dim(gexp)[2]])

library(HGNChelper)
gexp[,2]=findExcelGeneSymbols(gexp$V2, mog.map = read.csv(system.file("extdata/mog_map.csv", package =
                                                      "HGNChelper"), as.is = TRUE), regex = "impossibletomatch^")

rownames(gexp.filt)=make.unique(gexp[-1,2])
colnames(gexp.filt)=gexp[1,3:dim(gexp)[2]]

ihc.filt=ihc[match(colnames(gexp.filt), gsub("Sample_", "", ihc$Reddy_ID)),6:25]
rownames(ihc.filt)=colnames(gexp.filt)

#************************************************ Immunoscores ************************************************
dat_a3=gexp.filt[rownames(gexp.filt)%in%c("HLA-DMA",
                                "HLA-DMB",
                                "HLA-DPA1",
                                "HLA-DPB1",
                                "HLA-DRA",
                                "HLA-DRB1"),]

dat3=2^dat_a3+0.01
gm2=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm2)="HLAIIScore"

dat_a3=gexp.filt[rownames(gexp.filt)%in%c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM"),]

dat3=2^dat_a3+0.01
gm3=log2(t(apply(dat3, 2, gm_mean)))
rownames(gm3)="CytolyticScore"

classification1=data.frame(t(rep("medium", length(gm3))))
zscore=as.numeric(scale(t(gm3)))
classification1[zscore>=1]="high"
classification1[zscore<=(-1)]="low"
rownames(classification1)="CytolyticScore" 
colnames(classification1)=colnames(gexp.filt)

classification3=data.frame(t(rep("medium", length(gm2))))
zscore=as.numeric(scale(t(gm2)))
classification3[zscore>=1]="high"
classification3[zscore<=(-1)]="low"
rownames(classification3)="HLAIIScore" 
colnames(classification3)=colnames(gexp.filt)

classification=data.frame(t(rbind(classification1,classification3)), stringsAsFactors = F)

#************************************************ Immunoscores ************************************************
immunoscores=as.data.frame(t(rbind(gm2, gm3)),stringsAsFactors=F)
immunoscoresfm=make.features(df = immunoscores, datatype="SAMP", prefix="", make.pairwise = F)
classificationfm=make.features(df = classification, datatype="SAMP", prefix="", make.pairwise = F)
colnames(immunoscoresfm)=paste0("DLBCL_", colnames(immunoscoresfm))
colnames(classificationfm)=paste0("DLBCL_", colnames(classificationfm))

clindatfm=make.features(df = clin.filt, datatype="CLIN", prefix="", make.pairwise = F)
colnames(clindatfm)=paste0("DLBCL_", colnames(clindatfm))

mutfm=make.features(data.frame(t(mut.filt)), datatype="GNAB", prefix="")
colnames(mutfm)=paste0("DLBCL_", colnames(mutfm))

cnvfm=make.features(data.frame(t(cnvr.filt)), datatype="CNVR", prefix="")
colnames(cnvfm)=paste0("DLBCL_", colnames(cnvfm))

gexpfm=make.features(data.frame(t(gexp.filt)), datatype="GEXP", prefix="")
colnames(gexpfm)=paste0("DLBCL_", colnames(gexpfm))

ihcfm=make.features(ihc.filt, datatype="SAMP", prefix="")

l.fm=list(immunoscoresfm,classificationfm, clindatfm[,match(colnames(immunoscoresfm), colnames(clindatfm))], mutfm[,match(colnames(immunoscoresfm), colnames(mutfm))], cnvfm[,match(colnames(immunoscoresfm), colnames(cnvfm))], gexpfm) #ihcfm 

library(data.table)
fm=rbindlist(l.fm, use.names=T, fill=T)

fm=data.frame(fm, stringsAsFactors=F)
rownames(fm)=unlist(lapply(l.fm, rownames))

annot=data.frame(clin, stringsAsFactors = F)
annot$Sample..ID=paste0("DLBCL_", annot$Sample..ID)
annot=annot[match(colnames(fm), annot$Sample..ID),]

save(fm, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/DLBCL_CELL/REDDY_DLBCL_fm.Rdata")
save(annot, file="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/DLBCL_CELL/REDDY_DLBCL_annot.Rdata")


coordinates.subtype=annot$ABC.GCB..RNAseq.
annot2=annot[match(colnames(gexp.filt), annot$Sample..ID),]
gexp=gexp.filt
coordinates.subtype=annot2$ABC.GCB..RNAseq.

save(list = c("gexp", "coordinates.subtype"), file="Reddy_DLBCL_subtypes.Rdata")



