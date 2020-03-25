library(limma)
library(edgeR)
library(parallel)
library(gridExtra)

GIT_HOME="/research/users/ppolonen/git_home/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# all CGAss:
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)

# CCLE gexp data:
ccle=get(load("CCLE_RNASEQ_SYMBOL_COUNTS.Rdata"))
dge=DGEList(ccle)
# keep.exprs <- filterByExpr(dge) # these remove some of the CGA that are rarely expressed
# dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")
ccle=voom(dge, plot=T)$E

# CCLE RPKM
ccle2=get(load("CCLE_RNASEQ_SYMBOL_RPKM.Rdata"))

# CCLE meth data:
ccle.meth=get(load("CCLE_combined_methylation.Rdata"))

# CCLE annot:
annot=read.delim("CCLE_sample_info_file_2012-10-18.txt", header=T, stringsAsFactors = F)

find=intersect(intersect(colnames(ccle.meth), intersect(annot$CCLE.name, colnames(ccle))), colnames(ccle2))

ccle.meth=ccle.meth[,match(find, colnames(ccle.meth))]
ccle=ccle[,match(find, colnames(ccle))]
annot=annot[match(find, annot$CCLE.name),]
ccle2=ccle2[,match(find ,colnames(ccle2))]

# compute HLA-scores:
dat_a3=2^ccle[rownames(ccle)%in%c("HLA-DMA",
                                  "HLA-DMB",
                                  "HLA-DPB1",
                                  "HLA-DRA",
                                  "HLA-DRB1"),]+1

gm3=log2(t(apply(dat_a3, 2, gm_mean)))
rownames(gm3)="HLAII_SCORE"

dat2=rbind(ccle, gm3)

annotv=colnames(d)[colnames(d)%in%colnames(ccle.meth)[ccle.meth["CIITA@chr16_10971271",]>0.2]]


annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"plasma_cell_myeloma"]="MM"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"mantle_cell_lymphoma"]="MCL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"diffuse_large_B_cell_lymphoma"]="DLBCL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"chronic_lymphocytic_leukaemia-small_lymphocytic_lymphoma"]="CLL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"blast_phase_chronic_myeloid_leukaemia"]="CML"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"anaplastic_large_cell_lymphoma"]="ALCL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"adult_T_cell_lymphoma-leukaemia"]="TCL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"acute_myeloid_leukaemia"]="AML"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"acute_lymphoblastic_T_cell_leukaemia"]="T-ALL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"acute_lymphoblastic_B_cell_leukaemia"]="pre-B-ALL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"Hodgkin_lymphoma"]="CHL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%"Burkitt_lymphoma"]="BL"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%c("B_cell_lymphoma_unspecified")]="BCL, unspecified"
annot$Hist.Subtype_hem[annot$Hist.Subtype1%in%c("peripheral_T_cell_lymphoma_unspecified")]="TCL"

d=dat2
groups=factor(annot$Hist.Subtype_hem[grep("lymphoma|leukaemia|myeloma", annot$Hist.Subtype1)])

logicalVectors=lapply(unique(annot$Hist.Subtype_hem), function(g){annot$Hist.Subtype_hem%in%g})
names(logicalVectors)=unique(annot$Hist.Subtype_hem)
logicalVectors=logicalVectors[!is.na(names(logicalVectors))]

logicalVectors=logicalVectors[-6] # remove unspesified

a=names(sort(table(t.df[,3]), decreasing = T))

t.df=do.call(rbind, lapply(a, function(n){
  b=t.df[t.df[,3]%in%n,]
  b=b[order(b[,2], decreasing = T),]
}))

# data:
genelist=unique(t.df[,1])

ccle_antg=ccle[match(genelist,rownames(ccle)),]
antigens_expresssed=colSums(ccle2[match(genelist,rownames(ccle2)),]>0.5, na.rm = T)

disease=annot$Hist.Subtype_hem

go.thr=gsub("@.*", "", rownames(ccle.meth))
ccle.meth=ccle.meth[go.thr%in%genelist,]
go.thr=go.thr[go.thr%in%genelist]
  
# load("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/data/CCLE_METHYLATION/probelists.Rdata")

highest_anticor=unlist(lapply(unique(go.thr), function(g){
  
  if(!g%in%go.thr)return(NULL)
  if(!g%in%rownames(ccle))return(NULL)
  
  m=ccle.meth[go.thr%in%g,grep("lymphoma|leukaemia|myeloma", annot$Hist.Subtype1)]
  cc=as.numeric(ccle[rownames(ccle)%in%g,grep("lymphoma|leukaemia|myeloma", annot$Hist.Subtype1)])
  
  if(all(rowSums(is.na(m))>140))return(NULL)
  
  # remove probes that are not found in most datasets:
  m=m[!rowSums(!is.na(m))<20,]
  
  if(dim(m)[1]==0)return(NULL)
  
  res=cor(cc,t(m), use="pairwise.complete.obs", method = "spearman")
  
  # rm poor correlations:
  res=res[res<(-0.2)]
  
  if(!length(res))return(NULL)
  
  anticor=rownames(m[which.min(res),])
}))


# data frame
df_anno=data.frame("Disease"=disease, "Number_Antigens_expressed"=as.integer(antigens_expresssed))

# sort based on antigens expressed and disease
sample_order=order(match(annot$Hist.Subtype_hem, c("T-ALL", "pre-B-ALL", "AML", "CML", "CLL","MM","BL","CHL","DLBCL","MCL","BCL, unspecified","TCL","ALCL")), antigens_expresssed)

df_anno=df_anno[sample_order,]

ccle_antg=ccle_antg[,sample_order[!is.na(df_anno$Disease)]]
df_anno=df_anno[!is.na(df_anno$Disease),]

ccle_antg_meth=ccle.meth[,match(colnames(ccle_antg), colnames(ccle.meth))]

cname=colnames(ccle_antg)
colnames(ccle_antg)=gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","", colnames(ccle_antg))
colnames(ccle_antg_meth)=gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE","", colnames(ccle_antg_meth))

ccle_antg_meth=ccle_antg_meth[rownames(ccle_antg_meth)%in%highest_anticor,]
ccle_antg_meth=ccle_antg_meth[match(genelist, gsub("@.*", "", rownames(ccle_antg_meth))),]

ccle_antg_meth2=ccle_antg_meth
ccle_antg_meth2[ccle_antg_meth2<=0.2]=0 # this works best, tested
ccle_antg_meth2[ccle_antg_meth2>0.2]=1

ccle_antg_meth2=ccle_antg_meth2[!rowSums(is.na(ccle_antg_meth2))>100,]

sum_hypometh=colSums(ccle_antg_meth2==0, na.rm = T)

cor(df_anno$Number_Antigens_expressed, sum_hypometh, use="complete.obs")

df_anno=data.frame(df_anno, "Number_Antigens_methylated"=sum_hypometh)

library(ComplexHeatmap)
ha = HeatmapAnnotation(df=df_anno[,1,drop=F], Number.Antigens = anno_barplot(df_anno$Number_Antigens_expressed, axis = T, gp = gpar(fill = "indianred")), Number.Me.Antigens = anno_barplot(sum_hypometh, axis = T, gp = gpar(fill = "darkgoldenrod")), height = unit(4, "inch"))

ccle_antg=t(scale(t(ccle_antg)))
ccle_antg[ccle_antg<(-4)]=-4
ccle_antg[ccle_antg>(4)]=4

pdf("Fig6D_CCLE_ComplexHeatmap_testis_antigens.pdf", height = 15, width = 30)
Heatmap(ccle_antg, top_annotation = ha, name = "TestisAntigens", cluster_columns = F, cluster_rows = F, col=colorRamp2(c(-4,-2, 0,2, 4), c("#063061","#579ec9", "white", "#d45d4a", "#67001e")))
dev.off()

# pdf("CCLE_ComplexHeatmap_testis_antigens_methylation.pdf", height = 15, width = 30)
Heatmap(ccle_antg_meth, top_annotation = ha, name = "TestisAntigens", cluster_columns = F, cluster_rows = F)
# dev.off()



