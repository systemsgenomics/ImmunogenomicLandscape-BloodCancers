GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/pathway_analysis/functions.GSEA.R"))

library(data.table)
library(parallel)
library(GSVA)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

fm=get(load("GSE98588_fm.Rdata"))
annot=get(load("GSE98588_annot.Rdata"))
profile=get(load("GSE98588_DLBCL_mixtureM_profile.Rdata"))

# exclude testis dlbcl
profile=profile[,!colnames(fm)%in%"DLBCL_LS2208"]
annot=annot[!colnames(fm)%in%"DLBCL_LS2208",]
fm=fm[,!colnames(fm)%in%"DLBCL_LS2208"]

# CGAs
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)


# Choose GCB or ABC, or ""
name="ABC"
GSE="GSE98588_DLBCL"

if(name=="GCB"){
  GSE="GSE98588_DLBCL_GCB"
  annot=annot[fm["B:SAMP:COO_byGEP_GCB",]==1,]
  profile=profile[,fm["B:SAMP:COO_byGEP_GCB",]==1]
  fm=fm[,fm["B:SAMP:COO_byGEP_GCB",]==1]
}

if(name=="ABC"){
  GSE="GSE98588_DLBCL_ABC"
  annot=annot[fm["B:SAMP:COO_byGEP_ABC",]==1,]
  profile=profile[,fm["B:SAMP:COO_byGEP_ABC",]==1]
  fm=fm[,fm["B:SAMP:COO_byGEP_ABC",]==1]
}

gexp=fm[grepl("N:GEXP:", rownames(fm)),]
rownames(gexp)=gsub("N:GEXP:", "", rownames(gexp))

res=fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/GSE98588_DLBCL_antigen_correlations.tsv", data.table = F)

scores=c("N:SAMP:CytolyticScore",  "N:SAMP:HLAIScore",  "N:SAMP:HLAIIScore")

feats=c("N:SAMP:numberOfCNAs","N:SAMP:numberOfMutations", "N:SAMP:numberOfChromosomalRearrangements")

mut=c("B:GNAB:KLHL6", "B:GNAB:CD58","B:GNAB:SGK1", "B:GNAB:CD83", "B:GNAB:MYD88","B:GNAB:HIST1H1E", "B:GNAB:HIST1H2BK","B:CNVR:6Q:LOSS","B:GNAB:BTG1","B:GNAB:HLA-A", "B:GNAB:ETV6", "B:GNAB:UBE2A", "B:CNVR:1P13_1:LOSS", "B:GNAB:SPEN", "B:GNAB:NFKBIA", "B:GNAB:GNA13")

other.genes="N:GEXP:CD58"

#************************************** GSEA **************************************

GENESETS="Combined_pathway_signatures_2017_filtered_robust.gmt"
WD=file.path(getwd(), "GSEA")
OUTDIR=file.path(getwd(), "GSEA")

if(name=="GCB"){
  # run GSEA:
  # command=run.GSEA(data=gexp, cls.vector = as.numeric(fm["N:SAMP:nCGA",]), datatype = "N", GENESETS = GENESETS, dataname = GSE, clsname = "nCGA", WD=WD, OUTDIR=OUTDIR)
  # try(system(command))
  
  a=read.delim(file.path(WD, "GSE98588_DLBCL_GCB_nCGA_continuous_phenotype.Gsea.1552654676702/gsea_report_for_feat_pos_1552654676702.xls"), stringsAsFactors = F)
  b=read.delim(file.path(WD, "GSE98588_DLBCL_GCB_nCGA_continuous_phenotype.Gsea.1552654676702/gsea_report_for_feat_neg_1552654676702.xls"), stringsAsFactors = F)
}

if(name=="ABC"){
  # run GSEA:
  # command=run.GSEA(data=gexp, cls.vector = as.numeric(fm["N:SAMP:nCGA",]), datatype = "N", GENESETS = GENESETS, dataname = GSE, clsname = "nCGA", WD=WD, OUTDIR=OUTDIR)
  # try(system(command))
  
  a=read.delim(file.path(WD, "GSE98588_DLBCL_ABC_nCGA_continuous_phenotype.Gsea.1552654669505/gsea_report_for_feat_pos_1552654669505.xls"), stringsAsFactors = F)
  b=read.delim(file.path(WD, "GSE98588_DLBCL_ABC_nCGA_continuous_phenotype.Gsea.1552654669505/gsea_report_for_feat_neg_1552654669505.xls"), stringsAsFactors = F)
}
#*****************************************************************************************************

c=rbind(a, b)

# get GSVA visualization for the pathways
library(GSVA)
library(parallel)

# Geneset list
Onc.pathways=read.delim(GENESETS, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(GENESETS, sep = '\t'), na.rm = T)), fill = TRUE)

# Make list
listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
B=A[!A==""&!A=="NA"]}, mc.cores=6)

names(listA) <- Onc.pathways[,1]

# visualize using GSVA
viz_scores=gsva(expr = data.matrix(gexp), gset.idx.list = listA, parallel.sz=8, method="gsva", tau=0.25)

# make a complex heatmap
significant=c(c[c$FWER.p.val<0.001,1])

if(length(significant)>11)significant=significant[1:11]

dat_plot=viz_scores[match(significant, toupper(rownames(viz_scores))),]

profile[profile==-1] = 0
profile[data.matrix(gexp)<5]=0
f=sort(rowSums(profile[rownames(profile)%in%unique(t.df[,1]),]), decreasing = T)
f2=signif(f/dim(profile)[2],2)
names(f2)=f2
f2=f2[!f==0]
f=f[!f==0]

ag=c(paste0("N:GEXP:", names(f)))

ag=ag[ag%in%rownames(fm)]

add=data.frame(scale(t(pl2)))
add[add>2]=2
add[add<(-2)]=-2

pl1=data.matrix(fm[scores,order(fm["N:SAMP:nCGA",], decreasing = F)])
pl1=t(scale(t(pl1)))
pl1[pl1>2]=2
pl1[pl1<(-2)]=-2

pl2=data.matrix(fm[feats,order(fm["N:SAMP:nCGA",], decreasing = F)])
pl2[pl2>250]=250
# pl2=t(scale(t(pl2)))
# pl2[pl2>2]=2
# pl2[pl2<(-2)]=-2

pl3=data.matrix(fm[ag,order(fm["N:SAMP:nCGA",], decreasing = F)])
pl3=t(scale(t(pl3)))
pl3[pl3>2]=2
pl3[pl3<(-2)]=-2

pl4=data.matrix(fm[mut,order(fm["N:SAMP:nCGA",], decreasing = F)])

pl5=data.matrix(fm[other.genes,order(fm["N:SAMP:nCGA",], decreasing = F)])
pl5=t(scale(t(pl5)))
pl5[pl5>2]=2
pl5[pl5<(-2)]=-2

pl6=data.matrix(dat_plot[,order(fm["N:SAMP:nCGA",], decreasing = F)])

ann=data.frame("n.CGA"=t(fm["N:SAMP:nCGA",order(fm["N:SAMP:nCGA",], decreasing = F)]), stringsAsFactors = F)

library(ComplexHeatmap)
library(circlize)
library(multipanelfigure)

rownames(pl3)=gsub(".*.:", "", rownames(pl3))
rownames(pl4)=gsub(".*.:", "", rownames(pl4))
rownames(pl6)=gsub("-.*.|_HOMO_SAPIENS", "", rownames(pl6))

p1=Heatmap(pl1, top_annotation = HeatmapAnnotation("n.CGA"=anno_barplot(ann, height = unit(10, "mm")), df = annot$COO_byGEP), cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="cluster", col = colorRamp2(c(2, 0, -2), c("red", "white", "blue")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl1)[1]+10, "mm"))

p2=Heatmap(pl2, top_annotation = HeatmapAnnotation("nr.cnv"=anno_barplot(data.frame(pl2[1,]), height = unit(10, "mm")), "nr.mut"=anno_barplot(data.frame(pl2[2,]), height = unit(10, "mm"), ylim=c(0,250)),"nr.strrearr"=anno_barplot(data.frame(pl2[3,]), height = unit(10, "mm"))), cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="Samp-features", col=colorRamp2(c(-2, 0, 2), c("grey85", "white", "indianred")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl2)[1]+10, "mm"))

p3=Heatmap(pl3, cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="Scaled log2\nexpression", col=colorRamp2(c(-2,0,2), c("blue", "white", "red")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl3)[1], "mm")) + Heatmap(f2, width = unit(5, "mm"), row_names_gp = gpar(fontsize = 6))

p4=Heatmap(pl4, cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="Mutations", col=colorRamp2(c(0,1), c("grey85", "darkgreen")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl4)[1], "mm")) 

p5=Heatmap(pl5, cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="Scaled log2\nexpression", col=colorRamp2(c(-2,0,2), c("blue", "white", "red")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl5)[1], "mm"))

p6=Heatmap(pl6, cluster_columns = F, cluster_rows = F, row_names_side = "left", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="Scaled log2\nexpression", col=colorRamp2(c(-0.6, 0, 0.6), c("grey50", "white", "red")), show_column_names = F, width = unit(40, "mm"), height = unit(2*dim(pl6)[1], "mm"))

panels=list(p1, p2, p3, p4, p5, p6)
panels=panels[!is.null(panels)]
nrows=sum(c(dim(pl1)[1], dim(pl2)[1], dim(pl3)[1], dim(pl4)[1], dim(pl5)[1], dim(pl6)[1]))

figure <- multi_panel_figure(width = 200, height = 220+nrows*3, rows = length(panels), columns = 1, panel_label_type = "lower-alpha")

for(i in seq(panels))figure <- fill_panel(figure,panels[[i]], row = i, column = 1)

save_multi_panel_figure(figure, filename=paste0("Fig6H_FigS6H_", "_", GSE, "_CGA_heatmap.pdf"))



gexp=fm[grepl("N:GEXP:", rownames(fm)),]
rownames(gexp)=gsub("N:GEXP:", "", rownames(gexp))

# All pathways
GENESETS="/research/work/ppolonen/genesets/Combined_pathway_signatures_2017_filtered_robust.gmt"
WD="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/"
OUTDIR="/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/"

# get GSVA visualization for the pathways
# Geneset list
Onc.pathways=read.delim(GENESETS, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(GENESETS, sep = '\t'), na.rm = T)), fill = TRUE)

# Make list
listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
B=A[!A==""&!A=="NA"]}, mc.cores=6)

names(listA) <- Onc.pathways[,1]

viz_scores=gsva(expr = data.matrix(gexp), gset.idx.list = listA, parallel.sz=8, method="gsva", tau=0.25)

expressed_testis_num=colSums(profile[rownames(profile)%in%unique(t.df$gene),])
feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0 CGA"
feat_class[expressed_testis_num>=1&expressed_testis_num<=2]="1-2 CGA"
feat_class[expressed_testis_num>=3]=">3 CGA"

logicalVectors=lapply(unique(feat_class), function(cl)feat_class%in%cl)
names(logicalVectors)=paste("ABC", unique(feat_class))

annof=data.frame("HLAI"=annot$HLAIScore, "HLAII"=annot$HLAIIScore, "CytolyticScore"=annot$CytolyticScore, "TNFA_SIGNALING_VIA_NFKB"=viz_scores[rownames(viz_scores)%in%"TNFA_SIGNALING_VIA_NFKB-MsigDB_HALLMARKS",])

genelist=c("HLAI", "HLAII", "CytolyticScore", "TNFA_SIGNALING_VIA_NFKB")
p.all=lapply(genelist, plot.boxplot, logicalVectors = logicalVectors[c(1,3,2)], data = t(annof),order.bl = F,spread = T)

TestGeneWilcox("HLAI", data = t(annof), logicalVectors = logicalVectors[c(1,3,2)], logicalVector_normals = logicalVectors[c(1,3,2)], ALTERNATIVE = "less")
TestGeneWilcox("HLAII", data = t(annof), logicalVectors = logicalVectors[c(1,3,2)], logicalVector_normals = logicalVectors[c(1,3,2)], ALTERNATIVE = "less")
TestGeneWilcox("CytolyticScore", data = t(annof), logicalVectors = logicalVectors[c(1,3,2)], logicalVector_normals = logicalVectors[c(1,3,2)], ALTERNATIVE = "less")
TestGeneWilcox("TNFA_SIGNALING_VIA_NFKB", data = t(annof), logicalVectors = logicalVectors[c(1,3,2)], logicalVector_normals = logicalVectors[c(1,3,2)], ALTERNATIVE = "less")

ggsave(plot = p.all[[1]], filename = paste0(GSE, "_HLA.pdf"), width = unit(3.25, "cm"), height = unit(3, "cm"))
ggsave(plot = p.all[[2]], filename = paste0(GSE, "_HLAII.pdf"), width = unit(3.25, "cm"), height = unit(3, "cm"))
ggsave(plot = p.all[[3]], filename = paste0("FigS6I_", GSE, "_CytScore.pdf"), width = unit(3.25, "cm"), height = unit(3, "cm"))
ggsave(plot = p.all[[4]], filename = paste0(GSE, "_TNFA.pdf"), width = unit(3.25, "cm"), height = unit(3, "cm"))
