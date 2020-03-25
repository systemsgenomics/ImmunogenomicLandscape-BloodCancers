GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "pathway_analysis/functions.GSEA.R"))

library(GSVA)
library(parallel)
library(ComplexHeatmap)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")
t.df = read.delim("t.antigen_df.txt", stringsAsFactors=F, header=T)
t.df=t.df[order(t.df[,3]),]

annot = get(load("Hemap_immunology_Annotations.Rdata"))

# gexp data
data=t(get(load("data9544_with_gene_symbols.RData")))
data=data[,colnames(data)%in%annot$GSM.identifier..sample.]
  
# data:
genelist=t.df[t.df[,3]%in%"Cancer_Myeloma",1]

# take MM
GSE=c("GSE16716,GSE24080")
# GSE=c("GSE19784")

gexp=data[,annot$GSE.identifier..experiment.%in%GSE]
annot2=annot[annot$GSE.identifier..experiment.%in%GSE,]

# GSE="All_samples"
# gexp=data[,annot$colorClass=="MM"]
# annot2=annot[annot$colorClass=="MM",]


# rank patients by number of testis antigens expressed
# profiles
profile=get(load("mixtureM_profile.Rdata"))
profile[profile==-1] = 0
profile2=profile[,colnames(profile)%in%annot2$GSM.identifier..sample.]

# take only high expressed into account
profile2[data.matrix(gexp)<5]=0

expressed_testis_num=colSums(profile2[rownames(profile2)%in%unique(t.df$gene),])
feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0_Antigens"
feat_class[expressed_testis_num>=1&expressed_testis_num<=4]="1to4_Antigens"
feat_class[expressed_testis_num>=5&expressed_testis_num<=6]="5to6_Antigens"
feat_class[expressed_testis_num>=7]="over7_Antigens"

GENESETS="Combined_pathway_signatures_2017_filtered_robust.gmt"
WD=file.path(getwd(), "GSEA")
OUTDIR=file.path(getwd(), "GSEA")

# # run GSEA:
# command=run.GSEA(data=gexp, cls.vector = as.numeric(expressed_testis_num), datatype = "N", GENESETS = GENESETS, dataname = "MM", clsname = "Testis_antigen_groups", WD=WD, OUTDIR=OUTDIR)
# try(system(command))
# 
# command=run.GSEA(data=gexp, cls.vector = as.numeric(annot2$HLAIScore), datatype = "N", GENESETS = GENESETS, dataname = paste0("MM_", GSE[1]), clsname = "HLAI", WD=WD, OUTDIR=OUTDIR)
# try(system(command))
# 
# command=run.GSEA(data=gexp, cls.vector = as.numeric(annot2$HLAIIScore), datatype = "N", GENESETS = GENESETS, dataname = paste0("MM_", GSE[1]), clsname = "HLAII", WD=WD, OUTDIR=OUTDIR)
# try(system(command))

a=read.delim(file.path(WD, "MM_Testis_antigen_groups_continuous_phenotype.Gsea.1549532122926/gsea_report_for_feat_neg_1549532122926.xls"), stringsAsFactors = F)
a=read.delim(file.path(WD, "MM_Testis_antigen_groups_continuous_phenotype.Gsea.1549532122926/gsea_report_for_feat_pos_1549532122926.xls"), stringsAsFactors = F)

# # GSE19784
# a=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_Testis_antigen_groups_continuous_phenotype.Gsea.1549542936287/gsea_report_for_feat_pos_1549542936287.xls", stringsAsFactors = F)
# b=read.delim("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/MM_Testis_antigen_groups_continuous_phenotype.Gsea.1549542936287/gsea_report_for_feat_neg_1549542936287.xls", stringsAsFactors = F)
 
c=rbind(a,b)

# get GSVA visualization for the pathways
# Geneset list
Onc.pathways=read.delim(GENESETS, stringsAsFactors = FALSE, header=F, col.names = paste("V",1:max(count.fields(GENESETS, sep = '\t'), na.rm = T)), fill = TRUE)

# Make list
listA=mclapply(1:length(Onc.pathways[,1]), function(i){A=as.character(Onc.pathways[i,3:length(Onc.pathways),])
B=A[!A==""&!A=="NA"]}, mc.cores=6)

names(listA) <- Onc.pathways[,1]

viz_scores=gsva(expr = data.matrix(gexp), gset.idx.list = listA, parallel.sz=8, tau=0.25)

# make a complex heatmap
significant=c(c[c$FWER.p.val<0.001,1])[1:10]
dat_plot=viz_scores[match(significant, toupper(rownames(viz_scores))),]

annof=data.frame("HLAI"=annot2$HLAIScore, "HLAII"=annot2$HLAIIScore)

df_anno=data.frame("cytogeneticAbnormalities"=as.double(annot2$MM_CYTOGENETIC_ABNORMALITIES),"age"=annot2$AGE, "iss"=as.character(annot2$MM_ISS), "iss1"=annot2$MM_ISS==1,"iss2"=annot2$MM_ISS==2,"iss3"=annot2$MM_ISS==3, "OS"=annot2$OS_Time, "PFS"=annot2$PFS_Time, "number.T.Antigenes"=as.numeric(expressed_testis_num), "antigen.group"=feat_class)

expressed_testis_num2=expressed_testis_num[order(expressed_testis_num)]

df_anno=df_anno[match(names(expressed_testis_num2), rownames(df_anno)),]
dat_plot=dat_plot[,match(names(expressed_testis_num2), colnames(dat_plot))]
annof=annof[match(names(expressed_testis_num2), colnames(gexp)),]
feat_class=feat_class[match(names(expressed_testis_num2), names(feat_class))]
gexp=gexp[,match(names(expressed_testis_num2), colnames(gexp))]

ha = HeatmapAnnotation(df=df_anno,NUM_TESTIS = anno_barplot(expressed_testis_num2, axis = T, gp = gpar(fill = "indianred")), HLA = anno_barplot(annof$HLAI, axis = T, gp = gpar(fill = "indianred")), HLAII = anno_barplot(annof$HLAII, axis = T, gp = gpar(fill = "darkgoldenrod"),ylim = c(5,10)))
ha2 = HeatmapAnnotation(df=data.frame(as.character(feat_class)), height = unit(20, "mm"))

pl=data.frame(df_anno$antigen.group, annof$HLAI)

logicalVectors=lapply(unique(feat_class), function(cl)feat_class%in%cl)
names(logicalVectors)=unique(feat_class)

rownames(dat_plot)=gsub("-.*.|_IN_CANCER_HOMO_SAPIENS-WIKIPW", "", rownames(dat_plot))

pdf("FigS6G_CGA_GSVA_GSE16716_GSE24080.pdf", height = 8, width = 8)
Heatmap(dat_plot, bottom_annotation = ha,use_raster = F, top_annotation=HeatmapAnnotation("# CGA"=anno_barplot(data.frame(expressed_testis_num2), height = unit(10, "mm"))), cluster_columns = F, cluster_rows = F, row_names_side = "right", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="GSVA\nScore", col=colorRamp2(c(-0.5,0,0.5), c("grey50", "white", "red")), show_column_names = F, width = unit(40, "mm"), height = unit(1.25*20, "mm"))
Heatmap(t(annof$HLAII), cluster_columns = F, cluster_rows = F, row_names_side = "right", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="HLAII\nScore", show_column_names = F, width = unit(40, "mm"), height = unit(1.25*1, "mm"))

a=t(scale(annof$HLAII))
a[a>2]=2
a[a<(-2)]=-2

Heatmap(a, cluster_columns = F, cluster_rows = F, row_names_side = "right", column_names_gp = gpar(fontsize = 8), row_names_gp = gpar(fontsize = 6), name="HLAII\nScore", show_column_names = F, width = unit(40, "mm"), height = unit(1.25*dim(a)[1], "mm"))
dev.off()


# verify HLA association

GSE=c("GSE16716,GSE24080")
# GSE=c("GSE19784")
# GSE=c("GSE15695,GSE21349", "GSE16716,GSE24080", "GSE17306", "GSE19784", "GSE2134")

annot2=annot[annot$GSE.identifier..experiment.%in%GSE,]

# gexp data
gexp=data[,colnames(data)%in%annot2$GSM.identifier..sample.]

t.df=t.df[order(t.df[,3]),]

# data:
genelist=t.df[t.df[,3]%in%"Cancer_Myeloma",1]
genelist=t.df[,1]

# rank patients by number of testis antigens expressed
# profiles
profile[profile==-1] = 0
profile2=profile[,colnames(profile)%in%annot2$GSM.identifier..sample.]

# take only high expressed into account
profile2[data.matrix(gexp)<5]=0

expressed_testis_num=colSums(profile2[rownames(profile2)%in%unique(t.df$gene),])
feat_class=expressed_testis_num
feat_class[expressed_testis_num==0]="0 CGA"
feat_class[expressed_testis_num>=1&expressed_testis_num<=4]="1-4 CGA"
feat_class[expressed_testis_num>=5&expressed_testis_num<=6]="5-6 CGA"
feat_class[expressed_testis_num>=7]=">7 CGA"

logicalVectors=lapply(unique(feat_class), function(cl)feat_class%in%cl)
names(logicalVectors)=unique(feat_class)

annof=data.frame("HLAI"=annot2$HLAIScore, "HLAII"=annot2$HLAIIScore)

genelist=c("HLAI", "HLAII")
p.all=lapply(genelist, plot.boxplot, logicalVectors = logicalVectors, data = t(annof), order = T,spread = F)

TestGeneWilcox("HLAII", data = t(annof), logicalVectors = logicalVectors, logicalVector_normals = logicalVectors, ALTERNATIVE = "less")
TestGeneWilcox("HLAII", data = t(annof), logicalVectors = logicalVectors, ALTERNATIVE = "less")

ggsave(plot = p.all[[2]], filename = "FigureS6G_HLAII_MM_hemap.pdf", width = unit(3.25, "cm"), height = unit(3, "cm"))
