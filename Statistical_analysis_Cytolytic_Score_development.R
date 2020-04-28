GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/useful_functions.R"))
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))

library(Hmisc)
library("caTools")
library(parallel)
require(Hmisc)
require(graphics)
library(corrplot)
library(ggplot2)
library(stats)
library(reshape2)
library(RColorBrewer)
library(ComplexHeatmap)

# data loading
setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# load annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))

# profiles
profile=data.matrix(get(load("mixtureM_profile.Rdata")))
profile[profile==-1] = 0
profile=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]

# gexp
data=t(get(load("data9544_with_gene_symbols.RData")))

# get same samples
data=data[, colnames(data)%in%annot$GSM.identifier..sample.]

data_sub=data[,!annot$immunoNormals%in%c("")]
profile_sub=profile[,!annot$immunoNormals%in%c("")]
annot_sub=annot[!annot$immunoNormals%in%c(""),]

iNorm=get.logical(annovector = list(annot_sub$immunoNormals))
iNorm=iNorm[!names(iNorm)%in%c("")]

CYTOLYTIC=as.logical(unlist(iNorm[names(iNorm)%in%"CD8+Tcell"])|unlist(iNorm[names(iNorm)%in%"NKCell"]))

# genes DE and enriched in CD8/NK cells
hgt.mw.cd8=FUN_HGT_MW_TEST(genes = rownames(profile_sub), logicalVector = CYTOLYTIC, name = "Cytotoxic", logicalVector_normals = iNorm[!names(iNorm)%in%c("CD8+Tcell", "NKCell")], profile = t(profile_sub), data=t(data_sub), TEST_FAILS = T, HG_PVAL = 0.00001, MW_PVAL = 0.01)
hgt.mw.cd8.filt=hgt.mw.cd8[hgt.mw.cd8$FoldChange>1.5&hgt.mw.cd8$adj.pvalue>5&!grepl("Bcell|HematopoieticStemCell|Erythroid|Macrophage|Monocyte|DendriticCell", hgt.mw.cd8$fails),]
genes=hgt.mw.cd8.filt[order(hgt.mw.cd8.filt$adj.pvalue, decreasing = T),1]

# make boxplot to verify:
iNorm=get.logical(annovector = list(annot_sub$immunoNormals))
iNorm=iNorm[!names(iNorm)%in%c("")]

GENELIST="CD8_NK_normals"

p.all=lapply(genes, plot.boxplot, iNorm, data=data_sub, order=T, N=T, spread=T)
ggsave(paste0(gsub(".txt", "", GENELIST), ".pdf"), do.call(marrangeGrob, append(list(grobs=p.all, nrow=3, ncol=1),list(top=NULL))), width = 200 , height = 300, units = "mm", dpi=250, limitsize = FALSE)

# plot=t(scale(t(data_sub[rownames(data_sub)%in%genes,])))
plot=data_sub[rownames(data_sub)%in%genes,]

plot2=aggregate_by_median(data = plot, xs=annot_sub$immunoNormals)

plot2=plot2[,!colnames(plot2)%in%c("TcellActivatedMononuclear", "MyeloidProgenitor2", "MyeloidProgenitor1", "Mononuclear", "LymphNode", "LangerhansCell", "GerminalCentreCell",  "GerminCentre")]

pdf("Fig1B.pdf", width = 5, height=10)
Heatmap(plot2, show_column_names = T)
Heatmap(t(scale(t(plot2))), show_column_names = T)
dev.off()

# combinations
combinations=genes

GENELIST="Compare_NK_CD8"
p.all=lapply(combinations, plot.boxplot, iNorm[c(5,13)], data=data_sub, order=T, N=T, spread=T)
ggsave(paste0(gsub(".txt", "", GENELIST), ".pdf"), do.call(marrangeGrob, append(list(grobs=p.all, nrow=3, ncol=1),list(top=NULL))), width = 200 , height = 300, units = "mm", dpi=250, limitsize = FALSE)

# choose cutoff, take sorted and unsorted for disease, you know the label, see what is the cutoff for nr2 component?

#**************** STEP1: first confirm that these genes have correct profile ******************
# test where difference of pure - unpure populations is highest
library(mclust)

#**********************************************************************************************
# test where difference of pure - unpure populations is highest
library(mclust)
TEST=function(gene, subgroup, annovector, filter=NULL, name_filt=NULL, set=set){
  
  if(any(subgroup%in%annovector)){
    filt=annovector%in%subgroup
  }else{
    filt=annovector%in%grep(eval(subgroup), annovector, value=T)
  }
  
  A=as.numeric(set[rownames(set)%in%gene,])
  
  if(is.null(filter)){
    name=paste(gene, unique(annovector[filt])[1], "\n", subgroup)
    model=A[filt]

  }else{
    name=paste(gene, name_filt, "\n", subgroup)
    model=A[filt&filter]
  }
  return(model)
}

FUN_CLASSIF=function(gene, A, B, model){
  correct=c(rep("A", length(A)), rep("B", length(B)))
  
  missA=table(correct[model$classification==1])
  missB=table(correct[model$classification==2])
  missB[2]=ifelse(is.na(missB[2]), length(B), missB[2])
  
  # how many were classified correctly
  valA=signif(missA[1]/length(A), 3)*100
  valB=signif(missB[2]/length(B), 3)*100
  valC=signif(sum(model$uncertainty>0.1)/length(model$uncertainty),2)*100
  
  NAME=paste0(gene, " correct class: ", valA, "% and ", valB, "%, uncertain samples: ",valC, "%")
  return(NAME)
}

FUN_gene=function(gene, set){
  
  # comparisons pure population vs. mixed and blast %
  has_blast=na.omit(colnames(blast_per)[as.numeric(blast_per)>0])
  blast_pure=na.omit(colnames(blast_per)[as.numeric(blast_per)>=95])
  blast_mix=na.omit(colnames(blast_per)[as.numeric(blast_per)<95&as.numeric(blast_per)>25])
  blast_high=na.omit(colnames(blast_per)[as.numeric(blast_per)>=75&as.numeric(blast_per)<95])
  blast_low=na.omit(colnames(blast_per)[as.numeric(blast_per)<=25])
  blast_medium=na.omit(colnames(blast_per)[as.numeric(blast_per)<75&as.numeric(blast_per)>25])
  
  pure=table(annot$Category.specifying.lineage.tumor.origin[annot$GSM.identifier..sample.%in%has_blast])
  
  # no expression, pure, sorted samples
  pure1=TEST(gene, "BCL", annot$Category.specifying.lineage.tumor.origin,name_filt="CD19", annot$CELLS_SORTED==1, set=set)
  pure2=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin, name_filt="CD34|blast", annot$CELLS_SORTED==1, set=set)
  pure3=TEST(gene, "CellLine", annot$colorClass, name_filt="CellLine", !grepl("TCL", annot$Category.specifying.lineage.tumor.origin), set=set)
  
  # pure but positive
  pure_pos=TEST(gene, "CD8+Tcell|CD8+TcellActivated|NaturalKillerCell", annot$Category.specifying.lineage.tumor.origin,annot$GSM.identifier..sample.%in%annot$GSM.identifier..sample., set=set)
  
  # mix - test
  mix1=TEST(gene, "BCL", annot$Category.specifying.lineage.tumor.origin,name_filt="na|none", annot$CELLS_SORTED==0, set=set)
  mix2=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="na|none", annot$CELLS_SORTED==0, set=set)
  
  # different blast percentages
  mix5=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST_>95%", annot$GSM.identifier..sample.%in%blast_pure, set=set)
  mix6=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST_>75%", annot$GSM.identifier..sample.%in%blast_high, set=set)
  mix7=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST_medium", annot$GSM.identifier..sample.%in%blast_medium, set=set)
  mix8=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST_<25%", annot$GSM.identifier..sample.%in%blast_low, set=set)
  mix9=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST_<80%", annot$GSM.identifier..sample.%in%blast_mix, set=set)
  mix10=TEST(gene, "AML", annot$Category.specifying.lineage.tumor.origin,name_filt="BLAST%", annot$GSM.identifier..sample.%in%has_blast, set=set)
  perc=blast_per[annot$GSM.identifier..sample.%in%has_blast&annot$Category.specifying.lineage.tumor.origin%in%"AML"]
  
  model_pure1=Mclust(c(pure1,pure_pos),G=2)
  model_pure2=Mclust(c(pure2,pure_pos),G=2)
  model_pure3=Mclust(c(pure3,pure_pos),G=2)
  
  model2=Mclust(mix1, G=2)
  model3=Mclust(mix2, G=2)
  model5=Mclust(mix5, G=2)
  model6=Mclust(mix6, G=2)
  model7=Mclust(mix7, G=2)
  model8=Mclust(mix8, G=2)
  model9=Mclust(mix9, G=2)
  
  #*************** Check model performance *******************
  NAME_BCL=FUN_CLASSIF(gene, pure1, pure_pos, model_pure1)
  NAME_AML=FUN_CLASSIF(gene, pure2, pure_pos, model_pure2)
  NAME_CLINE=FUN_CLASSIF(gene, pure3, pure_pos, model_pure3)
  #************************************************************

  A=set[gene,annot$GSM.identifier..sample.%in%has_blast&annot$Category.specifying.lineage.tumor.origin%in%"AML"]
  
  B=annot[annot$GSM.identifier..sample.%in%names(A)[A>8],]
  B$Sample.isolation
  
  A=set[gene,annot$colorClass%in%"CellLine"]
  B=annot[annot$GSM.identifier..sample.%in%names(A)[A>8],]
  
  plot(model_pure1, what = "density", xlab="pure_BCL_vs_CD8_NK",sub=NAME_BCL, xlim=c(2, 15))
  # plot(model, what = "uncertainty", xlab="pure_BCL_vs_CD8_NK",sub=NAME_BCL, xlim=c(2, 15))

  plot(model_pure2, what = "density", xlab="pure_AML_vs_CD8_NK",sub=NAME_AML, xlim=c(2, 15))
  # plot(modelB, what = "uncertainty", xlab="pure_AML_vs_CD8_NK",sub=NAME_AML, xlim=c(2, 15))
  
  plot(model_pure3, what = "density", xlab="pure_CellLine_vs_CD8_NK",sub=NAME_CLINE, xlim=c(2, 15))
  # plot(modelB, what = "uncertainty", xlab="pure_AML_vs_CD8_NK",sub=NAME_AML, xlim=c(2, 15))
  
  plot(model2, what = "density", xlab="mix_BCL",sub=gene, xlim=c(2, 15))
  # plot(model2, what = "uncertainty", xlab="mix_BCL",sub=gene, xlim=c(2, 15))

  plot(model3, what = "density", xlab="mix_AML",sub=gene, xlim=c(2, 15))
  # plot(model3, what = "uncertainty", xlab="mix_AML",sub=gene, xlim=c(2, 15))

  # blast percentage, should decrease cytolytic activity
  boxplot(list(perc[scale(mix10)>1.25], perc[scale(mix10)<1.25&scale(mix10)>-1.25], perc[scale(mix10)<(-1.25)]), names = c("high", "medium", "low"), main=paste("blast %", gene))
  
  # plot(model5, what = "density", xlab="mix_AML_blast_>95%",sub=gene, xlim=c(2, 15))
  # # plot(model5, what = "uncertainty", xlab="mix_AML_blast_>95%",sub=gene, xlim=c(2, 15))
  # 
  # plot(model6, what = "density", xlab="mix_AML_blast_>80%",sub=gene, xlim=c(2, 15))
  # # plot(model6, what = "uncertainty", xlab="mix_AML_blast_>80%",sub=gene, xlim=c(2, 15))
  # 
  # plot(model7, what = "density", xlab="mix_AML_blast_medium",sub=gene, xlim=c(2, 15))
  # # plot(model7, what = "uncertainty", xlab="mix_AML_blast_medium",sub=gene, xlim=c(2, 15))
  # 
  # plot(model8, what = "density", xlab="mix_AML_blast_<25%",sub=gene, xlim=c(2, 15))
  # # plot(model8, what = "uncertainty", xlab="mix_AML_blast_<25%",sub=gene, xlim=c(2, 15))
  # 
}

# blast percentage
blast_per=gsub(">=80%","", annot$Purity.Tumor.cell.content)
blast_per=gsub("blast%: |>=|%|blast cell percentage: |t_cell_purity: |>|;|blast count, % of sample, -1=unavailable : ","", blast_per)
blast_per[blast_per=="high"]=95
blast_per[blast_per=="-1"]=NA
blast_per[as.numeric(blast_per)>100]=100
blast_per[blast_per%in%c("n./a.", "na")]=NA
blast_per=t(data.matrix(as.numeric(blast_per)))
rownames(blast_per)="N:SAMP:BLAST_PERCENTAGE:::::"
colnames(blast_per)=colnames(data)

# Gaussian mixture model densities
pdf("cytolytic_genes_best_separation.pdf")
par(mfrow=c(3,2))
lapply(genes, FUN_gene, set=data)
dev.off()

# these were promising:
combinations=c("GZMA", "GZMB", "PRF1", "GNLY", "GZMH", "GZMM")

library(ggridges)
library(dplyr)
library(forcats)

# no expression, pure, sorted samples
pure=(annot$Category.specifying.lineage.tumor.origin%in%c("AML", "BCL"))&annot$CELLS_SORTED==1
pure=(annot$Category.specifying.lineage.tumor.origin%in%c("AML", "BCL")|annot$colorClass%in%c("CellLine")&!grepl("TCL|NKTCL", annot$Category.specifying.lineage.tumor.origin))&annot$CELLS_SORTED==1

# pure but positive
pure_pos=annot$immunoNormals%in%c("CD8+Tcell", "NKCell")

# mixture of normal cells and tumor
mix=(annot$Category.specifying.lineage.tumor.origin%in%c("AML", "BCL")|annot$colorClass%in%c("CellLine")&!grepl("TCL", annot$Category.specifying.lineage.tumor.origin))&annot$CELLS_SORTED==0

annov=rep("", length(annot$Category.specifying.lineage.tumor.origin))
annov[pure]="pure"
annov[pure_pos]="CD8+/NK"
annov[mix]="mix"

d=data.frame(t(data[rownames(data)%in%combinations,]), "annov"=annov, annot$CLASS2)
d=d[!annov=="",]

plot_gexp=melt(d)

p=ggplot(aes(y = plot_gexp$variable), data = plot_gexp) +
  geom_density_ridges(aes(x = value, fill = annov), 
                      alpha = .8, color = "white", from = 2, to = 15) +
  labs(x = "Gene expression",
       y = "CD8+/NK expressed gene",
       title = "pure vs. mixed population and CD8+ T / NK cells",
       subtitle = "Gene expression density") +
       scale_y_discrete(expand = c(0.01, 0)) +
         scale_x_continuous(expand = c(0.01, 0), breaks=c(2,4,6,8,10,12,14)) +
         # scale_fill_cyclical(breaks = c("1980 Indy", "1980 Unionist"),
         #                     labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
         #                     values = c("#ff0000", "#0000ff", "#ff8080", "#8080ff"),
         #                     name = "Sample population", guide = "legend") +
         theme_ridges(grid = FALSE)

ggsave(p, filename = "Fig1B_S1C.pdf")

# GZMB exoression in pure, excluded
combinations=c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM")
dat_a=data[rownames(data)%in%combinations,]
dat=2^dat_a+0.01
gm=log2(t(apply(dat, 2, gm_mean)))

# MCP counter data
MCP=get(load("MCP_counter_data.Rdata"))
rownames(MCP)=paste0("N:SAMP:DUFVA_", rownames(MCP), ":::::")
MCP=MCP[,colnames(MCP)%in%colnames(data)]

# cibersort
results=read.delim("CIBERSORT-Results.txt", row.names = 1, header=T, stringsAsFactors = F)
cibersort=t(results)
cibersort=cibersort[,colnames(cibersort)%in%colnames(data)]
rownames(cibersort)=paste0("N:SAMP:DUFVA_CIBERSORT_", gsub(" |-|\\.", "_", rownames(cibersort)), ":::::")
results_fm=data.matrix(cibersort)

# GSVA
library(GSVA)
load("dufva_bindea_2013_geneset_listA_tempfile.Rdata")
gsva_es=gsva(data, listA, tau=0.25, parallel.sz=4)
save(gsva_es, file="bindea_scores.Rdata")

# bindea sets, computed above
load("bindea_scores.Rdata")

# plot data
cor(MCP["N:SAMP:DUFVA_MCP_Cytotoxic_lymphocytes:::::",], as.numeric(gm), use="pairwise.complete.obs")

add2=colSums(cibersort[grepl("DUFVA_CIBERSORT_T_cells_CD8|DUFVA_CIBERSORT_NK_cells", rownames(cibersort)),])
cor(add2, as.numeric(gm), use="pairwise.complete.obs")

cor(gsva_es["DUFVA_BINDEA_CYTOTOXIC_CELLS",], as.numeric(gm), use="pairwise.complete.obs")


library(ggplot2)

dat=data.frame(t(data.matrix(rbind(gm,add2, MCP["N:SAMP:DUFVA_MCP_Cytotoxic_lymphocytes:::::",], gsva_es["DUFVA_BINDEA_CYTOTOXIC_CELLS",]))))
colnames(dat)=c("CytolyticScore", "CibersortCD8NK", "MCPCytotoxiclymphocytes", "Bindeacytotoxiccells")

library(ggpubr)

pdf("FigS1GtoI.pdf", height = 2.5, width = 2.5)
ggscatter(dat, x = "CytolyticScore", y = "CibersortCD8NK",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Cytolytic Score") +
  ylab("Cibersort CD8+NK") +
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(limits = c(2,15), breaks = c(2, 5, 8, 11, 14))

ggscatter(dat, x = "CytolyticScore", y = "Bindeacytotoxiccells",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Cytolytic Score") +
  ylab("Bindea cytotoxic cells") +
  scale_y_continuous(limits = c(-1,1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_x_continuous(limits = c(2,15), breaks = c(2, 5, 8, 11, 14))

ggscatter(dat, x = "CytolyticScore", y = "MCPCytotoxiclymphocytes",
          size = 1.5,
          add = "reg.line",  # Add regression line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) +
  stat_cor(method = "spearman", label.sep = "\n") +
  xlab("Cytolytic Score") +
  ylab("MCP Cytotoxic lymphocytes") +
  scale_y_continuous(limits = c(2,15), breaks = c(2, 5, 8, 11, 14)) +
  scale_x_continuous(limits = c(2,15), breaks = c(2, 5, 8, 11, 14))
dev.off()




#*************************** Step 2, validation using other tools ****************************

combinations=c("GZMA", "GZMB", "PRF1", "GNLY", "GZMH", "GZMM")

FUN=function(n){
  A=combs(combinations, n)
  
  res=t(apply(A, 1, function(v){
    dat_a=data[rownames(data)%in%v,]
    dat=2^dat_a+0.01
    gm=log2(t(apply(dat, 2, gm_mean)))
  }))
  
  rownames(res)=t(apply(A, 1, function(v)paste0(v, collapse="_")))
  return(res)
}

cyt_mat=do.call(rbind, lapply(2:length(combinations), FUN))
colnames(cyt_mat)=colnames(data)

A=cyt_mat[,annot$immunoNormals!=""]

annot_sub$immunoNormals[annot_sub$immunoNormals%in%c("CD8+Tcell", "NKCell")]="Cytolytic Cells"

plot2=t(scale(t(aggregate_by_median(data = A, xs=annot_sub$immunoNormals))))

Heatmap(plot2)

# MCP counter data
MCP=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/feat_matrices/data9544/MCP_counter_data.R"))
rownames(MCP)=paste0("N:SAMP:DUFVA_", rownames(MCP), ":::::")
MCP=MCP[,colnames(MCP)%in%colnames(data)]

# cibersort
results=read.delim("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data9544/CIBERSORT-Results.txt", row.names = 1, header=T, stringsAsFactors = F)
cibersort=t(results)
cibersort=cibersort[,colnames(cibersort)%in%colnames(data)]
rownames(cibersort)=paste0("N:SAMP:DUFVA_CIBERSORT_", gsub(" |-|\\.", "_", rownames(cibersort)), ":::::")
results_fm=data.matrix(cibersort)

# GSVA
gsva=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data8238/GSVA/data8238_all_samples_dufva_immunological_genes_updated_2016_GSVA.Rdata"))
bindea=get(load("/research/groups/sysgen/PROJECTS/HEMAP/HEMAP/dat2figs/data/data8238/GSVA/data8238_all_samples_dufva_bindea_2013_geneset_GSVA.Rdata"))
gsva_es=rbind(gsva, bindea)
gsva_es=gsva_es[,colnames(gsva_es)%in%colnames(data)]

# match cols gsva
gsva_es2=gsva_es[,match(colnames(data), colnames(gsva_es))]
colnames(gsva_es2)=colnames(data)

add=t(data.frame(colSums(results_fm[grepl("DUFVA_CIBERSORT_T_cells_CD8|DUFVA_CIBERSORT_NK_cells_activated", rownames(results_fm)),])))
rownames(add)="CIBERSORT_CYTOLYTIC_SUM"

add2=t(data.frame(colSums(results_fm[grepl("DUFVA_CIBERSORT_T_cells_CD8|DUFVA_CIBERSORT_NK_cells", rownames(results_fm)),])))
rownames(add2)="CIBERSORT_CYTOLYTIC_SUM2"

data_test_cytolytic=rbind(data[gsub("N:GEXP:|:::::", "", rownames(data))%in%combinations,], cyt_mat, results_fm[grep("CD8|NK", rownames(results_fm)),,drop=F], gsva_es2[grep("DUFVA", rownames(gsva_es2)),])
data_test_cytolytic=rbind(cyt_mat, results_fm[grep("CD8|NK", rownames(results_fm)),,drop=F],add,add2,MCP, gsva_es2[grepl("CD8|NK|CYTOTOXIC", rownames(gsva_es2))&grepl("DUFVA", rownames(gsva_es2)),])

# remove non-meaningful samples:
exclude=annot$colorClass%in%"CellLine"|annot$Main.category%in%c("NonCancerImmortalized", "B-Lymphoid", "T-Lymphoid", "Erythroid", "Myeloid", " T-Lymphoid")
include=annot$Main.category%in%c("AML")

#coeff. + p-value matrix all reg. elements (from analyzeRepeats.pl)
rcorr_matrix_all = rcorr(t(as.matrix(data_test_cytolytic))) #row against row
rcorr_all_p = round(data.matrix(rcorr_matrix_all$P), digits=3)
rcorr_all_coef = round(data.matrix(rcorr_matrix_all$r), digits=3)

A=melt(rcorr_all_coef)
A=melt(rcorr_all_p)

#******************** Check score correlations as different combinations, choose the ones with overal good correlation ************************
# Therefore, we have checked that: 
# 1) genes and the cytolytic score work well in hematological context
# 2) genes are in agreement with published tools
# c("GZMA", "PRF1", "GNLY", "GZMH", "GZMM") were chosen in the end. GZMB excluded at it was expressed in pure cancer cells (mixture model).

find=c("N:SAMP:DUFVA_BINDEA_CYTOTOXIC_CELLS_GSVA:::::", "N:SAMP:DUFVA_MCP_Cytotoxic_lymphocytes:::::")
check=c("N:SAMP:DUFVA_BINDEA_CYTOTOXIC_CELLS_GSVA:::::", "N:SAMP:DUFVA_MCP_Cytotoxic_lymphocytes:::::", as.character(unique(rownames(cyt_mat))))
val=sapply(1:length(check), function(i){
  mean(rcorr_all_coef[rownames(rcorr_all_coef)%in%check[i],colnames(rcorr_all_coef)%in%find[!find%in%check[i]]])})

res=cbind(check, val)
res=res[order(res[,2], decreasing = T),]
        
#**********************************************************************************************************************************************        
        
       


# visualizations
#######################################################################

coef=rcorr_all_coef
pval=rcorr_all_p

# cluster data
data2 <- scale(t(coef)) # Z-score
data2 <- t(coef) # Z-score

ord <- hclust( dist(data2, method = "euclidean"), method = "ward.D")
ord=ord$labels[ord$order]

ord_final=match(ord, rownames(data2))

#petri's plot & stars
plot.data <- melt(coef, id="type")  # Convert to long
plot.data$stars <- cut(pval, breaks=c(-Inf, 0.001, 0.01, 0.05, 1), label=c("***", "**", "*", ""))  # Create column of significance labels
plot.data$p.value <- melt(pval, id="adj. P-value")  # Calculate p-values for t-values
plot.data$variable <- factor(plot.data$Var2, levels = rownames(data2)[ord_final])
plot.data$type <- factor(plot.data$Var1, levels = rownames(data2)[ord_final])

plot.data=plot.data[abs(plot.data$value)>0.5,]

#plot everything
myPalette <- colorRampPalette(rev(brewer.pal(8, "RdBu")))
sc <-scale_fill_gradientn(colours=myPalette(8), limits=c(-1,1))
p <- ggplot(aes(x=type, y=variable, fill=value), data=plot.data)
fig12 <- p + geom_tile() + sc +
  # geom_text(aes(label=stars), color="black", size=5) +
  labs(y=NULL, x=NULL, fill="Correlation Coefficient") +
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
fig12

plot.data=plot.data[!plot.data$value==1,]

# A=plot.data[grepl("DUFVA_BINDEA_CD8_T_CELLS_GSVA", plot.data[,1]),]
# A2=A[order(A[,3], decreasing = T),]
# A2=A2[1:50,]

# A=plot.data[grepl("DUFVA_BINDEA_NK_CELLS_GSVA", plot.data[,1]),]
# A3=A[order(A[,3], decreasing = T),]
# A3=A3[1:50,]

# A=plot.data[grepl("DUFVA_CIBERSORT_T_cells_CD8", plot.data[,1]),]
# A4=A[order(A[,3], decreasing = T),]
# A4=A4[1:50,]

A=plot.data[grepl("DUFVA_NK_CD8_COMBINED_GSVA", plot.data[,1]),]
A5=A[order(A[,3], decreasing = T),]
A5=A5[1:50,]

# A=plot.data[grepl("DUFVA_CIBERSORT_NK_cells_activated", plot.data[,1]),]
# A6=A[order(A[,3], decreasing = T),]
# A6=A6[1:50,]

A=plot.data[grepl("DUFVA_BINDEA_CYTOTOXIC_CELLS_GSVA", plot.data[,1]),]
A7=A[order(A[,3], decreasing = T),]
A7=A7[1:50,]

A=plot.data[plot.data[,1]%in%"CIBERSORT_CYTOLYTIC_SUM",]
A8=A[order(A[,3], decreasing = T),]
A8=A8[1:50,]

# A=plot.data[grepl("CIBERSORT_CYTOLYTIC_SUM2", plot.data[,1]),]
# A9=A[order(A[,3], decreasing = T),]
# A9=A9[1:10,]

A=plot.data[grepl("MCP_Cytotoxic_lymphocytes", plot.data[,1]),]
A10=A[order(A[,3], decreasing = T),]
A10=A10[1:50,]

# dd=rbind(A2, A3, A4, A5, A6, A7, A8, A9)
dd=rbind(A7, A10, A8)

cor_feats=plot.data[!plot.data$Var1%in%rownames(cyt_mat)&!plot.data$Var2%in%rownames(cyt_mat),]
cor_feats=cor_feats[order(cor_feats[,3], decreasing = T),]
test2=sort(table(cor_feats$Var1))
cor_feats=cor_feats[cor_feats$Var2%in%names(test2)[test2>=4],]
cor_feats$variable <- factor(cor_feats$Var2, levels = unique(cor_feats$variable[order(cor_feats$value, decreasing = T)]))


dd=dd[dd$Var2%in%rownames(cyt_mat),]

dd=dd[order(dd$value),]

dd$variable <- factor(dd$Var2, levels = unique(dd$variable[order(dd$value, decreasing = T)]))

A=plot.data[grepl("CIBERSORT_CYTOLYTIC_SUM2|DUFVA_BINDEA_CYTOTOXIC_CELLS_GSVA|MCP_Cytotoxic_lymphocytes", plot.data[,1]),]
A=A[order(A[,3], decreasing = T),]
ca=A[A[,2]%in%c("GZMA_PRF1", "GZMA_PRF1_GNLY_GZMH_GZMM", "MCP_Cytotoxic_lymphocytes"),1:4]
# this could be used as a table

test=sort(table(dd$Var2))
dd2=dd[dd$Var2%in%names(test)[test>=2],]

dd2$variable <- factor(dd2$Var2, levels = unique(dd2$variable[order(dd$value, decreasing = T)]))

#plot everything
myPalette <- colorRampPalette(rev(brewer.pal(8, "RdBu")))
sc <-scale_fill_gradientn(colours=myPalette(8), limits=c(-1,1))
p <- ggplot(aes(x=type, y=variable, fill=value), data=dd)
fig12 <- p + geom_tile() + sc +
  # geom_text(aes(label=stars), color="black", size=5) +
  labs(y=NULL, x=NULL, fill="Correlation Coefficient") +
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
fig12

#plot everything
myPalette <- colorRampPalette(brewer.pal(8, "Reds"))
sc <-scale_fill_gradientn(colours=myPalette(8), limits=c(0.5,1))
p <- ggplot(aes(x=type, y=variable, fill=value), data=dd2)
fig12 <- p + geom_tile() + sc +
  # geom_text(aes(label=stars), color="black", size=5) +
  labs(y=NULL, x=NULL, fill="Correlation Coefficient") +
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
fig12

pdf("cytolytic_correlations.pdf", width = 20, height = 15)
fig12
dev.off()

