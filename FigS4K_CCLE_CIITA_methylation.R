GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))

library(limma)
library(edgeR)
library(parallel)
library(gridExtra)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# CCLE gexp data:
ccle=get(load("CCLE_RNASEQ_SYMBOL_COUNTS.Rdata"))
dge=DGEList(ccle)
keep.exprs <- filterByExpr(dge)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method="TMM")
ccle=voom(dge, plot=T)$E

# CCLE meth data:
ccle.meth=get(load("CCLE_combined_methylation.Rdata"))

# CCLE annot:
annot=read.delim("CCLE_sample_info_file_2012-10-18.txt", header=T, stringsAsFactors = F)

find=intersect(colnames(ccle.meth), intersect(annot$CCLE.name, colnames(ccle)))

ccle.meth=ccle.meth[,match(find, colnames(ccle.meth))]
ccle=ccle[,match(find, colnames(ccle))]
annot=annot[match(find, annot$CCLE.name),]

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

genelist=c("HLAII_SCORE", c("HLA-DMA","HLA-DMB","HLA-DPB1","HLA-DRA","HLA-DRB1"))
p.all=lapply(genelist, plot.boxplot, logicalVectors = logicalVectors,data = d,order = T,spread = T, sample.annotation = annotv, outlier.size=1.5)

ggsave("FigS4K_CCLE_CIITA_methylation.pdf", do.call(marrangeGrob, append(list(grobs=p.all[1], nrow=3, ncol=2),list(top=NULL))), width = 205 , height = 250, units = "mm", dpi=150, limitsize = FALSE)
