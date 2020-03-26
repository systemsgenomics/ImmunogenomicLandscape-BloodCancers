GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "statistics/useful_functions.R"))
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))

library(parallel)
library(ComplexHeatmap)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# annotations
annot=get(load("Hemap_immunology_Annotations.Rdata"))

load("data9544_with_gene_symbols.RData")
data_org = data[rownames(data)%in%annot$GSM.identifier..sample.,]

profile=data.matrix(get(load("mixtureM_profile.Rdata")))
profile[profile==-1] = 0
profile=profile[,colnames(profile)%in%annot$GSM.identifier..sample.]

nonAPC=c("CD4+Tcell","CD8+Tcell","Tcell","Erythroid","NKCell","Neutrophil","RegulatoryTcell","Gamma-delta-Tcell")
APC=c("DendriticCell", "Bcell", "MemoryBcell", "NaiveBcell", "AlveolarMacrophage", "M1-Macrophage","M2-Macrophage", "DendriticCell")
genelist=rownames(profile)

lv=annot$immunoNormals%in%nonAPC|annot$immunoNormals%in%APC

data_sub=data_org[lv,]
annot_sub=annot[lv,]
profile_sub=profile[,lv]

# get annotations for normal cells non-APC
normals=annot_sub$immunoNormals%in%nonAPC
logicalVector_normals=get.logical(annovector = list(annot_sub$immunoNormals), filterv = normals)
normals_comb=Reduce("|", logicalVector_normals)

# APC cells
normals2=annot_sub$immunoNormals%in%APC
logicalVector=get.logical(annovector = list(annot_sub$immunoNormals), filterv = normals2)
logicalVector_APC=list(normals2)
names(logicalVector_APC)="APC"
lv_comb=Reduce("|", logicalVector)

# genelist, all genes
genelist=colnames(data_org)

#******************************** HLAII candidates APC vs non APC *******************************************

test_APC_nonAPC=do.call(rbind, mclapply(genelist, function(g){
  res=TestGeneWilcox(g, t(data_sub), logicalVectors = logicalVector_APC, logicalVector_normals = list("nonAPC"=normals_comb), ALTERNATIVE = "greater", prettynum = F)
  return(res)
}, mc.cores=8))
test_APC_nonAPC$adj.P.value=p.adjust(test_APC_nonAPC$p, method = "bonferroni")

save(test_APC_nonAPC, file="HLAIIScore_APC_nonAPC.Rdata")

#************************************************************************************************************

load("HLAIIScore_APC_nonAPC.Rdata")
test_APC_nonAPC_filt=test_APC_nonAPC[test_APC_nonAPC$adj.P.value<1e-10&test_APC_nonAPC$FC>4,]

test=do.call(rbind, mclapply(unique(test_APC_nonAPC_filt[,1]), function(g){
  res=TestGeneWilcox(g, t(data_sub), logicalVectors = logicalVector, logicalVector_normals = logicalVector_normals, ALTERNATIVE = "greater", prettynum = F)
  plot2=aggregate_by_median(data = t(res$FC), xs=res$Group1)
  rownames(plot2)=g
  return(plot2)
  }, mc.cores=10))

save(test, file="HLAIIScore_APC_nonAPC_Allgroups.Rdata")

# high difference for all cell types, FC over 2 for all celltypes
B=names(test[rowSums(test>2)==7,1])

tab=rbind(test_APC_nonAPC_filt[test_APC_nonAPC_filt[,1]%in%B,], test_APC_nonAPC_filt[!test_APC_nonAPC_filt[,1]%in%B,])
tab[,2]=prettyNum(tab[,2])
tab[,3]=prettyNum(tab[,3])
tab[,6]=prettyNum(tab[,6])

write.table(tab[tab[,1]%in%B,], "TableS4_HLAIIScore_APC_genes.txt", col.names=T,row.names=F, quote = F, sep="\t")

# test_APC_nonAPC_filt=fread("TableS4_HLAIIScore_APC_genes.txt",data.table=F)

dat=data_org[normals,colnames(data)%in%B]
dat2=data_org[normals2,colnames(data)%in%B]
dat3=data_org[,colnames(data_org)%in%B]

# cormat=cor(dat, method="spearman")
# cormat2=cor(dat2)
cormat3=cor(dat3)

# to save space final figure was to show only cluster on top:
pdf("FigureS4B_HLAII_corHeatmap2.pdf", width = 6.5, height = 6)
Heatmap(cormat3, column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6))
dev.off()

plot=t(data_sub[,colnames(data_sub)%in%B])
plot2=aggregate_by_median(data = plot, xs=annot_sub$immunoNormals)

pdf("FigureS4A_HLAII_Significant_genes.pdf", width = 3, height = 6)
Heatmap(plot2, column_names_gp = gpar(fontsize = 6), row_names_gp = gpar(fontsize = 6))
dev.off()

# This is what we used in the end (HLAII score): 
c("HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DRB1")
