library(Seurat)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)


GIT_HOME="/research/users/ppolonen/git_home/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/pathway_analysis/functions.GSEA.R"))

# plot CLL genes:
load("CLL_D0_scRNA.Rdata")

load("Hemap_microenvironment_summary_statistics.Rdata")
load("Hemap_cytolytic_correlated_genes_TableS2_onlysignif.Rdata")

dat=res_all.filt[res_all.filt$disease%in%"CLL",]
dat=dat[order(dat$category,-as.numeric(dat$adj.pval), decreasing = T),]

# only take genes expressed in cluster in CLL:
dat=dat[dat$gene%in%markers.all$gene[markers.all$p_val_adj<0.001],]

# take max20:
dat=do.call(rbind, lapply(unique(dat$category), function(v)head(dat[dat$category%in%v,], 20)))

DE.genes=dat$gene[dat$significant]

scmat.filt=scmat[,scmat[["SingleR.label"]][,1]%in%c("Memory B-cells", "Monocytes", "naive B-cells", "NK cells", "Class-switched memory B-cells", "CD8+ Tem", "CD8+ Tcm")]

scmat.filt[["SingleR.label"]][,1][scmat.filt[["SingleR.label"]][,1]%in%"Class-switched memory B-cells"]="Memory B-cells"

scmat.filt[["SingleR.label"]]=factor(scmat.filt[["SingleR.label"]][,1], levels=c( "Monocytes", "naive B-cells", "Memory B-cells", "CD8+ Tcm","CD8+ Tem", "NK cells"))

cor.hm=as.numeric(dat$Rho[match(DE.genes, dat$gene)])

pdf("FigureS2D.pdf", height = 8, width = 4.25)
plot.DotPlot(scmat.filt, features = DE.genes[DE.genes%in%rownames(scmat)], cols = c("white", "red"), group.by = "SingleR.label", dot.scale = 4)
Heatmap(cor.hm, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#1e08ff", "#764bfd", "white", "#f66b4b", "#f4060d")))
dev.off()


load("FIMM_AML_scRNA.Rdata")

dat=res_all.filt[res_all.filt$disease%in%"AML",]
dat=dat[order(dat$category,-as.numeric(dat$adj.pval), decreasing = T),]

# only take genes expressed in cluster in CLL:
dat=dat[dat$gene%in%markers.all$gene[markers.all$p_val_adj<0.001],]

# take max20:
dat=do.call(rbind, lapply(unique(dat$category), function(v)head(dat[dat$category%in%v,], 20)))

DE.genes=dat$gene[dat$significant]

a=table(scmat[["SingleR.label"]][,1])

scmat.filt=scmat[,scmat[["SingleR.label"]][,1]%in%names(a)[a>500]]

DE.genes=unique(c(DE.genes))

cor.hm=as.numeric(dat$Rho[match(DE.genes, dat$gene)])

pdf("FigureS2E.pdf", height = 5.5, width = 5)
plot.DotPlot(scmat.filt, features = DE.genes[DE.genes%in%rownames(scmat)], cols = c("white", "red"), group.by = "SingleR.label", dot.scale = 4)
Heatmap(cor.hm, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#1e08ff", "#764bfd", "white", "#f66b4b", "#f4060d")))
dev.off()
