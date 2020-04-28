GIT_HOME="/research/users/ppolonen/git_home/ImmunogenomicLandscape-BloodCancers/"
source(file.path(GIT_HOME, "common_scripts/featurematrix/functions_generate_fm.R"))
source(file.path(GIT_HOME, "common_scripts/visualisation/plotting_functions.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/functions_statistics.R"))
source(file.path(GIT_HOME, "common_scripts/statistics/statistics_wrappers.R"))
library(parallel)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(gridExtra)
library(survcomp)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# annotations
annot=get(load("Hemap_immunology_Annotations.Rdata"))

# data
load("data9544_with_gene_symbols.RData")
data = t(data[rownames(data)%in%annot$GSM.identifier..sample.,])
profile=data.matrix(get(load("mixtureM_profile.Rdata")))
profile[profile==-1] = 0
profile=profile[,colnames(profile)%in%colnames(data)]

# geometric mean:
gm=annot$CytolyticScore

# logical vectors for all diseases:
subclass2=get.logical(annovector = list(annot$subclasses), filterv = annot$CELLS_SORTED==0&annot$Sample.type%in%c("Prolif", "Cancer"))
tbly2=get.logical(annovector = list(annot$tbLY), filterv = annot$CELLS_SORTED==0&annot$Sample.type%in%c("Prolif", "Cancer"))
tbly2=tbly2[!names(tbly2)%in%"Lymphoma_BCL"]

list_cancers=lapply(c(subclass2, tbly2), list)
list_cancers=unlist(list_cancers, recursive=F)

l=c("BCL_DLBCL", "Lymphoma_BCL_CHL", "Lymphoma_BCL_MCL", "Lymphoma_BCL_FL","CLL", "CML", "AML","Prolif_Myeloproliferative_MDS", "pre-B-ALL",  "T-ALL")
logicalvectors=list_cancers[match(l, names(list_cancers))]
names(logicalvectors)=gsub("^BCL_|Lymphoma_BCL_|Prolif_Myeloproliferative_", "", names(logicalvectors))

#********************************** compute summary statistics: **********************************

# Then compute cor, cor.test and TNK FC for same cancers:
res=do.call(cbind, mclapply(seq(logicalvectors), function(i){
  logv=logicalvectors[[i]]
  NAME=names(logicalvectors)[i]
  cors=cor(t(data[,logv]), as.numeric(gm)[logv], method = "spearman", use = "complete.obs")
}, mc.cores=10))

NK_log=grepl("CD8|NaturalKiller", annot[,5])

fc_res=do.call(cbind, mclapply(seq(logicalvectors), function(i){
  logv=logicalvectors[[i]]
  NAME=names(logicalvectors)[i]

  FC=apply(data[,], 1, function(v){
    mean(v[NK_log])-mean(v[logv])
  })

}, mc.cores=10))

res_pval=do.call(cbind, mclapply(seq(logicalvectors), function(i){
  logv=logicalvectors[[i]]
  NAME=names(logicalvectors)[i]
  pvals=-log10(apply(data[,logv,drop=F], 1, cor.test.p, as.numeric(gm)[logv]))
  pvals[!is.finite(pvals)]=224
  return(pvals)
}, mc.cores=10))

# these matrices can be used for plotting:
colnames(res)=names(logicalvectors)
colnames(fc_res)=names(logicalvectors)
colnames(res_pval)=names(logicalvectors)

# sorted samples, see if genes are expressed in cancer
subclass2=get.logical(annovector = list(annot$subclasses), filterv = annot$CELLS_SORTED==1&annot$Sample.type%in%c("Prolif", "Cancer"))

# diseases used minus FL with 7 samples
l=c("AML", "BCL", "CLL", "CML", "pre-B-ALL", "Prolif_Myeloproliferative_MDS", "T-ALL")
logicalvectors.sorted=subclass2[names(subclass2)%in%l]
names(logicalvectors.sorted)=gsub("^BCL_|Lymphoma_BCL_|Prolif_Myeloproliferative_", "", names(logicalvectors.sorted))

pval.sorted=fisher.matrix.lv(logicalVector.list = logicalvectors.sorted, data = profile, to.log10 = T, adjust.method = "BH", fisher.alternative = "greater")
pval.sorted=pval.sorted[match(rownames(res), rownames(pval.sorted)),]

# All tests
save(list = c("pval.sorted", "fc_res", "res", "res_pval", "pval.sorted"), file = "Hemap_microenvironment_summary_statistics.Rdata")

#*****************************************************************************************************

# load from above R object:
load("Hemap_microenvironment_summary_statistics.Rdata")

res_all=do.call(rbind, mclapply(seq(logicalvectors), function(i){
  logv=logicalvectors[[i]]
  NAME=names(logicalvectors)[i]
  
  # correlate genes with cytolytic activity
  pvals=apply(data[,logv,drop=F], 1, cor.test.p, as.numeric(gm)[logv])
  FC=fc_res[match(names(pvals), rownames(fc_res)),NAME]
  cors=res[match(names(pvals), rownames(res)),NAME]
  
  mean.exp=apply(data[,logv,drop=F], 1, median)
  
  data_cor=data.frame("gene"=names(cors), "Rho"=as.numeric(cors), "pval"=pvals, "adj.pval"=p.adjust(pvals), "FCtoTNK"=as.numeric(FC), "disease"=NAME, "mean.expression"=mean.exp, stringsAsFactors = F)
  data_cor=data_cor[order(data_cor[,2], decreasing = T),]
  
  if(!dim(data_cor)[1])return(NULL)
  
  # significant
  data_cor$significant=data_cor$adj.pval<0.01&abs(data_cor$FCtoTNK)>2&abs(data_cor$Rho)>0.2
  
  # add category:
  data_cor$category=""
  data_cor$category[data_cor$FCtoTNK>0&data_cor$Rho>0&data_cor$significant]="CTL/NK gene"
  data_cor$category[data_cor$FCtoTNK<0&data_cor$Rho>0&data_cor$significant]="Stromal/cancer gene (Rho > 0)"
  data_cor$category[data_cor$FCtoTNK<0&data_cor$Rho<0&data_cor$significant]="Stromal/cancer gene (Rho < 0)"
  
  return(data_cor)
}, mc.cores=10))

res_all.filt=res_all[res_all$significant,]

res_all.filt=res_all.filt[order(res_all.filt$category, abs(res_all.filt$Rho), abs(res_all.filt$FCtoTNK), -res_all.filt$adj.pval, decreasing = T),]

write.table(res_all[res_all$gene%in%res_all.filt$gene,], "Hemap_cytolytic_correlated_genes_TableS2.txt", col.names=T,row.names=F, quote = F, sep="\t")

for(i in 2:5)res_all.filt[,i]=prettyNum(signif(res_all.filt[,i], 2))

write.table(res_all.filt, "TableS2.txt", col.names=T,row.names=F, quote = F, sep="\t")
save(res_all.filt, file="Hemap_cytolytic_correlated_genes_TableS2_onlysignif.Rdata")

#***************************************************************************************************************************
fun.plot=function(dat, main, xlab, ylab, aval=2, bval=0.2,cval=0.01,genelist=NULL, MAX=25, SIZE=3, T.SIZE=8, height=7.5, width=6){

  if(!is.null(genelist)){
    dat$significant=dat$gene%in%genelist
  }else{
    dat$significant=dat$adj.pval<cval&abs(dat$Rho)>bval&abs(dat$FCtoTNK)>aval
  }
  
  groups=unique(dat$category)
  dat$plot=F
  
  for(g in groups){
    a2=dat$FCtoTNK[dat$category%in%g]
    b2=dat$Rho[dat$category%in%g]
    gene=dat$gene[dat$category%in%g]
    test=tail(gene[order(abs(a2^2*b2^2))], MAX)
    dat$plot[dat$gene%in%test&dat$category%in%g&dat$significant]=T
  }
  
  dat=dat[order(dat$plot),]
  
  dat$size.point=ifelse(dat$plot, 2.4, 0.6)
  
  myColors <- data.frame(c("grey85","#E41A1C", "grey85","#377EB8"), c("", "CTL/NK gene", "Stromal/cancer gene (Rho < 0)", "Stromal/cancer gene (Rho > 0)"), stringsAsFactors = F)
  myColors=myColors[match(unique(dat$category), myColors[,2]),1]
  names(myColors) <- unique(dat$category)
  
  colScale <- scale_colour_manual(name = "category", values = myColors)
  
  library(ggrastr)
  
  p=ggplot(dat, aes(FCtoTNK, Rho, colour = category)) +
    geom_point_rast(aes(size = dat$size.point), raster.width = width, raster.height = height) +
    scale_size_continuous(range = c(0.6, 2.5))+
    colScale +
    theme_classic(base_size = T.SIZE) +
    ggtitle(main) +
    labs(x = xlab, y = ylab, fill = "") +
    scale_y_continuous(breaks = pretty(dat$Rho, n = 6))+
    
    geom_text_repel(
      data = subset(dat, dat$plot),
      aes(label = subset(dat, dat$plot)$gene),
      size = SIZE,
      color="black",
      force=1.5,
      segment.size = 0.1,
      box.padding = unit(0.01, "lines"),
      point.padding = unit(0.01, "lines"),
      family="Helvetica"
    )+
    theme(legend.position = "none", legend.direction = "horizontal",
          axis.line = element_line(size=0.5, colour = "black"),
          panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(),
          plot.title = element_text(size = T.SIZE, face = "bold", family="Helvetica"),
          text=element_text()
    )+
    theme(axis.text.x = element_text(colour="grey20",size=T.SIZE,face="plain", family="Helvetica"),
          axis.text.y = element_text(colour="grey20",size=T.SIZE,face="plain", family="Helvetica"),  
          axis.title.x = element_text(colour="grey20",size=T.SIZE,face="plain", family="Helvetica"),
          axis.title.y = element_text(colour="grey20",size=T.SIZE,face="plain", family="Helvetica"))  
  return(p)
}


# Figure2A:
res_comb=d13=rowMeans(res[,c("DLBCL", "MCL", "FL", "CHL")])
fc_res_comb=rowMeans(fc_res[,c("DLBCL", "MCL", "FL", "CHL")])

pval_comb=apply(10^-res_pval[,c("DLBCL", "MCL", "FL", "CHL")], 1, combine.test, method = "z.transform")

data_comb=data.frame("gene"=names(res_comb), "Rho"=as.numeric(res_comb), "adj.pval"=pval_comb, "FCtoTNK"=as.numeric(fc_res_comb), "disease"="BCL", stringsAsFactors = F)

# significant
data_comb$significant=data_comb$adj.pval<0.01&abs(data_comb$FCtoTNK)>1.5&abs(data_comb$Rho)>0.2

# add category:
data_comb$category=""
data_comb$category[data_comb$FCtoTNK>0&data_comb$Rho>0&data_comb$significant]="CTL/NK gene"
data_comb$category[data_comb$FCtoTNK<0&data_comb$Rho>0&data_comb$significant]="Stromal/cancer gene (Rho > 0)"
data_comb$category[data_comb$FCtoTNK<0&data_comb$Rho<0&data_comb$significant]="Stromal/cancer gene (Rho < 0)"

# these were picked from above analysis, as they were linked to immunological regulation and myeloid infiltrate
monoc=c(c("CD14", "LYZ", "S100A8", "CD68", "CD163","MS4A4A", "MS4A6A","TLR8","LILRB2","IDO1","IFI27", "CXCL9", "CXCL10", "CXCL11", "CCL8", "CCL18", "CCL19", "TNFSF13B","CXCL13", "C1QA", "C1QB", "C1QC", "HBB"), c("IFNG", "CCL5","TNFSF10", "CD226", "CD160","KLRF1","KLRG1", "KLRC3", "EOMES", "GZMA", "PRF1", "CD8A"))

lymphomagenes=unique(res_all.filt$gene[res_all.filt$disease%in%c("DLBCL", "MCL", "FL", "CHL")])

NAME="Figure2A.pdf"
dat=fun.plot(dat = data_comb,
             main="B cell lymphomas",
             genelist=unique(monoc[monoc%in%lymphomagenes]),
             # genelist2=unique(lymphomagenes),
             aval=1.5, 
             bval=0.35,
             cval = 0.01,
             MAX = 60,
             SIZE=2.5,
             T.SIZE = 10,
             ylab="Average correlation with cytolytic score",
             xlab= paste0("Average expression fold change\n(CTL/NK cells vs. ", "BCL", ")"))

ggsave(filename = NAME, plot = dat, width = 85, height = 120, units = "mm", dpi=150)


# FigureS2:
NAME="FigureS2C.pdf"
p.all=lapply(seq(logicalvectors)[-4], function(i){
  fun.plot(dat = res_all[res_all$disease%in%names(logicalvectors)[i],],
           main=names(logicalvectors)[i],
           aval=2, 
           bval=0.2,
           cval = 0.01,
           MAX = 20,
           SIZE=1.75,
           T.SIZE = 8,
           ylab="correlation with cytolytic activity",
           xlab= paste0("expression fold change\n(CTL/NK cells vs. ", names(logicalvectors)[i], ")"))
})
ggsave(NAME, ggpubr::ggarrange(plotlist = p.all, ncol = 3, nrow = 3),  width = 152, height = 200, units = "mm", dpi=300)

# plot dotplot:
data.plot=data.frame(res_all$gene, res_all$disease, res_all$Rho, abs(res_all$FCtoTNK), "category"=res_all$category)
data.plot=data.plot[data.plot$res_all.gene%in%res_all.filt$gene,]
colnames(data.plot)=c("features", "id", "variable.1", "variable.2", "category")
data.plot$id=factor(data.plot$id,  levels=c("DLBCL","CHL", "FL", "MCL", "CLL", "CML", "AML", "MDS", "pre-B-ALL", "T-ALL"))

#************************************ Figure 2B: ***********************************

TNK.genes=rownames(pval.sorted)[rowSums(pval.sorted>2)==0&rownames(pval.sorted)%in%res_all.filt$gene[res_all.filt$category=="CTL/NK gene"]]

monoc=c(c("IFNG", "CCL5","TNFSF10", "CD226", "CD160","KLRF1","KLRG1", "KLRC3", "EOMES"), c("CD14", "LYZ", "S100A8", "CD68", "CD163","MS4A4A", "MS4A6A","TLR8","LILRB2","IDO1","IFI27", "CXCL9", "CXCL10", "CXCL11", "CCL8", "CCL18", "CCL19", "TNFSF13B","CXCL13", "C1QA", "C1QB", "C1QC", "HBB"))

d=data.plot[data.plot$features%in%monoc,]
d=d[order(match(d$features, monoc)),]

d$features=factor(as.character(d$features), levels = monoc)

# pdf("Hemap_microenvironment_Figure2_stroma.pdf", height = 4, width = 3.25)
# plot.DotPlot.df(data.plot = d, name.variable.1 = "Rho", name.variable.2 = "FC", cols = c("blue", "white", "red"), dot.scale = 3, col.min = -0.5, col.max = 0.75)
# dev.off()
# 
# pdf("Hemap_microenvironment_TNK_genes.pdf", height = 7, width = 3.25)
# plot.DotPlot.df(data.plot = data.plot[data.plot$features%in%TNK.genes,], name.variable.1 = "Rho", name.variable.2 = "FC", cols = c("blue", "white", "red"), dot.scale = 3, col.min = -0.5, col.max = 0.75)
# dev.off()

# plot dotplot:
data.plot=data.frame(res_all$gene, res_all$disease, res_all$Rho, res_all$mean.expression, "category"=res_all$category)
data.plot=data.plot[data.plot$res_all.gene%in%res_all.filt$gene,]
colnames(data.plot)=c("features", "id", "variable.1", "variable.2", "category")
data.plot$id=factor(data.plot$id,  levels=c("DLBCL","CHL", "FL", "MCL", "CLL", "CML", "AML", "MDS", "pre-B-ALL", "T-ALL"))

d=data.plot[data.plot$features%in%monoc,]
d=d[order(match(d$features, monoc)),]
d$features=factor(as.character(d$features), levels = monoc)


pdf("Figure2B.pdf", height = 4, width = 4.25)
plot.DotPlot.df(data.plot = d, name.variable.1 = "Rho", name.variable.2 = "Median expression", cols = c("blue", "white", "red"), dot.scale = 3, col.min = -0.5, col.max = 0.75)
dev.off()


#*********************** Plot genes in normal cells ************************
normals2=get.logical(annovector = list(annot$immunoNormals),filterv = !annot$immunoNormals=="")
normals2=normals2[!names(normals2)%in%c("TcellActivatedMononuclear", "LangerhansCell", "Gamma-delta-Tcell", "AlveolarMacrophage")]

normals=get.logical(annovector = list(annot$plotNormals),filterv = !annot$plotNormals=="")

normals=append(normals, get.logical(annovector = list(annot$immunoNormals),filterv = grepl("Macrophage", annot$immunoNormals)))
normals=normals[!names(normals)%in%c("PBMC", "Langerhans Cell", "Langerhans cell",  "Gamma-delta-Tcell", "Macrophage", "AlveolarMacrophage")]

normals=append(normals, normals2[names(normals2)%in%"CD4+Tcell"])

median_res=do.call(cbind, mclapply(seq(normals), function(i){
  logv=normals[[i]]
  NAME=names(normals)[i]
  
  FC=apply(data[,], 1, function(v){
    median(v[logv])
  })
  
}, mc.cores=10))

colnames(median_res)=names(normals)

subclass2=get.logical(annovector = list(annot$subclasses), filterv = annot$CELLS_SORTED==1&annot$Sample.type%in%c("Prolif", "Cancer"))
l=c("AML", "BCL", "CLL", "CML", "pre-B-ALL", "Prolif_Myeloproliferative_MDS", "T-ALL")
logicalvectors.sorted=subclass2[names(subclass2)%in%l]
names(logicalvectors.sorted)=gsub("^BCL_|Lymphoma_BCL_|Prolif_Myeloproliferative_", "", names(logicalvectors.sorted))

median_res_sorted=do.call(cbind, mclapply(seq(logicalvectors.sorted), function(i){
  logv=logicalvectors.sorted[[i]]
  NAME=names(logicalvectors.sorted)[i]
  
  FC=apply(data[,], 1, function(v){
    median(v[logv])
  })
  
}, mc.cores=10))

colnames(median_res_sorted)=names(logicalvectors.sorted)

# picked, plot
d=t(scale(t(median_res[match(monoc, rownames(median_res)),])))
d[d>2]=2
d[d<(-2)]=2

pdf("Figre2B_NormalCells.pdf", width = 4, height = 7)
Heatmap(d, cluster_rows = F)
Heatmap(median_res[match(monoc, rownames(median_res)),], cluster_rows = F)
Heatmap(median_res_sorted[match(monoc, rownames(median_res_sorted)),], cluster_rows = F)
dev.off()

d=t(scale(t(median_res[match(TNK.genes, rownames(median_res_sorted)),])))
d[d>2]=2
d[d<(-2)]=2

# pdf("Fig2_normal_cells_allTNK_genes.pdf", width = 4, height = 15)
# Heatmap(d, cluster_rows = F)
# Heatmap(median_res[match(TNK.genes, rownames(median_res)),], cluster_rows = F)
# Heatmap(median_res_sorted[match(TNK.genes, rownames(median_res_sorted)),], cluster_rows = F)
# dev.off()
