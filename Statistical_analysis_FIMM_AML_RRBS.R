library(RnBeads)
# Set some analysis options
rnb.options(
  filtering.sex.chromosomes.removal=T,
  identifiers.column="SampleID",
  replicate.id.column="treatment",
  import.bed.style="bismarkCov",
  enforce.memory.management=T,
  assembly="hg38",
  differential.report.sites=FALSE,
  filtering.sex.chromosomes.removal = TRUE,
  import.table.separator="\t"
)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

options(fftempdir=file.path(getwd(), "temp"))

if(!dir.exists(file.path(getwd(), "temp")))dir.create(file.path(getwd(), "temp"))
rnb.set=load.rnb.set("analysis.zip", temp.dir=file.path(getwd(), "temp"))
diffmeth=load.rnb.diffmeth("analysis/")

# laad custom annot:
rnb.load.annotation("CpGregion_annotations.Rdata", "CpGregion")
rnb.load.annotation("CpGregion_CIITA_annotations.Rdata", "CpGregion_CIITA")
rnb.load.annotation("CpGexpanded_annotations.Rdata", "CpG.expanded")

#****************** volcanoplot to see differential methylation for each comparison: ******************
library(EnhancedVolcano)
# type="promoters"
# type="tiling"
# type="genes"
# type="cpgislands"
# type="CpG.expanded"

type="CpGregion"

comparison <- get.comparisons(diffmeth)[1]

tab.sites <- get.table(diffmeth, comparison, type, return.data.frame=TRUE)

which.diffmeth <- abs(tab.sites$mean.mean.diff)>0.1 & tab.sites$comb.p.adj.fdr<0.05

tab.sites$significant=which.diffmeth

lab = paste(tab.sites$Chromosome, tab.sites$Start, tab.sites$End, sep=":")

aa <- annotation(object = rnb.set, type = type)
annotated.dmrs <- data.frame(aa, tab.sites)

annotated.dmrs[grep("CpG:41.91|CpG:21.172|CpG:20.179",rownames(annotated.dmrs)),]

pdf("Fig4G_CIITA_CpG_Volcanoplot.pdf", width=5, height=5)
EnhancedVolcano(tab.sites, title = gsub(" \\(.*.", "", comparison), subtitle = type, xlab = "mean difference",
                     # lab = paste(tab.sites$Chromosome, tab.sites$Start, tab.sites$End, sep=":"), drawConnectors=T,
                     lab =rownames(annotated.dmrs), drawConnectors=F,widthConnectors = 1,
                     # selectLab = grep("CpG:41.91",rownames(annotated.dmrs), value=T),
                     x = 'mean.mean.diff', ylim=c(0,5),
                     y = 'comb.p.adj.fdr', col=c("grey75", "grey75","grey75","grey75"),
                     pointSize = 2,vline = -0.058,vlineCol = "red",vlineType = "solid",
                     pCutoff = 0, FCcutoff = 0)

EnhancedVolcano(tab.sites[grepl("CpG:41.91|CpG:21.172|CpG:20.179",rownames(annotated.dmrs)),], title = gsub(" \\(.*.", "", comparison), subtitle = type, xlab = "mean difference",
                     # lab = paste(tab.sites$Chromosome, tab.sites$Start, tab.sites$End, sep=":"), drawConnectors=T,
                     lab =rownames(annotated.dmrs)[grepl("CpG:41.91|CpG:21.172|CpG:20.179",rownames(annotated.dmrs))], drawConnectors=F, widthConnectors = 1,
                     selectLab = grep("CpG:41.91",rownames(annotated.dmrs), value=T),
                     x = 'mean.mean.diff', pointSize = 3,vline = -0.058,vlineCol = "red",vlineType = "solid",
                     y = 'comb.p.adj.fdr', col=c("red2", "red2","red2","red2"),
                     pCutoff = 0, FCcutoff = 0, xlim = c(-1,1), ylim=c(0,5))
dev.off()

# 
# # Check some global
# rnb.set.m <- meth(rnb.set, type="cpgislands")
# 
# rnb.set.m2=t(rnb.set.m[!rowSums(is.na(rnb.set.m))>20,])
# 
# d=do.call(rbind, lapply(1:24, function(i)prop.table(table(rnb.set.m2[i,]<0.1))))
# rownames(d)=rnb.set@pheno$treatment
# 
# d2=do.call(rbind, lapply(1:24, function(i)prop.table(table(rnb.set.m2[i,]>0.9))))
# rownames(d2)=rnb.set@pheno$treatment
# 
# res.cor <- cor(t(rnb.set.m2), method = "pearson", use = "complete.obs")
# 
# # Load required packages
# library(magrittr)
# library(dplyr)
# library(ggpubr)
# # Cmpute MDS
# mds <- rnb.set.m2 %>%
#   dist() %>%          
#   cmdscale() %>%
#   as_tibble()
# colnames(mds) <- c("Dim.1", "Dim.2")
# # Plot MDS
# 
# pdf("MDS.pdf", width = 8, height = 8)
# ggscatter(mds, x = "Dim.1", y = "Dim.2", 
#           label = rownames(rnb.set.m2),
#           size = 1,
#           repel = TRUE)
# dev.off()
# 
# res.cor <- 1-cor(t(rnb.set.m2), method = "pearson", use = "complete.obs")

dist.m=dist(rnb.set.m2, method = "euclidean")

fit <- hclust(dist.m, method="ward")

pdf("FigS4L_Euclidean.pdf")
plot(fit) # display dendogram
dev.off()

# CIITA CpG island
rnb.set.m <- meth(rnb.set, type="CpGregion_CIITA")
aa <- annotation(rnb.set, type="CpGregion_CIITA")

rownames(rnb.set.m)=rownames(aa)

comparison <- get.comparisons(diffmeth)[1]
tab.promoters <- get.table(diffmeth, comparison, "CpGregion_CIITA", return.data.frame=TRUE)

plots=lapply(c(2,1,3,5,4), function(i){
  group=gsub("_1|_2|_3","", colnames(rnb.set.m))
  group=ifelse(grepl("DC", colnames(rnb.set.m)),"DAC", "Rest")
  
  dat.list=lapply(unique(group), function(g)rnb.set.m[i,group%in%g])
  names(dat.list)=paste(unique(group), rownames(aa)[i])
  return(dat.list)
})
plots=unlist(plots, recursive = F)


a=plot.boxplot.list(plots, color.v = rep("grey75",length(plots)), ylab = "Methylation beta", spread = T, outlier.size = 1, y.range = c(0,1))

ggsave(plot = a, filename = "Fig4H_boxplot_CIITA_CpG.pdf",device = "pdf", height = 4, width = 4)

ind=4

# group=gsub("_1|_2|_3","", colnames(rnb.set.m))
# group=ifelse(grepl("DC", colnames(rnb.set.m)),"DAC", "Rest")
# 
# dat.list=lapply(unique(group), function(g)rnb.set.m[ind,group%in%g])
# names(dat.list)=unique(group)
# a=plot.boxplot.list(dat.list, color.v = rep("grey75",length(dat.list)), ylab = "Methylation beta", spread = T, outlier.size = 1.5, y.range = c(0,1))
# 
# ggsave(plot = a, filename = "boxplot_CIITA_CpG.pdf",device = "pdf", height = 3, width = 3)
# 
# write.table(rnb.set.m, "CIITA_CpGisland.txt", quote = F, row.names = T, col.names = T, sep="\t")
#
# # CIITA promoter
# rnb.set.m <- meth(rnb.set, type="promoters")
# aa <- annotation(rnb.set, type="promoters")
# 
# ind=aa$symbol%in%"CIITA"
# group=gsub("_1|_2|_3","", colnames(rnb.set.m))
# group=ifelse(grepl("DC", colnames(rnb.set.m)),"DAC", "Rest")
# 
# dat.list=lapply(unique(group), function(g)rnb.set.m[ind,group%in%g])
# names(dat.list)=unique(group)
# plot.boxplot.list(dat.list, color.v = rep("grey75",length(dat.list)), ylab = "Methylation beta")
# 
# 
# # region analysis
# type="CpGregion"
# rnb.set.m <- meth(rnb.set, type="CpGregion")
# comparison <- get.comparisons(diffmeth)[1]
# tab.sites <- get.table(diffmeth, comparison, type, return.data.frame=TRUE)
# which.diffmeth <- abs(tab.sites$mean.mean.diff)>0.05 & tab.sites$comb.p.adj.fdr<0.01
# 
# aa <- annotation(object = rnb.set, type = type)
# annotated.dmrs <- data.frame(aa, tab.sites)
# annotated.dmrs$significant=which.diffmeth
# 
# save(list = c("annotated.dmrs", "rnb.set.m"), file="annotated.dmrs.CpGregion.Rdata")
# 
# # region analysis
# type="promoters"
# comparison <- get.comparisons(diffmeth)[1]
# rnb.set.m <- meth(rnb.set, type="promoters")
# tab.sites <- get.table(diffmeth, comparison, type, return.data.frame=TRUE)
# which.diffmeth <- abs(tab.sites$mean.mean.diff)>0.05 & tab.sites$comb.p.adj.fdr<0.01
# 
# aa <- annotation(object = rnb.set, type = type)
# annotated.dmrs <- data.frame(aa, tab.sites)
# annotated.dmrs$significant=which.diffmeth
# 
# save(list = c("annotated.dmrs", "rnb.set.m"), file="annotated.dmrs.promoters.Rdata")
# 
# 
# region=cbind(as.character(annotated.dmrs$Chromosome), annotated.dmrs$Start, annotated.dmrs$End, rownames(annotated.dmrs))
# region=region[which.diffmeth,]
# region[,4]=gsub("@.*.", "", region[,4])
# 
# region=region[!duplicated(region[,4]),]
# 
# write.table(region[,4], "significant_CpG.bed", quote=F, sep="\t", row.names = F, col.names = F)
