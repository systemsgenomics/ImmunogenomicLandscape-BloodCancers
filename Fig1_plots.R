source("/research/users/ppolonen/git_home/common_scripts/featurematrix/functions_generate_fm.R")
source("/research/users/ppolonen/git_home/common_scripts/visualisation/plotting_functions.R")

# WD
setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

# gexp data
data=get(load("data9544_with_gene_symbols.RData"))

# annotations
annot = get(load("Hemap_immunology_Annotations.Rdata"))
data=data[rownames(data)%in%annot[,1],]

# make logical vectors
annot$tbLY[annot$tbLY%in%c("Lymphoma_BCL_DLBCL_GCB", "Lymphoma_BCL_DLBCL_ABC")]="Lymphoma_BCL_DLBCL"
tbly=get.logical(list(annot$tbLY), filterv = annot$CELLS_SORTED==0&!grepl("TCL", annot$tbLY)&!annot$CLASS2%in%"Cancer_Lymphoma_BCL_DLBCL_testicular")
testicular=list("DLBCL_testicular"=annot$CLASS2%in%"Cancer_Lymphoma_BCL_DLBCL_testicular")
tbly=tbly[sapply(tbly, sum)>5]

tbly=append(tbly, testicular)

lv=get.logical(list(annot$colorClass), filterv = annot$CELLS_SORTED==0&!grepl("TCL|T-ALL", annot$colorClass)&annot$Sample.type%in%c("Cancer", "Prolif"))
names(lv)[3]="MDS"

# barplots
v=sapply(lv, function(lv2){
  nr.high=sum(annot$CytolyticScore.1[lv2]=="high")
  nr_samples=sum(lv2)
  nr.high/nr_samples
})
v=sort(v, decreasing = T)
v=signif(v*100,2)

df2=data.frame(y=sort(v, decreasing = T), x=factor(names(v), levels = names(v)))

# barplots Figure1:

pdf("Fig1F.pdf",width = unit(4, "cm"),height = unit(3.5, "cm"))
plot.boxplot("CytolyticScore", logicalVectors = lv, data = t(data.frame("CytolyticScore"=annot$CytolyticScore)), order.bl = T, intercept.y = 7.76, spread = F,outlier.size = 0.5,  color.v = c("#c08e94","#4ec173", "#f9c155", "#cd7bc1", "#266c34","#8a4bac"))
dev.off()

pdf("Fig1G.pdf", width = unit(3, "cm"),height = unit(2, "cm"))
ggplot(data=df2, aes(x=x, y=v, fill=x)) +
  geom_bar(stat="identity", position=position_dodge())+scale_fill_manual(values = c("#8a4bac", "#c08e94","#f9c155", "#266c34", "#4ec173",  "#cd7bc1"))+
  geom_text(aes(label=df2$y), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# barplots
v=sapply(tbly, function(lv2){
  nr.high=sum(annot$CytolyticScore.1[lv2]=="high")
  nr_samples=sum(lv2)
  nr.high/nr_samples
})
v=sort(v, decreasing = T)

names(v)=gsub("Lymphoma_BCL_", "", names(v))
v=signif(v*100,2)
df2=data.frame(y=sort(v, decreasing = T), x=factor(names(v), levels = names(v)))

pdf("Fig1H.pdf", width = unit(3, "cm"),height = unit(2.5, "cm"))
ggplot(data=df2, aes(x=x, y=y, fill="#c08e94")) +scale_fill_manual(values = c("#c08e94"))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=df2$y), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()
