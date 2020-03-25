GIT_HOME="/research/users/ppolonen/git_home/common_scripts"
source(file.path(GIT_HOME, "scRNA/functions.scRNA.analysis.R"))
source(file.path(GIT_HOME, "visualisation/plotting_functions.R"))
library(Seurat)

setwd("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/HEMAP_IMMUNOLOGY/Published_data_figures")

load("FIMM_AML_scRNA.Rdata")
FIMM=scmat

load("HCA_scRNA.Rdata")
HCA=scmat

load("AML_Galen_scRNA.Rdata")
GALEN=scmat
GALEN[["batch"]][,1]=gsub("AML|.D0", "", GALEN[["batch"]][,1])

colors.group=data.table::fread("colors_lineage.txt", data.table = F, header = F)

# plot proportion of different T-cell types
prop1=data.frame(prop.table(table(FIMM[["SingleR.label"]][,1], FIMM[["batch"]][,1]), margin=2))
prop2=data.frame(prop.table(table(GALEN[["SingleR.label"]][,1], GALEN[["batch"]][,1]), margin=2))
prop3=data.frame(prop.table(table(HCA[["SingleR.label"]][,1], HCA[["batch"]][,1]), margin=2))
prop4=data.frame(prop.table(table(HCA[["SingleR.label"]][,1], rep("HCA", dim(HCA)[2])), margin=2))

# plot lineage to map:
s.anno=data.table::fread("/research/groups/sysgen/PROJECTS/HEMAP_IMMUNOLOGY/petri_work/scRNA/samples_scRNA_allData.txt", data.table=F)[c(1,3,4),]
s.anno[is.na(s.anno)]=""

# plot MDS-like samples to map:
FIMM[["MDS"]]=ifelse(grepl("5897|3667|5249", FIMM[["batch"]][,1]), "MDS-like", "other")

plot.scatter.seurat(scmat = FIMM, colors.group = data.frame("V1"=c("MDS-like", "other"), "V2"=c("darkgreen", "orange"), stringsAsFactors = F), seurat.feature = "MDS", name="Figure3G.pdf", cores=1, SIZE = 0.25, rasterize=F, width = 64*6, height = 74*6, text.size=10*6)
plot.scatter.seurat(scmat = FIMM, colors.group = colors.group[colors.group[,1]%in%prop1$Var1[prop1$Freq>0.001],],seurat.feature = "SingleR.label", name="Figure3F.pdf", cores=1, SIZE = 0.25, rasterize=F, width = 64*6, height = 74*6, text.size=10*6)

cells=gsub(":.*.", "", levels(Idents(FIMM)))

pdf("FigureS3K.pdf", height = 8, width = 9)
DimPlot(FIMM, cols = colors.group[match(cells, colors.group[,1]),2], label = T)
DimPlot(FIMM, cols = colors.group[match(cells, colors.group[,1]),2], label = F)
dev.off()

plot.proportion=function(proportions, celltype, sample.order=NULL, colors.group, max.y=1){
  proportions.cd=proportions[grepl(celltype, proportions$Var1),]
  
  if(!is.null(sample.order))proportions.cd$Var2=factor(proportions.cd$Var2, levels=sample.order)
  proportions.cd$Var1=factor(proportions.cd$Var1, levels=unique(as.character(proportions.cd$Var1)))
  myColors=colors.group[match(levels(proportions.cd$Var1), colors.group$V1),2]
  names(myColors) <- levels(proportions.cd$Var1)
  colScale2 <- scale_fill_manual(values = myColors)
  p2=ggplot(data=proportions.cd,
            aes(x=Var2, y=Freq, fill=Var1)) + colScale2 +
    geom_bar(stat = "identity", inherit.aes = T) +
    ggtitle(paste(" ", " ", sep="\n")) +
    theme(legend.position = "none", legend.direction = "horizontal",
          axis.line.y = element_line(size=1, colour = "black"),
          axis.line.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(colour="grey20",size=10,face="plain", family="Helvetica", angle = 90, hjust = 1),  
          axis.text.y = element_text(colour="grey20",size=10,face="plain", family="Helvetica"),  
          axis.ticks.length = unit(2,"mm"),
          panel.border = element_blank(), panel.background = element_blank(),
          plot.title = element_text(size = 10, face = "bold", family="Helvetica"),
          text=element_text())+
    theme(plot.margin=unit(c(1,1,1,1), "mm"))+ coord_cartesian(ylim=c(0,max.y))
  return(p2)
}

p1=plot.proportion(proportions = prop1, celltype = "HSC|GMP|CMP|MPP|MEP|Monocyte", 
                sample.order=c("5897", "3667","5249", "6333","5238","5750", "6386", "6187"),
                colors.group=colors.group)

p2=plot.proportion(proportions = prop2, celltype = "HSC|GMP|CMP|MPP|MEP|Monocyte", 
                   colors.group=colors.group)

p3=plot.proportion(proportions = prop3, celltype = "HSC|GMP|CMP|MPP|MEP|Monocyte", 
                   colors.group=colors.group)

p4=plot.proportion(proportions = prop4, celltype = "HSC|GMP|CMP|MPP|MEP|Monocyte", 
                   colors.group=colors.group)

ggsave(plot = p1, filename = "Figure3H_blast_proportion.pdf", width = 1.5, height = 2)
# ggsave(plot = p2, filename = "Blast_proportion_Galen_AML.pdf", width = 6, height = 2)
# ggsave(plot = p3, filename = "Blast_proportion_HCA_samples.pdf", width = 6, height = 2)
ggsave(plot = p4, filename = "Figure3H_blast_proportion_HCA.pdf", width = 2, height = 2)


p1=plot.proportion(proportions = prop1, celltype = "NK|CD8. Tem", 
                   sample.order=c("5897", "3667","5249", "6333","5238","5750", "6386", "6187"),
                   colors.group=colors.group, max.y = 1)

p2=plot.proportion(proportions = prop2, celltype = "NK|CD8. Tem", 
                   colors.group=colors.group, max.y = 1)

p3=plot.proportion(proportions = prop3, celltype = "NK|CD8. Tem", 
                   colors.group=colors.group, max.y = 1)

p4=plot.proportion(proportions = prop4, celltype = "NK|CD8. Tem", 
                   colors.group=colors.group, max.y = 1)

ggsave(plot = p1, filename = "Figure3H_CTL_NK_proportion.pdf", width = 1.5, height = 2)
# ggsave(plot = p2, filename = "TNK_proportion_Galen_AML.pdf", width = 6, height = 2)
# ggsave(plot = p3, filename = "TNK_proportion_HCA_samples.pdf", width = 6, height = 2)
ggsave(plot = p4, filename = "Figure3H_CTL_NK_proportion_HCA.pdf", width = 2, height = 2)

#************************* Integrated analysis, these were plotted in the Figure3K *************************
load("markers.all.filt.allgenes.Rdata") # from above

genes.up=fread("TableS3_Significant_genes_MDSlike_bulk.txt", data.table = F)

MDSup=markers.all.filt[markers.all.filt$gene%in%genes.up&markers.all.filt$avg_logFC>0.25,]

# interesting candidates from cellphonedb
feats=c("KLRF1", "CLEC2B", "NCR3", "BAG6", "HAVCR2", "LGALS9", "TGFBR3", "TGFB1", "IFNGR1", "IFNGR2", "IFNG", "LAG3")

# Check if significant in different category:
blasts=feats%in%markers.all.filt$gene[grepl("MDS-like",markers.all.filt$cluster)&grepl("HSC|MPP|MEP|GMP|CMP|Monocyte",markers.all.filt$cluster)]
NK=feats%in%markers.all.filt$gene[grepl("MDS-like",markers.all.filt$cluster)&grepl("CD8|NK",markers.all.filt$cluster)]
bulk=feats%in%genes.up$Gene[genes.up$FDR<0.001]
data.frame(feats, blasts, NK, bulk)

pdf("Figure3K_CellPhoneDB.pdf", width = 25, height = 6)
FeaturePlot(FIMM, c("KLRF1", "CLEC2B"), blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001, pt.size = 0.01, coord.fixed = F, sort.cell = T)
FeaturePlot(FIMM, c("NCR3", "BAG6"), blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001, pt.size = 0.01, coord.fixed = F, sort.cell = T)
FeaturePlot(FIMM, c("HAVCR2", "LGALS9"), blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001, pt.size = 0.01, coord.fixed = F, sort.cell = T)
FeaturePlot(FIMM, c("TGFBR3", "TGFB1"), blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001, pt.size = 0.01, coord.fixed = F, sort.cell = T)
dev.off()

#************************************** Analysis of IFNG *******************************************
load("Combined_pathway_signatures_2017_filtered_robust.Rdata")
add.scores=listA[grep("INTERFERON_GAMMA_SIGNALING", names(listA))]

add.scores=append(add.scores, list("IFNG_Receptor"=c("IFNGR1", "IFNGR2")))

fm.f=GetAssayData(FIMM)

gm.objects=do.call(rbind, lapply(seq(add.scores), function(i){
  dat3=fm.f[rownames(fm.f)%in%add.scores[[i]],]
  gm=log2(t(apply(dat3, 2, gm_mean))) # done to normalized values
  rownames(gm)=names(add.scores)[i]
  return(gm)
}))

# also add to seurat object:
for(i in seq(add.scores)){
  FIMM[[names(add.scores)[i]]] <- gm.objects[i,]
}

pdf("Fig3K_IFNG_analysis.pdf", width = 25, height = 6)
FeaturePlot(FIMM, c("IFNG", "IFNG_Receptor"), sort.cell = T, pt.size = 0.01, blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001)
FeaturePlot(FIMM, c("IFNG_Receptor", "IFNG"), sort.cell = T, pt.size = 0.01, blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001)
FeaturePlot(FIMM, c("INTERFERON_GAMMA_SIGNALING-REACTOME_MsigDB_c2", "IFNG"), sort.cell = T, pt.size = 0.01, blend = T, cols=c("red", "darkblue"), blend.threshold = 0.001)
plot.DotPlot(FIMM, features = unique(c("IFNGR1", "IFNGR2","IFNG", "IFNG_Receptor", "INTERFERON_GAMMA_SIGNALING-REACTOME_MsigDB_c2", "CytolyticScore")))
dev.off()

pdf("FigS3O_IFNG_analysis_vln.pdf", width = 5, height = 3)
FIMM[["MDSlike"]]=ifelse(grepl("5897|3667|5249", FIMM[["batch"]][,1]), "MDS-like", "other")
FIMM[["MDSlike"]][grepl("6333|6187", FIMM[["batch"]][,1]),1]="HLAII_CIITA_low"
FIMM[["batch2"]]=factor(FIMM[["batch"]][,1], levels=c("5897", "3667", "5249", "5238", "5750", "6187", "6333", "6386"))
VlnPlot(FIMM[,grepl("MPP|GMP|CMP|MEP|HSC", FIMM[["SingleR.label"]][,1])], c("IFNG_Receptor"), group.by = "MDSlike", log = F, pt.size = 0.1)
VlnPlot(FIMM[,grepl("MPP|GMP|CMP|MEP|HSC", FIMM[["SingleR.label"]][,1])], c("INTERFERON_GAMMA_SIGNALING-REACTOME_MsigDB_c2"), group.by = "MDSlike", log = F, pt.size = 0)
VlnPlot(FIMM[,grepl("MPP|GMP|CMP|MEP|HSC", FIMM[["SingleR.label"]][,1])], c("INTERFERON_GAMMA_SIGNALING-REACTOME_MsigDB_c2"), group.by = "batch2",cols = c(rep("darkgreen", 3), rep("orange", 5)), log = F, pt.size = 0)
dev.off()